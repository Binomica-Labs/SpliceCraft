# Hosting the SpliceCraft web demo (splicecraft.bio)

The web demo is **the real `splicecraft.py`**, run server-side and streamed to a
browser terminal via **[`textual-serve`](https://github.com/Textualize/textual-serve)**.
There is no separate web codebase, so the browser demo is **1-to-1 with the
terminal tool by construction** — same braille maps, same theme, same behaviour.
What makes it safe + public-ready is **demo mode** (`SPLICECRAFT_DEMO`) plus the
in-app web-tier lockdown.

> **Not `textual-web`.** Textualize's hosted `textual-web.io` relay was
> decommissioned (~2025 — the domain no longer resolves), so the old "tunnel out
> to a relay + `--signup`" route is a dead end. We self-host with `textual-serve`
> (its supported replacement) behind Caddy — cleaner anyway: the demo runs on our
> own domain with our own TLS and no third-party dependency.

## Architecture

```
splicecraft.bio / www  ──▶ Caddy (file_server) ─────────▶ /var/www/splicecraft/index.html   (static landing page)
demo.splicecraft.bio   ──▶ Caddy (reverse_proxy, TLS) ──▶ textual-serve :8000
                                                            └─▶ spawns `SPLICECRAFT_DEMO=web splicecraft` per session
```

- **Apex** (`splicecraft.bio` / `www`) serves the marketing **landing page**
  (`web/index.html` in this repo).
- **`demo.splicecraft.bio`** is the live demo. textual-serve must own a domain
  root, so the demo gets its own subdomain; the landing page's "Try the live
  demo" button links there.

## The two tiers

| `SPLICECRAFT_DEMO` | Use | Behaviour |
|---|---|---|
| `local` | previewing the demo build on your own machine | **full features**, sandboxed data dir + seeded examples |
| `web`   | the public deployment | **restrained** — see lockdown below |

Both force a **fresh ephemeral data dir per process** (`/tmp/splicecraft-demo-…`)
and ignore `$SPLICECRAFT_DATA_DIR` / `$XDG`, so a session can never read or write a
real library. Each launch seeds the FFE Golden-Braid / MoClo worked set (see
**Seed**) and shows a version banner.

## What the `web` tier locks (defence-in-depth, in-app)

Enforced in-process, regardless of UI path:

- **No network egress** — every fetch helper (`fetch_genbank`, PyPI update check,
  NCBI taxid, Kazusa, genome datasets, HMM-DB / arbitrary-URL opener) raises, so
  no SSRF / outbound abuse.
- **No host-FS reads** — open-file, load-part-from-file, entry-vector-from-file,
  image-attach, sequencing import, etc. are all gated.
- **No agent HTTP API** — `_start_agent_api` refuses even if `--agent` is passed.
- **No destructive / bulk-data ops** — data export/import (migrate), master
  delete, and the per-format exports (GenBank / FASTA / GFF / CommercialSaaS) show
  a polite "install SpliceCraft for the full tool" message.
- **Sequence-size cap** — `_apply_record` refuses sequences > 100 kb.

The **science tools stay fully usable** on the seed (Scrub, operon domestication,
primer design, Simulator, BLAST, codon optimisation) — pure compute, not gated.

## Deploying (DigitalOcean + textual-serve + Caddy)

A small always-on Linux box runs textual-serve; Caddy fronts it with automatic
HTTPS. An *idle* session costs ~0 % CPU and ~120 MB RAM, so 1 vCPU / 2 GB is
plenty for a demo (≈10 concurrent before memory is the ceiling).

### 1 — Droplet
DO console → **Ubuntu 24.04 LTS**, Basic **$12/mo (2 GB)**. Enable Monitoring +
IPv6 (free). Add your SSH key, note `<DROPLET_IP>`.

### 2 — Base + install
```bash
ssh root@<DROPLET_IP>
adduser --disabled-password --gecos "" demo
apt-get update && apt-get install -y pipx python3-venv
ufw allow OpenSSH && ufw allow 80/tcp && ufw allow 443/tcp && ufw --force enable
su - demo -c 'pipx install splicecraft'
python3 -m venv /home/demo/serve-venv && /home/demo/serve-venv/bin/pip install textual-serve
fallocate -l 2G /swapfile && chmod 600 /swapfile && mkswap /swapfile && swapon /swapfile
echo '/swapfile none swap sw 0 0' >> /etc/fstab
```
`splicecraft` and `textual-serve` live in separate envs — textual-serve only
*spawns* the `splicecraft` CLI, so their deps never entangle.

### 3 — Launcher: `/home/demo/serve_splicecraft.py`
```python
from textual_serve.server import Server

server = Server(
    # `exec` is load-bearing (see gotchas)
    # SPLICECRAFT_MAP_ASPECT pins the circular-map cell ratio — xterm.js
    # doesn't answer the CSI 16t self-query, so the web needs the explicit
    # pin or the map renders as a tall ellipse (see gotchas).
    command="exec env SPLICECRAFT_DEMO=web SPLICECRAFT_MAP_ASPECT=2.35 PATH=/home/demo/.local/bin:/usr/local/bin:/usr/bin:/bin /home/demo/.local/bin/splicecraft",
    host="127.0.0.1",                              # behind Caddy; never exposed directly
    port=8000,
    title="SpliceCraft — live demo",
    public_url="https://demo.splicecraft.bio",     # emit wss://demo…/ws, not 127.0.0.1
    templates_path="/home/demo/sc-templates",      # custom template → quit redirects (gotchas)
)
server.serve()
```

### 4 — systemd: `/etc/systemd/system/splicecraft-demo.service`
```ini
[Unit]
Description=SpliceCraft web demo (textual-serve)
After=network-online.target
Wants=network-online.target

[Service]
User=demo
WorkingDirectory=/home/demo
ExecStart=/home/demo/serve-venv/bin/python /home/demo/serve_splicecraft.py
Restart=always
RestartSec=3
MemoryMax=1400M
NoNewPrivileges=true
PrivateTmp=true

[Install]
WantedBy=multi-user.target
```
`systemctl enable --now splicecraft-demo`.

### 5 — Caddy: `/etc/caddy/Caddyfile`
```caddy
{
    email you@example.com          # REQUIRED — see gotchas
}

splicecraft.bio, www.splicecraft.bio {
    encode zstd gzip
    root * /var/www/splicecraft
    file_server
}

demo.splicecraft.bio {
    encode zstd gzip
    reverse_proxy 127.0.0.1:8000
}
```
Copy `web/index.html` → `/var/www/splicecraft/index.html`, then
`systemctl reload caddy`. Caddy auto-issues Let's Encrypt certs once DNS points
at the box and 80/443 are open.

### 6 — DNS (at the registrar, e.g. GoDaddy)
| Type | Name | Value |
|---|---|---|
| A | `@` | `<DROPLET_IP>` |
| A | `demo` | `<DROPLET_IP>` |
| CNAME | `www` | `splicecraft.bio` |

## Load-bearing gotchas (each was a real bug)

- **`exec` in the launch command.** textual-serve runs `sh -c "<cmd>"`; without
  `exec`, `splicecraft` is a *grandchild* of that shell. On disconnect the shell
  dies but splicecraft is orphaned (reparented to init) and **spins ~90 % CPU
  forever**. `exec` makes splicecraft replace the shell so it's killed cleanly.
- **Orphan reaper.** Belt-and-suspenders cron — `/etc/cron.d/sc-reap` running
  `*/2 * * * * root` a one-liner that `kill -9`s any `ppid==1` splicecraft (a hard
  network drop can still orphan one).
- **Terminal font size.** SpliceCraft requires ≥ 100×30. textual-serve's 16 px
  default web terminal can dip under 100 cols on small/split windows. Patch
  `textual_serve/server.py`'s `font_size` default 16 → 13 (or append
  `?fontsize=13`). A nightly `splicecraft` upgrade does **not** touch textual-serve,
  so the patch persists.
- **Valid ACME `email` in the Caddyfile.** With none, Caddy tries to register a
  Let's Encrypt account with the literal contact `default` → LE rejects it
  (`invalidContact`) and **new** certs fail (already-cached certs keep working,
  which hides the cause). Set a real email in the global block.
- **`public_url`.** Behind Caddy/HTTPS this makes textual-serve emit
  `wss://demo.splicecraft.bio/ws` + `https://…/static` instead of `127.0.0.1:8000`.
- **Quit → website.** `templates_path` points at a copy of textual-serve's
  templates (Jinja vars preserved) whose `app_index.html` adds a `MutationObserver`
  that runs `window.location.replace("https://splicecraft.bio")` when the session
  closes (`document.body` gains `-closed`).
- **Stale DNS after repointing.** Repointing the apex away from a parked page
  leaves local caches serving the old site — `resolvectl flush-caches` and clear
  the browser's `chrome://net-internals/#dns`.

## Seed

A demo session is seeded offline with the FFE Golden-Braid / MoClo worked set
from **`splicecraft_demo_plasmids.py`** (bundled in the wheel, loaded via
`json.loads`): L0 entry vectors (UPD + α/ω), a "Demo parts" bin of L0 parts, two
transcription units, and the cscA/cscB + pCambia plasmids — so Constructor /
Parts Bin / Scrub all have worked examples. Regenerate the sidecar from a source
library; never hand-edit. **Mind the private-name guard** — landing-page /
example copy must avoid `.private-names` tokens (a guarded reporter-gene name slipped into the landing-page session text and tripped `release.py` once).

## Staleguard — keeping the demo 1-to-1

The demo serves whatever `pipx` has, so "never behind the terminal build" has two
layers. **Primary: `release.py` refreshes it on every release** — after the tag is
pushed and PyPI is publishing, `_refresh_web_demo()` waits for the sdist to go
live, then `ssh root@147.182.141.165 'sc-update'` (the droplet IP — connecting by
the domain trips host-key verification under BatchMode) and verifies `version.txt` +
a 200 from the demo. It's **non-fatal** (the release is already shipped; a failure
just prints the manual command) and opt-out via `SPLICECRAFT_SKIP_DEMO_REFRESH=1`
(offline, or a laptop without the droplet's SSH key). Override the SSH target with
`SPLICECRAFT_DEMO_HOST`.

**Backstop: a nightly redeploy** catches anything the release-time refresh missed
(skipped, SSH down, PyPI slow):
```cron
30 4 * * * root /usr/local/bin/sc-update   # pipx upgrade splicecraft + restart demo + refresh /version.txt
```
The script's canonical copy is **[`scripts/sc-update`](../scripts/sc-update)** in
this repo — edit it there and redeploy with
`scp scripts/sc-update root@<droplet>:/usr/local/bin/sc-update`, so a droplet
rebuild can't silently lose a fix. It passes `--no-cache-dir` to pip for a
load-bearing reason: without it `pipx upgrade` reads pip's **cached** copy of the
PyPI index and reports "already at latest `<previous>`" for a version published
minutes earlier, even though the droplet can fetch the new release perfectly
well. That made the release-time refresh no-op on every release and defer to the
nightly cron — which appeared to "fix" it only because the cache had expired by
morning. It looked exactly like PyPI propagation lag and was misdiagnosed as such
until an A/B on the droplet showed the same `pipx upgrade` failing without the
flag and succeeding with it (2026-09-28).

After `release.py` ships X.Y.Z it pulls it automatically (or the next cron run /
a manual `sc-update` does);
the banner shows `v{__version__}` so drift is visible at a glance. `sc-update`
also writes the live version to `/var/www/splicecraft/version.txt`, which the
landing page's nav **version chip** fetches same-origin (allowed by the apex CSP
`connect-src 'self'`) so the chip always shows the deployed version. The in-app
"newer version?" check is off in the web tier (egress is blocked), so the
staleguard lives at the **deploy** layer.

## Live-status light

The landing nav's **"live" dot is a real end-to-end health check**, not a CSS
animation. `/etc/cron.d/sc-health` runs `/usr/local/bin/sc-health` every minute;
it probes `https://demo.splicecraft.bio/` (a genuine HTTP 200 means DNS + Caddy +
TLS + textual-serve are all answering) and writes
`/var/www/splicecraft/demo-status.json`. The landing page fetches that
same-origin (allowed by the apex CSP `connect-src 'self'`) and shows a green
pulse only on `"demo":"up"`, otherwise a red **demo offline** dot. A staleness
guard compares the response's `Date` and `Last-Modified` headers (both
server-provided, so no client-clock dependency); if the file is > 5 min stale —
i.e. the probe cron died — the light fails to **offline** rather than showing
false-live.

## Per-session isolation (if traffic grows)

textual-serve runs one subprocess per visitor (not a container). Today's
protection is **demo-mode sandbox + web lockdown + `MemoryMax` + the orphan
reaper**. If abuse / scale demands more, run each session in a minimal container
(no real data, no host mounts) + cgroup limits + a concurrency cap, behind the
same Caddy.

## Final human QA

`test_demo_mode.py` covers the lockdown contract, that the seed loads (plasmids +
parts bin + entry vectors), and that the core tools run. Still worth a manual
click-through of every modal in both tiers after a deploy.
