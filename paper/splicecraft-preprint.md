# SpliceCraft: a terminal-native plasmid workbench for molecular cloning that runs offline on commodity hardware

**Sebastian S. Cocioba**¹

¹ Binomica Labs, Long Island City, New York, USA
ORCID: [0000-0002-6821-2996](https://orcid.org/0000-0002-6821-2996)
Correspondence: sebastian@binomicalabs.org

**Running title:** SpliceCraft: an offline terminal plasmid workbench

**Keywords:** molecular cloning · plasmid · DNA sequence editor · terminal user interface · open-source software · Golden Gate · MoClo · GoldenBraid · primer design · codon optimization · research software sustainability · low-resource research · science education · local-first software · agentic AI

---

## Abstract

Almost every cloning experiment now passes through software, which makes design the cheapest place to catch a mistake and the costliest to miss one. The tools have split in two: open-source command-line software most bench biologists never reach, and commercial platforms sold as cloud subscriptions needing an account, a browser and an annual fee. Academic pricing is often generous, so money is often not the obstacle. The licence terms are. Two of the three platforms I surveyed state that when a subscription lapses the software reverts to a restricted viewer with most functions disabled: you rent the capability and keep little more than the ability to look.

I present SpliceCraft, a free MIT-licensed plasmid workbench that runs in a terminal, offline. It draws braille plasmid maps, edits sequences, designs primers, simulates Golden Gate, Gibson and restriction assemblies, optimises codons, verifies constructs against sequencing reads, and keeps a lab notebook. Version 1.3.19 is 222,000 lines of pure Python under 12,277 tests, 202 MB installed and 87 MB resident. On a terminal without braille glyphs the map falls back to ASCII rather than refusing to start.

Four commitments shape it. Data integrity outranks every feature: each write passes an atomic write, a four-layer backup chain that once recovered an accidental overwrite of my own 160 MB library, and an authorisation chokepoint added so that recovery is never needed again. Established engines are integrated, not rewritten; the two exceptions were benchmarked first. Every score ships its measured operating point, and models that failed to beat chance are absent, with the numbers that killed them. And every operation works without the interface: 303 endpoints with discoverable schemas let someone describe a construct in plain language to a model of their choosing, with no sequence leaving the machine.

Two questions the software raises get their own sections: why abandonment is the normal fate of scientific software, and what building a tool to be inherited requires. On method, SpliceCraft was built by one bench biologist with continuous large-language-model assistance. The useful part is what went wrong: nine failure modes found in this codebase, including a module that shipped with 38 passing tests and no caller. Generated code fails by looking right, which is the one failure a scientific user cannot catch unaided.

Available at [github.com/Binomica-Labs/SpliceCraft](https://github.com/Binomica-Labs/SpliceCraft), PyPI `splicecraft`, Zenodo concept DOI [10.5281/zenodo.21960400](https://doi.org/10.5281/zenodo.21960400). Runs in a browser at [splicecraft.bio](https://splicecraft.bio).
---

## 1. Introduction

### 1.1 The design step is software

A modern cloning project is a sequence of decisions made in software and then executed at the bench: which enzyme, which overhang, which primer, which codon table, which assembly order. The physical work takes days and costs money; the design work takes minutes and costs nothing. This asymmetry means the software is the cheapest place in the whole workflow to be careful, and the most expensive place to be wrong. A primer that anneals 4 bp from where the display claims it does will consume a week and a round of sequencing before anyone suspects the drawing.

The software that mediates this step has polarised into two populations that serve different people.

The first is open-source, command-line, and library-level: Biopython [1], Primer3 [2], HMMER via PyHMMER [3], BLAST [4], and the long tail of EMBOSS-era utilities. It is powerful, scriptable, inspectable, free, and — as a 2023 review written specifically for wet-lab scientists observes — frequently not reached by the people who would benefit most, who "are not formally trained to apply bioinformatic tools and may therefore assume that they lack the necessary experience to do so themselves" [5]. The barrier is not intelligence; it is that no one hands a bench biologist an integrated workbench made of these parts.

The second population does hand them one, and charges for it. Commercial sequence-design platforms are good at integration and good at interface. The engineering in them is real and this paper does not dispute it. Over the last decade most of them have moved to cloud hosting, per-seat licensing and annual renewal. Table 1 gives the list prices their own pricing pages carried on the day this manuscript was written.

**Table 1.** Published list pricing for three widely used commercial sequence-design platforms, as stated on each vendor's own pricing page. All figures in US dollars per year unless noted. Accessed 2026-11-11: `snapgene.com/pricing`, `geneious.com/pricing`, `benchling.com/pricing` and `benchling.com/academic`. Figures are quoted, not converted or estimated; a fourth vendor's academic page refused automated access and is therefore omitted rather than estimated. The Geneious five-seat figures are not five times the single-seat ones because they are a different product: single seats are a Personal subscription, while any group is a Team subscription priced strictly per user ($685 academic, $2,500 corporate) with no single-seat tier, so a Team seat costs slightly more than a Personal one rather than less. SnapGene's academic ladder is linear to five seats and discounts at ten. Both vendors label the non-academic tier "Corporate".

| | Student | Academic, 1 seat | Academic, 5 seats | Commercial, 1 seat | Commercial, 5 seats |
|---|---|---|---|---|---|
| SnapGene | $149 | $360 | $1,800 | $1,945 | $9,725 |
| Geneious Prime | $200 | $660 | $3,425 | $2,350 | $12,500 |
| Benchling | free academic plan (Notebook + Molecular Biology) | — | — | not published; "Contact us" / "Request a Demo" | not published |

Three observations about this table matter more than the numbers in it, and the first two are concessions.

**The academic terms are often generous.** Benchling's academic plan is free, and has been offered to a stated "more than 200,000 scientists across 7,500+ universities and research institutions." Both SnapGene and Geneious provide course licences for teaching **at no charge to the students**, on condition that the instructor holds an active paid subscription. For a funded graduate student at a degree-granting institution, cost is frequently not the binding constraint at all, and any argument pretending otherwise should be disregarded.

**For commercial users the numbers are, however, simply large.** A five-person company wanting both SnapGene and Geneious Prime pays $22,225 per year, recurring, before a single reagent — and Benchling's commercial pricing is not published, so a startup cannot even budget against it without entering a sales conversation. For the lean-startup reader of §1.3 this is a material fraction of a first-year budget.

**The structurally interesting finding is in the terms, not the prices.** Two things the vendors state plainly:

1. **Expiry reverts the tool to a restricted viewer.** SnapGene "reverts to Viewer mode when your subscription expires. In Viewer mode you can view, annotate and share your sequence files, however most other functions are disabled." Geneious Prime "reverts to restricted mode when the subscription has expired," where you "can view your data and import and export your files, however most other functions are disabled." Both vendors say *most*, not all, and both leave more than reading: annotating and sharing in one case, import and export in the other. What lapses is the analytical capability — the simulation, the design, the search — which is the part that was rented. When a grant ends, the designs remain readable and the ability to do the work stops — which is the §1.4 continuity problem arriving by a different route. The tool has not disappeared; your access to it has.
2. **Buying outright is available but deliberately terminal.** SnapGene's permanent licence "only allows you to access the version of SnapGene that you purchased, with no option for upgrades to future versions," and is "only available in single-seat form." Ownership is offered on the condition that it does not keep up.

And a third thing, which the pricing pages state as eligibility criteria rather than as policy: **the academic rate is a gate somebody else controls.** SnapGene's academic pricing requires that "your group must be a degree-granting institution or nonprofit registered charity and use SnapGene exclusively for academic research"; its student tier "can only be owned by an individual as opposed to a University department" and requires enrolment "in an active degree-granting program." These are ordinary, reasonable commercial terms. They also exclude, by construction, a community laboratory, a biohacker space, an independent researcher, a secondary-school teacher outside a degree-granting institution, a two-person company, and a self-funded amateur. The present author was in the last of those categories when this project began, which is a disclosed bias and also the point: the populations least able to pay commercial rates are frequently the same populations ineligible for the academic ones.

### 1.2 What the subscription actually buys

It is worth being precise about the value proposition, because the argument of this paper depends on it.

Underneath the interface, most of what these platforms compute is open source or long published. Melting temperatures come from nearest-neighbour thermodynamics, implemented in Primer3 or something equivalent. Homology search is BLAST. Profile search is HMMER. Enumerating restriction sites is a regular-expression scan against a published catalogue. Simulating an assembly is sticky-end bookkeeping, worked out in the literature of the 1970s and 1980s. Reading a file is the INSDC specification.

The hard engineering in a commercial platform lies almost entirely in two places. **Integration**: making those pieces agree with each other about coordinates, strands, topology and provenance. And **interface**: making the answer legible to a biologist in a hurry.

Integration and interface are real, difficult, expensive work. They deserve to be paid for. But they are **capital** costs, not **marginal** ones: once the coordinate conventions agree and the map draws correctly, keeping them agreeing costs maintenance, not rent. A subscription converts a one-time engineering cost into a recurring charge, and attaches to that charge a set of consequences that have nothing to do with the engineering: an account requirement, a network requirement, a browser requirement, custody of the user's sequences on someone else's infrastructure, terms that can change, features that can move behind a higher tier, eligibility criteria the user does not set, and — per §1.1, in the vendors' own words — a working tool that becomes a read-only one on the day the money stops.

None of that is a criticism of the engineering, and none of it is remedied by the price being low or even zero. A free account is still an account; a hosted notebook still needs the link to be up. The question this paper asks is narrower: given that the computational core is open and the integration cost has already been paid once by somebody, does the recurring charge and its attached dependencies have to exist at all? I take the position that it does not, and §§2–8 are the existence proof.

For a well-funded laboratory these consequences are an annoyance. For other users they are disqualifying.

### 1.3 The users for whom this matters

**The bench scientist** who works over SSH on a shared institutional machine, or who wants the design tool in the same terminal as the rest of the analysis, and does not want a browser tab holding the lab's unpublished constructs.

**The high-school or undergraduate educator** who has a room of ageing machines or Chromebooks, a filtered network, no budget line for software, and a legal obligation not to create accounts for minors on third-party services.

**The lean startup** with three people and a freezer, for whom per-seat licensing on four tools is a meaningful fraction of the first year's reagent budget, and for whom data custody is an investor-diligence question.

**And, most consequentially, researchers where resources and access are limited.** A 2026 review in *Nature Communications*, reporting a survey of African high-performance-computing custodians (74 responses, 36 complete, spanning 22 countries), lists the operating conditions directly: "slow internet speeds, system downtime due to loadshedding (electricity supply outages), limited storage space… and lack of computer resources" [10]. Of the respondents, 28 relied on a locally hosted cluster and nine on another institution's, and **none on cloud HPC at all** (the review's own counts; it reports them against 36 complete responses, though the two figures and the percentages printed beside them do not reconcile with that denominator, so I quote the counts and not the proportion). The authors are careful that the survey targeted infrastructure custodians and is used "only for contextual illustration" rather than as a prevalence estimate — but as a description of where the compute actually sits, it is pointed. A separate narrative review of sub-Saharan Africa adds the capital barrier upstream of any software question: many laboratories "still lack in-country next-generation sequencers and must send samples abroad," and establishing a sequencing facility "can be expensive, requiring US$100,000–US$700,000 in instrument capital alone (before bioinformatics compute and maintenance)" [11] — a long-standing estimate that reaches that review through a 2016 survey from a 2011 one, so it should be read as an order of magnitude rather than a current quotation. A cloud-subscription design tool fails in each of those conditions independently. It fails when the power goes out, it fails when the link saturates, it fails when the card is declined, and it fails in a way that the user cannot repair, inspect, or route around.

A tool that must work in those conditions has unusually strict requirements, and meeting them happens to produce a better tool for everyone else: it must be small, it must be fast on old silicon, it must work with no network, it must not lose data on a power cut, and it must not require trust in a vendor's continued existence.

### 1.4 Abandonware is the normal outcome, and it is measured

There is a second reason to care about where the design tool lives, and it is not about money. Scientific software decays, the decay has been quantified repeatedly, and the numbers are much worse than the field's working assumptions.

Mangul and colleagues web-mined the **36,702** published URLs of omics software resources appearing between 2005 and 2017 and found that **almost 28% were no longer accessible at the address published in the paper that introduced them** (11.9% unreachable by timeout, 15.9% broken) [6]. The dependence on age is the part worth sitting with: of software referenced before 2012 (n = 15,439), **41.9% was unavailable**; for 2012 onwards (n = 21,263), 17.5%. Kern and colleagues screened **2,396** bioinformatics web tools continuously over 133 days and found availability declining almost linearly with age — about 90% for tools published in 2019–2020 and **about 50% for tools published in 2010** [7]. Within their observation window, only 31% of tools were *always* reachable; 48.4% were intermittent and **20.6% never responded at all**. A decade earlier, Schultheiss and colleagues had tracked **927** web servers published in one journal's annual Web Server issues and could positively confirm working functionality for only **45%** of them [8].

Installability is no better than availability. In the same study, Mangul and colleagues attempted to install **98** tools: **28% could not be installed at all** because of defects in the implementation, and **57.1% failed when the documented instructions were followed exactly**. Expert manual intervention — permitted to install missing dependencies and to modify source — **failed for 66% of the tools that had failed automatic installation** (a figure that sits awkwardly against the same paper's 28%, and is quoted here only against its own denominator). Installing 98 tools consumed **71 hours**. On average eight commands were required where the manual documented 3.9 [6]. One finding in that paper is directly actionable, and SpliceCraft is built around it: **tools distributed through a well-maintained package manager were always installable, while tools not shipped through one failed to install in 32% of cases** [6]. The always-installable arm rests on about nine tools, so it is a thin result — but it is the only measurement of the question anyone has published.

This is not a tooling inconvenience. A method whose software cannot be obtained or installed is a method that cannot be reproduced, and the published paper describing it quietly becomes a historical document.

### 1.5 Why it happens: the incentive gap, and whose problem it is

The cause is structural, and it has also been measured. In Schultheiss and colleagues' survey of corresponding authors — 872 contacted, 274 replying — **78% of respondents reported that their service had been developed solely by students and researchers without a permanent position** (the figure given in that paper's abstract; its Results text and Figure 7 both report 58% for the same proposition, and I quote the abstract while noting the discrepancy), and for **24% of services there was no plan for maintenance once the original developers left** [8]. The authors state the consequence plainly: these services "are in danger of falling into disrepair after the original developers move to another institution."

That is the whole mechanism. Scientific software is overwhelmingly a by-product of a thesis or a fixed-term contract, written by someone whose career is advanced by the paper and not by the program, funded by a grant with an end date, and hosted on infrastructure tied to an institutional affiliation they are about to lose. Nothing in the reward system pays for year six. Maintenance generates no publication, no authorship, and no line on a CV; a graduate student who built a widely-used tool and then left for industry has done nothing wrong by any incentive they were given. The tool simply stops, and nobody is accountable, because accountability was never assigned.

The incentive demonstrably exists in aggregate. It just never reaches the person who would have to act on it. Kern and colleagues find a tested association between a tool's reachability and its citation count, where Schultheiss had earlier observed the same trend without testing it [7,8]. And a 2023 analysis across three areas of protein function and disorder prediction found that methods with a working web server drew, on average, **three times as many annual citations** as methods available only as source code, lacking both, or gone altogether — a median of 17.9 against 5.9 for source-only and 5.2 for no longer available [9]. The same study found that differences in *predictive performance* were unlikely to move citation rates, on the one of its three areas where a robust independent assessment existed to measure performance against. Being available beat being good. Its authors conclude that "a relatively low-cost investment into the provision and long-term support of web servers would substantially increase the impact of bioinformatics tools" [9].

So maintenance pays about threefold in the currency the field actually spends, and it pays more reliably than being better does. But it pays the *original paper*, years later, to somebody who has moved on, in a form no hiring committee will ever credit to the act of maintenance. The return is real, it is measurable, and it lands on nobody in a position to act on it.

Two responses are available. The first is to hand the problem to a vendor, which relocates the dependency instead of removing it. A company can be acquired, discontinue a product, raise a price, or move a feature into a tier the lab cannot afford. When it does, the lab has less recourse than it had with the abandoned academic tool, because it never held the source.

The second response is **software independence**: a lab, a classroom or a whole field able to run its own methods, on its own hardware, from source it can read, over data in a format it can repair, from an artefact archived somewhere that outlasts both me and my institution. That is not an ideological preference. On a ten-year horizon it is the only arrangement in which "reproducible" means anything.

I would go further. Maintaining a widely-depended-upon scientific tool is a **civic obligation**, of the kind that attaches to infrastructure and not to products. Publish a tool that other people's experiments rely on and you have taken on a duty: keep it running, or hand it over cleanly, in a state where a stranger can pick it up. The field asks for the first duty, funds neither, and has almost no convention for the second. §10 sets out what discharging the second one looks like in practice. §9 describes the method that, for the first time, makes the first one affordable for one person.

### 1.6 ApE: the existence proof, and the inspiration

SpliceCraft is not the first argument of this kind, and it is important to say whose argument it continues.

**ApE (A plasmid Editor)**, written by M. Wayne Davis and published with Erik M. Jorgensen, is a free, cross-platform plasmid editor "written over the last 17 years by a molecular biologist for molecular biologists" [13]. It visualises and edits plasmids, reads a remarkable breadth of formats — including FASTA, GenBank, EMBL, GCG, pDraw, GFF3, DNAStrider and Serial Cloner, plus Sanger chromatograms — and simulates PCR, Gibson, restriction-ligation, Golden Gate and recombinase assembly. Its authors have since written a Methods in Molecular Biology chapter on using it for *in silico* Golden Gate design [14]. It is distributed free of charge and runs on Windows, macOS and Linux; it needs no account and no subscription, and nothing in the core workflow requires a connection, though it does offer optional links out to NCBI BLAST and Entrez and an upload path to a LabArchives account.

ApE is also, against the base rates in §1.4, a **counterexample to the abandonment statistics**: it has been maintained and distributed free for over two decades by one practising scientist. That is precisely why it is the right thing to build on rather than around.

ApE settled the question of whether a single practitioner can build and sustain a plasmid editor that working scientists actually use. SpliceCraft is built in its lineage and in deliberate deference to it: ApE stores its formatting in feature qualifiers beginning `/ApEinfo_`, and SpliceCraft reads and writes the two colour qualifiers of that family, `/ApEinfo_fwdcolor` and `/ApEinfo_revcolor`, so a plasmid annotated in one program keeps its colours in the other. That interoperability is not incidental; it is the point. The goal is not to replace ApE but to extend the free, local, no-account model into the parts of the workflow that come *after* the design — assembly simulation with honest provenance, construct verification against real reads, a laboratory notebook, scripting, and automation — and into the environments where even a desktop GUI is the wrong shape.

### 1.7 Why a terminal

The terminal is an unfashionable choice that turns out to be the correct engineering answer to "must run anywhere."

A text user interface has no GPU requirement and no browser engine; the whole program, interface mounted and a plasmid on the canvas, occupies about 87 MB of resident memory (Table 2), and that figure is dominated by the Python interpreter and its imports rather than by the data — it is essentially unchanged between a 1 kb record and a 50 kb one with a hundred features. It is remotable by construction: `ssh` into a machine and the interface is already there, at a bandwidth a graphical remote desktop cannot touch — which matters on a metered or saturated link, and matters again when the capable machine is a shared server in a different building. It composes with everything else a scientist already does in a shell. It runs on a Raspberry Pi, on a decade-old laptop, on a Chromebook in developer mode, and inside the Windows Subsystem for Linux. And one server-side process can place the identical interface in a web browser with no installation at all, which is how https://splicecraft.bio works — a classroom can therefore use the tool on machines where nothing can be installed, from a single Raspberry Pi on the local network, with no student accounts and no data leaving the room.

Most importantly, a terminal degrades instead of refusing. SpliceCraft draws its plasmid maps with Unicode braille glyphs (U+2800–U+28FF), which gives 2×4 addressable dots per character cell. Where the terminal cannot encode UTF-8, the launcher reconfigures the stream; where that genuinely fails, or where the font lacks braille glyphs, the map falls back to a 7-bit ASCII density ramp that renders on anything that can display ANSI. There is no configuration in which the program declines to start because the display is too humble.

---

## 2. Design principles

Six commitments shape the implementation. They are stated here because they explain the choices in §3–§6 better than a feature list would, and because several of them are unusual.

### 2.1 Local-first, offline by default

No account exists. No telemetry is collected. Nothing leaves the machine unless the user asks for it by name.

The features that do need a network — fetching a record from NCBI, submitting an online BLAST or EBI HMMER search, pulling a codon-usage table from Kazusa, downloading a profile database — are opt-in one at a time, and they all pass through a single hardened fetch layer whose egress gate fails closed. Programmatic callers are gated twice. An automated agent cannot submit a sequence to an online search service until the user has ticked the setting that allows it.

Keeping sequences out of the diagnostic log works by two mechanisms with different strengths, and the distinction is the whole of the honest claim. The structured event path — the one that records user-visible state changes — **enforces** a bound on every field it is given: a string over 200 characters is truncated to a prefix and a count, a byte blob is reduced to a type and a length without being rendered at all, and a Biopython `Seq` or `SeqRecord` is replaced by an opaque tag carrying only its identifier and length, precisely because the standard `repr` of those classes embeds the first several dozen bases. Nothing a caller passes to that path can take sequence content with it. Plain log calls are the weaker half: they do not pass sequences, but that is a convention rather than an enforcement, and nothing truncates them. So the accurate statement is that sequence content cannot reach a structured event and is merely unlikely to reach an ordinary log line — a partial guarantee, not a total one. The project's own architecture notes claimed the total version; §9.3(g) follows what that overstatement did across three drafts of this paper, including this one. The bug-report bundle ships logs, sanitised settings and system information, and inherits exactly the same split.

### 2.2 The data is the product

A plasmid library is months to years of irreplaceable work, and the user frequently has no institutional backup. Data integrity is therefore treated as the system's primary requirement, outranking features, performance and convenience. §4 describes the mechanism.

### 2.3 Simulate the steps; never fabricate the result

When SpliceCraft reports an assembled construct, it has done the assembly. Real amplicons out of a real primed-amplicon simulation. Real digestion by the enzymes you picked. Real sticky-end compatibility checks, real ligation, real circularisation. The product is whatever those steps produce, and it is never hand-built backwards from the answer you wanted.

That is not fastidiousness for its own sake. Three things fall out of it.

Provenance comes free. Every plasmid carries a lineage tree and a numbered protocol that reads like a bench recipe, because one was generated by the other.

The viewer can check the record against the molecule. If a stored step claims a restriction site the DNA does not contain, or a primer that no longer anneals, the display says so instead of repeating the claim as fact. Sequences and histories are never quietly edited to make a warning go away.

And the simulation can find failures that an asserted answer would hide. A backbone whose two ends happen to be compatible with each other will close with no insert at all, and that empty vector is where most empty colonies come from. So the traditional-cloning simulation names that route, says whether a diagnostic digest could even tell it apart from the construct you wanted, and suggests the fix.

### 2.4 Pure Python, by commitment rather than by convenience

SpliceCraft contains no compiled component of its own. Every line of it is Python, and that is a load-bearing decision rather than an artefact of how it was started.

It means there is no build step, no ABI to break on a Python upgrade, no platform-specific binary to produce, and no architecture on which the program itself fails to exist. It means a user can read the function behind the button they just pressed — which matters for trust, for teaching (§8.2), and for the project outliving its author. And, least obviously but most importantly, **it means identical results on every platform.** A mixed-language implementation that selects a fast native path where one is available and a portable path elsewhere is not one program; it is two programs that agree most of the time, and the disagreements surface as irreproducibility on exactly the hardware the maintainer does not own. §3.3 describes the case where that happened here, and what was removed to fix it.

Where pure Python is assumed to be too slow, that assumption is worth testing rather than inheriting: the built-in aligner is pure Python and ~10× faster than the Biopython dynamic-programming path it replaced on an 18 kb read, and the replacement GenBank reader is ~1.1× faster than Biopython's (§3.3, where both figures are measured with their method stated). Neither comparison is Python-versus-C: Biopython's GenBank parser is itself pure Python, and Biopython's compiled extensions serve its alignment, clustering and structure code rather than its flat-file readers. The gains came from algorithm and I/O structure, not from a language change — which is the honest version of the claim, and a weaker one than "pure Python beat C." The compiled code in the installation belongs to upstream dependencies whose algorithms genuinely need it, and every one of those is *integrated*, not reimplemented.

### 2.5 Integrate established tools; rewrite only on a measured defect; refuse the dependency the standard library covers

The default is to use the accepted implementation, in-process, rather than to rewrite it. Melting temperatures and oligo thermodynamics come from the real Primer3 C engine via `primer3-py` [2]. Profile and sequence search run HMMER 3 [22] compiled into the process via PyHMMER [3], which removes the need for the user to install an external BLAST+ or HMMER toolchain at all. The record model, the writer, and the EMBL/FASTA/AB1 readers are Biopython's [1]. Agarose mobility follows the Helling–Goodman–Boyer empirical relation [15]. Secondary structure uses nearest-neighbour parameters. Codon usage tables are fetched from Kazusa [16] with each species name checked against its taxonomic identifier and each table's basis recorded.

This is also the direct lesson of §1.4: a tool that is hard to install is a tool that will be unavailable, and the single strongest predictor of installability in the only study to measure it was distribution through a maintained package manager [6]. SpliceCraft is therefore `pipx install splicecraft`, with every dependency resolving to a prebuilt wheel, and a release-gate script that fails the release if a dependency loses wheel coverage on a supported platform. Two current version pins exist for no reason other than this: Biopython is capped below the release that dropped Intel-macOS wheels, so that no Mac user is asked to compile a large C extension, and one pure-Python dependency's floor is held below its latest release because raising it would make the conda-forge solver fail. Packaging decisions are accessibility decisions.

Two components *were* rewritten, and both are instructive because the rewrite was triggered by a defect in integration rather than by a preference, and both were benchmarked against the thing they replaced before shipping. §3.3 gives the numbers.

The third form of the discipline is turning down a dependency the standard library already covers. The command-line client imports nothing outside stdlib, so it runs under any Python. The robot driver talks to the OT-2 over `urllib` instead of taking on `requests`. The local-model assistant speaks HTTP to Ollama directly, rather than pulling in an SDK and a vector database to do it. SVG export is string concatenation. Eight runtime dependencies remain, one of them platform-conditional, and the number has gone down at least once.

### 2.6 Report the operating point; refuse to bluff

Every quantitative claim ships the measurement that justifies it, and the measurement is reported *in the output*, not only in the documentation. A scanner's threshold states the false-positive rate it was anchored at and the sensitivity that buys. An empty result set says how hard it looked, because an empty list reads as "clean" when it means "not checked." A search over the plasmid in front of you reports the number of base pairs searched, so "no off-targets" can never be misread as genome-wide. Where a model was built, measured and found not to work, the feature is absent and the number that killed it is published. §5 treats this at length; it is, I think, the single most important property of a tool used by people who have no core facility to check it against.

### 2.7 Degrade, never refuse

Braille map unavailable → ASCII map. `primer3-py` unavailable → documented, less accurate GC rule, with the degradation stated. PyHMMER unavailable (native Windows) → pure-Python search engine and an explanatory message. No C compiler → everything except one dependency installs from wheels. Clipboard tool absent → a four-tier cascade ending in a file the user can open. Corrupt data file → restore from the backup chain. Non-UTF-8 locale → reconfigure, then fall back. The program's refusals are reserved for cases where proceeding would be wrong, not for cases where the environment is poor.

---

## 3. Implementation

### 3.1 Architecture

SpliceCraft 1.3.19 is pure Python, declaring support for 3.10 and above and developed on 3.14, built on the Textual terminal-application framework [23], in a **hub-and-layered-siblings** layout: an application hub holding the controller, the three main panels and the application-coupled dialogs, plus 38 flat sibling modules holding every cleanly separable layer — biology primitives and the restriction scanner, the persistence engine, the data-access layer, GenBank serialisation, the INSDC reader, file I/O, the construction simulator, sequence analysis, the codon optimiser, primer design, the regulatory and splice and initiation-context models, gels, the notebook, the backup engine, the map-image exporter, CRISPR, the plant cis-regulatory engine, the robot driver, the chat engine, and the agent API.

The siblings are strictly layered L0→L7 with no upward imports, and that property is enforced by a test rather than by convention: `tests/test_import_layers.py` fails the suite on an import cycle, an upward import, or a sibling missing from the packaging manifest. Shared mutable state lives in one module accessed by attribute, because a by-value import binds a stale copy that silently desynchronises from runtime writes — a class of bug worth a dedicated guard test. The public API surface is pinned byte-for-byte in a baseline file, so an accidental rename of anything a script might call fails the suite.

Totals: ~222,000 lines of implementation across 39 modules; ~167,500 lines and 12,277 tests across 146 test files — a test-to-implementation ratio of roughly 0.75:1.

### 3.2 Footprint and performance

Table 2 reports measured figures rather than estimates. "Runs on a potato" is a design target, and a design target has to be measurable.

**Table 2.** Measured footprint and representative operation costs. Install size measured on Linux / Python 3.14 as a clean virtual environment with the package and every dependency installed from PyPI, counting `pip` itself and the NumPy that Biopython pulls in — i.e. the whole directory a user would have, not the package alone. Memory and map figures were measured through a *mounted* application instance driven by the test harness, because a bare render call short-circuits and reports a meaningless zero. Per-operation ceilings are the enforced regression budgets rather than single observations; the figures given as single values are measured best-of-N.

| Quantity | Measured value |
|---|---|
| Complete end-user install, every dependency, incl. `pip` and NumPy | 202 MB on disk |
| Runtime dependencies | 8 (one platform-conditional); wheels on x86-64, one ARM source build |
| Peak resident memory, import + full-catalogue restriction scan | 73 MiB |
| Resident memory, interface mounted, empty canvas | 85 MiB |
| Resident memory, interface mounted, 50 kb / 100-feature plasmid loaded | 87 MiB |
| Import time (cold) | 619 ms |
| Full 202-enzyme restriction scan, 1 kb | 2.1 ms |
| Restriction scan, 10 kb (regression ceiling) | < 150 ms |
| Digest with 6 enzymes, 50 kb circular, 100 features carried | < 10 ms |
| GenBank parse, 100 kb / 500 features | < 20 ms |
| Codon optimisation, 500 residues | < 2 ms |
| Six-frame ORF scan, 50 kb circular | < 75 ms |
| Map repaint, 20 kb / 2,000 features | < 30 ms |
| Map per-frame ring cost, 10 Mbp / 10,000 features | 0.0013 ms † |
| Map per-frame ring cost, 10 Mbp / 200,000 features | 0.0014 ms † |
| Map per-frame ring cost, 10 Mbp / 1,000,000 features | 0.0012 ms † (flat; the aggregate path is O(buckets)) |
| Full map repaint through the harness, 10 Mbp / 200,000 features | 86 ms † |
| Nanopore-scale read alignment, 18 kb at ~99% identity with indels | 0.48 s vs 5.05 s (10.6×) ‡ |
| Full test suite | ~5 min on 8 cores |

† Measured for this manuscript, not gated. The per-frame figures use the same method as the regression test `test_aggregate_cost_does_not_grow_with_feature_count` — a warm aggregate ring, mean of 20 — which gates the *shape* (cost must not grow with feature count) but not an absolute value. The flatness from 10,000 to 1,000,000 features is the claim that matters and it reproduces; the absolute figure is machine-specific. The full-repaint figure drives the mounted widget through the test harness and is therefore an upper bound that includes the framework's own render pipeline. Initial `load_record` at 200,000 features took 8.3 s on this machine, which is a real cost and is *not* interactive — the constant-time property applies to redraw, not to first load.

‡ Measured for this manuscript against Biopython's `PairwiseAligner` on identical input (18 kb reference, ~0.8% substitutions plus six indels of 1–30 bp, full traceback materialised in both engines, same process, same machine). The scoring is the application's own default and the figure does not reproduce without it: match 2.0, mismatch −1.0, gap open −2.0, gap extend −0.5. Substituting unit-cost linear gaps roughly halves Biopython's time and the ratio with it, so a reader checking this must use the affine parameters above. Earlier project notes record 0.7 s against 8.6 s (~12×) on a real sequencing read; that measurement is not reproducible from the repository and is superseded here by one whose method is stated. The conclusion is unchanged and the margin is slightly smaller.

**Nineteen** operations are enforced as **regression gates**: `tests/perf-baseline.json` pairs each with a ceiling, and the suite fails if the best of N samples exceeds it. Loosening a ceiling requires a written justification in the same change. This exists because processing speed on modest hardware is a feature that silently erodes; six operations were gated before a 2026 audit, and thirteen more were added when the audit observed that the digest, the codon optimiser, the scanners, the readers and the writers could all have regressed unnoticed. The rows above marked † and ‡ are **not** among the nineteen, and the rows given as `< N ms` are the enforced ceilings rather than observations.

A note on honesty in benchmarking: each gated operation explicitly clears the cache it would otherwise hit. An uncleared enzyme-cut cache made the digest measure 0.01 ms, which is to say it measured nothing at all.

### 3.3 The two rewrites, and what they cost to justify

**The pairwise aligner.** Checking a construct against a long sequencing read needs global alignment. The fast implementation everyone reaches for is a C library, and it ships wheels for some platforms and not others. On the platforms without one, which means ARM Linux, recent Python versions, native Windows and a large share of fresh installations in practice, the fallback was Biopython's dynamic programming. On an 18 kb read that costs seconds, not a fraction of one: 5.05 s under the measurement in Table 2. Nobody waits five seconds per read on a laptop, let alone on a Pi.

A built-in bit-parallel Myers aligner [17] with Hirschberg linear-space recursion [18] was written in pure Python. Measured for this manuscript on identical input, an 18 kb reference against a read carrying ~0.8% substitutions and six indels of 1 to 30 bp, with the full traceback materialised by both engines in one process, it completes in **0.48 s against Biopython's 5.05 s, a ratio of 10.6×** (Table 2). Pure Python is viable here for a specific reason: Python's arbitrary-precision integers give the Myers bit vector for free, so the whole pattern is one integer and the 64-bit word-blocking a C implementation needs disappears. An IUPAC-aware common-prefix/suffix trim applied at every recursion level peels the long exact-match runs between scattered errors, and is the component the implementation treats as the principal lever; its *correctness* is gated by a dedicated test on scattered-error inputs, though the speedup attributed to it in project notes is a one-off measurement the repository does not reproduce. Correctness of the engine as a whole is established against a naive O(nm) dynamic-programming reference, using the same IUPAC-compatibility predicate, over **400** random pairs, asserting round-trip reconstruction *and* minimum-edit-distance optimality rather than mere agreement, plus 200 further pairs against the small-block base case, and verified to agree exactly with Biopython — same matches, mismatches and gaps — on a uniquely-optimal near-identical pair.

The instructive part is what happened next: the C library was **removed entirely**. It was faster where it existed, but it selected different co-optimal alignments, which produced phantom variants on reads missing an end. And because the release machine could not install it, the suite that gates releases had never once exercised that path. One engine on every platform, with identical results everywhere, was worth more than a 35× speedup available to some users and invisible to the gate. A test now fails if the dependency returns.

**The GenBank reader.** Every GenBank text the application reads now passes through SpliceCraft's own INSDC reader. The trigger was not performance but a class of silent defects: Biopython's reader dropped roughly ten bases per *unnumbered* sequence line, returned zero records for a file with a UTF-8 byte-order mark, raised out of its scanner on a tab-indented feature table, and produced `location=None` for an origin-spanning feature on a record whose LOCUS line was not marked circular, after which the next library save crashed. Files like these come from real commercial exports and real public repositories.

The replacement builds Biopython's own model objects, the same `SeqRecord`, `SeqFeature`, location and position classes, but it reads by what each line *is* instead of by which column it starts in. The contract is written down. **Produce the same record as Biopython wherever Biopython reads correctly.** Where Biopython refuses, crashes or reads wrong, produce what the file means, with a warning and a counted repair.

All of the differential evidence below was re-run against the current reader for this manuscript, not quoted from project notes. **313 of 313** NCBI reference records are identical. In my production library, **1,272 of 1,625** records are identical and **353 differ** — every one of them the same deliberate decode rather than a defect: SpliceCraft consumes the `SpliceCraft_strand` qualifier and restores strand 0 for a feature drawn without an arrow, where Biopython keeps the qualifier and defaults the feature to the plus strand. Of **214** local GenBank files, SpliceCraft reads all of them and Biopython refuses two; of the 212 both read, **171 are identical** and 41 differ, again by that decode.

Then the adversarial set: Biopython's own 50 GenBank test fixtures, which exist precisely to exercise edge cases. **40 are identical and 2 differ by design**, those being the deliberate repairs. The remaining **8 cannot be compared, and the direction is the part worth reporting: six are CONTIG or scaffold stubs that SpliceCraft's File > Open refuses outright**, because a record carrying no sequence of its own is nothing a plasmid editor can open, while Biopython returns a record with a fabricated sequence of the stated length. One file only Biopython refuses; one both refuse. A committed test round-trips a further 180 records written by Biopython's own writer.

Throughput, both readers in one process over the 313-record corpus (16.1 MB of text in memory), best of three, parsing from text already in memory: **35.8 MB/s against Biopython's 32.4, a ratio of 1.11×.** Repeated on the production library through the same harness it is 1.06×, so the defensible statement is a margin of roughly a tenth, not a third.

Two notes on that paragraph, because the first draft of it was wrong in instructive ways. Project notes recorded a larger differential corpus (1,303 library records, 196 local files) and a wider throughput gap (37–46 against 27–31 MB/s); the corpus figures were simply **stale**, since the library has grown. The throughput gap was something worse: an artefact of **comparing different things**, a reader-level measurement for one engine against a measurement that charged the other for file I/O and for a double-parse the application door used to perform. Measured identically, the margin is about 1.1× rather than ~1.4×, and that is the number above. The hardening pass that followed mutation-fuzzed damaged files with Biopython as co-oracle on every case. It found and fixed, among others, a quadratic string concatenation on a 200,000-line qualifier value and an unbounded recursion on a deeply nested location. The committed fuzzer runs **400** mutants from one fixed seed and the format scoreboard's default fuzz segment 480; the much larger ad-hoc run recorded in project notes is not reproducible from the repository, so the defensible claim is the one that is: every finding it produced is now a named regression test, and 24 of them fail against the pre-hardening reader.

Biopython remains the writer, and the test oracle. Both rewrites illustrate the same rule, which governs the whole dependency policy. *Do not respin an established tool. If the integration is itself the defect, replace it only with the differential numbers already in hand, and keep the original as your oracle.*

### 3.4 Correctness regime

Ten rules are marked **sacred**. They are the biology-correctness facts that are easy to get subtly wrong and expensive to get wrong at the bench. Palindromic enzymes are scanned forward only, because scanning both strands double-counts every site. A reverse-strand site keeps its forward coordinate. Complementation covers the full IUPAC alphabet, not just ACGT. There is one formula for the midpoint of an arc that crosses the origin, and one two-piece emission for a restriction site that straddles it. Feature length is computed wrap-aware. An edit shifts each part of a compound location separately, because asking Biopython for the start of a wrapped feature silently returns the smallest part and flattens the rest. Undo snapshots are deep-copied. Every save goes through the chokepoint.

Each rule has at least one regression test, and touching the function it governs trips that test in under two seconds. The point of naming them is not ceremony. It is that a future contributor — or a future model, working from a fresh context with no memory of why — gets a loud, fast, specific failure instead of a plausible-looking patch.

Above that, the suite applies:

- **Property-based testing** (Hypothesis, 21 properties at 100 to 300 examples each) on the biology layer's circular-coordinate and reverse-complement invariants, the primer-Tm fallback and wrap-aware CDS translation: that reverse-complement is involutive and agrees with Biopython on unambiguous bases, that wrap-aware feature length matches the closed-form wrap formula, that a computed label midpoint lies on the arc it labels, and that a palindromic enzyme yields exactly one site record. The codon optimiser's properties (ACGT-only output, in-frame, no internal stop, no codon→amino-acid mismatch) are asserted as conventional tests rather than generatively, which is a weaker guarantee and worth saying so.
- **Differential testing against independent references.** The CRISPR PAM scanner is validated against a brute-force first-principles matcher written separately from the engine and living in the test file rather than as a fixture: **24 of 24 exact on linear and 18 of 18 on circular** molecules, across four and three length classes respectively. Export conformance is checked twice — once against spec validators written from the published INSDC, ENA EMBL and GFF3 1.26 documents, and again by re-reading the output with parsers that did not come from this project.
- **Crash injection** on the save path: a forked child writing through the atomic-save path is `SIGKILL`ed inside a ~50 ms window, and recovery must yield one of the two intact payloads and **never a torn mix**. And **concurrency fuzzing** on the data layer: eight threads for 1.5 s of random save/load against a 20-entry library, asserting every entry survives as a dict with a unique non-empty id and the on-disk JSON still parses.
- **Mutation fuzzing** of damaged inputs: the committed GenBank fuzzer runs 400 mutants from one fixed seed under six mutation operators, asserting that nothing but `ValueError` is ever raised, no location is ever `None`, and no impossible feature key is ever produced; the format scoreboard's fuzz segment adds 480 more. A GFF3 round-trip fuzz covers 28 cases across seven feature shapes × both strands × coding and non-coding, with **the GenBank path as co-oracle** — a difference counts as a loss only if the GenBank round trip preserved it, which is what stops the test flagging shapes no format can express. Larger one-off fuzz runs are recorded in project notes; this paragraph deliberately cites only what a reader can run.
- **Performance gates** (§3.2), **static type checking** (pyright) and **linting** (ruff) in CI.
- A methodological rule learned the hard way: **a regression test counts only once it has been seen failing on a pre-fix snapshot.** A source-string assertion and a bulk-ticked checklist once concealed twelve of twenty-six items that had not actually been fixed.

Every test runs offline against synthetic records and monkeypatched paths, and an autouse fixture makes it impossible for a test to write to real user data.

---

## 4. Data integrity: the interlocks

This section is the centre of the paper's claim to be a tool you can trust with the only copy of something.

### 4.1 Four layers on every write

All user data persists as human-readable JSON. That is a choice, not laziness. Anyone who outlives this project should be able to recover their work without it, and plain text is what makes that possible. Every write goes through one sanctioned function, and that function does four things.

1. **Atomic write with read-back validation.** `tempfile.mkstemp` → `os.fsync` → `os.replace`, plus a parent-directory fsync so the rename is journalled, with the prior version copied aside first. Before the swap, the freshly written temporary file is **read back and validated**: a full parse plus entry-count check below 8 MB, and a cheap non-empty-and-well-terminated check above it. A silently truncated write is therefore caught, and the live file left intact instead of overwritten with garbage.
2. **Rotating generational backups**, ten retained per file, ordered by an explicit generation counter rather than by wall-clock stamp. (Order by the stamp and two backups written in the same second, or a clock that steps backwards, file the newest as the oldest. Retention then deletes it immediately.) Older generations are gzipped at launch and the newest kept plain, so the fast recovery path never has to decompress. An aggregate per-file byte cap and a byte-identical dedup stop a 274 MB library from accreting 2.7 GB of backups.
3. **Daily snapshots**, thirty days retained, written once per calendar day at launch.
4. **A suspicious-shrink guard.** A save that would discard more than half the entries spills them to a recovery directory *before* the overwrite proceeds. A save that would discard more than 90% is **refused outright**, unless an explicit bypass has been armed. Legitimate shrinks take a separate, exempt path: switching the active collection rewrites a mirror file from one collection's contents to another's, which drops entries without losing them. That exemption exists because a guard which fires on correct behaviour is a guard somebody eventually disables.

Recovery runs the chain backwards. A corrupt main file is restored from the immediate backup. If both are corrupt, the loader walks the generational rotation from newest to oldest, decompressing as it goes, and restores from the first file that parses, setting the corrupt one aside for forensics. Restoring a spill *merges* it into the current file instead of replacing it, since a spill holds only what one shrinking save dropped. Plasmids left with no collection are gathered into a recovered collection at the next launch, never overwritten.

### 4.2 Content-addressed sequence storage

The bulk of every entry is its full GenBank text. Each unique text is stored exactly once as an immutable blob named by its own SHA-256, with entries carrying a small reference. A blob's filename *is* the hash of its content, so a present blob is correct by construction; writes are re-hashed after landing and reads re-hash and refuse a mismatch, so a truncated blob fails loudly instead of being returned as a silently wrong or empty sequence. A save that cannot persist and verify a blob **aborts before writing any metadata that would reference it**. A missing blob leaves that one entry's text empty, keeps the reference, and resolves itself when the blob is restored — the entry is never dropped. Garbage collection *quarantines* rather than deletes, refuses to run if any metadata file is unreadable, treats everything referenced by backups as live, and skips blobs younger than an hour to protect an in-flight save.

### 4.3 Interlocks against the tool's own maintainers

Two mechanisms exist specifically because the most dangerous writer is an authorised one.

**A write chokepoint.** Every domain save function refuses to run unless the caller is one of four sanctioned contexts (the application's own main loop, the test suite's data-protection fixture, the API server, or an explicitly sandboxed verification harness). An ad-hoc script that imports the library and calls a save function raises `RuntimeError`. This is an unusual thing to build — it defends the user against the developer — and it exists because of a documented incident. On 2026-05-22 an unsandboxed diagnostic script, written to set up a test, called the collections save function and wrote directly over my real 160 MB library, rotating the previous good state to `.bak`. The four-layer net and the spill directory recovered the data in full. The chokepoint was added so that the recovery is never needed again, and the incident is recorded verbatim in the contributor handover document rather than quietly fixed.

**Pre-update snapshots.** `splicecraft update` snapshots the entire library, collections, parts, primers, features, grammars, codon tables, settings, autosaves and format sidecars *before* invoking the package manager, and aborts the upgrade if the snapshot cannot be taken. Snapshots live in a *sibling* directory of the data directory, so that a hypothetical recursive-deletion bug in a new version cannot reach them. Restore verifies the schema version, checks the target against a whitelist, rejects path separators and traversal, and re-verifies a SHA-256 against the manifest before any file is replaced; a manifest with a missing or empty hash is refused. The restore is itself reversible.

### 4.4 Bounded inputs, bounded blast radius

Every external input has a documented ceiling: a 1 GB load cap, a 50 MB default on API file loads, 64 MB on a single record, a 500 MB/50 MB-per-member/2000-member cap on sequencing archives, streaming decompression with explicit length limits on compressed history, and per-fetch response caps on every network source. Destructive operations are previewed before they run and undo as a unit: a bulk deletion restores all N entries on one keypress, and if it cannot restore all of them it restores none and says why — partial recovery reads as complete recovery, which is worse than an honest failure. Name collisions on import always ask, with three outcomes (skip, copy, overwrite) and never a silent choice. The full-wipe operation is triple-gated and, by deliberate policy, **permanently unreachable from the API**: an automated agent may create, edit, delete and undo, but it may never wipe a library.

Interactive edits debounce a three-second autosave into a per-record crash file, so a power cut mid-edit surfaces the survivor at the next launch. The data directory takes a per-process lock with stale-PID detection. A session killed on a shared filesystem therefore releases its lock when the program next starts, instead of locking the user out of their own library. That last one was a real field report before it was a feature.

---

## 5. Honest calibration: reporting what a score means, and refusing what does not work

A feature list is a weak argument for a scientific tool. I think the strongest argument SpliceCraft can make is about epistemics, and it applies most to exactly the users this paper is addressed to: someone without a core facility, a staff bioinformatician, or a second software licence to cross-check against is at the mercy of the confidence their tool projects.

So the rule is: **a number is reported together with the measurement that justifies it, a threshold is reported together with the error rate it was anchored at, and a model that does not beat chance does not ship — but its numbers do.**

### 5.1 Scores that state their operating point

The **plant splice-site model** scores donors and acceptors as log-odds against a decoy built from real-but-unused GT/AG dinucleotides sitting inside verified introns. Its thresholds are anchored at a **1% false-positive rate**. The sensitivity that buys, measured, is **0.459 for donors and 0.333 for acceptors**, and both figures ship with every scan.

So the output says what follows from that, and so does the user documentation: **a clean scan is not a clean sequence.** About half of real donor sites fall below the bar. A real intron whose own acceptor scores under threshold is the model doing exactly what was measured, not a bug. An untrained or uncalibrated model raises instead of returning an empty list, because a silent empty scan reads as "clean" when what it means is "not checked". (The module docstring used to claim a sensitivity-anchored threshold, which overstated the scan's reach two- to threefold. §9.3(g) has that story.)

The **plant translation-initiation model** is a log-odds position-weight matrix over a −12..+4 window. It was trained to tell real start codons apart from the upstream ATGs a uORF screen actually has to rank, not from base composition. It reports a percentile against real plant starts, and it ships beside the −3-purine/+4-G consensus rule instead of replacing it, with both halves in one result and each labelled with its own basis. They can disagree, and the disagreement is informative: the rule calls a poly-T leader with a +4 G "adequate", while the matrix puts the same site at the 0th percentile.

Held-out AUC is **0.766 against the rule's 0.712** (ΔAUC +0.054, 95% CI [+0.035, +0.074], gene-resampled), and it wins in six of seven species held out one at a time and ties in the seventh. The capability that actually justifies building it is narrower than any of that: inside each of the rule's three buckets, it still separates real starts from upstream ATGs (AUC 0.62, 0.68, 0.66). A three-level rule cannot order anything within its own levels. "Which of my four upstream ATGs matters most" had no answer before and has one now.

The **bacterial promoter and terminator scanners** take their thresholds from real elements measured against a matched random background, with the numbers written beside the constants. The promoter default works out to roughly one chance hit per kilobase. Terminators are reported at ΔG ≤ −8 kcal/mol and called strong at ≤ −12, about one per 25 kb. That asymmetry is on purpose: a spurious terminator makes the transcription map declare a gene insulated when it is not. Neither scanner returns a rate. The promoter score is a composite of named components, labelled uncalibrated, and good for ranking candidates against each other.

The **plant cis-element catalogue** carries, for each of its 34 entries, a measured rate in real promoters and in a composition-matched background, with enrichment derived from those two numbers so that a label cannot drift away from its evidence. The uninformative entries are kept and labelled as such rather than quietly dropped, because the useful answer to "does my promoter have a GT1-motif" is not a hit. It is "that motif turns up in 94.8% of random 400 bp windows". TATA-box scanning refuses to run without a transcription-start anchor: `TATAWA` occurs about once per 429 bp of real plant DNA, roughly 47 times on a 10 kb plasmid, so presence carries no information and position carries all of it. It peaks at −40 to −31, 3.17-fold over background.

### 5.2 Capabilities deliberately absent, with the numbers that killed them

**Table 3.** Models that were built, measured, and refused. Each number is retrievable at runtime from the corresponding calibration endpoint.

| Refused capability | Measurement that killed it |
|---|---|
| Plant poly(A) / cleavage-site prediction | A position-specific model trained on 600 real plant poly(A) sites located them **worse than chance** (±20 nt in 27.2% vs 29.1% expected). No arrow is drawn and no field exists to hold one. |
| Plant promoter-strength score | Best composite reached **AUC 0.646**; 5.0% sensitivity at a usable background rate. |
| "Plant-specific" Y-patch (pyrimidine patch) scan | Measured **1.05×** over background. |
| Initiation-context window beyond +4 | A −12..+9 window scores **0.791**, but all of the gain is the coding-versus-UTR composition shift — it measures "is there coding downstream," not initiation. |
| Initiation rate or protein yield | Not produced. The matrix ranks contexts against each other; it is not a trained predictor of output. |
| Scoring an ATG with fewer than 12 nt upstream | Returns `None` with a note, never `0.0` — because 0.0 is a *real* score meaning "as likely a start as an upstream ATG." Costs 4.2% of real upstream ATGs. |
| Genome-wide CRISPR off-target search | No index exists and one is not faked. The search covers the molecule supplied and reports `searched_bp` in its own output so a caller cannot overstate it. |
| CRISPR on-target efficiency score | Not produced. Guides are triaged on named mechanistic properties (a `TTTT` that terminates Pol III, a GC window out of range, a self-folding spacer) and a coarse tier, because a number that looked like a published efficiency score and was not would be worse than no number. |

A related integrity rule: during reference-set construction for the plant models, **92% of GenBank single-gene "TSS" annotations turned out to be start codons**, so the pipeline requires an annotated 5′UTR. That finding is published with the calibration rather than absorbed silently.

### 5.3 Silence is not safety

One failure shape has turned up in this codebase's own audits more often than any other: **a component that goes quiet instead of refusing.**

The full chain, from a real case. A digest cannot proceed, so it returns nothing. A wrapper turns that nothing into `None`, which looks identical to "no vector was configured". A check further downstream sees no input, finds no problem, and says so. At the end of it, a working acceptor plasmid scored 20 out of 100, and so did a broken one. Four components behaved reasonably and the answer was worthless.

The rule adopted throughout is that an operation which cannot answer has to say it cannot answer. Scanners raise instead of returning an empty list. An empty filter result is distinguishable from a question nobody asked. A verdict that covers only part of the evidence names the part it covers. A sequencing read with nothing legible in it comes back "inconclusive", never "mismatch", because no information is not evidence against a clone.

---

## 6. Capabilities

Table 4 maps the feature set onto the workflow rather than onto the menu bar. Full reference documentation ships with the source.

**Table 4.** Capability map by workflow stage.

| Stage | What SpliceCraft does |
|---|---|
| **View / edit** | Circular and linear braille maps with lane-stacked features; per-base sequence panel with live translation; wrap-aware six-frame ORF finder; end-for-end flip and origin re-cut as real undoable edits; fuzzy subsequence search on both strands; publication-quality PNG/SVG map export |
| **Restriction** | 202-enzyme catalogue plus user-defined enzymes; named enzyme collections scoping every scan; live superscript cut counts that decrement as a site is edited out; supplier-synonym resolution (Eco31I→BsaI, LguI→SapI, AarI→PaqCI); origin-spanning and Type IIS handling pinned by invariant tests |
| **Annotation** | User feature library plus a 117-element, 20-category preset catalogue in which every sequence was extracted from public GenBank records by multi-accession consensus — never typed from memory — with the accession and independent-submission count cited per entry; one-click annotation of markers, origins, promoters and polylinkers with preview and single-keystroke undo |
| **Primers** | Primer3-backed designer for detection, cloning, modular-grammar and generic primers; Designed→Ordered→Validated lifecycle; measured pair QC (hairpin flagged when it holds a 3′ end at the annealing temperature, dimer when stable above 47 °C); *in-silico* PCR across the whole library with identity, strand, position, amplicon and feature; rotation-invariant re-derivation of binding sites so a primer displays exactly where it anneals; CSV and IDT bulk-upload export |
| **Mutagenesis** | Residue-level editing through strand, frame offset, introns and the origin, with the resolved codon shown before anything changes; SOE-PCR design with a two-primer fallback near the ends, offered only when the primer genuinely carries the change; whole-plasmid restriction-site cure by minimal point changes, silent across every overlapping frame, reported when a site cannot be cured silently |
| **Synthesis** | Linear composer with feature stripes and live translation; protein→CDS with frequency-matched codon optimisation, CAI/GC3/rare-codon analysis, and forbidden-site scrubbing; codon tables built from an NCBI genome, a local CDS file, Kazusa or TSV; plant designs additionally screened with the trained splice matrix, with repair on request; eukaryotic expression-cassette assembler; long direct and inverted repeat detection |
| **Assembly** | Traditional restriction, Gibson [21], Golden Gate, GoldenBraid [20], MoClo [19] and user-defined grammars; every step simulated (§2.3); parent features carried forward so a finished construct traces back to its L0 parts; regenerated sites read as intact; each junction labelled with the bases that actually anneal; empty-vector, double-insert, short-chain and inverted-fragment outcomes enumerated, with a statement when no such route exists |
| **Verification** | Sanger `.ab1` and nanopore/Plasmidsaurus ingestion; partial reads anchored on the reference they actually span; inverted regions recognised as inversions rather than written off; basecall quality consulted so "1 real change, 1 in unreliable basecalls" replaces a flat count; per-base allele fractions from raw FASTQ, because a consensus is a single read and a population of distinct escapers consenses back to wild type; platform-matched noise floors; graded verification report that never rounds a true sub-100% identity up |
| **Analysis** | σ70 promoter and intrinsic terminator scans; whole-circle transcription mapping that separates "nothing reaches this gene" from "nothing the *host* polymerase reads reaches it"; mature-mRNA reconstruction from annotated introns with ranked upstream ATGs and an explicit list of which ones splicing removes; plant cis-regulatory annotation; CRISPR guide design for three Cas variants with annealed cloning oligos; agarose gel simulation on a published mobility curve; multi-sequence alignment lanes under the map |
| **Record-keeping** | Markdown lab notebook with projects, live cross-references between entries, plasmids, actions and gels, arbitrary attachments, full-text search, and reverse lookup from a plasmid to what was written about it; protocol generation from a plasmid's recorded build steps, which writes only steps that were actually recorded; Markdown and standalone-HTML export |
| **Automation** | 303-endpoint local HTTP API with self-describing request schemas; stdlib-only CLI client; read-only attachment alongside a running GUI (reads answer, writes return 409 naming the lock holder); Opentrons OT-2 liquid-handler control with deck layout, protocol compilation, on-robot simulation, and interlocked execution |

---

## 7. AI integration, and why calibration is its precondition

SpliceCraft is AI-enabled in two distinct senses, and both are designed around the assumption that the model will sometimes be wrong.

### 7.1 A machine-drivable instrument

The 303-endpoint local API exists so that every operation a human can perform is also performable by a script or an agent — a standing design requirement, applied in the same change as each new user-facing feature rather than retrofitted. Endpoints self-describe their request schemas at `/tools`. Subtle operations are paired preview-and-commit. A published Claude Skill teaches a language model the parts that are easy to get wrong: how to attach read-only to a session the user already has open, the canvas-versus-library distinction, the coordinate convention, which endpoints must be paired, and — importantly — **the failure modes that look like success**, such as an assembly that cannot close answering HTTP 200 with `ok: false`.

That last item is the whole design thesis for agentic use. An agent cannot tell a plausible answer from a correct one; it can only propagate what the tool told it. An instrument that returns an empty list for "I could not check" will cause a model to report a clean construct. The calibration discipline of §5 is therefore not a separate virtue from the AI integration — it is its precondition.

### 7.2 Natural language as a first-class entry point

The consequence of full parity is worth stating as its own claim, because it changes who can use the program: **the terminal interface is not the only way in, and nobody has to learn it.**

The coverage is complete by policy, not by accident. Every operation a user can invoke is also an endpoint, written in the same change as the feature rather than bolted on later. That is a standing rule of the project.

Run `splicecraft --headless` and there is no interface at all: just a local HTTP server, a readiness probe, 303 endpoints, and `/tools`, which hands back any endpoint's request schema on demand. So a frontier model needs no prior knowledge of SpliceCraft to drive it correctly. It can discover the whole surface while it works, and the published skill supplies the conventions a schema cannot express.

In practice this means a scientist can describe the construct they want in plain language, to a model of their choosing, and have the design performed, simulated, verified and filed — without opening the program, learning a keybinding, or reading documentation. The operations are the same ones the interface performs, with the same interlocks, the same refusals and the same provenance, because they are literally the same code paths.

Several consequences matter specifically for the audiences of §1.3:

- **The learning curve becomes optional rather than mandatory.** The barrier the literature identifies for wet-lab scientists is not capability but the absence of an on-ramp [5]. A conversational on-ramp to a rigorous engine is a different proposition from a command-line tool with a manual.
- **It sidesteps the language barrier.** A 2026 review of African bioinformatics capacity names English dominance in "training, online resources, and professional networks" as a neglected barrier to equitable access, observing that English-only delivery "increased uptake in Anglophone countries but limited participation from Francophone Africa," and that remedies "require deliberate strategies such as translation and localization of core materials" [12]. SpliceCraft's own interface strings and documentation are in English and this paper makes no claim to have solved that; what it offers is a second door. A frontier model is multilingual even when the tool is not, so a researcher can describe what they want in their own language and operate an engine that has no language at all.
- **It is an accessibility layer.** A text interface driven by keyboard shortcuts is not equally usable by everyone. A conversational interface to the same endpoints is usable by a screen-reader user, by someone for whom precise keyboard navigation is difficult, and by anyone who would simply rather describe the outcome than operate the tool.
- **The hardware bill stays trivial.** A headless instance on a single-board computer plus a model the user already has access to is a complete design workstation. The expensive component — inference — is the one the user can choose to run in the cloud, locally (§7.3), or not at all, and in every case the sequences stay on the machine unless an online lookup is explicitly enabled.

The honest caveat is that this inherits whatever the model brings. A model will occasionally ask for the wrong operation, and §7.3's interlocks and §5's refusals are what stand between that and a damaged library. The design assumption throughout is not that the model is reliable; it is that the instrument must be.

### 7.3 A local assistant, with interlocks

The in-application assistant talks to a **local Ollama model**: no cloud service, no API key, no sequence leaving the machine, and no additional heavyweight dependency (the client is standard-library HTTP). It streams Markdown answers, exposes a context budget, and can optionally ground answers in a research corpus with citations, grown by a built-in paper scraper or by indexing the user's own library into a pack that is private by construction.

With agent mode enabled, it drives the same endpoints the API exposes, with live visible effect in the application. The interlocks:

- **It asks before every write** by default.
- **Destructive whole-library operations are not reachable at all** — not gated, not confirmed, absent.
- **Physical robot motion always asks**, regardless of mode.
- **Online database lookups are off until the user enables them**, and even then only the query string is transmitted, never the sequence.
- **Anything read from the web is treated as data, never as instructions**, and once the assistant has read external content, even hands-off mode reverts to asking before it changes user data.
- **The turn ends by stating what actually happened.** If a requested change failed or never ran, the turn says so regardless of what the model's prose claims — the report comes from the application, not from the model.
- A model without native tool calling is given a JSON protocol the runtime enforces, so capability degrades rather than failing.

A benchmark script grades local models on everyday tasks, because "which small model can actually do this" is an empirical question that depends on the user's hardware.

---

## 8. Educational and capacity-building value

### 8.1 Deployment properties that matter in a classroom

No cost, no accounts, no network. One Raspberry Pi on the local network can serve the same interface to a whole class of browsers from a single process, which turns machines nobody is allowed to install software on into workstations: school-managed laptops, Chromebooks, a locked-down lab image. No student signs up with a third party. No student data leaves the room. The installation is 202 MB and the program holds about 87 MB while running, so the oldest machine in the room is still a candidate. For a teacher without reliable connectivity, none of the core workflow needs a connection at all.

### 8.2 Pedagogical properties that are not accidents

Several design decisions are better teaching devices than they are features.

**The mechanism is visible.** The superscript cut count beside each multi-cutting enzyme ticks down as a student edits a site out of the sequence. That is a causal loop a student can close themselves, in seconds, without an instructor — and it teaches what a recognition site *is* far more durably than a diagram does.

**The provenance is a protocol.** Since the program simulates each step instead of asserting the result (§2.3), a construct's History reads as a numbered bench recipe for what the software just did. A student designs the construct on screen, reads the protocol the program wrote, and goes and does it. The two cannot drift apart, because one produced the other. The generator writes only steps that were actually recorded, so there is no invented transformation hiding in it.

**The refusals teach calibration literacy.** "This motif occurs in 94.8% of random 400 bp windows, so finding one means nothing" is a statistical idea that most undergraduate curricula never make concrete, and here it arrives attached to the student's own sequence. So does "this scan's threshold is set at a 1% false-positive rate, which means it finds about 46% of real donor sites" — i.e. that a negative result has a sensitivity, and that tools have operating points. A tool that modelled confident silence instead would be teaching the opposite lesson.

**It is readable.** The implementation is pure Python with a flat module layout and a handover document that explains the invariants and the pitfalls. A motivated student can read the function that computes the thing they just clicked. For capacity building this is the difference between a tool and a curriculum: the barrier reported for wet-lab scientists in the literature [5] is not capability but the absence of an on-ramp, and an inspectable integrated tool is an on-ramp.

### 8.3 Sustainability

Every release is archived on Zenodo with its own version DOI, so a methods section can cite the exact version that produced a design, and the archive survives whatever happens to the repository, to me, or to my institution. Data sits in human-readable JSON at a documented path. The licence is MIT. Formats go both ways. So the project's failure mode is "unmaintained", not "inaccessible". Whoever outlives it keeps their library, their files, and the ability to read the code that made them.

---

## 9. How this was built: LLM-assisted development, and the discipline that makes it legitimate

### 9.1 Disclosure

SpliceCraft was written by one practising bench biologist, who is not a professional software engineer, with heavy and continuous help from large language models, over about two years. Nothing of this size — 222,000 lines, 39 modules, 12,277 tests, a 303-endpoint API — comes out of one self-taught developer at this quality in this time without that help.

I say so plainly and early, for two reasons. Hiding it would be dishonest in a paper whose argument is about calibrated honesty. The more useful reason: thousands of people are answering the question "can scientific software be built this way?" in practice right now, mostly without publishing an account of what goes wrong. I happen to have an unusually complete record of what went wrong here, since every audit sweep is numbered and written down, including the sweeps that caught the previous sweep's mistakes. That record is the most transferable thing this project has.

My position is narrow and conditional. **LLM assistance is a large, real capability multiplier for scientific software. It is safe for scientific software only alongside edge-casing, hardening and auditing considerably more aggressive than normal practice.** Not "tests". Audits that start from the assumption that the code is wrong and set out to prove it.

### 9.2 Generated code fails differently, and that is the whole argument

A tired human writes code that is wrong in visible ways. It is half-finished, or clumsy, or it does not compile, and the variable names give away where the confusion set in.

A language model writes code that is wrong in plausible ways. It compiles. It reads idiomatically. It carries a docstring describing what it was supposed to do. It handles the three cases you mentioned and quietly mishandles the fourth. It reaches for the right-looking helper from the wrong layer. Asked for an answer, it returns a number in the place where the correct output was a refusal.

In most software, plausible-but-wrong is a bug. In scientific software it is worse than that, because the person reading the output cannot check it. That inability is why they opened the tool. A map that draws a primer four bases from where it anneals does not look wrong. A splice-site scan that returns an empty list because its model failed to load looks exactly like a clean sequence. The user's only real defence is the tool's own willingness to say no — which is the last property generated code is likely to have, because saying no is not what anyone prompted for.

So §5 and this section are one argument told twice. Calibrated honesty is what the software owes its user. Adversarial auditing is what the developer owes the software.

### 9.3 A failure taxonomy, with instances from this codebase

Each of the following was found in this project, is recorded in its audit log, and recurred more than once. They are listed because I have not seen them enumerated elsewhere with specifics, and because each one implies a specific countermeasure.

**(a) Code that exists, passes tests, and is reached by nothing.** The most striking instance: an entire eukaryotic-cassette module was packaged into the wheel, registered at its architectural layer, covered by **38 passing tests**, and reachable from **no production caller and no endpoint**. The only non-test mentions of it anywhere in the tree were its own docstring and two prose comments. It had never run for a user. (The project's own notes recorded "100+ tests" for this episode; the committed file at that commit holds 38 test functions and no parametrisation.) The tests passed because the tests imported it directly. Separately, the laboratory notebook maintained three search indexes on every save that **nothing ever read**. A third instance surfaced during the final fact-check of this paper: a bounded-`repr` helper sitting in the hub, covered by its own tests, duplicating a function that already lives in the logging module and is wired into the event path — the live one does the work, the hub copy is called by nothing, and the two have different default limits. That near-duplicate is also what made the §9.3(g) story below take four passes to land, because an auditor who greps for the orphan finds a dead function and concludes the bound does not exist. This shape is near-invisible: coverage is high, tests are green, the code is good, and the feature does not exist. Generated code is unusually prone to it because a model asked to build a subsystem will build the subsystem, and wiring it in is a different request. *Countermeasure:* test reachability explicitly, by asserting that a production entry point calls the thing. Then periodically enumerate the public surface and ask what calls each member of it.

**(b) Tests that assert the fix instead of the behaviour.** A review round once ticked 26 items as fixed. Re-verified item by item, **7 were complete, 7 were partial, and 12 were untouched.** So 19 of 26 were not actually done. Four distinct causes were identified, and only one of them is the obvious one: ticking in bulk; assertions on *source strings* (one checked that a method name appeared in a function's source, and stayed green after the logic moved and the bug returned); several "fixes" that were themselves regressions; and one *behavioural* test that passed vacuously because its fixture did not reproduce the shape the auditor had reported. *Countermeasure, now a rule in the contributor document:* **a regression test counts only once it has been seen failing on a pre-fix snapshot of the code.** Checking out the parent commit and watching the test go red is the only evidence that it tests anything — and the fourth cause above is why even a behavioural assertion needs that check.

**(c) Harnesses that stub the mechanism under test.** A hand-rolled test rig built the library widget via `__new__`, wired it to a hand-written stand-in for the application, and replaced the very undo-push method under test with a one-line append. Seven tests passed. They were proving nothing, and they concealed a **real data-loss bug in bulk deletion**: the undo stack was capped at 25 entries, so deleting 30 rows and pressing undo restored 25 of them *and reported success*. The fix raised the cap, made a batch push atomic, and discarded any batch the cap had broken, because a partial restore that announces completion is worse than an honest refusal. *Countermeasure:* stub I/O only — notifications, disk, the screen stack — never the logic being asserted.

**(d) Bugs inside the fixes.** This is the one I was least prepared for. An audit of the data-safety layer returned four findings, and **three of the four were inside the previous audit's data-safety fixes.** Every identifier involved had been introduced by that earlier repair commit and left untouched since. The sharpest: a concurrency barrier meant to prevent a race inspected a queue slot that the worker pops *before* it takes the lock, so for essentially the whole life of every deferred write the barrier answered "nothing in flight." In a different audit, two consecutive hardening rounds found mostly defects in the preceding round's repairs. (The fourth finding is *not* attributable to the earlier fixes. It is a gap in the delete-side write guard that predated them by about five months. The project's own notes credit all four to the prior round; they are wrong, and this paper inherited that error until the claim was checked against the commit dates.) *Countermeasure:* the fix is new code and inherits the full suspicion of new code. Audit the patch, not just the bug.

**(e) Reimplementing what already exists under another name.** In one session, four of the five capabilities proposed as missing **already existed**. They were spelled under internal names that no search for the feature's own name would ever surface. *Countermeasure, now a rule:* grep for the mechanism, meaning the verb, the data structure, or the endpoint list. Never the feature name. Every time, the real gap turned out to be narrower and more useful than the feature nearly duplicated on top of it.

**(f) Confident, wrong accounts of the cause.** This one recurs, and it comes from humans as readily as from models. Across several field reports and generated diagnoses, the symptom was real and the stated cause was wrong. Act on the stated cause and you get a plausible fix that fixes nothing. *Countermeasure:* reproduce before accepting any diagnosis, your own included, and most of all one that arrives already explained. The fact-check of this manuscript produced a fresh example, recorded in the drafting notes.

**(g) Prose drifting into fiction.** A module docstring described its scoring thresholds as anchored at 95% sensitivity. They were anchored at a 1% false-positive rate, which is a different and much less reassuring thing. The documentation overstated the scan's reach two- to threefold, and said so for months. A second instance was found while fact-checking *this manuscript*, and it took four passes to reach the bottom of — three of them wrong, in different directions. The architecture notes stated that "sequence content is never logged." The first correction weakened that to a leak bound rather than an absence, and §2.1 of the second draft said so. An adversarial pass then found a bounded-`repr` helper with **no production caller**, concluded there was no bound at all, and §2.1 of the third draft said *that* — the overcorrection, written with more confidence than either claim it replaced. The fourth pass found the actual mechanism: a second, near-identically-named function in the logging module, wired into the structured event path, which does bound every field and renders a Biopython sequence object as an opaque tag. The truth was a split guarantee that none of the three previous statements described. Each wrong version was produced by reading the previous one, or by grepping for a name taken from the previous one, instead of asking what the logging path actually calls. Nothing checks prose, including a preprint's, including a preprint's correction to its own correction — and an auditor armed with the wrong function name is no safer than the prose that gave it to them. *Countermeasure:* where a number or a guarantee matters, make it **executable**. Ship the calibration as machine-readable metadata beside the model. Have a test assert the properties that prove the coordinates are right, and have the calibration script exit non-zero when they are not. The number in the paper and the number in the program should come from the same place, and the claim in the paper should be read off the function, not off the summary of the function.

**(h) Intermittent failures dismissed as flakes.** A test that failed under parallel execution and passed eight times out of eight in isolation was assumed flaky. Reading what the assertion actually received, an empty string where a typed value belonged, revealed a **real race**, and in all four of the application's search bars: tearing the widget down changed its value, which posted a fresh change event that armed a *new* debounce timer with nothing left to cancel it. The teardown re-armed the timer it had just cancelled. *Countermeasure:* "passes in isolation" is not evidence of a flake. Read the value the assertion got.

**(i) Audits that trust the previous audit.** Thirty-odd numbered sweeps into the project, one pass was run differently: re-derive every behaviour from the code, and refuse to treat the accumulated invariants document as authoritative. It found real bugs that all the earlier sweeps had missed. They had missed them by inheriting each other's conclusions. *Countermeasure:* audit from the source now and then, and treat your own documentation as a hypothesis.

### 9.4 What actually worked

The countermeasures that earned their keep, in rough order of value returned:

1. **Differential oracles: a second implementation, not a second opinion.** The highest-yield technique here by a wide margin. Ask a model to review its own output and you find style problems. Ask whether two independently derived implementations agree and you find correctness problems. The CRISPR scanner is checked against a brute-force reference written separately from the engine, 24 of 24 linear and 18 of 18 circular. Exports are validated against spec-derived validators and then re-read by parsers from outside this project. The replacement GenBank reader had to produce records identical to Biopython's across 313 reference records, the whole production library and 180 writer round-trips before it was allowed anywhere near a user (§3.3). An oracle you did not write yourself is the only cheap ground truth available.

2. **Generated edge cases, not imagined ones.** Property-based testing on the biology layer and mutation fuzzing on the parsers turned up defects no amount of careful thinking produced, among them a quadratic string concatenation and an unbounded recursion. Ask anyone to list edge cases, human or model, and you get the memorable ones. A fuzzer gives you the ones that matter.

3. **Measurement instead of assertion.** Every threshold in §5 is anchored on a measured error rate. Every gated operation in §3.2 has a ceiling and a stated procedure. Every refusal in Table 3 carries the number that caused it. That is what turns "the model said this works" into something a reader can check for themselves. It also catches you deceiving yourself with your own benchmark: one cache left warm made the digest measure 0.01 ms, which is to say it measured nothing at all.

4. **Executable invariants with stable names.** Ten sacred rules plus a tagged catalogue of pitfalls, each backed by a regression test. The tag does more work than it sounds like it would. It gives a model something greppable to retrieve: a specific, written reason why the obvious-looking change is wrong. That is precisely the knowledge that does not survive into a fresh context window.

5. **A handover document addressed to whoever comes next.** The repository's contributor file is written for the next arrival, human or model. It carries the rules, the traps that keep recurring, and the incidents behind them, including the day an unsandboxed diagnostic script overwrote my real library (§4.3). Keeping the *reason* for a constraint attached to the constraint is what stops someone deleting it as redundant.

6. **Static typing as a separate oracle.** A type checker catches null-dereference defects the test suite does not, because the two fail in different directions. The auditable instance here: the gate's own file list had drifted until it covered four of roughly twenty-seven modules, and **24 latent optional-type errors had piled up behind the gap**. Correcting the list found them. No test ever would have. The same gate has aborted two releases outright. (Project notes carry a sharper anecdote, five defects missed by a stated test count, but neither number is reproducible from the repository, so the checkable one appears above instead. That substitution is itself an instance of §9.3(g).)

7. **Scoped adversarial audits run in parallel, with written scope boundaries** — followed by a pass over the audit's own fixes, per §9.3(d).

None of items 1 through 7 is novel. They are ordinary practice in safety-critical software. The honest summary of this section is not that LLM-assisted development needs some special discipline. It is that **generating code removes the rate limit that used to keep undisciplined projects small.** Writing 200,000 lines by hand takes long enough that you learn the domain's traps on the way through. Generating them does not. Same discipline as always; what changed is the volume of plausible output, and that is what makes skipping it newly expensive.

### 9.5 So is it legitimate?

Yes, on a condition that cannot be met with enthusiasm.

The capability gain is real. This program would not exist otherwise, and the honest counterfactual is not a more carefully hand-written version of it. The counterfactual is nothing at all, which is also the counterfactual for every tool this field does not have.

But the defects arrive looking correct, and the people using scientific software cannot spot a correct-looking defect in the one domain they came for help with. That puts the whole weight on the developer's willingness to go after their own work. Build this way without differential oracles, generated edge cases, measured thresholds, executable invariants and audits of the audits, and you have not built a fast version of a good project. You have built a confident version of a wrong one. In biology the cost of that is somebody's experiment, somebody's month, and somebody's reagent budget.

The line I would draw: **if you cannot say what measurement would prove your tool wrong, you have not finished building it.** That holds for the code, the thresholds, the documentation, and this paper.

## 10. Maintenance as an obligation, and what LLM assistance changes about it

§1.5 argued that abandonment is a structural outcome of an incentive system, not a moral failing of the people caught in it. This section is what I think follows.

### 10.1 Discharging the second duty: building to be inherited

If a widely-used tool carries an obligation either to maintain it or to hand it over cleanly, then "hand over cleanly" needs an operational definition. Mine, implemented:

- **A licence that permits a successor.** MIT, not "source available".
- **Independent archival of every release.** Each release is deposited to Zenodo with its own version DOI, so the artefact that produced a published design outlives the repository, me, my institution and my hosting bill. A methods section can cite the exact version, and that citation keeps resolving after the project stops.
- **Distribution through a maintained package manager.** In the only study to measure it — on a thin sample, but the only one there is — that is the difference between installable and not [6] (§2.5).
- **User data in human-readable JSON at a documented path.** Nothing about recovering a library should require the program that wrote it. Whoever outlives this project keeps theirs.
- **Bidirectional format support.** Data can leave. Nothing is held in place by export friction.
- **A test suite a stranger can change the code with.** 12,277 tests are not a quality badge. They are the mechanism by which somebody who has never met me makes a change and finds out whether they broke it. Of everything on this list, that is the one that matters most.
- **Binding compatibility promises.** Since 1.0 the project commits to semantic versioning, backward compatibility at the file-format, API and CLI levels, a migrator for every data-schema change, and a two-version support window. Written down as a gate, not an intention.
- **Written instructions for the successor**, carrying the invariants, the traps and the incidents behind them (§9.4, item 5).

Put together, that changes the project's worst case from "inaccessible" to "unmaintained". For a user those are very different outcomes, and only one of them loses the science.

### 10.2 What LLM assistance changes

This is where the two halves of the paper meet, and it is the most important claim in either.

The constraint behind the numbers in §1.4 was never that maintenance is intellectually hard. It is that **maintenance costs maintainer-hours nobody will fund**, and the people who hold the domain knowledge are the same people whose careers have moved on. A fixed cost per hour, multiplied by a reward of zero, produces abandonment as regularly as the studies measure it.

LLM assistance changes the first term. It does nothing about the incentive problem, and neither does anything else in this paper. What it changes is the cost to one domain expert, working alone in the evenings with no engineering staff and no grant line, of keeping a tool alive, correct, and answering bug reports from people at a bench.

That cost is what has been killing these tools, and it has just dropped. It is also the strongest case I can make for §9's method. Not that the method writes software faster. That it may put *sustained* scientific software within reach of the people who actually use it, who are the same people the current funding model has no way to pay.

The condition from §9.5 carries over undiminished. A field full of tools generated fast and audited never would be worse than the field we have now, not better. It would swap software that is honestly absent for software that is confidently wrong. The claim here is narrower: that generated leverage, plus adversarial discipline, plus building to be inherited, is together a workable answer to a problem the literature has been measuring for fifteen years without solving.

### 10.3 Why the usual abandonment mechanism does not apply here

The mechanism those studies describe is specific. The tool was a means to an end: a thesis, a paper, a contract. When the end arrived, the means stopped being maintained [8]. Nobody lost interest in software. They lost the reason the software existed for them.

That mechanism is structurally absent in this case, for a reason worth stating because it is the most credible sustainability claim this paper can make. **SpliceCraft is my sole plasmid design tool.** Every construct designed at this bench is designed in it. As a direct consequence of the program reaching that standard, I have **divested entirely from paid design software** — there is no subscription running in parallel, no licensed alternative held in reserve, and no fallback to retreat to if maintenance lapses.

So the maintainer is also the most dependent user, and what I depend on is daily bench work, not a deliverable with an end date. That is a structural fact, not a promise of diligence. The bug reports come out of my own cloning, successful and otherwise, which is why they tend to be specific and why they tend to get fixed fast. The person who has to live with the bug is the person who can fix it, and I have to live with it tomorrow.

The visible evidence of that alignment, as of this writing: **344 tagged releases since the first public version on 2026-04-12, 39 of them in the last 90 days**, with the version described here released the day before this manuscript was prepared. Releases typically ship in the same week a problem surfaces at the bench.

Two honest qualifications. First, **this is a work in progress and is presented as one.** It is actively developed, it is not finished, and the roadmap is driven by what fails at the bench rather than by a plan; §11 lists what is missing and §5.2 lists what was tried and refused. The maintenance commitment is indefinite, not the feature set. Second, structural alignment is not immortality: a single-maintainer project still carries single-maintainer risk, and this argument changes the *mechanism* of that risk from incentive expiry — which is predictable, common, and arrives on a schedule — to the ordinary contingencies of one person's life, which is a different and in practice much longer horizon. §10.1 is the hedge against the part that remains, and it is deliberately built so that the hedge does not depend on anyone's good intentions, mine included.

## 11. Limitations

Stated plainly, because a paper arguing for honest calibration cannot be coy about its own.

**Scale.** SpliceCraft is a plasmid and small-genome workbench. Visualisation scales further than expected, and for an instructive reason: the binding constraint was feature *count*, not sequence length. The redraw path is now flat in feature count — measured at 0.0013, 0.0014 and 0.0012 ms per frame for 10,000, 200,000 and 1,000,000 features on a 10 Mbp molecule (Table 2) — and that flatness, not the absolute figure, is what a regression test enforces. **First load is not flat:** building the ring is linear in features, and opening a 10 Mbp / 200,000-feature record took 8.3 s on the test machine. So a chromosome can be *explored* smoothly once open, and opening it is a wait.

**Storage is the open wall, and my figure for it is weak.** A 1 Gbp GenBank file is on the order of 1.3 GB of text, with parse and write times on the order of tens of seconds — inherent to a text format storing 60 bases per line with position numbers. I give these as orders of magnitude on purpose: they come from a one-off measurement recorded in project notes, nothing in the repository reproduces them, and I did not re-measure at that scale for this paper. Routine chromosome-scale work needs a binary or indexed store; that is identified and deliberately not yet built, and anyone depending on the figure should measure it themselves.

**Continuous integration covers one Python version, not the five that are supported.** The package declares support for Python 3.10 through 3.14 and development happens on 3.14, but CI runs the full suite on **3.12 only**; a four-version matrix (3.10–3.13) ran previously and was narrowed for turnaround time. Declared support therefore rests on language-level compatibility and on two versions being exercised in practice, not on continuous evidence across all five. This is a real gap and not a theoretical one: the project has already been bitten by precisely this shape, when a dependency that installed on some platforms and not on the release machine produced incorrect alignment results on the platforms the release gate never exercised — the defect that caused that dependency to be removed outright (§3.3). A matrix would have caught it earlier. Readers on 3.10, 3.11 or 3.13 should treat those versions as supported-by-declaration.

**Native Windows is unverified on hardware.** The UTF-8 console and VT setup is exercised only through mocked system calls on Linux CI, because no Windows runner is available. The platform matrix says "should work," not "verified," and will keep saying so until someone confirms it on real hardware. The Windows Subsystem for Linux path is fully supported.

**One ARM dependency needs a compiler.** `primer3-py` ships no `aarch64` wheel, so on ARM Linux and Apple Silicon that one dependency compiles at install and a C toolchain must be present first. Without it primer design degrades rather than failing — melting temperatures fall back to the 2(A+T)+4(G+C) rule and the pair QC states that the oligos are not secondary-structure checked — but that degradation is real enough that `primer3-py` stays a hard requirement. Everything else arrives as a wheel. A release-gate script tracks exactly which dependencies lack wheels on which platforms, so this statement stays current rather than aspirational.

**No genome-wide off-target search, no efficiency score, no poly(A) prediction, no promoter-strength score.** See §5.2; these are refusals, not roadmap items, and each is published with its number.

**Single-maintainer bus factor.** One developer, which is the risk ApE carried through seventeen productive years [13]. §10.3 argues that the usual abandonment mechanism cannot operate here, since my own bench work runs on the tool and the paid alternatives are cancelled. But that argument covers incentive expiry, not contingency. For what is left, the mitigations are built not to depend on anyone's good intentions, mine included: an MIT licence, every release archived on Zenodo beyond the reach of my institution, user data as plain text at a documented path, binding compatibility promises, 12,277 tests so an unfamiliar contributor can change something and find out at once whether they broke it, and instructions written for the next maintainer rather than the current one (§10.1).

**Work in progress.** The feature set is not stable and is not presented as stable; it is driven by what fails at this bench. Readers should take §6 as a description of version 1.3.19 and §5.2 as the current state of what has been tried and refused, both of which will move.

**Not peer-reviewed in parts.** The trained models in §5 are calibrated and their calibrations are published and retrievable at runtime, but they have not been externally validated by an independent group. The refusals are the honest version of that caveat; the positive claims should be read as measured-by-the-author until someone else measures them.

---

## 12. Conclusion

The hard computational parts of plasmid design were solved and published long ago. What was missing was never an algorithm. It was somewhere to put them all together that did not also require an account, a subscription, a browser, a connection, a recent machine, and faith in a company's continued interest. ApE showed that one working scientist can supply that much for the design step and then keep it alive for twenty years [13,14]. SpliceCraft tries to cover the rest of the workflow — the assembly, the verification, the notebook, the automation — in the most portable interface anyone has, built so that the data outranks the features and no number ships without the measurement behind it. For the people who would rather not learn an interface at all, every operation also answers in plain language.

I am not claiming SpliceCraft beats the commercial platforms at what they do. The claim is narrower. Recurring charges, and the dependencies that ride along with them, are not an inherent property of integrated tooling. For someone on a shared server, a teacher with a room of ten-year-old laptops, a company of two, or a lab where the power and the bandwidth both come and go, a program that is small, local, repairable and candid about its own limits is the better instrument. Not the compromise option. The better one.

Two larger things follow from having built it.

The field's software layer is eroding, at a rate several groups have now measured, for reasons nobody involved chose. Better intentions will not fix that. Tools built to be inherited might: licensed so a successor is allowed, archived somewhere an institution cannot withdraw it from, storing data in a form that can be read and mended by hand, and carrying enough tests that someone who has never met me can make a change and know whether they broke it.

And for the first time, a tool of this size is maintainable by the person who actually uses it at the bench. That matters because the bench is precisely the population the funding model excludes, and so it is precisely where the abandonment numbers come from. The capability is real. It is also conditional. Generated code fails by looking right, scientific users cannot catch code that merely looks right in the domain they came for help with, and so this way of working is only honest if the developer attacks their own output harder than a reviewer would. If you cannot say what measurement would prove your tool wrong, you are not finished building it.

This one is offered freely, in the state described above, to a community that taught an outsider how to do the work and never once asked what I was qualified for.

---

## Availability and requirements

- **Project home:** https://github.com/Binomica-Labs/SpliceCraft
- **Project site and browser demo (no installation):** https://splicecraft.bio
- **Package:** PyPI, `pipx install splicecraft`
- **Archive:** Zenodo, concept DOI [10.5281/zenodo.21960400](https://doi.org/10.5281/zenodo.21960400) [24] (resolves to the newest release; each release also receives its own version DOI)
- **Version described here:** 1.3.19 (released 2026-11-10)
- **Licence:** MIT
- **Operating systems:** Linux, macOS, Windows Subsystem for Linux (fully supported); native Windows (implemented, unverified on hardware); ARM Linux including Raspberry Pi OS 64-bit (one-time C toolchain for one dependency)
- **Requirements:** Python ≥3.10 (declared support 3.10–3.14; continuous integration runs the suite on 3.12, development on 3.14 — see Limitations); a UTF-8-capable terminal with a braille-glyph font for the graphical map, with automatic 7-bit ASCII fallback otherwise; 256-colour ANSI
- **Citation:** `splicecraft --citation` emits an APA reference and BibTeX pinned to the installed version; the repository ships `CITATION.cff`

## Data availability

No new experimental data are reported. The software, its test suite, its calibration scripts and its documentation are in the repository under the MIT licence. Model calibration figures quoted in §5 are retrievable at runtime from the corresponding endpoints (`regulettuce-calibration`, `splice-model`, `initiation-model`) and are embedded in the shipped model metadata. Reference sets for the trained models are regenerable from public sequence databases with the included calibration scripts. Performance figures in Table 2 are reproducible with `tests/test_perf_regression.py` and the measurement procedure stated in the caption.

## Competing interests

The author declares no competing interests. SpliceCraft is distributed free of charge under an MIT licence; the author derives no revenue from it, holds no financial interest in it, and holds no interest in any competing product. Two disclosures are made in the interest of completeness rather than because they constitute a conflict: development is supported in part by small voluntary user donations (see Funding), and the author is himself the software's primary user, which is the basis of the sustainability argument in §10.3 and should be read as the interested testimony it is.

## Funding

This work received no institutional, governmental, or commercial funding. Development has been supported in part by small voluntary donations from users via Ko-fi (https://ko-fi.com/atinygreencell). No funder had any role in the design, implementation, evaluation, or reporting of the software or of this manuscript.

## Acknowledgements

This program is a thank-you note.

I came to molecular biology from outside it — an amateur, self-taught, without a degree in the field or a laboratory to belong to — and the research community let me in anyway. Scientists answered emails they had no obligation to answer. People sent me plasmids, protocols, reprints, and corrections. Nobody asked for my credentials before helping me, and a great many people helped.

My largest debt is to the **#biotwitter** community: the loose, sprawling, substantially unpaid network of bench scientists, principal investigators, graduate students, technicians, hobbyists and biohackers who answer questions from strangers at two in the morning for no reason other than that someone asked. They explained Type IIS overhang logic to me patiently and more than once. They told me when I was wrong, which is the more valuable service. They sent me constructs, reagents, reading lists and encouragement, and several of them have since filed bug reports against this software with the care of a peer reviewer and none of the credit. That community has scattered across several platforms over the last few years and is not what it was, which makes acknowledging it feel less like a formality and more like getting it on the record while the record still means something. It taught an outsider how to do this work. Everything good in this program came from there first.

SpliceCraft exists because that generosity turned out to be transitive: the only sensible thing to do with knowledge you were handed for free is to build something and hand it on. If this program saves a graduate student a failed cloning week, or lets a high-school class design a real construct on a laptop their school already owns, or lets a laboratory somewhere with unreliable power and metered bandwidth do its own design work on its own hardware without asking anyone's permission — then it has paid back the part of the debt that can be paid.

Particular thanks to **M. Wayne Davis**, whose ApE proved that one person who actually does the bench work can build the tool the bench needs, and who has been giving it away for over two decades. This program is built in that tradition and tries to be worthy of it.

Finally, and with some sincerity, I thank the software-as-a-service industry, whose pricing department supplied the activation energy for this entire project. No funding panel has ever been as motivating as a renewal notice, and no reviewer has ever requested a feature as persuasively as moving it behind a higher tier. The subscriptions have since been cancelled (§10.3). Any errors remaining in this work are my own; only the motivation was outsourced.

---

## References

> **Verification status.** References **1–9 and 12–22** were confirmed field-by-field during drafting — authors, title, journal, volume, issue, pages, year and DOI exactly as printed below. References 1–9 and 12–15, 18–21 were checked against PubMed records; 16 and 17, being computer-science papers not indexed in PubMed, were checked against DBLP and their ACM DOIs. Every quantitative figure quoted from references 6, 7, 8 and 9 in §1.4–§1.5 was read from those articles' own abstracts or results text, never from a secondary summary or a search snippet. The list prices and licence terms in Table 1 were read from the vendors' own pricing pages on the access date given in that table's caption, not from third-party comparison sites, which disagreed with each other and with the vendors. Only references 10 and 11 remain marked ⟨check⟩. Attribution note: bibliographic metadata for the PubMed-confirmed references was retrieved from PubMed.

1. Cock PJA, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, Friedberg I, Hamelryck T, Kauff F, Wilczynski B, de Hoon MJL. Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics*. 2009;25(11):1422–3. doi:[10.1093/bioinformatics/btp163](https://doi.org/10.1093/bioinformatics/btp163)

2. Untergasser A, Cutcutache I, Koressaar T, Ye J, Faircloth BC, Remm M, Rozen SG. Primer3 — new capabilities and interfaces. *Nucleic Acids Research*. 2012;40(15):e115. doi:[10.1093/nar/gks596](https://doi.org/10.1093/nar/gks596)

3. Larralde M, Zeller G. PyHMMER: a Python library binding to HMMER for efficient sequence analysis. *Bioinformatics*. 2023;39(5):btad214. doi:[10.1093/bioinformatics/btad214](https://doi.org/10.1093/bioinformatics/btad214)

4. Altschul SF, Gish W, Miller W, Myers EW, Lipman DJ. Basic local alignment search tool. *Journal of Molecular Biology*. 1990;215(3):403–10. doi:[10.1016/S0022-2836(05)80360-2](https://doi.org/10.1016/S0022-2836(05)80360-2)

5. Laub V, Devraj K, Elias L, Schulte D. Bioinformatics for wet-lab scientists: practical application in sequencing analysis. *BMC Genomics*. 2023;24(1):382. doi:[10.1186/s12864-023-09454-7](https://doi.org/10.1186/s12864-023-09454-7)

6. Mangul S, Mosqueiro T, Abdill RJ, Duong D, Mitchell K, Sarwal V, Hill B, Brito J, Littman RJ, Statz B, Lam AK-M, Dayama G, Grieneisen L, Martin LS, Flint J, Eskin E, Blekhman R. Challenges and recommendations to improve the installability and archival stability of omics computational tools. *PLoS Biology*. 2019;17(6):e3000333. doi:[10.1371/journal.pbio.3000333](https://doi.org/10.1371/journal.pbio.3000333)

7. Kern F, Fehlmann T, Keller A. On the lifetime of bioinformatics web services. *Nucleic Acids Research*. 2020;48(22):12523–33. doi:[10.1093/nar/gkaa1125](https://doi.org/10.1093/nar/gkaa1125)

8. Schultheiss SJ, Münch M-C, Andreeva GD, Rätsch G. Persistence and availability of Web services in computational biology. *PLoS ONE*. 2011;6(9):e24914. doi:[10.1371/journal.pone.0024914](https://doi.org/10.1371/journal.pone.0024914)

9. Song J, Kurgan L. Availability of web servers significantly boosts citations rates of bioinformatics methods for protein function and disorder prediction. *Bioinformatics Advances*. 2023;3(1):vbad184. doi:[10.1093/bioadv/vbad184](https://doi.org/10.1093/bioadv/vbad184)

10. Panji S, Ras V, Gichoya J, Julius R, Maslamoney S, Motshegwa T, Sharaf A, Skelton M, Tunga M, Mulder N. Unleashing potential: assessing Africa's readiness for the data science revolution to impact health. *Nature Communications*. 2026;17(1):5138. doi:[10.1038/s41467-026-71454-4](https://doi.org/10.1038/s41467-026-71454-4)

11. Ofosu W, Kotey FCN, Odoom A, Ntim OK, Donkor ES. The potential of bioinformatics in revolutionizing disease diagnosis in sub-Saharan Africa: a narrative review. *Health Science Reports*. 2026;9(4):e71912. doi:[10.1002/hsr2.71912](https://doi.org/10.1002/hsr2.71912)

12. Hussein J, Myovela H. Reframing bioinformatics capacity development in Africa: from training supply to research-driven demand. *Journal of Community Genetics*. 2026;17(2):36. doi:[10.1007/s12687-026-00877-9](https://doi.org/10.1007/s12687-026-00877-9)

13. Davis MW, Jorgensen EM. ApE, A Plasmid Editor: a freely available DNA manipulation and visualization program. *Frontiers in Bioinformatics*. 2022;2:818619. doi:[10.3389/fbinf.2022.818619](https://doi.org/10.3389/fbinf.2022.818619)

14. Davis MW, Jorgensen EM. Using ApE for in silico Golden Gate cloning. *Methods in Molecular Biology*. 2025;2850:79–87. doi:[10.1007/978-1-0716-4220-7_5](https://doi.org/10.1007/978-1-0716-4220-7_5)

15. Helling RB, Goodman HM, Boyer HW. Analysis of endonuclease R-EcoRI fragments of DNA from lambdoid bacteriophages and other viruses by agarose-gel electrophoresis. *Journal of Virology*. 1974;14(5):1235–44. doi:[10.1128/JVI.14.5.1235-1244.1974](https://doi.org/10.1128/JVI.14.5.1235-1244.1974)

16. Nakamura Y, Gojobori T, Ikemura T. Codon usage tabulated from international DNA sequence databases: status for the year 2000. *Nucleic Acids Research*. 2000;28(1):292. doi:[10.1093/nar/28.1.292](https://doi.org/10.1093/nar/28.1.292)

17. Myers G. A fast bit-vector algorithm for approximate string matching based on dynamic programming. *Journal of the ACM*. 1999;46(3):395–415. doi:[10.1145/316542.316550](https://doi.org/10.1145/316542.316550)

18. Hirschberg DS. A linear space algorithm for computing maximal common subsequences. *Communications of the ACM*. 1975;18(6):341–3. doi:[10.1145/360825.360861](https://doi.org/10.1145/360825.360861)

19. Weber E, Engler C, Gruetzner R, Werner S, Marillonnet S. A modular cloning system for standardized assembly of multigene constructs. *PLoS ONE*. 2011;6(2):e16765. doi:[10.1371/journal.pone.0016765](https://doi.org/10.1371/journal.pone.0016765)

20. Sarrion-Perdigones A, Vazquez-Vilar M, Palací J, Castelijns B, Forment J, Ziarsolo P, Blanca J, Granell A, Orzaez D. GoldenBraid 2.0: a comprehensive DNA assembly framework for plant synthetic biology. *Plant Physiology*. 2013;162(3):1618–31. doi:[10.1104/pp.113.217661](https://doi.org/10.1104/pp.113.217661)

21. Gibson DG, Young L, Chuang R-Y, Venter JC, Hutchison CA 3rd, Smith HO. Enzymatic assembly of DNA molecules up to several hundred kilobases. *Nature Methods*. 2009;6(5):343–5. doi:[10.1038/nmeth.1318](https://doi.org/10.1038/nmeth.1318)

22. Eddy SR. Accelerated profile HMM searches. *PLoS Computational Biology*. 2011;7(10):e1002195. doi:[10.1371/journal.pcbi.1002195](https://doi.org/10.1371/journal.pcbi.1002195)

23. Textualize. Textual: a Rapid Application Development framework for Python. https://textual.textualize.io (accessed 2026-11-11)

24. Cocioba S. SpliceCraft: a terminal-native plasmid workbench for molecular cloning. Zenodo. doi:[10.5281/zenodo.21960400](https://doi.org/10.5281/zenodo.21960400)

---

## Figures and tables

Tables 1–4 are in-line above. Figure placeholders with draft captions follow; all are producible from the existing codebase.

**Figure 1. The workbench.** Full-window screenshot: braille circular map with lane-stacked features and restriction overlay (left), per-base sequence panel with live translation and feature stripes (right), library panel and menu bar. Suggested subject: pUC19 (GenBank L09137) with the restriction overlay set to unique cutters, so the superscript cut counts are visible. *Produce with:* `splicecraft L09137`, then the built-in SVG screenshot export; convert to PNG for submission.

**Figure 2. Architecture.** Layer diagram, L0 (pure biology, utilities, state, logging, splice/initiation models) through L7 (agent API), with the application hub above and arrows showing that imports only ever point downward. Annotate the three guard tests that enforce it (import layers, public-surface baseline, state single-source). Purpose: show that "one big Python program" is a deliberate, enforced structure rather than an absence of one.

**Figure 3. The write path.** Flow diagram of a single save: caller → authorisation chokepoint (refuse if unsanctioned) → backup rotation → atomic temp write → fsync → read-back validation → shrink guard (spill / refuse) → `os.replace` → parent fsync, with the recovery path (corrupt main → `.bak` → generational rotation → snapshot) drawn as a return arrow. Annotate the 2026-05-22 incident at the point the chokepoint now sits.

**Figure 4. Calibration.** Two panels. **(a)** Donor and acceptor log-odds score distributions for real versus decoy sites, with the 1%-FPR threshold drawn and the resulting sensitivities (0.459 / 0.333) labelled on the curves — the figure that makes "a clean scan is not a clean sequence" visually obvious. **(b)** ROC curves for the trained initiation-context matrix (AUC 0.766) against the −3/+4 consensus rule (0.712), with an inset showing within-rule-bucket discrimination (0.62 / 0.68 / 0.66) — the capability the rule structurally lacks.

**Figure 5. Simulate the steps.** A worked GoldenBraid assembly: four L0 parts → the digest fragments with overhangs labelled by the bases that anneal → the assembled transcriptional unit → the generated numbered protocol and lineage tree beside it. Purpose: show that provenance is a by-product of honest simulation, and double as the educational argument of §8.2.

**Figure 6. Deployment.** Three panels showing the same session: a native terminal, a browser at splicecraft.bio with no local installation, and an SSH session from a low-power client to a Raspberry Pi — with the measured footprint figures from Table 2 annotated. Purpose: make the portability claim concrete rather than asserted.

**Figure 7. The abandonware baseline.** Single panel, tool availability against years since publication, plotting the measured series from refs [6] and [7] on one axis (Mangul: 41.9% of pre-2012 software unavailable vs 17.5% from 2012 onward; Kern: ~90% availability for 2019–2020 publications falling to ~50% for 2010). Annotate the Schultheiss survey findings as text on the figure — 78% built solely by people without permanent positions, 24% with no maintenance plan — because the decay curve and its cause belong in one image. Purpose: make §1.4–§1.5 a figure rather than a paragraph of numbers; it is the most quotable thing in the paper.

**Figure 8 (optional). Headless natural-language operation.** A transcript-style figure: a plain-language request on the left, the sequence of endpoint calls the model made in the middle, and the resulting construct and its generated protocol on the right. Purpose: show §7.2 rather than assert it, and make concrete that the natural-language path and the interface path are the same code.

**Table 5 (recommended for a journal version).** The §9.3 failure taxonomy as a three-column table — failure mode · instance observed in this codebase · countermeasure adopted. §9.3 reads well as prose in a preprint, but a table is what people will want to cite and reuse, and it is the section most likely to be read on its own.

---

## Author's drafting notes ⟨delete before submission⟩

### Style pass (2026-11-11)

Measured before and after, on prose with table rows excluded. Em-dashes 169 → 80 (9.0 → 4.3 per 1,000 words), semicolons 135 → 86, "rather than" 62 → 39, the reassurance adverbs *deliberately / explicitly / genuinely / precisely* 29 → 18. Mean sentence length 24.1 → 22.2 words with the share under twelve words rising from ~19% to ~26%, which is the figure that matters: the tell is uniform rhythm, not any single word.

Four repeated signature phrases were removed. "A format a person can read is a format a person can repair" appeared **four times** verbatim; "a handover document addressed to", "a stranger can change the code" and "earns its keep" each appeared two or three times. Reusing a formulation instead of rethinking it is one of the clearer marks of generated prose.

Also cut: stagey deictic openers ("Here is the perverse part", "Here is the full chain"), awkward imperatives ("Read that as a structural fact"), the meta-comment "that is said without irony", and two overwrought figures ("the bill for that arrives as somebody's experiment", "§9.3(g) in miniature"). Sentence-initial *And / So / But* was trimmed from 31 to 22 — an overcorrection from the em-dash diet, which had traded one tic for a breezier one.

**Abstract cut to 375 words** from 726 (and 760 in the first draft). What went, and what it costs: the opening framing is one sentence rather than two; the feature list names representative assemblies instead of all of them; "39 modules" and the dependency count moved to §3; "verifies constructs against Sanger and nanopore reads" became "against sequencing reads"; and the abandonment statistics came out entirely. That last one was a deliberate choice rather than a space saving: all five figures now live only in §1.4–§1.5, and the abstract names abandonment as a topic the paper treats without arguing it there. The trade is that a reader skimming the abstract no longer meets the numbers, which are the most quotable thing in the paper; the gain is an abstract that makes one argument at a time. Everything else load-bearing survives: the read-only-on-expiry finding, all four commitments, the method disclosure and the failure-mode count. Getting below ~300 would now mean cutting the feature list or a commitment.

**Voice unified to first person singular.** The draft mixed four registers: "we" (27), "the author" (32), "he/his" (4) and "I" (3). Everything authorial is now "I", which removes constructions like "the author has made it his only design tool". Two exceptions are deliberate: the Competing Interests declaration keeps the conventional third person, and one inclusive "we" survives in §10.2 ("the field we have now"), where it means the reader and me together.

### Adversarial fact-check log (2026-11-11)

Every claim in the manuscript was re-checked against code, shipped metadata, tests, git history or the source article — on the assumption that it was wrong. Findings, including the embarrassing ones:

**Corrected — errors introduced in drafting:**

1. **"the pure-Python GenBank reader is ~1.4× faster than the C-assisted one."** FALSE. Biopython's GenBank parser is *itself* pure Python (`Bio/GenBank/Scanner.py`); its 13 compiled extensions serve alignment, clustering, PWM, Nexus, 2bit and structure code, not flat-file reading. The speed comparison is Python-vs-Python. §2.4 now says so and explicitly notes it is the weaker claim.
2. **"eight runtime dependencies (all distributed as prebuilt wheels)"** in the abstract contradicted §11's note that `primer3-py` has no ARM wheel and compiles at install. Abstract now scopes the wheel claim to x86-64 Linux / Intel macOS / Windows.
3. **"the whole program's resident footprint is tens of megabytes"** was asserted twice and never measured. Measured: **85 MiB** with the interface mounted and an empty canvas, **87 MiB** with a 50 kb / 100-feature plasmid loaded — essentially flat in plasmid size, because the figure is interpreter and imports. Both sites now give the number; Table 2 gained two rows.
4. **"a 2023 analysis found a median of 17.9 citations against 5.2"** — the article does not state that. It states a threefold average difference, plus the better finding that predictive performance did *not* affect citation rates. Traced to a search summary, not the paper. Corrected in §1.5.
5. **§1.3 cited two clauses no source supports.** "Consumables and maintenance costs that can exceed national budgets" appears only in ref [11], where that sentence is attributed to a development-news website; replaced with that paper's checkable figure (US$100,000–700,000 NGS instrument capital). "Data generated locally is analysed elsewhere" is stated by **none** of the four reviews examined — the nearest support is *physical samples* shipped abroad for sequencing, which is what §1.3 now says.
6. **§7.2 overstated its source.** The review names English dominance in "training, online resources, and professional networks"; the draft said "tooling, documentation and training". Two of four were unsupported. Reworded to the source's own scope, with SpliceCraft's English-only interface now conceded explicitly rather than implied away.
7. **"tested on 3.10 through 3.14."** CI runs **3.12 only** (`test.yml` says so, and records that a 3.10–3.13 matrix was narrowed for turnaround). Reworded throughout, and added to §11 as a genuine limitation with the precedent that already bit this project.
8. **"The diagnostic log never records sequence content."** Overstated, and inherited from the project's own architecture notes. The cited function is a bounded-`repr` guard: a long string is truncated to a ~50-character prefix plus a length. That is a leak *bound*, not an absence — and the guard's own docstring says so accurately while the prose two files away did not. §2.1 now states the mechanism; §9.3(g) uses it as a second worked example, since it is that failure mode occurring inside a paper about that failure mode.

**Verified against shipped metadata or runtime state — no change needed:** every §5 figure. Splice `fpr_target` 0.01, donor sensitivity 0.4593 / acceptor 0.3325, donor threshold 3.766 (plus AUC 0.962 / 0.9393, n=2,776 real and 50,000 decoy test sites). Kozak AUC 0.766 vs rule 0.712, ΔAUC 0.054, CI [0.035, 0.074], within-bucket 0.624 / 0.683 / 0.658, window −12..+4, `uatg_coverage` 0.958 (= the 4.2%), the three ATG columns reading 0.0 bits, `real_start_consensus` `AAAAAAGCAACAATGG`, A+T 0.518 vs 0.593. Regulettuce: GT1-motif `null_rate` **0.948** exactly, 34 elements, 3,089 promoters / 12,000 null windows / 600 poly(A) sites, TATAWA 1 per 429 bp with a −40..−31 peak at 3.17×, poly(A) 27.2% vs 29.1%, promoter composite AUC 0.646 at 5.0% sensitivity, Y-patch 40.6% vs 38.6% = 1.05×. Endpoint count **303** confirmed from the live registry, not a grep. No master-delete-shaped endpoint is registered. No telemetry or analytics dependency exists. Enzyme aliases Eco31I→BsaI, LguI→SapI, AarI→PaqCI all resolve (39 aliases). Three Cas variants. Agarose 0.5–4.0%. `--read-only` write refusal and its naming of the holding PID are covered by dedicated tests. The L2 write chokepoint fired on the author's own measurement probe during this audit, which is §4.3 demonstrating itself.

**Second wave — two independent adversarial audits of the prose-sourced claims.** Every number in §3.3, §3.4, §9.3 and §9.4 that came from project documentation rather than from code was re-checked against code, tests, stored baselines or git history, and several were re-measured from scratch. Further corrections:

9. **"32 MB" read-back threshold → 8 MB.** The constant is `8 * 1024 * 1024`, and has been since June 2026 — the comment beside it explains why it was lowered. `docs/data-safety.md` still says 32 MB and is where the paper got it. **Repo fix needed.**
10. **"Twenty operations are enforced as regression gates" → nineteen.** `perf-baseline.json` holds 19, and its own notes do the arithmetic (6 + 13). Worse, three Table 2 rows were implied to be gated and are not; they are now marked † and ‡ and the text says so.
11. **CRISPR "40/40 linear, 30/30 circular" → 24/24 and 18/18.** The parametrisation is 4 sizes × 6 and 3 sizes × 6. The brute-force oracle is real and genuinely independent — only the counts were wrong, and they appear to have been borrowed from an unrelated subsystem's inversion tests. **Repo fix needed** (`docs/invariants.md` [INV-198]).
12. **"12,000 damaged files" and "1,280 GFF3 cases" → 400 and 28.** Both larger figures are ad-hoc audit-run numbers in project notes with no committed harness. The 1,280 never multiplied out either — the paper's own stated factorisation (8 × 2 × 2) gives 32, which a reader with a calculator would have caught. The committed figures and the co-oracle design are now what is cited.
13. **"thousands of random IUPAC pairs" → 400** (plus 200 against the base case). The repo's own notes are internally inconsistent here, asserting "thousands" in the body and "over 400" in the same entry's test line.
14. **GenBank differential corpus and throughput, all re-measured.** 1,303 → **1,625** library records and 196 → **214** local files were simply stale (the library grew). "44 of 50 Biopython fixtures" → **40 of 50** on re-run. "4,500 synthetic round-trips" → the committed test does **180**. The throughput claim "37–46 vs 27–31 MB/s" was comparing *different things* — a reader-level figure against one that charged Biopython for file I/O and a since-removed double parse. Measured identically, in one process: **1.11×** (see item 23 — the first corrected pair, 37.8 vs 32.9, was itself basis-mismatched).
15. **Aligner "0.7 s vs 8.6 s (~12×)" → 0.48 s vs 5.05 s (10.6×)**, re-measured on identical input with the method stated, because the original pair existed only as a source comment. The "1.95 s → 0.70 s" trim-lever figure appears nowhere in the repository at all and has been replaced by a statement about what the trim's tests actually gate.
16. **Chromosome figures.** "843 ms load / 8 ms rotation at 200,000 features" are prose-only and did **not** reproduce: measured 8.3 s to load and 86 ms for a full harness repaint. But the claim the auditor flagged as unsupported — **0.001 ms per frame at 1,000,000 features** — *does* reproduce (0.0012 ms), and now has a measurement behind it for the first time. The gated property is the flatness, which held across 10,000 → 1,000,000.
17. **"100+ passing tests" on the orphaned module → 38.** The committed file holds 38 test functions and no parametrisation. **Repo fixed.** One footnote on this item, because it is the sharpest thing the audit produced: an auditor reported the error was also live in the published changelog, and an earlier version of this log repeated that. It is not — the changelog's only "100+" refers to megabytes. A confident, specific, wrong attribution, accepted because it arrived already evidenced, and caught only on the way to acting on it. That is §9.3(f), occurring inside the fact-check of §9.3(f), with the author as the one who propagated it.
18. **"All four worst findings were inside the prior fixes" → three of four.** The fourth, a delete-side guard gap, predates those fixes by five months; the commit that made them never touched it. **Repo fix needed** (`docs/invariants.md` [INV-234]).
19. **"12 of 26 items not fixed" → 7 complete, 7 partial, 12 untouched** (19 of 26 incomplete) — the draft *understated* the failure — and the single-cause attribution was wrong: four causes are documented, including a behavioural test that passed vacuously.
20. **"pyright caught 5 defects that 8,400 tests missed" → unverifiable and removed.** Both numbers exist in one prose note, the release bundles feature and fix in one commit so the five cannot be separated, and the suite was probably nearer 9,000 at that commit. Replaced with the auditable figure: 24 latent optional-type errors accumulated behind a gate whose file list covered 4 of ~27 modules.
21. **Hypothesis scope narrowed.** It is one file, 21 tests, covering reverse-complement and circular-coordinate invariants — not the codon-optimiser properties the draft attributed to it, which are conventional tests. Now stated as the weaker guarantee it is.
22. **Mechanism wording corrected** in two places: the stubbed rig built a *widget* via `__new__` with a stand-in application (not the application itself), and the debounce race was a teardown *re-arming* a timer it had just cancelled (not one surviving a missed cancel). Both are tighter and more interesting than the draft.

**Verified in the second wave, unchanged:** 313/313 NCBI records identical (re-run); the naive-DP oracle's existence and its optimality — not merely agreement — assertion; the CRISPR brute-force oracle's independence and location in the test file; the GFF3 GenBank co-oracle; crash-injection by `SIGKILL` mid-write with a never-torn assertion; eight-thread concurrency fuzzing; all five data-safety mechanics in §4.1 quoted from code, including fsync failures being re-raised rather than swallowed, generation-ordered backups with the clock explicitly not used for ordering, the 50%/90% spill-and-refuse thresholds, and SHA-256 blob verification refusing a mismatch on read; the orphaned module's genuine unreachability and the date it was wired up; the four already-existing features under other names; and the splice-docstring episode, which is the best-evidenced item in §9 — the old text, the change date, and a model file three months older carrying `fpr_target: 0.01` all check out, and 0.95/0.4593 = 2.07 with 0.95/0.3325 = 2.86 lands exactly on the stated "two- to threefold."

**Two claims refined rather than corrected:** the "better in 6 of 7 species" Kozak result is real but the *rule's* per-species AUCs are computed by the calibration script and not persisted into the shipped metadata, so it is reproducible-with-the-script rather than runtime-retrievable; and it is actually "better in six, tied in the seventh," which is stronger than the draft said. The "92% of GenBank single-gene TSS annotations are start codons" figure lives in the module docstring and calibration script rather than the metadata dict, though the *filter* it justifies is in the metadata.

**Verified in this draft.** Every number in §3, §4, §5 and Table 2 was taken from the v1.3.19 source, its calibration metadata, its test suite, or a measurement made during drafting, rather than from documentation prose. Specifically measured during drafting: install footprint (202 MB, clean venv from PyPI), peak RSS (73 MiB), import time (619 ms), 1 kb full-catalogue scan (2.1 ms), collected test count (**12,277**), endpoint count (**303** unique), enzyme catalogue (**202**), presets (**117** in **20** categories), cis-elements (**34**), sibling modules (**38** plus hub), implementation lines (~222,000), test lines (~167,500).

Release cadence in §10.3 measured from git tags: **344 tags, first `v0.2.1` on 2026-04-12, 39 tags in the trailing 90 days**, `v1.3.19` on 2026-11-10.

Quantitative claims quoted from the literature in §1.4–§1.5 were read from the source articles' own abstracts and results sections, not from secondary summaries. One figure was caught and corrected during drafting precisely because of that rule: a search summary attributed "median 17.9 vs 5.2 citations" to reference [9], which the article does not state — what it states is a threefold average difference, plus the more interesting finding that differences in *predictive performance* did not affect citation rates. The paper now says the latter. Mentioned because it is the §9.3(f) failure mode occurring inside the writing of §9.3(f).

**Discrepancy found while drafting:** `README.md` advertised "4,000+ tests" against a suite collecting **12,277**. Both this and the 32 MB → 8 MB figure in `docs/data-safety.md` have since been fixed in the repository; the two "repo fix needed" flags above are discharged.

### Third-pass adversarial audit (2026-11-11, after the log above)

A fresh pass re-derived every checkable claim from source, on the explicit instruction to treat the log above as a hypothesis rather than as settled. It confirmed roughly 140 figures exactly — every §5 calibration number, every §4 constant, all 19 perf gates and their six Table 2 ceilings, the CRISPR 24/24 and 18/18, the GFF3 28 cases, the fuzzer's 400 and 6, the aligner oracle's 400 + 200, the crash-injection and concurrency parameters, 12,277 tests, 303 endpoints from the live registry, 202 enzymes, 344 tags, the cassette module's 38 tests and zero parametrisation at the commit named, the pyright gate's 4-of-27 coverage, and §9.3(b)'s 7/7/12 with all four causes. It also found the following, each now fixed in the body:

23. **§3.3's three differential identity numbers were acceptance counts.** Re-running the project's own `format_scoreboard.py` against v1.3.19: the library is **1,272 of 1,625** identical, not 1,625 of 1,625; local files are **171 of 212**, not 213 of 215. The corpus sizes matched exactly, so the run had happened — the *acceptance* row was read as the *agreement* row, twice. Every one of the 353 library differences turned out to be a single intentional decode (`SpliceCraft_strand` → strand 0), which is a better paragraph than the one it replaced.
24. **§3.3 inverted the adversarial set.** "8 are files Biopython refuses or returns empty" — Biopython refuses **2** of the 50; SpliceCraft refuses **7**. Six of the eight uncomparable fixtures are CONTIG/scaffold stubs that SpliceCraft deliberately declines and Biopython returns with a fabricated sequence. The arithmetic 40+2+8=50 is what let the false attribution survive two previous checks.
25. **§3.3 claimed a test that cannot exist:** "verified identical to the C library it replaced." That library is removed and two assertions forbid its return; the committed oracle is Biopython.
26. **§2.1 corrected twice, the first correction wrong.** The audit found a bounded-`repr` helper with no production caller and concluded there was no leak bound; §2.1 was rewritten to say so. Checking that against `docs/architecture.md` before committing revealed a **second, near-identically-named function** in the logging module that *is* wired into the structured event path and does bound every field — capping strings at 200 characters, reducing byte blobs to a size tag, and rendering a Biopython `Seq`/`SeqRecord` as an opaque tag because its `repr` embeds bases. The truth is a split guarantee: enforced on structured events, convention-only on plain log calls. §2.1 now says that; the orphaned hub duplicate is a real §9.3(a) instance in its own right and is noted there. This is the single most-corrected sentence in the paper and the reason §9.3(g) now runs to four passes.
27. **§3.4's "300 examples per property"** is 100 for six properties, 200 for seven and 300 for eight; the same sentence also understated the file's scope.
28. **§11 contradicted §2.7** on the primer3 fallback. A Tm fallback does exist; what is absent is structure checking. The false version originates in a comment in `scripts/check_dep_wheels.py`.
29. **Literature.** Ref [8]'s 78% is the abstract's figure and its own Results and Figure 7 say 58%; the tested reachability-citation association is Kern's, not Schultheiss's; ref [6]'s 66% does not reconcile with its own 28%; the package-manager finding rests on about nine tools; ref [10]'s 28 + 9 of 36 sums to 37; ref [11]'s capital figure originates in a 2011 survey; ApE offers optional NCBI and LabArchives network features, so "no network" was wrong; and ref [13] names no individual `/ApEinfo_` qualifier, so the colour-qualifier spellings are observed interoperability, not a sourced claim.
30. **Item 4 of the log above removed a true figure.** Song & Kurgan *does* state the medians (17.9 working server, 5.9 source-only, 5.2 no longer available). They are restored in §1.5.
31. **Table 1's terms were harder than the vendors'.** Both say "most other functions are disabled", and both leave annotate-and-share or import-and-export; the course licences are free to students but require the instructor to hold a paid seat. Both softened.
32. **Two measurements were slightly optimistic.** Throughput re-measures at 1.11×, not 1.15×. The aligner's 10.6× reproduces exactly (10.7×) but only under the application's affine gap costs, now stated in the Table 2 footnote because unit-cost gaps halve it.

One methodological note worth keeping, because it is §9.3(f) with me as the author again: the audit's first two "findings" were both its own bugs — a re-seeded RNG inside a generator expression that produced a homopolymer with no restriction sites at all, and unit-cost gap scoring in the aligner comparison. Both would have been publishable-looking corrections to correct numbers. The rule that caught them was checking that the benchmark found anything before believing what it timed.

**Repo fixes made in the same batch as this manuscript:** four stale `splicecraft._DATA_DIR` references that have raised `AttributeError` since the attribute moved to `_state` — in the write chokepoint's own remediation message, in the sandbox docstring that message points at, and in two `scripts/` sandbox guards, the last of which means the guard crashes instead of guarding; the "no fallback" comment in `scripts/check_dep_wheels.py` that §11 inherited; and the logging paragraph in `docs/architecture.md`, which now distinguishes the enforced path from the convention and warns against the orphaned near-duplicate that misled this audit.

**One finding from this pass was withdrawn.** The audit reported that `--citation` emits the concept DOI labelled as a version, so following the tool's output would not pin a release. It does not: the command already prints a closing note saying the concept DOI resolves to the newest release and telling the reader to pick the exact version from the record's version list. The audit had read the first twelve lines of a fifteen-line output. Together with the two measurement bugs recorded above — the homopolymer and the gap costs — that makes three occasions in one audit where the checker was wrong and the thing being checked was right, each caught only by looking a second time at something that already looked settled. §9.3(f) does not stop applying to the person doing the auditing.

**Decisions to make before submission.**

1. **Title.** The draft title is descriptive and keyword-rich. Alternatives, in decreasing formality: *"SpliceCraft: a terminal-native plasmid workbench for molecular cloning"* (matches `CITATION.cff` and the Zenodo record exactly, which keeps citations consistent); *"SpliceCraft: a free, offline plasmid workbench for laboratories without a software budget"*; *"SpliceCraft: the whole cloning workflow, offline, on the computer you already have."*
2. **Naming commercial platforms — done, and check Table 1 before submission.** §1.1 now names SnapGene, Geneious Prime and Benchling with list prices read off their own pricing pages on 2026-11-11, plus verbatim quotes of the expiry, perpetual-licence and academic-eligibility terms. Two deliberate choices there: (a) the argument leads with **two concessions** (Benchling's academic plan is free; SnapGene and Geneious both give free course licences) before making its point, because a reviewer would otherwise raise both in the first paragraph of their report and the paper would look like it was arguing in bad faith; (b) the load-bearing claim is now about **terms rather than prices** — the read-only-on-expiry clause and the no-upgrades perpetual licence — which is both more interesting and immune to a price change. **DNASTAR was deliberately omitted**: its academic pricing page returned HTTP 403 to an automated request, and an estimated number would have been worse than an absent one. Add it by hand from a browser if you want a fourth row. Re-verify the whole table shortly before you submit; these pages change without notice, and a stale price is the single easiest thing for a hostile reader to discredit.
3. **bioRxiv category.** Suggest *Synthetic Biology* primary, *Bioinformatics* secondary.
4. **Consider splitting §9 out as its own paper.** §9 (the LLM-assisted development failure taxonomy) is arguably the most novel and most broadly useful content here, and it is currently buried at position nine in a software paper, where the audience that most needs it — people building scientific software this way right now — will never find it. It stands alone as a short methods/perspective piece, with SpliceCraft as the case study and §1.4–§1.5 as the motivation. Keeping both is also defensible: it is the evidence that the software's quality claims are not self-assessed. My recommendation: **keep it here for the preprint**, since it is what makes the rest credible, and extract an expanded version later if it draws attention.
5. **Length.** Now ~10,000 words of body text. Comfortable for bioRxiv, long for a journal: *Bioinformatics* would need ~2,000 words plus supplement, *BMC Bioinformatics* could take most of it, *PLoS Computational Biology* (Software/Education section) is probably the best fit for the paper as written, and *Journal of Open Source Software* wants ~1,000 words — a very different, much shorter paper, but peer-reviewed and fast, and it would give the software a reviewed citation to sit beside this preprint.
6. **Ko-fi URL verified 2026-11-11.** `ko-fi.com/atinygreencell` resolves, is live, and carries 28 supporters. It could not be checked automatically (Ko-fi returns HTTP 403 to automated clients and the handle is not search-indexed), so it was confirmed visually in a browser. The Funding statement stands as written. The page also describes you as "Amateur plant biologist … Creator and maintainer of SpliceCraft", which is consistent with the framing in §1 and the Acknowledgements.

7. **Name consistency — now one-sided.** The paper is bylined *Sebastian S. Cocioba*, and the Ko-fi profile uses that same form, so the middle initial is your public-facing name. `CITATION.cff` and `.zenodo.json` are the odd ones out at *Cocioba, Sebastian*. Aligning them (`given-names: Sebastian S.`) would make the preprint, the software citation and the Zenodo record resolve to one author string in indexes. It edits released metadata, so it remains your call, but the evidence now points one way.

8. **The Acknowledgements closing paragraph is a joke, deliberately.** It is deadpan rather than bitter and it lands the paper's thesis one last time, but it is the one paragraph a conservative reviewer might ask you to cut. It is also, in my view, the most memorable thing in the manuscript and the reason people will send it to each other. Keep it; just know it is a choice.
7. **Figures 4 and 7 need plotting.** Figure 4's panels are derivable from the existing calibration scripts, though the splice panel needs the score arrays retained rather than summarised — likely a small flag on the calibrator. Figure 7 is two published series plus three annotations and can be drawn from the numbers already in §1.4.
8. **Decide how hard to lean on §10.3.** "I cancelled my paid subscriptions and this is now my only design tool" is the single most credible sustainability statement in the paper, and it is unusual to see one stated so plainly. It also commits you in public. It is written as a structural argument rather than a promise, which is the right framing, but read it once more with that in mind.
