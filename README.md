# Takivisor

Quick QC reviewer for **bacterial** isolate sequencing: give it an assembly (and, ideally, the reads and the QC
reports your pipeline already made) and it tells you whether the result looks trustworthy, what is wrong if not,
and what to do first.

- **Research use only, not for clinical decisions.**
- Bacteria only (species table: 4,741 NCBI bacterial species).
- Python 3.8+, standard library only. Takivisor never runs other tools; it reads their reports and does its own
  checks on the FASTA/FASTQ files. When a report would help and is missing, it says which tool to run.
- Works with any pipeline's output (Illumina, ONT, PacBio HiFi/CLR; SPAdes, SKESA, Unicycler, Flye, Raven, Canu,
  Autocycler, Hybracter, Dragonflye, myloasm...).

## Install

```
pip install /path/to/takivisor        # gives the `takivisor` command
```
or run it in place: `python3 takivisor.py ...` (keep `species_genome_size.tsv` and `demo/` next to it).

## Quick start

```
takivisor --demo                       # makes ./takivisor_demo and reviews it (real genome, synthetic reads)
takivisor RUN_FOLDER                   # a pipeline output folder: assembly, reads and reports are found in it
takivisor --assembly asm.fasta --reads R1.fastq.gz R2.fastq.gz --species "Klebsiella pneumoniae" \
          --reports qc_folder/
```

Always give `--species` (or an identity report: Kraken2, sylph-tax or GTDB-Tk); without it, genome size and GC
cannot be checked against the species. Typos get a "did you mean" suggestion.

**Folder mode** picks the final/polished assembly over intermediate ones, filtered/trimmed reads over raw ones,
skips assembler work folders, and prints what it used. It refuses folders with more than 8 assembly candidates
(probably several runs; point it at one run). Any option overrides what it found.

## What you get

| Output | Where |
|---|---|
| Text report (Markdown) | standard output |
| HTML report | `takivisor_report.html` in the run folder, or `<assembly>.takivisor.html`; `--html PATH` / `--no-html` |
| JSON | `--json PATH` (see *Output for pipelines*) |

The verdict is **PASS**, **REVIEW** or **DETECTED ISSUE** (stored as `FAIL` in the JSON and exit code 4). The HTML page is written for people new to QC: a one-line
verdict, six status tiles (size, contigs, depth, completeness, contamination, accuracy), a **Start here** list
in the order to fix things (inputs, then data, then contamination, then assembly), split into **In the lab** and
**In bioinformatics** tabs that give each step's action for that side, a card per finding with plain
words and *What to do*, and charts (cumulative length vs the species' range, contig coverage vs length, read
lengths, read GC vs the assembly). Technical details, all metrics and a glossary are collapsed below.

## What the data is good for

Every report rates six uses of the data: species identification, MLST/serotype, AMR/virulence gene presence,
claims that a gene is absent, relatedness (cgMLST, SNP distances, outbreaks) and genetic context (plasmid vs
chromosome, complete genome), as **OK**, **use with caution** or **not recommended**, or **not checked** when the check that would rate it is
unreliable on this run (e.g. ONT lambda DCS reads that Kraken2/Bracken count as *E. coli*, which inflate the second
species; the read-GC chart marks the DCS spike). Ratings follow the EU
reference laboratories' dry-lab QC guideline (2025); where a peer-reviewed study found data below it still
usable, the rating is "use with caution" and the report names both (e.g. Illumina at 25x: below the EU 30x, but
assembly-based cgMLST was ~100% correct from 20x in Palma et al. 2022). Long reads follow nanopore studies
(e.g. cgMLST accurate from 35x with Dorado sup; Genome Research 2025); PacBio HiFi depth is not rated (no
peer-reviewed bacterial benchmark). Nanopore data basecalled with a high-accuracy or
fast model is "use with caution" for typing and relatedness, and the note adds that Pasteur BIGSdb only accepts
R10.4.1 + Dorado sup at >= 50x. A ConFindr report rates same-species mixing; incomplete AMRFinderPlus/RGI calls
make gene presence and gene-absence claims "use with caution". Checks that could not run (no CheckM2, no identity
report, no ConFindr, no mosdepth, unknown basecalling model) are listed as "not checked". A second species in an identity report from 1% up to the EU limit
makes gene, typing and context uses "use with caution" without changing the verdict (AllTheBacteria requires
>= 99% one species; Raven 2020). Each finding also gives a lab action and a
bioinformatics action. Rules and sources: `fitness()` and `ACTIONS` in `takivisor.py`.

## What it checks

47 checks, each with a stable id (shown in every report). All cut-offs are in the `T` dict at the top of
`takivisor.py`, each with its source and validation status.

- **Inputs**: reports belong to this assembly (length/name/contig-name match), species known, files are what they
  claim, contig names unique (repeated names usually mean two assemblies joined into one file).
- **Reads**: sequencing technology (detected from read names/qualities), depth (WARN below 30x, FAIL below 15x;
  100x for Trycycler), Illumina Q30 by read length (80% up to 2x150, 75% up to 2x250, 70% up to 2x300), ONT read
  length, SRA-lite placeholder qualities, ONT lambda DCS control reads, read GC vs assembly GC and per-read GC
  outliers. Where another body sets a different depth, the report notes it without changing the verdict: looser
  minimums are added to the depth WARN (*M. tuberculosis* 15x, WHO TB catalogue 2023; *Listeria* and
  *Campylobacter* 20x, PulseNet 2026), stricter ones are an INFO note (`depth_below_agency`: PulseNet 40x for
  Illumina *E. coli*/*Shigella*, *Vibrio*, *Cronobacter*, *Yersinia*, *C. botulinum*; Pasteur BIGSdb 50x for
  nanopore-only data).
- **Depth breadth** (mosdepth report): WARN `low_breadth` when < 95% of the genome has >= 10x (EU guideline 2025);
  for *M. tuberculosis* an INFO note below 95% at >= 20x (WHO TB catalogue 2023). The report is used when its
  length matches the assembly, or when it is the only mosdepth report and was mapped to a reference (no shared
  contig names, e.g. H37Rv).
- **Species and contamination**: genome size and GC vs the species' NCBI/RefSeq range, and size vs 99.8% of the
  species' GenBank drafts; identity reports (Kraken2/Bracken, sylph-tax, GTDB-Tk; Bracken is used over the raw Kraken2 report of the same sample, as in the EU guideline), where a second species at >= 5% is a
  WARN (>= 2% for *Listeria monocytogenes* and coagulase-positive staphylococci, EU per-pathogen limit) and >= 20% a FAIL (close relatives in the same genus, and Escherichia/Shigella, stay INFO); BUSCO/CheckM
  lineage placement; CheckM2 contamination (FAIL above 5%); Kraken2 per-contig output; low-coverage contigs; ConFindr
  (same-species mixing: WARN when it reports contamination, noting Pasteur BIGSdb's 5% limit).
- **AMR gene calls**: AMRFinderPlus (NDARO) and RGI (CARD) calls cut by a contig end, partial (< 90% of the
  reference), with an internal stop codon, or RGI Loose hits are listed in a note (`amr_calls_incomplete`) with
  what to check before reporting those genes.
- **Assembly**: contig count and N50 for short reads (species limit = 99.5th percentile of contigs >= 200 bp in
  the species' GenBank drafts for the 61 species with >= 1,000 drafts, otherwise EnteroBase's 600; N50 20 kb),
  replicon-aware limits for long reads, N content, circular chromosome, suspect circles, Canu repeat contigs,
  lambda DCS in the assembly.
- **Accuracy**: indel errors in unpolished long-read assemblies (CheckM2 gene length, BUSCO fragmented, or
  Takivisor's own ORF length check when no report is given).
- **Completeness**: CheckM2 (WARN below 95%, FAIL below 50%; BUSCO is reported but does not judge completeness:
  bacteria_odb10 under-scores divergent genomes), QUAST genome fraction.

Limits follow the EU reference laboratories' 2025 dry-lab QC guideline, ECDC proficiency tests, bactscout and
QualiBact where they agree, and were checked for false alarms on 13,599 NCBI-accepted drafts and the validation
sets.

Reports read: BUSCO `short_summary*.txt`, CheckM2 `quality_report.tsv`, legacy CheckM, QUAST `report.tsv`,
fastp `*.json`, NanoPlot `NanoStats.txt`, Kraken2 report and per-contig output (+ DB folder for
`ktaxonomy.tsv`), sylph-tax `.sylphmpa`, GTDB-Tk summary, Flye `assembly_info.txt`, ConFindr `confindr_report.csv`, AMRFinderPlus `*amrfinder*.tsv`,
RGI `*rgi*.txt`, mosdepth `<prefix>.mosdepth.summary.txt` (+ `<prefix>.mosdepth.global.dist.txt`
next to it). Reports with no genome length (Kraken2, legacy CheckM, ConFindr) are matched by sample name
(a one-sample ConFindr report is used as is; ConFindr names a sample by its file name before `_R1`); use `--sample NAME` when file names are generic. AMRFinderPlus/RGI
reports are matched by contig names.

## Output for pipelines

**Exit codes** (with `--exit-code`; without it Takivisor exits 0 after a successful review):

| Code | Meaning |
|---|---|
| 0 | PASS |
| 3 | REVIEW |
| 4 | FAIL |
| 1 | input error (unreadable, empty or wrong-type file, ambiguous run folder) |
| 2 | argument error |

**JSON, schema 1.** These keys are stable; a breaking change bumps `schema`:

| Key | Type |
|---|---|
| `schema` | integer, currently 1 |
| `version` | Takivisor version string |
| `sample` | sample name used |
| `verdict` | `"PASS"`, `"REVIEW"` or `"FAIL"` (reports show FAIL as "DETECTED ISSUE") |
| `findings` | list of `{"severity": "INFO"/"WARN"/"FAIL", "id": check id, "message": text, "action": text}` |

Check ids are stable; new ids may be added. Informational keys added 2026-10-05: `fitness` (per use: `rating`, `reasons` with `basis`, `not_checked`) and `actions` (per finding id: `lab`, `bioinformatics`). `message` and `action` are human text and may be reworded. Everything
else in the JSON (`assembly`, `reads`, `reports`, `species`, `inputs`, `species_match`) is informational and may
change between versions. Verdict rule: any FAIL → FAIL; any WARN except the data-only ones (`illumina_q30`,
`fastp_low_pass`, `ont_short_reads`) → REVIEW; otherwise PASS.

## Validation

Validation results (the data sets are not included in this repository):

| Set | What | Result |
|---|---|---|
| False-positive sweep | 34 complete RefSeq genomes, GC 25-72%, 0.58-12 Mb | all PASS after fixes |
| Sensitivity | 17 deliberately broken sets (contaminant mixes, low depth, cut reads) | caught; limits documented |
| GenBank held-out | 20,796 genomes with NCBI atypical/ANI labels | 0.02% flagged by size/GC |
| ONT read-GC | 29 public ONT runs + the RefSeq genome built from each | 0 false read-GC alarms |
| Illumina drafts | 24 public paired Illumina runs + the draft built from each | 0 false alarms; full-run depth 45-707x by ENA base counts |

## Limitations

- Bacteria only. Not for metagenomes or mixed samples.
- Low-level contamination by a species of similar GC is only caught with a CheckM2 or identity report; a mix of
  strains of the same species only with a ConFindr report (checked on real ConFindr 0.8.2 output; without its
  restricted rMLST database ConFindr only covers Escherichia, Salmonella and Listeria, other genera show as "not checked").
- FASTQ must be standard 4-line records; wrapped or cut-off files stop with an error rather than giving wrong
  numbers.
- RGI writes `<name>.txt`: name the file with "rgi" (e.g. `sample.rgi.txt`) so it is found; AMRFinderPlus output
  needs "amrfinder" in its name.
- Depth breadth is genome-wide only (mosdepth `global.dist.txt`, 1% resolution); per-gene depth (WHO: a gene is
  missing if any base < 10x) and `samtools coverage` output are not read.
- Read-based checks assume the FASTQ is a random sample of the run. A file exported from a position-sorted
  BAM/CRAM, or only its first reads, can give false read-GC warnings.

## Licence

Not licensed yet: all rights reserved until a licence is chosen. The bundled species table is derived from public
NCBI data.
