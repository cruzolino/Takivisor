#!/usr/bin/env python3
"""Takivisor v1.0 - quick QC reviewer for bacterial reads + assemblies.

Detects sequencing tech (per read file) and assembler (from contig headers),
computes cheap metrics, reads existing QC tool reports, checks everything
against context-aware thresholds and prints a markdown report. Stdlib only.
"""
import argparse, datetime, difflib, functools, gzip, json, math, os, re, sys, time, zlib
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
# Bundled data sits next to this file, or under <prefix>/share/takivisor after `pip install`.
_data = lambda name: next((d / name for d in (HERE, Path(sys.prefix) / "share" / "takivisor") if (d / name).exists()),
                          HERE / name)
SPECIES_TSV = _data("species_genome_size.tsv")
LONG_READ_ASM = {"flye", "canu", "raven", "trycycler", "unicycler", "hifiasm", "autocycler", "hybracter"}
SHORT_READ_ASM = {"spades", "skesa", "shovill"}
ILLUMINA_POLISHERS = {"pilon", "polypolish", "polca"}
SEV_ORDER = {"FAIL": 0, "WARN": 1, "INFO": 2}
LONG_TECHS = {"ont", "hifi", "clr"}  # clr = PacBio continuous long reads (subreads)
# --exit-code: verdict as exit status. 1 (bad input) and 2 (bad arguments) stay reserved for errors.
EXIT_CODES = {"PASS": 0, "REVIEW": 3, "FAIL": 4}
# What reports show for FAIL (verdict and finding severity); JSON and exit codes keep "FAIL" (schema 1).
SHOWN = {"FAIL": "DETECTED ISSUE"}
shown = lambda x: SHOWN.get(x, x)

# Every cut-off in one place. Status: "sweep" = no false alarms on 34 complete RefSeq genomes
# (2026-09-30); "caught" = fired on real/real-derived defect data; "untested" = never fired on bad data.
T = {
    "gc_slack": 1.0,              # NCBI gc_percent rounded to 0.5 | caught (mislabelled E. coli); 0.02% of
                                  # 20,796 held-out GenBank genomes flagged (all 4 look genuinely wrong)
    "gc_min_genomes_fail": 20,
    # Size outside the NCBI species range by at most this fraction of the nearest bound: WARN, beyond: FAIL
    # (#16 decision 2026-10-01: the range is a statistical cut). Design choice, untested on bad data.
    # Real contamination is far bigger: 3% S. aureus reads in Kpne gave +2.5 Mb (6.7 -> 8.1 Mb, 21% over).
    "size_warn_frac": 0.10,    # fewer RefSeq genomes than this: GC outlier is WARN, not FAIL
    "n_pct": 3,                   # EnteroBase | caught (real GenBank scaffolds 5.9%, 7.7% N)
    "long_contigs_min": 10,       # floor for long-read contig count | caught (Flye 12, Canu 14)
    "replicon_slack": 3,          # + species RefSeq replicon p90 | sweep (Borrelia 13)
    "circ_min_len": 1e6,          # "chromosome-sized" circular contig
    "near_genome_frac": 0.8,      # largest contig this close to genome size = probably complete
    "suspect_circ_frac": 0.7,     # circular >=1 Mb but < this x genome size, in a fragmented asm
    "big_linear_len": 100_000,    # linear contig counted as "other large piece"
    "short_contigs_max": 600,     # EnteroBase | caught (3% contam, 1841 contigs); 8x Kpne gave 518
                                  # Fallback only: species with >= 1,000 GenBank drafts use their own limit
                                  # (contigs >= 200 bp, p99.5; build_species_table.py). 2026-10-05 test on
                                  # 13,599 NCBI-accepted Illumina drafts: QualiBact-style species limits added
                                  # 0.9 pts flagged, +6 pts of NCBI-"contaminated" drafts caught (54 -> 60%).
    "short_n50_min": 20_000,      # EnteroBase | caught (real GenBank drafts N50 1.6-13.8 kb); 8x Kpne 25.8 kb passed
    "lowcov_rel": 0.1,            # contig cov < this x median ... | caught at 3% contam, not 15%
    "lowcov_frac": 0.05,          # ... covering more than this fraction of the assembly
    # % reads outside the assembly's 99% window-GC range. Set from clean sets only (2026-09-30):
    # Illumina DRR125067 2.23%; ONT CFT073 55x/100x 0.28/0.30%, DRR204858 0.73%. ~2-3x the clean max.
    "read_gc_outlier_pct": {"illumina": 5.0, "long": 2.0},
    # Held-out: clean ONT S. aureus 1.44, E. faecium 0.87, A. baumannii 0.36 (no alarm); P. aeruginosa
    # 19.5% was REAL: 18% of its reads are ONT lambda control (DCS), 2.7 kb of which Raven built into the
    # chromosome. Real P. aeruginosa reads in S. aureus: 1% missed, 3% caught; close-GC E. faecium (38 vs
    # 33) invisible even at 20%. Illumina + simulated S. aureus: caught at 10%, not 5% (same as mean GC).
    "read_gc_outlier_clean": {"illumina": "2%", "long": "0.3-1.4%"},
    "read_gc_min_reads": 2000,
    "dcs_read_step": 50, "dcs_read_min_hits": 3,   # 15-mers every 50 bp ...
    "dcs_len": (3000, 4200), "dcs_span": 0.6,      # ... of ~3.56 kb reads, matching >=60% of their length
    "dcs_asm_step": 5, "dcs_asm_min_hits": 50,     # exact 21-mers every 5 bp; 50 hits ~ 250 bp of DCS
    "dcs_reads_pct": 0.5,         # % of long reads that are DCS: PA run 19.6%; E. coli runs 0%
    "read_gc_diff": 2,            # points | 6 matched sets <=0.54, wrong species 6.5; caught 3%/15% contam
    "depth_min": 30,              # EU EURL 2025, ECDC EQA, bactscout, GenomeTrakr ONT (was 20, TakiLine) |
                                  # caught (Canu 18.5x, Kpne 8x, ONT 10x); full real runs: Illumina 45-707x
                                  # (24/24 pass), ONT 21-1,482x (1/29 WARN: 21x hybrid). FAIL below half.
    # Other bodies' depth minimums: notes only, the EU limit keeps the rating (user policy 2026-10-05).
    # (species or genus or None = any, tech or None = any, min x, who). Looser -> note on low_depth;
    # stricter -> INFO depth_below_agency. PulseNet PNQ09 v10 App. 3 (2026-03-24, Illumina); Pasteur
    # BIGSdb accepted technologies (ONT R10.4.1 + Dorado sup, >= 50x); WHO TB catalogue 2023 (median > 15x).
    "depth_other": (
        ("Mycobacterium tuberculosis", None, 15, "the WHO TB catalogue 2023 (median depth)"),
        ("Listeria monocytogenes", "illumina", 20, "PulseNet 2026"),
        ("Campylobacter", "illumina", 20, "PulseNet 2026"),
        *((g, "illumina", 40, "PulseNet 2026") for g in (
            "Escherichia", "Shigella", "Vibrio cholerae", "Vibrio parahaemolyticus", "Vibrio vulnificus",
            "Cronobacter", "Yersinia", "Clostridium botulinum")),
        (None, "ont", 50, "Pasteur BIGSdb (ONT R10.4.1 + Dorado sup only)"),
    ),
    "depth_min_trycycler": 100,   # TakiLine
    # Depth breadth from a mosdepth report: % of the mapped-to sequence at >= x. EU EURL 2025: >= 95% at 10x
    # (WARN; no FAIL tier until real runs show where one belongs). WHO TB catalogue 2023: 95% at 20x (note only).
    "breadth_min": (10, 95),
    "breadth_other": (("Mycobacterium tuberculosis", 20, 95, "the WHO TB catalogue 2023"),),
    # % bases >= Q30 by read length (EU EURL 2025: >=80% up to 2x150, 75% up to 2x250, 70% up to 2x300;
    # bactscout WARN < 80 / FAIL < 70; Illumina specs) | caught (real MiSeq runs 50.9%, 73.3%); Illumina set
    # 84.5-100% (24/24 pass).
    "q30_min_by_len": ((151, 80), (251, 75), (10**9, 70)),   # max read length incl. the extra cycle
    "fastp_pass_min": 80,         # heuristic | caught (same real runs: 56.7%, 75.5% passed fastp)
    "ont_n50_fail": 1000,         # Flye's 1 kb minimum overlap | caught (real ONT run ERR12031453, N50 518 bp)
    "ont_n50_warn": 5000,         # rRNA operons 5-7 kb | caught (reads cut <1.5 kb)
    "unicycler_long_max_contigs": 50,  # without reads: Unicycler this closed = long-read mode
    "busco_frag_indel": 20,       # % | sweep max 8.9; unpolished 43-53 (first set), held-out 6 species 8-50
    "gene_len_indel": 220,        # aa | sweep min 265.5; unpolished 136-164, held-out 60-278 (5/6 caught)
    "cut_at_ends_n50": 50_000,
    "assembler_min_agreement": 0.5,  # share of contig names that must match one assembler (#2 review)
    # Identity reports: the top species must hold this % of bacterial species-level assignments to be
    # compared with --species (or to fill it in); a second species above id_second_share is noted.
    # Second species: WARN >= 5% (EU EURL 2025 critical < 5%, GHRU ConFindr > 5% fail, KlebNET 5%, Raven 2020
    # 4%), FAIL >= 20% (bactscout). Same genus or Escherichia/Shigella stays INFO: classifier splits between
    # close relatives are not contamination (EURL: C. jejuni/C. coli cross-hits; QualiBact merges Shigella).
    "id_min_share": 50, "id_second_share": 5, "id_second_fail": 20,
    # EU EURL 2025 per-pathogen tables: Kraken2 <= 2% for Listeria monocytogenes and coagulase-positive staphylococci.
    "id_second_share_species": {"Listeria monocytogenes": 2, **{f"Staphylococcus {x}": 2 for x in (
        "aureus", "argenteus", "schweitzeri", "intermedius", "pseudintermedius", "delphini", "coagulans", "lutrae")}},
    # Fitness only (never a finding): a second species from this % is "detected" -> CAUTION. AllTheBacteria HQ
    # needs sylph >= 99% one species; Raven 2020 controls < 0.4% Kraken reads (= < 1% contamination).
    # CheckM2 is not used for this tier (contamination MAE 1.2 +- 1.3%, Chklovski 2023).
    "id_second_caution": 1,
    # Kraken2 (standard DB) assigns ONT lambda DCS control reads to these species. Real runs: abby teste_local
    # (K. pneumoniae ONT, 2026-10-06) 92% of 'E. coli' and 95% of 'A. baumannii' reads were DCS (Takivisor's
    # detector; 96%/100% align to DCS with minimap2); virus_pipe ABL runs 2026-09-29 the same two species.
    "dcs_miscalled_as": ("Escherichia coli", "Acinetobacter baumannii"),
    # Per-contig Kraken2 output (#25): a contig's species = most k-mer votes at species level or below.
    # Real abby contigs (standard-8 DB, all left 'U' by Kraken2): chromosome 97% K. pneumoniae of 21,204
    # votes, plasmids 80% (228) and 72% (61). Cut-offs below are design choices, untested on contamination.
    "contig_min_votes": 20, "contig_min_share": 50,
    "contig_foreign_pct": 1.0,    # % of assembly bp on contigs called as another species -> WARN
    # Foreign contig at >= this x median contig coverage looks like a shared plasmid/mobile element, not
    # contamination (INFO). Real plasmids >= 80 kb: 0.40-2.46x chromosome (Autocycler archive, 7 plasmids
    # in ~165 assemblies; abby 0.42-0.5x); a contaminant at <= 20% abundance sits at <= 0.25x.
    "contig_lowcov_rel": 0.25,
    # All contigs not flagged on their own, votes pooled: a second species above this % of species-level
    # votes -> INFO. Clean abby run: E. coli 0.35%, K. quasipneumoniae 0.5%. Untested on contamination.
    "contig_pool_share": 3.0,
    "id_min_reads": 100,          # Kraken2 only: fewer species-level reads is noise (real sepsis report S1: 4)
    # ORF N50 (aa), code-appropriate. Set on calibration data only (2026-09-30): 34 clean genomes >= 216
    # (M. tuberculosis; Mycoplasma 298+ under code 4), 3 original unpolished 162-181.
    "orf_n50_indel": 200, "orf_n50_clean_min": 216,    # short-read asm below this: BUSCO F is contig-end cutting
    "busco_complete": 95,         # INFO only: 13/34 complete genomes score 57-93%
    "busco_dup": 5,               # % | sweep max 3.2; caught (mixes 15-44%)
    # CheckM2 for isolates (EU EURL 2025 critical: completeness >= 95, contamination <= 5; ECDC EQA 96/4;
    # was MIMAG 90/5-10) | sweep completeness >= 97.5, contamination max 3.81; caught (ONT 10x 80.9%, 11/12
    # mixes). NCBI-data test: +0.2 pts flagged, frameshift-excluded drafts caught 69 -> 91%.
    "ck_compl_warn": 95, "ck_compl_fail": 50,
    "ck_cont_warn": 5, "ck_cont_fail": 5,
    "quast_gf_min": 95,           # % of reference | caught (CFT073 vs K-12 82.9%; same strain 99.9%)
}


def _open(path):
    with open(path, "rb") as f:
        gz = f.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rb") if gz else open(path, "rb")


def _lines(path, handle=None):
    """Line iterator that turns a truncated/corrupt gzip into a clear exit, not a traceback.
    `handle` (a list) receives the open file, so callers can ask how far into it they are."""
    try:
        with _open(path) as fh:
            if handle is not None:
                handle.append(fh)
            yield from fh
    except (EOFError, gzip.BadGzipFile, zlib.error) as e:
        sys.exit(f"{path}: truncated or corrupt gzip ({e}); re-download or re-compress it")


class _Progress:
    """Live progress on stderr for one file: an updating line in a terminal, only start/end lines in logs."""
    def __init__(self, path, what):
        self.path, self.name, self.t0, self.last = path, Path(path).name, time.monotonic(), 0
        self.size, self.fh, self.tty = os.path.getsize(path), [], sys.stderr.isatty()
        print(f"takivisor: {what} {self.name} ({self.size / 1e6:,.0f} MB)...", file=sys.stderr, flush=True)

    def tick(self, n, bases):
        now = time.monotonic()
        if not self.tty or now - self.last < 0.5 or not self.fh:
            return
        self.last = now
        f = self.fh[0]
        frac = (f.fileobj.tell() if hasattr(f, "fileobj") else f.tell()) / max(1, self.size)
        left = (now - self.t0) * (1 - frac) / frac if frac > 0.02 else None
        print(f"\r  {frac:4.0%} · {n:,} reads · {_fmt_len(bases)}"
              + (f" · about {left:.0f} s left" if left is not None else "") + " " * 8,
              end="", file=sys.stderr, flush=True)

    def done(self, n, bases):
        if self.tty:
            print("\r" + " " * 72 + "\r", end="", file=sys.stderr)
        print(f"takivisor: read {self.name}: {n:,} reads, {_fmt_len(bases)} in {time.monotonic() - self.t0:.0f} s",
              file=sys.stderr, flush=True)


def nx(hist, frac=0.5):
    """Nx from a {length: count} histogram."""
    total, run = sum(l * n for l, n in hist.items()), 0
    for l in sorted(hist, reverse=True):
        run += l * hist[l]
        if run >= total * frac:
            return l
    return 0


# ---------------- assembly ----------------

# (assembler, header regex). Order matters: most specific first.
# Polishers often append text ("contig_1 polypolish", "NODE_..._pilon"), so no "$" anchors.
ASM_PATTERNS = [
    ("spades", re.compile(r"^NODE_\d+_length_\d+_cov_[\d.]+")),
    ("skesa", re.compile(r"^Contig_\d+_[\d.]+(_Circ)?(\s|$)")),  # real: "Contig_1_80.4955_Circ [topology=circular]"
    ("unicycler", re.compile(r"^\d+ length=\d+ depth=[\d.]+x")),
    # real: Autocycler 0.5 "1 length=5369556 circular=true" (no depth=, unlike Unicycler);
    # Hybracter 0.11 "chromosome00001 len=5369555 circular=true" / "plasmid00001 len=..."
    ("autocycler", re.compile(r"^\d+ length=\d+(\s+circular=\w+)?\s*$")),
    ("hybracter", re.compile(r"^(chromosome|plasmid)\d{5} len=\d+")),
    ("canu", re.compile(r"^tig\d{8}")),
    ("raven", re.compile(r"^Utg\d+ .*LN:i:")),
    ("hifiasm", re.compile(r"^[hp]tg\d{6}[lc](\s|$)")),
    ("trycycler", re.compile(r"^cluster_\d+_consensus")),
    ("flye", re.compile(r"^contig_\d+(\s|$)")),  # medaka/polypolish keep this name too
]
WRAPPER = re.compile(r"sw=(shovill|dragonflye)-?(\w*)")  # renamed contigs, real assembler in sw=


def detect_assembler(headers):
    votes = Counter()
    for h in headers:
        w = WRAPPER.search(h)
        if w:
            votes[w[2] or w[1]] += 1
            continue
        for name, pat in ASM_PATTERNS:
            if pat.search(h):
                votes[name] += 1
                break
    if not votes:
        return "unknown", 0.0
    name, n = votes.most_common(1)[0]
    # A few contigs keeping an assembler's naming style (1 of 41 = 2.4% in review) must not decide the
    # rules for the whole assembly; below a majority it's unknown, which has its own safe handling.
    return (name if n / len(headers) >= T["assembler_min_agreement"] else "unknown"), n / len(headers)


def detect_polishing(headers):
    joined = " ".join(headers).lower()
    return [p for p in ("pilon", "polypolish", "medaka", "racon", "polca") if p in joined]


def contig_info(header):
    """Circularity and coverage when the header carries them."""
    cov = re.search(r"_cov_([\d.]+)|depth=([\d.]+)x|^Contig_\d+_([\d.]+)(?:_Circ)?(?:\s|$)|cov=([\d.]+)", header)
    cov = float(next(g for g in cov.groups() if g).rstrip(".")) if cov else None
    circ = None  # circular=Y/N: Dragonflye; circular-yes/no: myloasm (real headers, Autocycler paper archive)
    if re.search(r"\bcircular[=-](true|yes|y)\b|suggestCircular=yes|XO:i:1|^[hp]tg\d{6}c|topology=circular", header, re.I):
        circ = True
    elif re.search(r"\bcircular[=-](false|no|n)\b|suggestCircular=no|XO:i:0|^[hp]tg\d{6}l", header, re.I):
        circ = False
    return cov, circ


def read_flye_info(path):
    info = {}
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        info[f[0]] = (float(f[2]), f[3] == "Y")
    return info


def assembly_metrics(path, flye_info=None):
    contigs, seq, header = [], [], None

    def flush():
        if header is not None:
            contigs.append((header, b"".join(seq)))

    for line in _lines(path):
        line = line.rstrip()
        if line.startswith(b">"):
            flush()
            header, seq = line[1:].decode(errors="replace"), []
        elif header is None and line:
            sys.exit(f"{path}: this is a FASTQ reads file, not an assembly. Pass it with --reads, and give the "
                     "assembly FASTA (lines starting with '>') with --assembly." if line.startswith(b"@") else
                     f"{path}: not a FASTA file; it starts with '{line[:30].decode(errors='replace')}'. "
                     "Takivisor needs the assembly as FASTA (header lines starting with '>'); "
                     "export it as FASTA from your assembler or annotation tool.")
        else:
            seq.append(line.upper())
    flush()
    lengths = [len(s) for _, s in contigs]
    total = sum(lengths)
    if not total:
        sys.exit(f"{path}: the file has no sequence in it (empty or headers only); check the assembly step finished.")

    headers = [h for h, _ in contigs]
    counts = Counter()
    for _, s in contigs:
        counts.update(s)  # C-level byte counting
    gc = counts[ord("G")] + counts[ord("C")]
    acgt = gc + counts[ord("A")] + counts[ord("T")]
    asm, conf = detect_assembler(headers)

    per_contig = []
    fi = read_flye_info(flye_info) if flye_info else {}
    for i, (h, l) in enumerate(zip(headers, lengths)):
        cov, circ = contig_info(h)
        name = h.split()[0] if h.split() else f"unnamed_{i + 1}"
        if name in fi:
            cov, circ = fi[name]
        per_contig.append({"name": name, "length": l, "cov": cov, "circular": circ,
                           "repeat": "suggestRepeat=yes" in h})

    orf = {}
    for code in (11, 4):
        L = sorted((x for _, seq in contigs for x in orf_lengths(seq, code)), reverse=True)
        run, half = 0, sum(L) / 2
        orf[code] = next((x // 3 for x in L if (run := run + x) >= half), 0)
    dcs = [{"name": pc["name"], "dcs_bp": h * T["dcs_asm_step"]}
           for pc, (_, seq) in zip(per_contig, contigs)
           if (h := dcs_hits(seq, _DCS_ASM, DCS_ASM_K, T["dcs_asm_step"])) >= T["dcs_asm_min_hits"]]
    hist = Counter(lengths)
    # 99% GC range of read-length windows (step W/4), per GC_CLASSES size: reads outside it look like
    # sequence the assembly doesn't contain. Windows with N are skipped.
    win = {}
    for W in GC_CLASSES:
        vals = sorted(int(100 * (w.count(b"G") + w.count(b"C")) / W)
                      for _, seq in contigs for i in range(0, len(seq) - W + 1, max(1, W // 4))
                      if b"N" not in (w := seq[i:i + W]))
        if len(vals) >= 200:
            win[W] = (vals[len(vals) // 200], vals[-1 - len(vals) // 200])
    return {
        "dcs_contigs": dcs,
        "orf_n50_aa": orf,
        "gc_window_range": win,
        "path": str(path),
        "assembler": asm,
        "assembler_confidence": round(conf, 2),
        "polishing_hints": detect_polishing(headers),
        "n_contigs": len(lengths),
        "n_contigs_200": sum(1 for x in lengths if x >= 200),
        "total_length": total,
        "largest": max(lengths),
        "n50": nx(hist, 0.5),
        "n90": nx(hist, 0.9),
        "gc_pct": round(100 * gc / acgt, 2) if acgt else 0.0,
        "n_pct": round(100 * counts[ord("N")] / total, 3),
        "contigs": per_contig,
        "dup_names": sorted(n for n, k in Counter(c["name"] for c in per_contig).items() if k > 1),
    }


# ---------------- reads ----------------

UUID = re.compile(r"^@[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}")
ILLUMINA_NAME = re.compile(r"^@\S*\s*\S*:\d+:[\w-]+:\d+:\d+:\d+:\d+")


# Read-GC is compared per length class against same-length windows of the assembly, so ONT reads of
# very different lengths are each judged against a matching spread. Classes are window sizes (bp).
GC_CLASSES = (150, 300, 1000, 3000, 10000)
# Illumina TruSeq / Nextera adapter starts: GC is counted only up to the first one (read-through
# adapters lowered raw GC from 56.9% to 54.8% on real run ERR10992850).
ADAPTERS = (b"AGATCGGAAGAGC", b"CTGTCTCTTATACACATCT")


# ONT DNA control strand (DCS): the 3.6 kb 3' end of phage lambda (NC_001416.1, 44,943-48,502; public
# domain). Real run ERR8958862 was 18% DCS reads and Raven built 2.7 kb of it into the P. aeruginosa
# chromosome, invisible to size/GC/CheckM2. Reads: sampled 15-mers tolerate ONT errors; assembly: exact
# 21-mers (a random 21-mer hit is ~1e-9 per lookup). Lambda-like prophages match too (29 of 51 pre-2014
# enterobacterial genomes: 250-4,515 bp), so an assembly hit only counts as DCS when the reads contain DCS.
DCS = (
    "CAAACTGAGCCGTAGCCACTGTCTGTCCTGAATTCATTAGTAATAGTTACGCTGCGGCCTTTTACACATGACCTTCGTGAAAGCGGGTGGCAGGAG"
    "GTCGCGCTAACAACCTCCTGCCGTTTTGCCCGTGCATATCGGTCACGAACAAATCTGATTACTAAACACAGTAGCCTGGATTTGTTCTATCAGTAA"
    "TCGACCTTATTCCTAATTAAATAGAGCAAATCCCCTTATTGGGGGTAAGACATGAAGATGCCAGAAAAACATGACCTGTTGGCCGCCATTCTCGCG"
    "GCAAAGGAACAAGGCATCGGGGCAATCCTTGCGTTTGCAATGGCGTACCTTCGCGGCAGATATAATGGCGGTGCGTTTACAAAAACAGTAATCGAC"
    "GCAACGATGTGCGCCATTATCGCCTGGTTCATTCGTGACCTTCTCGACTTCGCCGGACTAAGTAGCAATCTCGCTTATATAACGAGCGTGTTTATC"
    "GGCTACATCGGTACTGACTCGATTGGTTCGCTTATCAAACGCTTCGCTGCTAAAAAAGCCGGAGTAGAAGATGGTAGAAATCAATAATCAACGTAA"
    "GGCGTTCCTCGATATGCTGGCGTGGTCGGAGGGAACTGATAACGGACGTCAGAAAACCAGAAATCATGGTTATGACGTCATTGTAGGCGGAGAGCT"
    "ATTTACTGATTACTCCGATCACCCTCGCAAACTTGTCACGCTAAACCCAAAACTCAAATCAACAGGCGCCGGACGCTACCAGCTTCTTTCCCGTTG"
    "GTGGGATGCCTACCGCAAGCAGCTTGGCCTGAAAGACTTCTCTCCGAAAAGTCAGGACGCTGTGGCATTGCAGCAGATTAAGGAGCGTGGCGCTTT"
    "ACCTATGATTGATCGTGGTGATATCCGTCAGGCAATCGACCGTTGCAGCAATATCTGGGCTTCACTGCCGGGCGCTGGTTATGGTCAGTTCGAGCA"
    "TAAGGCTGACAGCCTGATTGCAAAATTCAAAGAAGCGGGCGGAACGGTCAGAGAGATTGATGTATGAGCAGAGTCACCGCGATTATCTCCGCTCTG"
    "GTTATCTGCATCATCGTCTGCCTGTCATGGGCTGTTAATCATTACCGTGATAACGCCATTACCTACAAAGCCCAGCGCGACAAAAATGCCAGAGAA"
    "CTGAAGCTGGCGAACGCGGCAATTACTGACATGCAGATGCGTCAGCGTGATGTTGCTGCGCTCGATGCAAAATACACGAAGGAGTTAGCTGATGCT"
    "AAAGCTGAAAATGATGCTCTGCGTGATGATGTTGCCGCTGGTCGTCGTCGGTTGCACATCAAAGCAGTCTGTCAGTCAGTGCGTGAAGCCACCACC"
    "GCCTCCGGCGTGGATAATGCAGCCTCCCCCCGACTGGCAGACACCGCTGAACGGGATTATTTCACCCTCAGAGAGAGGCTGATCACTATGCAAAAA"
    "CAACTGGAAGGAACCCAGAAGTATATTAATGAGCAGTGCAGATAGAGTTGCCCATATCGATGGGCAACTCATGCAATTATTGTGAGCAATACACAC"
    "GCGCTTCCAGCGGAGTATAAATGCCTAAAGTAATAAAACCGAGCAATCCATTTACGAATGTTTGCTGGGTTTCTGTTTTAACAACATTTTCTGCGC"
    "CGCCACAAATTTTGGCTGCATCGACAGTTTTCTTCTGCCCAATTCCAGAAACGAAGAAATGATGGGTGATGGTTTCCTTTGGTGCTACTGCTGCCG"
    "GTTTGTTTTGAACAGTAAACGTCTGTTGAGCACATCCTGTAATAAGCAGGGCCAGCGCAGTAGCGAGTAGCATTTTTTTCATGGTGTTATTCCCGA"
    "TGCTTTTTGAAGTTCGCAGAATCGTATGTGTAGAAAATTAAACAAACCCTAAACAATGAGTTGAAATTTCATATTGTTAATATTTATTAATGTATG"
    "TCAGGTGCGATGAATCGTCATTGTATTCCCGGATTAACTATGTCCACAGCCCTGACGGGGAACTTCTCTGCGGGAGTGTCCGGGAATAATTAAAAC"
    "GATGCACACAGGGTTTAGCGCGTACACGTATTGCATTATGCCAACGCCCCGGTGCTGACACGGAAGAAACCGGACGTTATGATTTAGCGTGGAAAG"
    "ATTTGTGTAGTGTTCTGAATGCTCTCAGTAAATAGTAATGAATTATCAAAGGTATAGTAATATCTTTTATGTTCATGGATATTTGTAACCCATCGG"
    "AAAACTCCTGCTTTAGCAAGATTTTCCCTGTATTGCTGAAATGTGATTTCTCTTGATTTCAACCTATCATAGGACGTTTCTATAAGATGCGTGTTT"
    "CTTGAGAATTTAACATTTACAACCTTTTTAAGTCCTTTTATTAACACGGTGTTATCGTTTTCTAACACGATGTGAATATTATCTGTGGCTAGATAG"
    "TAAATATAATGTGAGACGTTGTGACGTTTTAGTTCAGAATAAAACAATTCACAGTCTAAATCTTTTCGCACTTGATCGAATATTTCTTTAAAAATG"
    "GCAACCTGAGCCATTGGTAAAACCTTCCATGTGATACGAGGGCGCGTAGTTTGCATTATCGTTTTTATCGTTTCAATCTGGTCTGACCTCCTTGTG"
    "TTTTGTTGATGATTTATGTCAAATATTAGGAATGTTTTCACTTAATAGTATTGGTTGCGTAACAAAGTGCGGTCCTGCTGGCATTCTGGAGGGAAA"
    "TACAACCGACAGATGTATGTAAGGCCAACGTGCTCAAATCTTCATACAGAAAGATTTGAAGTAATATTTTAACCGCTAGATGAAGAGCAAGCGCAT"
    "GGAGCGACAAAATGAATAAAGAACAATCTGCTGATGATCCCTCCGTGGATCTGATTCGTGTAAAAAATATGCTTAATAGCACCATTTCTATGAGTT"
    "ACCCTGATGTTGTAATTGCATGTATAGAACATAAGGTGTCTCTGGAAGCATTCAGAGCAATTGAGGCAGCGTTGGTGAAGCACGATAATAATATGA"
    "AGGATTATTCCCTGGTGGTTGACTGATCACCATAACTGCTAATCATTCAAACTATTTAGTCTGTGACAGAGCCAACACGCAGTCTGTCACTGTCAG"
    "GAAAGTGGTAAAACTGCAACTCAATTACTGCAATGCCCTCGTAATTAAGTGAATTTACAATATCGTCCTGTTCGGAGGGAAGAACGCGGGATGTTC"
    "ATTCTTCATCACTTTTAATTGATGTATATGCTCTCTTTTCTGACGTTAGTCTCCGACGGCAGGCTTCAATGACCCAGGCTGAGAAATTCCCGGACC"
    "CTTTTTGCTCAAGAGCGATGTTAATTTGTTCAATCATTTGGTTAGGAAAGCGGATGTTGCGGGTTGTTGTTCTGCGGGTTCTGTTCTTCGTTGACA"
    "TGAGGTTGCCCCGTATTCAGTGTCGCTGATTTGTATTGTCTGAAGTTGTTTTTACGTTAAGTTGATGCAGATCAATTAATACGATACCTGCGTCAT"
    "AATTGATTATTTGACGTGGTTTGATGGCCTCCACGCACGTTGTGATATGTAGATGATAATCATTATCACTTTACGGGTCCTTTCCGGTGATCCGAC"
    "AGGTTACG"
)
_COMP = bytes.maketrans(b"ACGT", b"TGCA")


def _kmer_set(seq, k):
    seq = seq.encode()
    rc = seq.translate(_COMP)[::-1]
    return {x[i:i + k] for x in (seq, rc) for i in range(len(x) - k + 1)}


DCS_READ_K, DCS_ASM_K = 15, 21
DCS_GC = round(100 * sum(c in "GC" for c in DCS) / len(DCS), 1)   # 44.1
_DCS_READ, _DCS_ASM = _kmer_set(DCS, DCS_READ_K), _kmer_set(DCS, DCS_ASM_K)


def dcs_hits(seq, kset, k, step):
    return sum(seq[i:i + k] in kset for i in range(0, len(seq) - k + 1, step))


def is_dcs_read(s):
    """DCS reads are whole ~3.56 kb amplicons matching along most of their length. Lambda-like prophage
    reads (common in E. coli) are random-length genomic reads with a partial match: on real runs
    97% of DCS-hit reads in ERR8958862 had the amplicon shape, 0 of 195 in two E. coli runs did."""
    if not T["dcs_len"][0] <= len(s) <= T["dcs_len"][1]:
        return False
    k, st = DCS_READ_K, T["dcs_read_step"]
    pos = [i for i in range(0, len(s) - k + 1, st) if s[i:i + k] in _DCS_READ]
    return len(pos) >= T["dcs_read_min_hits"] and (pos[-1] - pos[0] + k) >= T["dcs_span"] * len(s)


# Stop-to-stop open reading frames: frameshift (indel) errors break genes into short ORFs, so ORF N50
# is a built-in stand-in for CheckM2's gene length when no report is given. Genera using genetic code
# 4 (TGA = Trp) get TAA/TAG only, or their ORFs look falsely short (Mycoplasma 197 aa -> 298 aa).
_STOPS = {11: re.compile(rb"(?=(TAA|TAG|TGA))"), 4: re.compile(rb"(?=(TAA|TAG))")}
CODE4_GENERA = {"mycoplasma", "mycoplasmoides", "mycoplasmopsis", "metamycoplasma", "mesomycoplasma",
                "malacoplasma", "ureaplasma", "spiroplasma", "mesoplasma", "entomoplasma"}


def orf_lengths(seq, code=11, min_nt=300):
    """Lengths (nt) of complete stop-to-stop open frames >= min_nt, both strands, no N."""
    out = []
    for s in (seq, seq.translate(_COMP)[::-1]):
        frames = [[], [], []]
        for m in _STOPS[code].finditer(s):
            frames[m.start() % 3].append(m.start())
        for pos in frames:
            for a, b in zip(pos, pos[1:]):
                if b - a - 3 >= min_nt and b"N" not in s[a + 3:b]:
                    out.append(b - a - 3)
    return out


def _gc_class(n):
    """Nearest class on a log scale; None below 100 bp (too short for a stable GC value)."""
    return None if n < 100 else min(GC_CLASSES, key=lambda c: abs(math.log(n / c)))


def fastq_stats(path):
    """Accumulated histograms: memory stays flat however many reads there are."""
    acc = {"lens": Counter(), "qhist": Counter(), "first": None, "gc": 0, "gc_bases": 0,
           "adapter_reads": 0, "read_gc": Counter(), "long_reads": 0, "dcs_reads": 0}
    prog = _Progress(path, "reading")
    it = _lines(path, prog.fh)
    n_reads = bases = 0
    for h in it:
        if not h.strip():
            continue  # blank lines (e.g. at the end of the file) are not reads
        if acc["first"] is None:
            if not h.startswith(b"@"):
                sys.exit(f"{path}: this is a FASTA file, but --reads needs FASTQ reads (with quality lines). "
                         "If it is your assembly, pass it with --assembly." if h.startswith(b">") else
                         f"{path}: not a FASTQ reads file (first line should start with '@').")
            acc["first"] = h.decode(errors="replace").rstrip()
        s = next(it, b"").rstrip().upper()
        plus = next(it, b"")
        q = next(it, b"").rstrip()
        # 4-line records only: wrapped or cut-off files would otherwise give wrong depth/quality silently
        if not h.startswith(b"@") or not plus.startswith(b"+") or len(q) != len(s):
            sys.exit(f"{path}: malformed FASTQ at read {n_reads + 1} (expected 4 lines per read: '@' header, "
                     "sequence, '+', quality of the same length). The file may be cut off (incomplete copy or "
                     "download) or wrapped over several lines; re-copy it or rewrite it as 4-line FASTQ.")
        acc["qhist"].update(q)
        acc["lens"][len(s)] += 1
        n_reads += 1
        bases += len(s)
        if n_reads % 5000 == 0:
            prog.tick(n_reads, bases)
        cut = min((i for i in (s.find(a) for a in ADAPTERS) if i >= 0), default=len(s))
        acc["adapter_reads"] += cut < len(s)
        t = s[:cut]
        g = t.count(b"G") + t.count(b"C")
        acc["gc"] += g
        acc["gc_bases"] += len(t)
        c = _gc_class(len(t))
        if c:
            acc["read_gc"][(c, int(100 * g / len(t)))] += 1
        if len(s) >= 1000:  # DCS is an ONT control; only long reads are scanned (keeps Illumina fast)
            acc["long_reads"] += 1
            acc["dcs_reads"] += is_dcs_read(s)
    if acc["first"] is None:
        sys.exit(f"{path}: no reads")
    prog.done(n_reads, bases)
    return acc


def summarize_reads(files, acc):
    lens, qhist = acc["lens"], acc["qhist"]
    nbases = sum(qhist.values())
    if not nbases:
        sys.exit(f"{files}: reads have no bases")
    err = sum(n * 10 ** (-(c - 33) / 10) for c, n in qhist.items())
    q = lambda t: 100 * sum(n for c, n in qhist.items() if c - 33 >= t) / nbases
    n_reads, run, median = sum(lens.values()), 0, 0
    for l in sorted(lens):
        run += lens[l]
        if run >= n_reads / 2:
            median = l
            break
    len_hist = Counter()
    for l, n in lens.items():
        if l:
            len_hist[int(10 * math.log10(l))] += l * n
    return {
        "files": [str(f) for f in files],
        "n_reads": n_reads,
        "bases": sum(l * n for l, n in lens.items()),
        "min_len": min(lens),
        "median_len": median,
        "max_len": max(lens),
        "read_n50": nx(lens),
        # probability-averaged: the honest overall error rate, lower than averaging Q scores
        "mean_q": round(-10 * math.log10(err / nbases), 1) if err else 60.0,
        "pct_q20": round(q(20), 1),
        "pct_q30": round(q(30), 1),
        # NCBI SRA-lite (and some re-exported data) replaces qualities with 1-2 placeholder values
        # (real run SRR11620701: only '?' Q30 and '$' Q3), so quality-based stats mean nothing
        "qual_simplified": len(qhist) <= 2,
        # tech-detection signatures (see detect_tech)
        "n_qual_values": len(qhist),
        "max_q": max(qhist) - 33,
        "pct_q60": round(100 * sum(n for c, n in qhist.items() if c - 33 >= 60) / nbases, 2),
        "frac_at_max_len": round(lens[max(lens)] / n_reads, 3),
        "near_mode_frac": round(sum(n for l, n in lens.items()
                                    if abs(l - lens.most_common(1)[0][0]) <= max(2, lens.most_common(1)[0][0] // 100))
                                / n_reads, 3),
        "gc_pct": round(100 * acc["gc"] / max(1, acc["gc_bases"]), 2),
        "adapter_pct": round(100 * acc["adapter_reads"] / n_reads, 1),
        **({"dcs_pct": round(100 * acc["dcs_reads"] / acc["long_reads"], 2)} if acc["long_reads"] else {}),
        "read_gc": acc["read_gc"],  # {(length class, GC%): reads}; used by the read-GC outlier check
        # bases per log10 length bin (10 per decade, key = floor(10*log10 len)); for the HTML length chart
        "len_hist": len_hist,
        "first_header": acc["first"][:200],
    }


PACBIO_NAME = re.compile(r"\bm\d{5}[a-z]?_\d{6}_\d{6}(?:_s\d+)?/\d+(/ccs|/\d+_\d+)?")
# Standard Illumina cycle counts (+1 for the extra base some runs keep).
ILLUMINA_LENGTHS = {c + d for c in (35, 50, 75, 100, 125, 150, 250, 300) for d in (0, 1)}


def detect_tech(r):
    """Read names first (ENA keeps them, SRA/NCBI often strips them), then signatures from how each
    platform reports quality and length. Evidence (2026-09-30, first 2-5k reads of real runs):
      HiFi, unbinned (Sequel/Sequel II CCS): 92-95% of bases >= Q60
      HiFi, Revio: long reads with only 7 binned quality values, max Q40
      CLR (Sequel, Sequel II subreads): every base Q0
      ONT (R9-R10, Guppy/Dorado, incl. amplicons): 30-69 quality values, <= 2.3% bases >= Q60
      Illumina (MiSeq/NextSeq/HiSeq/NovaSeq): longest read = a standard cycle count, 82-100% of reads
        at the modal length (>= 5% exactly at it even after hard trimming); ONT short/amplicon runs 2-33%"""
    h = r["first_header"]
    if "basecall_model_version_id=" in h or UUID.match(h) or re.search(r"\s[0-9a-f]{8}-[0-9a-f]{4}-", h):
        m = (re.search(r"basecall_model_version_id=(\S+)", h)     # MinKNOW / Dorado < 1.2
             or re.search(r"RG:Z:[0-9a-f-]{36}_((?:dna|rna)[\w.]*?@v[\d.]+)", h))   # Dorado >= 1.2 (ONT HTS spec)
        return "ont", "header" + (f" (model {m.group(1)})" if m else " (ONT read UUID)")
    pb = PACBIO_NAME.search(h)
    if pb:
        if (pb.group(1) or "").startswith("/") and pb.group(1) != "/ccs" or r["max_q"] == 0:
            return "clr", "header (PacBio subread name" + (", all bases Q0)" if r["max_q"] == 0 else ")")
        return "hifi", "header (PacBio read name, real quality values)"
    if ILLUMINA_NAME.match(h):
        return "illumina", "header (Illumina read name)"
    if r["pct_q60"] >= 50:
        return "hifi", f"stats ({r['pct_q60']}% of bases >= Q60)"
    long = r["max_len"] > 600
    if long and r["max_q"] == 0:
        return "clr", "stats (long reads, every base Q0: PacBio subreads)"
    if long and 3 <= r["n_qual_values"] <= 8 and r["max_q"] <= 41:
        return "hifi", f"stats (long reads with {r['n_qual_values']} binned quality values, as PacBio Revio)"
    if not long:
        # Untrimmed Illumina reads end exactly at the cycle count; even hard quality trimming leaves a
        # share there (real MiSeq cut at the first base < Q20: 16%), while ONT has ~1 read per length.
        if r["max_len"] in ILLUMINA_LENGTHS and (r["near_mode_frac"] >= 0.5 or r.get("frac_at_max_len", 0) >= 0.05):
            return "illumina", f"stats ({r.get('frac_at_max_len', r['near_mode_frac']):.0%} of reads end exactly at a {r['max_len']} bp cycle length)"
        return "ont", "stats (short reads of variable length: ONT amplicons or similar; pass --tech if wrong)"
    if r.get("qual_simplified"):
        return "ont", "stats (long reads; quality scores simplified, so HiFi can't be told apart - pass --tech hifi if PacBio)"
    return "ont", "stats (long, variable length)"


def read_sets(paths, force_tech=None):
    """Each file is classified on its own, then files of the same tech are merged,
    so R1+R2+ONT gives one Illumina set and one ONT set (hybrid)."""
    groups = {}
    for p in paths:
        acc = fastq_stats(p)
        tech, basis = (force_tech, "user") if force_tech else detect_tech(summarize_reads([p], acc))
        g = groups.get(tech)
        if g is None:
            groups[tech] = {**acc, "files": [p], "basis": basis}
            continue
        g["files"].append(p)
        for k in ("lens", "qhist", "read_gc"):
            g[k].update(acc[k])
        for k in ("gc", "gc_bases", "adapter_reads", "long_reads", "dcs_reads"):
            g[k] += acc[k]
    return {t: {**summarize_reads(g["files"], g), "source": "FASTQ", "tech_basis": g["basis"]}
            for t, g in groups.items()}


# ---------------- existing tool reports ----------------

def _num(s):
    return float(s.replace(",", "").rstrip("%"))


def parse_busco(p):
    txt = Path(p).read_text()
    m = re.search(r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%,n:(\d+)", txt)
    lin = re.search(r"lineage dataset is: (\S+)", txt)
    tot = re.search(r"(\d+)\s+Total length", txt)
    src = re.search(r"notation for file (\S+)", txt)
    return {"name": Path(src[1]).stem if src else None,"complete": float(m[1]), "single": float(m[2]), "duplicated": float(m[3]),
            "fragmented": float(m[4]), "missing": float(m[5]), "n": int(m[6]),
            "lineage": lin[1] if lin else None, "total_length": int(tot[1]) if tot else None}


def parse_checkm(p):
    """CheckM2 quality_report.tsv, legacy CheckM --tab_table, or legacy CheckM's default printed
    table (dashed rules, columns padded with spaces; real file from MultiQC test-data)."""
    lines = [l for l in Path(p).read_text().splitlines() if l.strip() and not l.startswith("---")]
    split = (lambda l: l.split("\t")) if "\t" in lines[0] else (lambda l: re.split(r"\s{2,}", l.strip()))
    cols = split(lines[0])
    if "Completeness" not in cols:
        raise KeyError("Completeness")
    rows = []
    for line in lines[1:]:
        r = dict(zip(cols, split(line)))
        get = lambda k, f=float: f(r[k]) if r.get(k) not in (None, "", "None") else None
        rows.append({"name": r.get("Name") or r.get("Bin Id"),
                     "tool": "CheckM2" if "Coding_Density" in cols else "CheckM",
                     "completeness": get("Completeness"), "contamination": get("Contamination"),
                     "coding_density": get("Coding_Density"),
                     "avg_gene_length_aa": get("Average_Gene_Length", lambda v: round(float(v), 1)),
                     "strain_heterogeneity": get("Strain heterogeneity"),
                     "marker_lineage": r.get("Marker lineage"),
                     "total_length": get("Genome_Size", int)})
    return rows


def parse_quast(p):
    """One row per assembly column (several assemblies, or --scaffolds' extra '_broken' column;
    real files from MultiQC test-data); load_reports picks the column matching --assembly."""
    rows = [l.split("\t") for l in Path(p).read_text().splitlines() if "\t" in l]
    head, by = rows[0], {r[0]: r[1:] for r in rows[1:]}
    if head[0] != "Assembly":
        raise ValueError("not a QUAST report.tsv")
    out = []
    for i, name in enumerate(head[1:]):
        get = lambda k: _num(by[k][i].split()[0]) if k in by and i < len(by[k]) and by[k][i].strip() not in ("", "-") else None
        # plain "Total length" skips contigs <500 bp, which breaks matching on fragmented assemblies
        total = get("Total length (>= 0 bp)") or get("Total length")
        out.append({"name": name or None, "total_length": int(total), "n_contigs": int(get("# contigs")),
                    "n50": int(get("N50")), "rrna_genes": get("# predicted rRNA genes"),
                    "genome_fraction": get("Genome fraction (%)"), "misassemblies": get("# misassemblies")})
    return out


def parse_fastp(p):
    d = json.loads(Path(p).read_text())
    b, a = d["summary"]["before_filtering"], d["summary"]["after_filtering"]
    if not b["total_reads"]:
        raise ValueError("fastp report has no reads")
    return {"n_reads": b["total_reads"], "bases": b["total_bases"], "bases_after": a["total_bases"],
            "pct_q20": round(b["q20_rate"] * 100, 1), "pct_q30": round(b["q30_rate"] * 100, 1),
            # GC after trimming: adapter read-through skews raw GC (real run ERR10992850: 54.8% raw vs
            # 56.9% trimmed, half the reads adapter-trimmed), and assemblies are built from trimmed reads
            "gc_pct": round(a.get("gc_content", b["gc_content"]) * 100, 2),
            "pct_passed": round(100 * a["total_reads"] / b["total_reads"], 1),
            "dup_pct": round(d.get("duplication", {}).get("rate", 0) * 100, 1)}


def parse_nanostats(p):
    """NanoPlot/NanoStat NanoStats.txt, both the classic 'Key: value' and the --tsv_stats layout
    (real files: MultiQC test-data, incl. NanoPlot 1.41.6 TSV and FASTA input without qualities)."""
    kv, lines = {}, Path(p).read_text().splitlines()
    if lines and len(lines[0].replace("General summary:", "x").split()) > 2:
        raise ValueError("multi-sample NanoStats (one column per sample); pass a per-sample file")
    for l in lines:
        k, sep, v = l.partition(":") if ":" in l else l.partition("\t")
        if sep and v.split():
            kv[k.strip()] = v.split()[0]
    pick = lambda *ks: _num(next(kv[k] for k in ks if k in kv))
    return {"n_reads": int(pick("Number of reads", "number_of_reads")),
            "bases": int(pick("Total bases", "number_of_bases")),
            "median_len": int(pick("Median read length", "median_read_length")),
            "read_n50": int(pick("Read length N50", "n50")),
            # NanoPlot averages per-read Q values, which overstates accuracy vs error-averaging.
            # Absent when NanoPlot was run on FASTA.
            **({"mean_read_q": pick("Mean read quality", "mean_qual")}
               if {"Mean read quality", "mean_qual"} & kv.keys() else {})}


def _sample_name(p):
    """Sample name from a report file name: 'teste_local_k2report.txt' -> 'teste_local'."""
    n = Path(p).name
    for suf in ("_nohuman", "_k2report", ".k2report", "_kraken2_report", "_kraken_report", ".kreport",
                ".sylphmpa", ".txt", ".tsv"):
        n = n.split(suf)[0] if suf in n else n
    return n


def _identity(tool, counts, p, n_reads=None):
    """Top species as a share of bacterial/archaeal species-level assignments (host, unclassified and
    reads stopping above species are left out: real UR03 report had 2,240 reads stopping at family)."""
    tot = sum(counts.values())
    if not tot:
        raise ValueError("no bacterial species-level assignments")
    (top, n1), *rest = sorted(counts.items(), key=lambda x: -x[1])
    second = rest[0] if rest else (None, 0)
    return {"tool": tool, "top": top, "top_share": round(100 * n1 / tot, 1), "second": second[0],
            "second_share": round(100 * second[1] / tot, 1), "name": _sample_name(p), "total_length": None,
            "n_assigned": n_reads}


def parse_identity(p):
    """Kraken2 report, sylph-tax (MetaPhlAn-style .sylphmpa) or GTDB-Tk summary; sniffed by content."""
    lines = Path(p).read_text().splitlines()
    if lines and lines[0].startswith("user_genome") and "classification" in lines[0]:  # GTDB-Tk
        cols = lines[0].split("\t")
        rows = [dict(zip(cols, l.split("\t"))) for l in lines[1:] if l.strip()]
        sps = [r["classification"].split("s__")[-1] for r in rows if "s__" in r.get("classification", "")]
        return _identity("GTDB-Tk", Counter(x for x in sps if x), p)
    if any(l.startswith("clade_name") for l in lines[:5]):  # sylph-tax / MetaPhlAn profile
        c = Counter()
        for l in lines:
            f = l.split("\t")
            if "|s__" in f[0] and "|t__" not in f[0] and f[0].startswith(("d__Bacteria", "d__Archaea", "k__Bacteria", "k__Archaea")):
                nm = f[0].split("|s__")[-1]   # GTDB 'Enterococcus_B faecium' keeps its suffix; MetaPhlAn 'Escherichia_coli'
                c[nm if " " in nm else nm.replace("_", " ")] += float(f[1])
        return _identity("sylph", c, p)
    # Kraken2 report: pct, clade reads, direct reads, [minimizers, distinct minimizers with
    # --report-minimizer-data,] rank, taxid, indented name. Kraken2 < 2.17 with 2025+ NCBI taxonomy could print
    # the domain at a root rank (R1/R2) instead of D, so domains are also recognised by name.
    c, domain = Counter(), None
    for l in lines:
        f = l.split("\t")
        if len(f) not in (6, 8) or not f[1].strip().isdigit():
            raise ValueError("not a Kraken2 report")
        rank, name = f[-3].strip(), f[-1].strip()
        if rank == "D" or (rank[:1] == "R" and name in ("Bacteria", "Archaea", "Eukaryota", "Viruses")):
            domain = name
        elif rank == "S" and domain in ("Bacteria", "Archaea"):
            c[name] += int(f[1])
    # Bracken writes the same layout (-w); known by file name (*bracken*, *.breport)
    tool = "Bracken" if "bracken" in Path(p).name.lower() or ".breport" in Path(p).name else "Kraken2"
    return _identity(tool, c, p, n_reads=sum(c.values()))


def parse_contig_kraken(p):
    """Kraken2 per-sequence output run on the contigs (C/U, name, taxid, length, k-mer trace).
    Kraken2's own call is ignored: with a confidence threshold long contigs come out 'U', so the
    trace's k-mer votes per taxid are kept and resolved to species later (needs ktaxonomy.tsv)."""
    rows = []
    for l in _lines(p):
        f = l.decode().rstrip("\n").split("\t")
        if len(f) < 5 or f[0] not in ("C", "U") or not f[3].isdigit():
            raise ValueError("not Kraken2 per-sequence output")
        votes = Counter()
        for x in f[4].split():
            t, _, n = x.rpartition(":")
            if t not in ("0", "A", "|"):
                votes[t] += int(n)
        rows.append((f[1], int(f[3]), votes))
    if not rows:
        raise ValueError("empty")
    return {"tool": "Kraken2 (contigs)", "rows": rows, "name": _sample_name(p), "total_length": None}


def load_ktaxonomy(dirs):
    """taxid -> species name for every taxon at or below species, from the Kraken2 DB's ktaxonomy.tsv
    (prebuilt DBs ship it). None if not given."""
    path = next((q for d in dirs for q in sorted(Path(d).rglob("ktaxonomy.tsv"))), None)
    if not path:
        return None
    node = {}
    for l in _lines(path):
        f = [x.strip() for x in l.decode().split("|")]
        if len(f) >= 5:
            node[f[0]] = (f[1], f[2], f[4])
    sp = {}
    for t in node:
        u, seen = t, set()
        while u in node and u not in seen and node[u][1][:1] == "S" and node[u][1] != "S":
            seen.add(u); u = node[u][0]   # S1/S2 (strain, subsp.) -> parent
        if u in node and node[u][1] == "S":
            sp[t] = node[u][2]
    return sp


def parse_confindr(p):
    """ConFindr confindr_report.csv, one row per sample. 0.8.2 writes Sample,Genus,NumContamSNVs,ContamStatus,
    BasesExamined,DatabaseDownloadDate (PercentContam only in some versions). ContamStatus False with 0 bases
    examined means nothing was checked: genus without a database (no rMLST) or Genus 'Error processing sample'
    (real runs, 2026-10-08), so it is read as not determined. Multi-genus 'A:B' is True from mash, kept."""
    lines = [l for l in Path(p).read_text().splitlines() if l.strip()]
    head = lines[0].split(",")
    if "ContamStatus" not in head:
        raise ValueError("not a ConFindr report")
    num = lambda v, f=float: f(float(v)) if v not in (None, "", "ND", "NA") else None
    out = []
    for line in lines[1:]:
        r = dict(zip(head, line.split(",")))
        st = r["ContamStatus"].strip().lower()
        bases = num(r.get("BasesExamined"), int)
        out.append({"name": r.get("Sample"), "tool": "ConFindr", "genus": r.get("Genus"),
                    "contaminated": True if st == "true" else False if st == "false" and bases != 0 else None,
                    "snvs": num(r.get("NumContamSNVs"), int), "pct": num(r.get("PercentContam")), "bases": bases})
    return out


# AMR gene calls that are incomplete. AMRFinderPlus Method codes (ncbi/amr test outputs): PARTIAL_CONTIG_END* =
# gene runs off a contig end, PARTIAL* = < 90% of the reference, INTERNAL_STOP. RGI: Prodigal partial=10/01/11 in
# ORF_ID = gene at a contig edge; < 90% of reference length (AMRFinderPlus's partial cut-off); Loose = discovery.
def parse_amrfinder(p):
    rows = [l.rstrip("\n").split("\t") for l in Path(p).read_text().splitlines() if l.strip()]
    head = rows[0]
    if "Method" not in head or not {"Element symbol", "Gene symbol"} & set(head):
        raise ValueError("not an AMRFinderPlus report")
    hits = []
    for r in rows[1:]:
        d = dict(zip(head, r))
        m = d["Method"]
        hits.append({"gene": d.get("Element symbol") or d.get("Gene symbol"), "contig": d.get("Contig id"),
                     "issue": "edge" if m.startswith("PARTIAL_CONTIG_END") else "partial" if m.startswith("PARTIAL")
                     else "stop" if m == "INTERNAL_STOP" else None})
    return {"name": None, "tool": "AMRFinderPlus", "hits": hits}


def parse_rgi(p):
    rows = [l.rstrip("\n").split("\t") for l in Path(p).read_text().splitlines() if l.strip()]
    head = rows[0]
    if "Best_Hit_ARO" not in head or "Cut_Off" not in head:
        raise ValueError("not an RGI report")
    hits = []
    for r in rows[1:]:
        d = dict(zip(head, r))
        orf = d.get("ORF_ID", "")
        part = re.search(r"partial=(\d\d)", orf)
        pct = d.get("Percentage Length of Reference Sequence", "").strip()
        pct = _num(pct) if pct else None
        # Contig column is '<contig>_<ORF number>' for Prodigal ORFs; rRNA-model rows hold coordinates there
        # and protein input leaves it empty (real RGI output, hAMRonization test data), so only ORFs name contigs.
        contig = d.get("Contig", "").rsplit("_", 1)[0] if " # " in orf else None
        hits.append({"gene": d["Best_Hit_ARO"], "contig": contig,
                     "issue": "edge" if part and part[1] != "00" else "partial" if pct is not None and pct < 90
                     else "loose" if d["Cut_Off"] == "Loose" else None})
    return {"name": None, "tool": "RGI", "hits": hits}


def parse_mosdepth(p):
    """mosdepth <prefix>.mosdepth.summary.txt (chrom length bases mean min max; 'total' and '*_region' rows skipped)
    plus <prefix>.mosdepth.global.dist.txt next to it: per chrom, cumulative fraction of bases at >= depth."""
    rows = [l.split("\t") for l in Path(p).read_text().splitlines() if l.strip()]
    if rows[0][:4] != ["chrom", "length", "bases", "mean"]:
        raise ValueError("not a mosdepth summary")
    contigs = {r[0]: int(r[1]) for r in rows[1:] if r[0] != "total" and not r[0].endswith("_region")}
    prefix = Path(p).name[:-len(".mosdepth.summary.txt")]
    out = {"name": prefix, "tool": "mosdepth", "total_length": sum(contigs.values()), "contigs": contigs,
           "mean": round(sum(float(r[2]) for r in rows[1:] if r[0] in contigs) / max(sum(contigs.values()), 1), 1),
           "at": {}}
    dist = Path(p).with_name(prefix + ".mosdepth.global.dist.txt")
    if dist.is_file():
        tot = {int(d): float(f) for c, d, f in (l.split("\t") for l in dist.read_text().splitlines() if l.strip())
               if c == "total"}
        # no row for a depth above the maximum: no base reaches it
        out["at"] = {x: round(100 * tot.get(x, 0.0), 1) for x in (1, 10, 20, 30)}
    return out


def mosdepth_row(rep, asm):
    """(report, what it was mapped to) if usable: this assembly (same length), or a reference with none of this
    assembly's contig names when it is the only mosdepth report given (WHO TB maps to H37Rv). Else None."""
    r = rep.get("mosdepth")
    if not r:
        return None
    if r.get("matches_assembly"):
        return r, "this assembly"
    if r["n_candidates"] == 1 and not r["contigs"].keys() & {c["name"] for c in asm["contigs"]}:
        return r, f"a {r['total_length']:,} bp reference"
    return None


def confindr_row(rep):
    """The ConFindr row for this sample: a name match, or the only row of a one-sample report."""
    cf = rep.get("confindr")
    return cf if cf and (cf.get("matches_assembly") or cf["n_candidates"] == 1) else None


def amr_hits(r, asm):
    """The report's hits if its contigs are in this assembly (None if they aren't: another assembly). A hit cut by
    the end of a contig the assembler marked circular spans the circle's start/end: issue 'wrap', a complete gene
    (#30: seen in real Autocycler-paper assemblies, e.g. sul2 on a 6.8 kb circular plasmid)."""
    names = {c["name"] for c in asm["contigs"]}
    used = {h["contig"] for h in r["hits"] if h["contig"] not in (None, "", "NA")}
    if used - names:
        return None
    circ = {c["name"] for c in asm["contigs"] if c.get("circular")}
    return [{**h, "issue": "wrap"} if h["issue"] == "edge" and h["contig"] in circ else h for h in r["hits"]]


def contig_species(rows, tax):
    """Per contig: (name, length, species-level votes, top species, its % share)."""
    out = []
    for name, length, votes in rows:
        c = Counter()
        for t, n in votes.items():
            if t in tax:
                c[tax[t]] += n
        n = sum(c.values())
        top, k = c.most_common(1)[0] if c else (None, 0)
        out.append({"name": name, "length": length, "votes": n, "top": top,
                    "share": round(100 * k / n, 1) if n else 0, "sp_votes": dict(c.most_common(5))})
    return out


REPORTS = {  # key: (filename globs, parser)
    "busco": (["short_summary*.txt"], parse_busco),
    "checkm": (["quality_report.tsv", "checkm*.tsv"], parse_checkm),
    "quast": (["report.tsv"], parse_quast),
    "fastp": (["*fastp*.json"], parse_fastp),
    "nanostats": (["*NanoStats*.txt"], parse_nanostats),
    "contig_kraken": (["*kraken*output*", "*.kraken", "*.kraken2", "*kraken*.out"], parse_contig_kraken),
    "confindr": (["*confindr_report*.csv"], parse_confindr),
    "mosdepth": (["*.mosdepth.summary.txt"], parse_mosdepth),
    "amrfinder": (["*amrfinder*.tsv", "*amrfinder*.txt"], parse_amrfinder),
    "rgi": (["*rgi*.txt", "*card*.txt"], parse_rgi),
    "identity": (["*k2report*.txt", "*kraken*report*", "*.kreport", "*.breport*", "*bracken*species*", "*.sylphmpa", "*gtdbtk*summary*.tsv"], parse_identity),
}


# Default output names shared by many runs and tools: matching on them proves nothing (#1 review:
# a legacy CheckM report from another run was "verified" because both assemblies were final_assembly).
GENERIC_NAMES = {"final_assembly", "assembly", "assembly_filtered", "draft_assembly", "contigs",
                 "scaffolds", "consensus", "polished", "before_rr", "genome", "final", "output",
                 "contig", "scaffold", "flye", "raven", "canu", "unicycler", "spades", "skesa",
                 "medaka", "pilon", "polypolish", "trycycler", "report", "sample", "quality_report"}


def identifying(name):
    """A name counts as evidence only if it isn't a generic default."""
    return name if name and name.lower() not in GENERIC_NAMES and len(name) >= 3 else None


MATCH_TOL = 0.0005  # same assembly => same length. Polishing changed Raven's CFT073 assembly by
                   # 0.15%, so a pre-polish report must NOT count as a match for the polished file.


def load_reports(dirs, asm_len, asm_name=None):
    """Every matching file is parsed. Assembly reports are matched to --assembly by total
    length (then file name), so a folder holding several runs still gives the right one."""
    found = {}
    for key, (globs, parser) in REPORTS.items():
        cands, failed = [], []
        for d in dirs:
            for g in globs:
                for p in sorted(Path(d).rglob(g)):
                    if not p.is_file():
                        continue
                    try:
                        r = parser(p)
                    except (KeyError, ValueError, IndexError, TypeError, StopIteration, AttributeError):
                        failed.append(str(p))  # same file name, different tool: expected
                        continue
                    cands += [{**x, "file": str(p)} for x in (r if isinstance(r, list) else [r])]
        if not cands:
            if failed:
                print(f"warning: found but could not parse as {key}: {', '.join(failed)}", file=sys.stderr)
            continue
        # identical files (e.g. BUSCO's two copies) collapse; rows of one multi-genome file stay distinct
        cands = list({(c["file"], c.get("total_length"), c.get("name")): c for c in cands}.values())
        for c in cands:
            c["name"] = identifying(c.get("name"))
        sized = [c for c in cands if c.get("total_length")]
        if sized:  # length within tolerance, then same name, then closest; BUSCO 'specific' over 'generic'
            diff = lambda c: abs(c["total_length"] - asm_len)
            best = min(sized, key=lambda c: (diff(c) > MATCH_TOL * asm_len, c.get("name") != asm_name,
                                             diff(c), "specific" not in c["file"]))
            best["matches_assembly"] = diff(best) <= MATCH_TOL * asm_len
        else:  # no length (legacy CheckM, old BUSCO): only a matching sample name verifies it
            # identity: Bracken over raw Kraken2 for the same sample (EU EURL 2025 uses Kraken2 + Bracken)
            best = min(cands, key=lambda c: (c.get("name") != asm_name, c.get("tool") != "Bracken"))
            best["matches_assembly"] = True if asm_name and best.get("name") == asm_name else None
        best["n_candidates"] = len(cands)
        found[key] = best
    return found


# ---------------- species ----------------

def _binomial(name):
    return " ".join(name.replace("'", "").replace("[", "").replace("]", "").split()[:2]).lower()


def species_range(name):
    """Exact name, then NCBI synonym, then genus+species only (drops subsp./strain),
    then 'E. coli'-style abbreviation if unambiguous. Returns (record, how matched)."""
    rows = []
    for line in open(SPECIES_TSV):
        if line.startswith(("#", "taxid")):
            continue
        f = line.rstrip("\n").split("\t")
        gc = {"gc_median": float(f[6]), "gc_lo": float(f[7]), "gc_hi": float(f[8]), "gc_n": int(f[9])} if f[6] else {}
        rep_p90 = {"replicons_p90": int(f[11])} if len(f) > 11 and f[11] else {}
        draft = ({"draft_n": int(f[14]), "draft_contigs": int(f[15]), "draft_size_lo": int(f[16]),
                  "draft_size_hi": int(f[17])} if len(f) > 17 and f[15] else {})
        rows.append(({"species": f[1], "taxid": f[0], "min": int(f[2]), "max": int(f[3]),
                      "expected": int(f[4]), "n_genomes": int(f[5]), **gc, **rep_p90, **draft,
                      "lineage": set(f[13].split(";")) if len(f) > 13 and f[13] else set()},
                     [a for a in f[10].split(";") if a] if len(f) > 10 else []))
    q = name.strip().lower()
    for rec, _ in rows:
        if rec["species"].lower() == q:
            return rec, "exact"
    b = _binomial(name)
    for rec, _ in rows:
        if _binomial(rec["species"]) == b:
            return rec, f"genus+species of '{name}'"
    for rec, aliases in rows:
        if any(_binomial(a) == b for a in aliases):
            return rec, f"NCBI synonym '{name}'"
    m = re.match(r"^([A-Za-z])\.\s*(\S+)$", name.strip())
    if m:
        hits = [rec for rec, _ in rows if rec["species"].lower().startswith(m[1].lower())
                and _binomial(rec["species"]).split()[1:] == [m[2].lower()]]
        if len(hits) == 1:
            return hits[0], f"abbreviation '{name}'"
        if hits:
            return None, "ambiguous abbreviation: " + ", ".join(h["species"] for h in hits[:5])
    return None, "not in NCBI size table"


def species_suggestions(name):
    """Closest table names for a likely typo (#15). Only suggested, never applied: on 12,503 real species
    absent from the table, auto-matching picked a wrong one 19-93 times (Cyclobacterium marinum ->
    Mycobacterium marinum), while 1-letter typos of table names were found 96-97% of the time."""
    names = {_binomial(l.split("\t")[1]): l.split("\t")[1] for l in open(SPECIES_TSV) if not l.startswith(("#", "taxid"))}
    return [names[m] for m in difflib.get_close_matches(_binomial(name), names, n=3, cutoff=0.85)]


@functools.lru_cache(maxsize=1)
def lineage_vocab():
    """Every taxon name in the table's lineages: a placement outside it can't be judged."""
    v = set()
    for line in open(SPECIES_TSV):
        f = line.rstrip("\n").split("\t")
        if len(f) > 13 and not line.startswith(("#", "taxid")):
            v.update(f[13].split(";"))
    return v


def placement_name(s):
    """BUSCO 'enterobacterales_odb10' / 'enterobacterales_odb12.2' (BUSCO 6) / 'actinobacteria_class_odb10', CheckM 'f__Enterobacteriaceae (UID5124)'
    -> NCBI-style lineage name as stored in the table (lower, '/'->'-', ' '->'_')."""
    if not s:
        return None
    s = re.sub(r"_odb\d+(\.\d+)*$", "", s.split(" (")[0].split("__")[-1].strip().lower())   # odb10, odb12, odb12.2
    s = re.sub(r"_(class|phylum)$", "", s)
    return s.replace("/", "-").replace(" ", "_")


def second_limit(species):
    """EU limit (%) for reads from a second species: 2% for Listeria monocytogenes and CPS, else 5%."""
    return T["id_second_share_species"].get(" ".join((species or "").split()[:2]), T["id_second_share"])


def dcs_suspect(idr, reads):
    """% of long reads that are lambda DCS when the classifier's second species is one that Kraken2/Bracken
    assign DCS reads to (so the share may be mostly DCS); else 0. Both inputs are Takivisor's own DCS check
    and the identity report's species name."""
    dcs = max((r.get("dcs_pct", 0) for r in (reads or {}).values()), default=0)
    ok = (idr and dcs >= T["dcs_reads_pct"] and idr.get("second") in T["dcs_miscalled_as"]
          and idr.get("tool") in ("Kraken2", "Bracken"))
    return dcs if ok else 0


def close_relatives(a, b):
    """Same genus or Escherichia/Shigella: classifier splits, not contamination."""
    g1, g2 = a.split()[0], b.split()[0]
    return g1 == g2 or {g1, g2} == {"Escherichia", "Shigella"}


def same_species(name, sp):
    """Compare via the table (so NCBI synonyms match), else by genus+species with GTDB suffixes stripped."""
    rec, _ = species_range(name)
    if rec:
        return rec["taxid"] == sp["taxid"]
    strip = lambda x: _binomial(re.sub(r"_[A-Z]+\b", "", x))
    return strip(name) == strip(sp["species"])


# ---------------- checks ----------------
# Cut-offs live in T at the top. Rule logic stays in code: the indel override, Unicycler's mode and
# report matching depend on context a flat rules table can't express.

def run_checks(asm, reads, sp, rep=None):
    rep = rep or {}
    reads = reads or {}
    out = []
    add = lambda sev, id_, msg, action: out.append(
        {"severity": sev, "id": id_, "message": msg, "action": action})
    if asm.get("dup_names"):
        d = asm["dup_names"]
        add("WARN", "duplicate_contig_names",
            f"{len(d)} contig name(s) used more than once (e.g. {', '.join(d[:3])}).",
            "Check the file holds one assembly (two concatenated assemblies repeat names); per-contig reports "
            "can't be matched to contigs while names repeat.")
    has_long =bool(reads.keys() & LONG_TECHS)
    mode = "long" if asm["assembler"] in LONG_READ_ASM else "short" if asm["assembler"] in SHORT_READ_ASM else None
    if asm["assembler"] == "unicycler":
        # Illumina-only Unicycler is a SPAdes assembly (real Kpne run: 148 contigs), so it only
        # counts as long-read when long reads are present, or, without reads, when it's near-closed.
        mode = ("long" if has_long else "short") if reads else (
            "long" if asm["n_contigs"] <= T["unicycler_long_max_contigs"] else "short")
    if mode is None:
        # Renamed contigs (user, NCBI): infer from the reads. With none, apply the short-read limits,
        # which are the most lenient (failing them is bad for any assembler: real GenBank drafts with
        # 3,567 contigs / N50 1.6 kb passed when they were skipped), and say the stricter long-read
        # checks (contig count, circularity) were not run.
        mode = ("long" if has_long else "short") if reads else "short"
        add("INFO", "assembler_unknown",
            "Assembler not recognised from contig names; "
            + (f"judged as a {mode}-read assembly from the reads given." if reads else
               "only the general N50/contig-count limits were applied; long-read checks (contig count, "
               "circularity) were skipped."),
            "Pass --assembler (e.g. flye, spades) for assembler-specific checks.")
    long_asm = mode == "long"
    total = asm["total_length"]

    if sp:
        if not sp["min"] <= total <= sp["max"]:
            off = (sp["min"] - total) / sp["min"] if total < sp["min"] else (total - sp["max"]) / sp["max"]
            add("WARN" if off <= T["size_warn_frac"] else "FAIL", "size_range",
                f"Assembly {total/1e6:.2f} Mb outside NCBI range for {sp['species']} "
                f"({sp['min']/1e6:.2f}-{sp['max']/1e6:.2f} Mb), {off:.0%} beyond it.",
                ("Check for contamination or a mixed sample, and for duplicated contig ends (too big)."
                 if total > sp["max"] else "Check the species label and whether the assembly is complete (too small)."))
        elif "draft_size_lo" in sp and not sp["draft_size_lo"] <= total <= sp["draft_size_hi"]:
            add("WARN", "size_range",
                f"Assembly {total/1e6:.2f} Mb outside 99.8% of {sp['draft_n']:,} GenBank {sp['species']} drafts "
                f"({sp['draft_size_lo']/1e6:.2f}-{sp['draft_size_hi']/1e6:.2f} Mb).",
                ("Check for contamination or a mixed sample, and for duplicated contig ends (too big)."
                 if total > sp["draft_size_hi"] else "Check the species label and whether the assembly is complete (too small)."))
        # NCBI gc_percent is rounded to 0.5, so allow 1 point of slack
        if "gc_lo" in sp and not sp["gc_lo"] - T["gc_slack"] <= asm["gc_pct"] <= sp["gc_hi"] + T["gc_slack"]:
            # A range from few genomes is weak evidence: GenBank held-out test (2026-09-30) found
            # small species whose real genomes fall outside it, so it can warn but not fail.
            thin = sp.get("gc_n", 0) < T["gc_min_genomes_fail"]
            add("WARN" if thin else "FAIL", "gc_range",
                f"GC {asm['gc_pct']}% vs {sp['species']} RefSeq {sp['gc_lo']}-{sp['gc_hi']}% "
                f"(median {sp['gc_median']}" + (f"; range from only {sp.get('gc_n', 0)} genomes)." if thin else ")."),
                "Wrong species label or a contaminated/mixed sample; confirm species (Kraken2/sylph/ANI).")
    if asm["n_pct"] > T["n_pct"]:
        add("WARN", "n_content", f"{asm['n_pct']}% N bases (>{T['n_pct']}%, EnteroBase).",
            "Use contigs instead of scaffolds, or reassemble.")

    if long_asm:
        # Complete genomes can have many replicons (Borrelia: up to 19), so the limit follows the
        # species' RefSeq Complete Genome replicon count (p90) + 3 slack, never below 10.
        limit = (max(T["long_contigs_min"], sp["replicons_p90"] + T["replicon_slack"])
                 if sp and "replicons_p90" in sp else T["long_contigs_min"])
        if asm["n_contigs"] > limit:
            add("WARN", "long_read_fragmented",
                f"{asm['n_contigs']} contigs from a long-read assembler (expected <= {limit}"
                + (f"; complete {sp['species']} genomes have up to ~{sp['replicons_p90']} replicons)." if limit > T["long_contigs_min"] else ")."),
                "Check read N50 and coverage; low-coverage long reads fragment the assembly.")
        circ = [c for c in asm["contigs"] if c["circular"] is not None]
        ref = sp["expected"] if sp else total
        # "chromosome-sized" scales down for small genomes (Mycoplasma 0.58 Mb: a closed circle is not a problem)
        chrom_min = min(T["circ_min_len"], T["near_genome_frac"] * ref)
        if circ and not any(c["circular"] and c["length"] > chrom_min for c in circ):
            # Raven marked a complete 5.24 Mb CFT073 chromosome XO:i:0, so a genome-sized
            # contig without the flag is only a note, not a warning.
            near = asm["largest"] >= T["near_genome_frac"] * ref
            add("INFO" if near else "WARN", "no_circular_chromosome",
                f"No chromosome-sized contig (>{chrom_min/1e6:.2g} Mb) marked circular by {asm['assembler']}"
                + (f" (largest {asm['largest']/1e6:.2f} Mb is near genome size)." if near else "."),
                "Probably complete; if you need it closed, check overlapping ends and rotate (dnaapler)."
                if near else "Chromosome not closed: more coverage or longer reads.")
        # A chromosome-sized circle that is far smaller than the genome, next to other large linear
        # pieces (raw Flye on E. coli: 2.17 Mb circle among 12 contigs, 5.3 Mb) is suspect.
        big_linear = [c for c in circ if not c["circular"] and c["length"] >= T["big_linear_len"]]
        for c in circ:
            if c["circular"] and T["circ_min_len"] <= c["length"] < T["suspect_circ_frac"] * ref and big_linear:
                add("INFO", "suspect_circular",
                    f"Contig {c['name']} ({c['length']/1e6:.2f} Mb) is circular but the genome is "
                    f"~{ref/1e6:.2f} Mb and {len(big_linear)} other contigs >=100 kb are linear.",
                    "Could be a false circularisation or collapsed repeat; secondary chromosomes are real in "
                    "some genera (Vibrio, Burkholderia). Check against a reference or with more coverage.")
    elif mode == "short":
        if sp and "draft_contigs" in sp:
            n, limit = asm.get("n_contigs_200", asm["n_contigs"]), sp["draft_contigs"]
            why = f"contigs >= 200 bp; more than 99.5% of {sp['draft_n']:,} GenBank {sp['species']} drafts"
        else:
            n, limit, why = asm["n_contigs"], T["short_contigs_max"], "EnteroBase"
        if n > limit:
            add("WARN", "short_read_fragmented", f"{n} contigs (>{limit}, {why}).",
                "Check depth and read quality, and look for contamination (too many contigs often means one of these).")
        if asm["n50"] < T["short_n50_min"]:
            add("WARN", "low_n50", f"N50 {asm['n50']:,} bp (<{T['short_n50_min'] // 1000} kb, EnteroBase).",
                "Check coverage and read quality.")

    # Lambda-like sequence in the assembly is only DCS contamination if the reads contain DCS: 29 of 51
    # pre-2014 enterobacterial genomes (before DCS existed) carry 250-4,515 bp of it as prophage.
    lam = asm.get("dcs_contigs")
    long_sets = [r for t_, r in reads.items() if t_ != "illumina" and "dcs_pct" in r]
    if lam:
        where = ", ".join(f"{c['name']} (~{c['dcs_bp']:,} bp)" for c in lam)
        if any(r["dcs_pct"] >= T["dcs_reads_pct"] for r in long_sets):
            add("WARN", "dcs_in_assembly",
                f"ONT lambda control (DCS) was assembled into the genome: {where}; the reads contain DCS.",
                "Remove DCS reads (NanoLyse, or chopper/minimap2 against lambda) and reassemble.")
        elif not long_sets:
            add("INFO", "lambda_like_sequence",
                f"Lambda-like sequence in {where}: either ONT control (DCS) contamination or a lambda-like "
                "prophage (common in E. coli).",
                "Pass the ONT reads to tell them apart, or check whether DCS reads were filtered before assembly.")
    reps = [c for c in asm["contigs"] if c.get("repeat")]
    if reps:
        add("INFO", "canu_repeat_contigs",
            f"{len(reps)} contig(s), {sum(c['length'] for c in reps)/1000:.0f} kb, flagged suggestRepeat=yes by Canu.",
            "Usually collapsed repeats; they inflate the contig count and are commonly dropped before submission.")

    covs = [(c["cov"], c["length"]) for c in asm["contigs"] if c["cov"] is not None]
    if covs:
        # length-weighted median contig coverage
        run, half, med = 0, sum(l for _, l in covs) / 2, 0
        for cv, l in sorted(covs):
            run += l
            if run >= half:
                med = cv
                break
        low = sum(l for cv, l in covs if cv < T["lowcov_rel"] * med)
        asm["median_contig_cov"] = med
        asm["low_cov_bp"] = low
        if low > T["lowcov_frac"] * total:
            add("WARN", "low_cov_contigs",
                f"{low/total*100:.1f}% of the assembly is in contigs <{T['lowcov_rel']:.0%} of median coverage ({med:.1f}).",
                "Likely contamination or junk: inspect/remove low-coverage contigs.")

    # An incomplete assembly understates genome size and so inflates depth (ONT 10x read as 13x on
    # a 4.1 Mb Raven assembly), so divide by the species' expected size when it is larger.
    genome = max(total, sp["expected"]) if sp else total
    for tech, r in reads.items():
        depth = r["bases"] / genome
        r["est_depth"] = round(depth, 1)
        if r.get("qual_simplified"):
            add("INFO", "quality_simplified",
                f"{tech} quality scores use only placeholder values (e.g. NCBI SRA-lite download), so Q20/Q30 "
                "and mean quality are not real and quality checks were skipped.",
                "Download full-quality reads (ENA, or SRA normalized format) if read quality matters.")
        # Reads from another organism (wrong file, mixed sample, host DNA) shift read GC away
        # from the assembly's. 6 correctly paired real sets (Illumina + ONT, 2 species) differed by
        # <=0.54 points; Kpne assembly vs E. coli ONT reads differed by 6.5. ponytail: 2-point cut-off
        # only catches whole-file mix-ups or large contaminant fractions (5% of reads 10 points off
        # moves read GC by 0.5); per-read GC bimodality would catch smaller ones.
        if r.get("dcs_pct", 0) >= T["dcs_reads_pct"]:
            add("INFO", "dcs_reads", f"{r['dcs_pct']}% of {tech} reads (>=1 kb) are ONT lambda control (DCS).",
                "Filter them before assembly (NanoLyse, or chopper/minimap2 against lambda); they can be "
                "assembled into the genome.")
        rng, rg = asm.get("gc_window_range"), r.get("read_gc")
        if rng and rg:
            n = sum(k for (c, g), k in rg.items() if c in rng)
            n_out = sum(k for (c, g), k in rg.items() if c in rng and not rng[c][0] <= g <= rng[c][1])
            if n >= T["read_gc_min_reads"]:
                r["gc_outlier_pct"] = round(100 * n_out / n, 2)
                kind = "illumina" if tech == "illumina" else "long"
                if r["gc_outlier_pct"] > T["read_gc_outlier_pct"][kind]:
                    add("WARN", "read_gc_outliers",
                        f"{r['gc_outlier_pct']}% of {tech} reads have a GC outside the 99% range of same-length "
                        f"stretches of the assembly (clean samples: ~{T['read_gc_outlier_clean'][kind]}).",
                        "Some reads look like another organism (contamination, host DNA or wrong file) that the "
                        "assembly doesn't contain; run Kraken2/sylph on the reads.")
        if "gc_pct" in r and abs(r["gc_pct"] - asm["gc_pct"]) > T["read_gc_diff"]:
            add("WARN", "read_gc_mismatch",
                f"{tech} reads GC {r['gc_pct']}% vs assembly GC {asm['gc_pct']}%.",
                "These reads may not belong to this assembly, or carry a contaminant/host fraction; "
                "check file pairing and run Kraken2/sylph on the reads.")
        need = T["depth_min_trycycler"] if asm["assembler"] == "trycycler" and tech != "illumina" else T["depth_min"]
        # Other bodies' minimums for this species/tech (EU limit only; Trycycler's 100x is its own rule).
        name = sp["species"] if sp else ""
        other = [(m, who) for scope, t, m, who in T["depth_other"] if need == T["depth_min"]
                 and t in (None, tech) and (scope is None or name == scope or name.startswith(scope + " "))]
        if depth < need:
            looser = [f"{who} accepts from {m}x" for m, who in other if m <= depth < need]
            add("FAIL" if depth < need / 2 else "WARN", "low_depth",
                f"Estimated {tech} depth {depth:.1f}x (min {need}x for {asm['assembler']})."
                + (f" Below EU {need}x; " + "; ".join(looser) + "." if looser else ""),
                "Sequence more.")
        for m, who in other:
            # Pasteur's ONT rule is for nanopore-only assemblies; with Illumina reads too, skip it.
            # ponytail: any Illumina depth counts; check it against Pasteur's 40x if hybrids get thin Illumina
            if need <= depth < m and not (tech == "ont" and "illumina" in reads):
                add("INFO", "depth_below_agency",
                    f"{TECH_NAME.get(tech, tech)} depth {depth:.1f}x meets the EU {need}x but is below the "
                    f"{m}x asked by {who}.",
                    f"Top up to >= {m}x before submitting to {who.split(' (')[0]}.")
        if tech == "illumina":
            q30_min = next(m for L, m in T["q30_min_by_len"] if r.get("max_len", 150) <= L)
            if not r.get("qual_simplified") and r.get("pct_q30", 100) < q30_min:
                add("WARN", "illumina_q30", f"{r['pct_q30']}% bases >=Q30 (<{q30_min}% for {r.get('max_len', '?')} bp reads).",
                    "Review the run and library quality with the sequencing facility; trimming won't fully fix it.")
            if r.get("pct_passed", 100) < T["fastp_pass_min"]:
                after = r["bases_after"] / genome if "bases_after" in r else None
                ok = after is not None and after >= need
                add("INFO" if ok else "WARN", "fastp_low_pass",
                    f"Only {r['pct_passed']}% of reads passed fastp filtering"
                    + (f"; {after:.1f}x left after filtering." if after is not None else "."),
                    "Enough data remains; no action needed." if ok else
                    "Low-quality run or library; check the fastp report for where reads are lost.")
        if tech == "ont" and "read_n50" in r:
            if r["read_n50"] < T["ont_n50_fail"]:
                add("FAIL", "ont_short_reads", f"Read N50 {r['read_n50']} bp (<1 kb).",
                    "Re-extract DNA for longer fragments (gentle extraction, size selection); with reads this short "
                    "long-read assemblers fail or fragment (Flye needs >= 1 kb overlaps).")
            elif r["read_n50"] < T["ont_n50_warn"]:
                add("WARN", "ont_short_reads", f"Read N50 {r['read_n50']:,} bp (<5 kb).",
                    "Aim for read N50 >= 5 kb (gentle DNA extraction, size selection) so repeats such as rRNA "
                    "operons (~5-7 kb) can be resolved.")

    techs = set(reads)
    if techs == {"illumina"} and long_asm:
        if ILLUMINA_POLISHERS & set(asm["polishing_hints"]):  # long-read assembly polished with Illumina
            add("INFO", "long_reads_not_given", f"{asm['assembler']} assembly but only Illumina reads/reports given.",
                "Long-read input was not assessed; pass the ONT/HiFi reads or NanoPlot report too.")
        else:
            add("WARN", "tech_assembler_mismatch", f"Illumina reads with {asm['assembler']} (long-read assembler).",
                "Check the right files were given.")
    if techs and techs <= LONG_TECHS and asm["assembler"] in SHORT_READ_ASM:
        add("WARN", "tech_assembler_mismatch", f"Long reads with {asm['assembler']} (short-read assembler).",
            "Check the right files were given.")

    # Reports that don't match this assembly are shown but not interpreted.
    for key in ("busco", "checkm", "quast"):
        r = rep.get(key)
        if r and r.get("matches_assembly") is False:
            add("WARN", "report_mismatch",
                f"{key} report is for a {r['total_length']:,} bp assembly, this one is {total:,} bp; ignored.",
                f"Point --reports at the outputs for this assembly ({r['file']} is the closest match).")
    usable = lambda k: rep.get(k) if rep.get(k) and rep[k].get("matches_assembly") is not False else None
    busco, ck, quast = usable("busco"), usable("checkm"), usable("quast")
    # A report that can't be tied to this assembly (no length, name differs) is still read, but it
    # can't fail the assembly on its own and every finding from it says so.
    unverified = {k for k, r in (("busco", busco), ("checkm", ck), ("quast", quast))
                  if r and r.get("matches_assembly") is None}
    for k in sorted(unverified):
        add("INFO", "report_unverified",
            f"{k} report has no genome length and no identifying sample name matching this assembly "
            f"(generic names like final_assembly don't count), so it may belong to another assembly.",
            f"Check {rep[k]['file']} is for this assembly (or pass --sample NAME); its findings are capped at WARN.")
    unv = lambda k: " [report not verified for this assembly]" if k in unverified else ""
    # Species identity (#22): an identity report from the pipeline (Kraken2 / sylph / GTDB-Tk) catches
    # close-relative mislabels that size/GC can't (7% of NCBI's ANI-failed genomes were GC-visible).
    idr = rep.get("identity") if rep.get("identity") and rep["identity"].get("matches_assembly") is not False else None
    if idr:
        tag = "" if idr.get("matches_assembly") else " [report not verified for this sample]"
        if idr.get("n_assigned") is not None and idr["n_assigned"] < T["id_min_reads"]:
            add("INFO", "identity_unclear",
                f"{idr['tool']} report has only {idr['n_assigned']} bacterial species-level reads{tag}; too few "
                "to confirm the species.", "Not an isolate-style read set (host-dominated or low biomass?).")
        elif idr["top_share"] < T["id_min_share"]:
            add("INFO", "identity_unclear",
                f"{idr['tool']} report has no dominant species (top: {idr['top']}, {idr['top_share']}% of "
                f"species-level assignments){tag}.", "Can't confirm the species; check for a mixed sample.")
        elif sp and same_species(idr["top"], sp) is False:
            sp["identity_conflict"] = f"{idr['tool']} says {idr['top']}"
            add("WARN", "species_mismatch",
                f"{idr['tool']} report says {idr['top']} ({idr['top_share']}% of bacterial species-level "
                f"assignments), not {sp['species']}{tag}.",
                "Check the sample label; if the label is right, the sample may be swapped or contaminated.")
        elif sp:
            sp["verified_by"] = f"{idr['tool']} ({idr['top_share']}%)"
        lim = second_limit(sp["species"] if sp else idr["top"])
        if idr["second"] and idr["second_share"] >= lim and idr.get("n_assigned", T["id_min_reads"]) >= T["id_min_reads"]:
            close = close_relatives(idr["top"], idr["second"])
            sev = ("INFO" if close else "FAIL" if idr["second_share"] >= T["id_second_fail"] and not tag
                   else "WARN")
            dcs = dcs_suspect(idr, reads)
            dcs_note = bool(dcs)
            add(sev, "identity_mixed",
                f"{idr['tool']} report: second species {idr['second']} at {idr['second_share']}%{tag}."
                + (f" {dcs:g}% of the long reads are ONT lambda DCS control, which Kraken2 assigns to "
                   f"{' or '.join(T['dcs_miscalled_as'])}: this share may be mostly DCS, not contamination." if dcs_note else ""),
                "Classifiers often split reads between close relatives; usually not contamination." if close else
                "Filter the lambda DCS reads (NanoLyse, or chopper/minimap2 against lambda) and re-run the classifier "
                f"before acting; if {idr['second']} stays above {lim:g}%, re-streak a single colony and re-sequence." if dcs_note else
                "Mixed culture or contamination: re-streak a single colony, re-extract and re-sequence "
                f"(EU limit: < {lim:g}% reads from another species).")
            if dcs_note and not close:   # this finding's own lab/bioinformatics actions replace ACTIONS
                out[-1]["acts"] = ("None yet: the second species may be ONT's lambda DCS control spike-in; re-streak "
                                   "only if it persists after DCS filtering.",
                                   "Remove DCS reads (NanoLyse or chopper) and re-run Kraken2 and Bracken; act on the new share.")
    # Contig identity (#25): close relatives with similar GC are invisible to size/GC and CheckM2 misses
    # ~5% mixes; a per-contig classification names the foreign contigs. Takivisor never runs the
    # classifier itself, it only reads its output, and advises running it when absent.
    ckr = rep.get("contig_kraken")
    cs = ckr and ckr.get("contigs")
    if cs:
        mine = {(c["name"], c["length"]) for c in asm["contigs"]}
        hit = sum((c["name"], c["length"]) in mine for c in cs) / len(cs)
        if hit < 0.9:
            add("WARN", "report_mismatch",
                f"Kraken2 contig output: only {hit:.0%} of its contigs match this assembly by name and length; ignored.",
                f"Point --reports at the Kraken2 output for this assembly ({ckr['file']}).")
            cs = None
    if cs:
        called = [c for c in cs if c["votes"] >= T["contig_min_votes"] and c["share"] >= T["contig_min_share"]]
        if sp:
            ref = sp["species"]
            other = [c for c in called if same_species(c["top"], sp) is False]
        else:  # no species: the species holding most classified bp is the reference
            bp = Counter()
            for c in called:
                bp[c["top"]] += c["length"]
            ref = bp.most_common(1)[0][0] if bp else None
            other = [c for c in called if c["top"] != ref]
        cov = 100 * sum(c["length"] for c in cs) / total
        foreign = (lambda name: same_species(name, sp) is False) if sp else (lambda name: name != ref)
        # Coverage cross-check (point 1): contaminants sit far below the main coverage, shared plasmids don't.
        med, ccov = asm.get("median_contig_cov"), {c["name"]: c["cov"] for c in asm["contigs"]}
        rel = lambda c: ccov.get(c["name"]) / med if med and ccov.get(c["name"]) is not None else None
        plasmid_like = [c for c in other if rel(c) is not None and rel(c) >= T["contig_lowcov_rel"]]
        fbp = sum(c["length"] for c in other)
        sbp = fbp - sum(c["length"] for c in plasmid_like)   # low or unknown coverage
        if other:
            pct = 100 * fbp / total
            mostly = pct >= 50   # most of the assembly is another species: a label problem, coverage can't tell
            sev = "WARN" if mostly or 100 * sbp / total >= T["contig_foreign_pct"] else "INFO"
            other.sort(key=lambda c: -c["length"])
            add(sev, "contig_contamination",
                f"{len(other)} contig(s), {fbp:,} bp ({pct:.2f}% of the assembly), classified as another species "
                f"than {ref} by Kraken2 k-mer votes: "
                + ", ".join(f"{c['name']} ({c['length']:,} bp) {c['top']} {c['share']}%"
                            + (f", {rel(c):.2f}x median cov" if rel(c) is not None else "") for c in other[:5])
                + (" ..." if len(other) > 5 else "") + "."
                + (f" {len(plasmid_like)} of them at >= {T['contig_lowcov_rel']}x median coverage, as expected for a "
                   "plasmid or mobile element shared with a relative." if plasmid_like and not mostly else ""),
                "Check the sample label; the assembly is mostly another species." if mostly else
                "Contamination or mixed sample if at low coverage; plasmids and mobile elements shared across "
                "species can also vote for a relative. Inspect these contigs (BLAST) before removing them.")
        elif called:
            add("INFO", "contig_identity",
                f"Kraken2 contig output: {len(called)} of {len(cs)} contigs (covering {cov:.0f}% of the assembly) "
                f"classified, all as {ref}.", "No contig from another species found.")
        # Pooled votes (point 2): contaminant contigs are often too small to reach contig_min_votes alone
        # (~0.4% of k-mers voted at species level with the 8 GB DB), but together they add up.
        flagged = {c["name"] for c in other}
        pool = Counter()
        for c in cs:
            if c["name"] not in flagged:
                pool.update(c.get("sp_votes", {}))
        n_pool = sum(pool.values())
        second = [(k, v) for k, v in pool.most_common() if foreign(k)][:1] if ref else []
        if second and n_pool and second[0][1] >= T["contig_min_votes"] \
                and 100 * second[0][1] / n_pool >= T["contig_pool_share"]:
            k, v = second[0]
            add("INFO", "contig_second_species",
                f"Kraken2 votes pooled over the contigs not flagged on their own: {k} has {100 * v / n_pool:.1f}% "
                f"({v:,} of {n_pool:,} species-level votes).",
                "Possible low-level contamination spread over small contigs; check with CheckM2 or the read-level "
                "identity report (second species).")
    elif ckr and ckr.get("contigs") is None:
        add("INFO", "contig_identity_unchecked",
            "Kraken2 contig output found but no ktaxonomy.tsv, so contigs could not be resolved to species.",
            "Add the Kraken2 database folder (it holds ktaxonomy.tsv) to --reports.")
    elif not ckr:
        add("INFO", "contig_identity_unchecked",
            "No per-contig classification; contaminant contigs from a species with similar GC can be missed "
            "(CheckM2 misses ~5% mixes).",
            "Run `kraken2 --db DB --output contigs.kraken.output ASSEMBLY.fasta` and pass the output and the "
            "DB folder (for ktaxonomy.tsv) via --reports.")
    # Coarse placement (#22 B): BUSCO auto-lineage / legacy CheckM marker lineage vs the species' lineage.
    if sp and sp.get("lineage"):
        for tool, place in (("BUSCO", busco and busco.get("lineage")), ("CheckM", ck and ck.get("marker_lineage"))):
            nm = placement_name(place)
            if nm and nm in lineage_vocab() and nm not in sp["lineage"]:
                add("WARN", "taxonomy_mismatch",
                    f"{tool} placed the genome in {nm} ({place}), which doesn't include {sp['species']}.",
                    "Wrong species label, or a genome from another group; check the sample and an identity report.")
    # ConFindr (same-species mixing, which Kraken2/CheckM2 miss). Sample-named rows: a row from another sample
    # is only used when it is the report's only row. Pasteur BIGSdb rejects >= 5% contamination.
    cf = confindr_row(rep)
    if rep.get("confindr") and not cf:
        add("INFO", "report_unverified",
            f"ConFindr report has {rep['confindr']['n_candidates']} samples and none is named like this assembly; "
            "not used.", f"Pass --sample NAME with the ConFindr sample name (its file name before _R1) "
            f"({rep['confindr']['file']}).")
    if cf and cf["contaminated"]:
        tag = "" if cf.get("matches_assembly") else " [report not verified for this sample]"
        several = ":" in (cf.get("genus") or "")
        add("WARN", "confindr_contamination",
            f"ConFindr found contamination{tag}: {cf['snvs'] if cf['snvs'] is not None else '?'} contaminating SNVs"
            + (f", ~{cf['pct']:g}% of reads" + (" (above Pasteur BIGSdb's 5%)" if cf["pct"] >= 5 else "")
               if cf["pct"] is not None else "")
            + (f"; several genera ({cf['genus']})" if several else "; same genus, so likely a mix of strains") + ".",
            "Re-streak a single colony and re-sequence; a mix of strains of one species is invisible to Kraken2 "
            "and CheckM2.")
    # Depth breadth (#32, mosdepth): a mean depth can hide regions with too few reads to call.
    md = mosdepth_row(rep, asm)
    if rep.get("mosdepth") and not md:
        add("WARN", "report_mismatch",
            f"mosdepth report is for a {rep['mosdepth']['total_length']:,} bp sequence that doesn't match this "
            f"assembly; ignored.", f"Point --reports at the mosdepth output for this assembly ({rep['mosdepth']['file']}).")
    elif md and not md[0]["at"]:
        add("INFO", "breadth_unchecked", f"mosdepth summary found but no {md[0]['name']}.mosdepth.global.dist.txt "
            "next to it, so depth breadth could not be read.", "Keep the global.dist.txt file mosdepth writes.")
    elif md:
        r, against = md
        x, need = T["breadth_min"]
        if r["at"][x] < need:
            add("WARN", "low_breadth",
                f"mosdepth: {r['at'][x]:g}% of {against} has >= {x}x (EU guideline >= {need}%; mean {r['mean']:g}x).",
                "Parts of the genome have too few reads to call variants or genes there; sequence more, or check "
                "for regions the reads don't cover (plasmid lost, wrong reference).")
        name = sp["species"] if sp else ""
        for scope, x2, need2, who in T["breadth_other"]:
            if (name == scope or name.startswith(scope + " ")) and r["at"].get(x2, 0) < need2:
                add("INFO", "depth_below_agency",
                    f"mosdepth: {r['at'].get(x2, 0):g}% of {against} has >= {x2}x; {who} asks for {need2}%.",
                    f"Top up the reads before using this genome for {who}.")
    # AMR gene calls (NDARO AMRFinderPlus / CARD RGI): incomplete calls tell which genes the assembly cut or broke.
    for key in ("amrfinder", "rgi"):
        r = rep.get(key)
        if not r:
            continue
        hits = amr_hits(r, asm)
        if hits is None:
            add("WARN", "report_mismatch", f"{r['tool']} report names contigs that are not in this assembly; ignored.",
                f"Point --reports at the {r['tool']} output for this assembly ({r['file']}).")
            continue
        bad = Counter(h["issue"] for h in hits if h["issue"] and h["issue"] != "wrap")
        wrap = sorted({h["gene"] for h in hits if h["issue"] == "wrap"})
        if bad:
            words = {"edge": "cut by a contig end", "partial": "partial (< 90% of the reference)",
                     "stop": "with an internal stop codon", "loose": "Loose (discovery-level) hits"}
            genes = sorted({h["gene"] for h in hits if h["issue"] and h["issue"] != "wrap"})
            # #30 real data (E. coli CFT073, Raven unpolished, 22k indels): 9 PARTIAL + 1 INTERNAL_STOP inside one
            # contig, and 30 of 44 genes not called at all; polished assemblies: 0 incomplete, 0 missing.
            acts = (["genes at contig ends may be split or missed: more depth or long reads closes them"]
                    if bad["edge"] else []) \
                + (["genes cut short or stopped inside a contig are often indel errors in nanopore assemblies, which "
                    "also hide other genes completely: polish (Medaka) or basecall with sup; in a polished genome "
                    "they can be real truncated genes"] if bad["partial"] or bad["stop"] else []) \
                + (["CARD advises Loose hits for discovery, not surveillance"] if bad["loose"] else [])
            add("INFO", "amr_calls_incomplete",
                f"{r['tool']}: {sum(bad.values())} of {len(hits)} gene calls are incomplete: "
                + ", ".join(f"{n} {words[k]}" for k, n in bad.items())
                + f" ({', '.join(genes[:5])}{', ...' if len(genes) > 5 else ''}).",
                "Before reporting these genes as present, partial or absent: " + "; ".join(acts) + ".")
        if wrap:
            add("INFO", "amr_calls_incomplete",
                f"{r['tool']}: {len(wrap)} gene(s) span the start/end of a circular contig ({', '.join(wrap[:5])}"
                f"{', ...' if len(wrap) > 5 else ''}); reported as cut by a contig end but complete in the genome.",
                "Rotate the circular contigs (e.g. dnaapler) and re-run the gene caller to see them whole.")
    if not ck:
        # Sensitivity tests: size + GC caught >= 10-20% of a different-GC contaminant, never Salmonella
        # below 50%; BUSCO duplication caught 20% Salmonella; read-GC outliers caught 3% of a real
        # different-GC contaminant in ONT reads; CheckM2 caught 11/12 mixes down to 5%. The note lists
        # what actually ran instead of assuming nothing did (#3 review).
        ran = (["genome size/GC vs species"] if sp else []) + (["BUSCO duplication"] if busco else []) \
            + (["identity report (second species)"] if idr else []) \
            + (["Kraken2 contig votes"] if cs else []) + (["ConFindr (same-species mixing)"] if cf else []) \
            + [f"{t_} read-GC outliers" for t_, r in reads.items() if "gc_outlier_pct" in r] \
            + (["low-coverage contigs"] if covs else [])
        add("INFO", "contamination_unchecked",
            "No CheckM2 report; contamination checked only by "
            + (", ".join(ran) if ran else "nothing (no species, reports or reads given)")
            + ". Low-level contamination (below ~10-20%, or ~3% for different-GC ONT reads) and close "
              "relatives with similar GC can still be missed.",
            "Run CheckM2 and pass its quality_report.tsv via --reports.")
    cap = lambda sev, k: "WARN" if sev == "FAIL" and k in unverified else sev

    # Frameshift (indel) errors split genes: BUSCOs turn 'fragmented' and genes get short.
    # Watson & Warr 2019 Nat Biotech 37:124. Sweep 2026-09-30: 34 complete RefSeq genomes (GC 25-72%,
    # 0.58-12 Mb) had gene length 265-376 aa and BUSCO F <=8.9%; 3 unpolished ONT assemblies had
    # 136-164 aa and F 43-53%. Cut-offs sit between the groups with margin on both sides.
    # A fragmented short-read assembly also cuts genes at contig ends, so BUSCO F only
    # counts as an indel sign when contigs are long.
    indel = []
    cut_at_ends = not long_asm and asm["n50"] < T["cut_at_ends_n50"]
    if busco and busco["fragmented"] > T["busco_frag_indel"]:
        if cut_at_ends:
            add("WARN", "busco_fragmented", f"BUSCO fragmented {busco['fragmented']}% with N50 {asm['n50']:,} bp.",
                "Genes cut at contig ends: improve contiguity (coverage, read length) before judging completeness.")
        else:
            indel.append(f"BUSCO fragmented {busco['fragmented']}%")
    if ck and ck["avg_gene_length_aa"] is not None and ck["avg_gene_length_aa"] < T["gene_len_indel"]:
        indel.append(f"CheckM2 mean gene length {ck['avg_gene_length_aa']} aa (typical ~300)")
    # Built-in fallback when CheckM2 isn't given; skipped where contig ends cut genes (like BUSCO F).
    code4 = bool(sp) and sp["species"].split()[0].lower().strip("[]") in CODE4_GENERA
    orf = (asm.get("orf_n50_aa") or {}).get(4 if code4 else 11)
    if orf and not ck and not cut_at_ends:
        asm["orf_n50_used"] = orf
        if orf < T["orf_n50_indel"]:
            indel.append(f"ORF N50 {orf} aa (clean genomes >= {T['orf_n50_clean_min']})"
                         + ("" if sp else "; if this is a Mycoplasma/Spiroplasma (genetic code 4), pass --species"))
    if indel:
        add("WARN", "indel_errors", "Signs of frameshift errors: " + "; ".join(indel) + ".",
            "Polish (Medaka for ONT, Polypolish/Pilon with Illumina) or re-basecall with a more accurate "
            "model; re-check completeness/contamination after that, they are unreliable until then.")
    caveat = " (likely caused by the indel errors above)" if indel else ""

    if busco:
        # BUSCO completeness reflects the organism as much as the assembly: 13 of 34 complete
        # RefSeq genomes scored 57-93% on bacteria_odb10 while CheckM2 gave 97.5-100%.
        # So CheckM2 judges completeness; BUSCO C is informational.
        if busco["complete"] < T["busco_complete"]:
            add("WARN" if indel else "INFO", "busco_low",
                f"BUSCO complete {busco['complete']}% (<{T['busco_complete']}%){caveat}{unv('busco')}.",
                "Fix the sequence errors first (polish), then re-check." if indel else
                "Complete genomes of reduced/divergent lineages (Mycoplasma, Chlamydia, Borrelia...) score "
                "57-93% on bacteria_odb10; " + ("rely on CheckM2 completeness." if ck else
                "run CheckM2 for a completeness verdict."))
        if busco["duplicated"] > T["busco_dup"]:
            add("WARN", "busco_duplicated", f"BUSCO duplicated {busco['duplicated']}% (>{T['busco_dup']}%){unv('busco')}.",
                "Mixed sample (another strain or species: E. coli + 20% Salmonella gave 15%, Kpne + 3% S. aureus "
                "reads 44%), duplicated contigs, or untrimmed overlapping contig ends (Canu/Raven).")
    if ck:
        c, k, tool = ck["completeness"], ck["contamination"], ck["tool"]
        if c is not None and c < T["ck_compl_warn"]:
            add(cap("FAIL" if c < T["ck_compl_fail"] and not indel else "WARN", "checkm"), "checkm_completeness",
                f"{tool} completeness {c}% (<{T['ck_compl_warn']}%, EU EURL isolate limit){caveat}{unv('checkm')}.",
                "Fix the sequence errors first (polish), then re-check." if indel else
                "Check depth and read length; more or longer reads usually complete the assembly.")
        if k is not None and k > T["ck_cont_warn"]:
            add(cap("FAIL" if k > T["ck_cont_fail"] and not indel else "WARN", "checkm"), "checkm_contamination",
                f"{tool} contamination {k}% (>{T['ck_cont_warn']}%, EU EURL isolate limit){caveat}{unv('checkm')}.",
                "Fix the sequence errors first (polish), then re-check." if indel else
                "Mixed sample or contaminant contigs: check species (Kraken2/sylph) and low-coverage contigs.")
    if quast and quast.get("genome_fraction") is not None and quast["genome_fraction"] < T["quast_gf_min"]:
        add("WARN", "quast_genome_fraction", f"QUAST genome fraction {quast['genome_fraction']}% of the reference{unv('quast')}.",
            "Check that the reference is the same strain; if it is, more or longer reads may fill the missing regions.")

    return sorted(out, key=lambda f: SEV_ORDER[f["severity"]])


# ---------------- beginner guide: plain words + "Start here" order ----------------
# id: (stage, step, plain). Stage orders the to-do list: 1 inputs, 2 sequencing data, 3 contamination,
# 4 assembly fixes, 5 usually a consequence of the earlier stages (listed under "re-run").
GUIDE = {
    "species_unknown": (1, "Fix the species name (--species)", "The species name was not recognised, so genome size and GC could not be checked."),
    "species_mismatch": (1, "Confirm which organism this sample really is", "The identity report says this is a different species than the label."),
    "taxonomy_mismatch": (1, "Confirm which organism this sample really is", "BUSCO/CheckM placed the genome in a group that does not contain the labelled species."),
    "gc_range": (1, "Confirm which organism this sample really is", "The genome's GC content does not fit the labelled species: wrong label or a mixed sample."),
    "report_mismatch": (1, "Give the reports that belong to this assembly", "A tool report was made for a different assembly, so it was ignored."),
    "report_unverified": (1, "Give the reports that belong to this assembly", "A report could not be tied to this assembly by name or length, so it counts less."),
    "tech_assembler_mismatch": (1, "Check that the right files were given", "The reads do not fit the assembler used (e.g. long reads with a short-read assembler)."),
    "read_gc_mismatch": (1, "Check that the right files were given", "The reads' GC content differs from the assembly's: they may be from another sample."),
    "long_reads_not_given": (1, "Add the long reads (or NanoPlot report)", "The assembly was made from long reads, but only short-read data was given to check."),
    "low_depth": (2, "Sequence more", "Too few reads for this genome size; assemblies get fragmented and less accurate."),
    "illumina_q30": (2, "Check the Illumina run quality", "Fewer high-quality bases than a typical good Illumina run."),
    "fastp_low_pass": (2, "Check the Illumina run quality", "Many reads were thrown away by quality filtering."),
    "ont_short_reads": (2, "Get longer nanopore reads (DNA extraction, size selection)", "The nanopore reads are short, which limits how complete the assembly can be."),
    "quality_simplified": (2, "Download the reads with original quality scores", "The read qualities are placeholders (NCBI SRA-lite), so quality could not be judged."),
    "dcs_reads": (2, "Remove the lambda control reads before assembly", "Part of the reads is ONT's lambda control DNA, not your sample."),
    "checkm_contamination": (3, "Check for contamination", "CheckM found marker genes from more than one organism."),
    "contig_contamination": (3, "Check for contamination", "Some contigs look like a different species."),
    "contig_second_species": (3, "Check for contamination", "Small contigs together point to a second species."),
    "read_gc_outliers": (3, "Check for contamination", "Some reads have a GC content the assembly does not contain, e.g. another organism."),
    "identity_mixed": (3, "Check for contamination", "The identity report sees a second species in the reads."),
    "busco_duplicated": (3, "Check for contamination", "Many single-copy genes appear twice: a mixed sample or duplicated contigs."),
    "low_cov_contigs": (3, "Check for contamination", "Part of the assembly has far fewer reads than the rest, typical of contamination."),
    "dcs_in_assembly": (3, "Remove the lambda control sequence from the assembly", "ONT lambda control DNA was assembled into the genome."),
    "lambda_like_sequence": (3, "Check for contamination", "Lambda-like sequence was found in the assembly."),
    "indel_errors": (4, "Polish the assembly", "Small insertion/deletion errors break genes; typical of unpolished nanopore assemblies."),
    "n_content": (4, "Use contigs instead of scaffolds", "The assembly has many N's (gaps of unknown sequence)."),
    "busco_low": (5, "", "Many expected core genes were not found complete."),
    "busco_fragmented": (5, "", "Many expected core genes are broken."),
    "checkm_completeness": (5, "", "CheckM estimates part of the genome is missing."),
    "long_read_fragmented": (5, "", "The genome is in more pieces than expected for long reads."),
    "short_read_fragmented": (5, "", "The genome is in very many pieces."),
    "low_n50": (5, "", "Contigs are short, so genes near breaks may be incomplete."),
    "no_circular_chromosome": (5, "", "The chromosome was not closed into a circle (it may still be complete)."),
    "quast_genome_fraction": (5, "", "Part of the reference genome is missing from the assembly."),
    "size_range": (5, "", "The genome is larger or smaller than usual for this species."),
    "suspect_circular": (5, "", "A circular contig is much smaller than the genome; it may be wrongly joined."),
    "canu_repeat_contigs": (5, "", "Canu flagged some contigs as repeats; they are usually extra copies."),
    "assembler_unknown": (5, "", "The assembler could not be recognised, so only general checks were used."),
    "identity_unclear": (5, "", "The identity report could not confirm the species."),
    "contig_identity": (5, "", "The contigs were checked against a species database."),
    "contig_identity_unchecked": (5, "", "Contigs were not checked against a species database."),
    "duplicate_contig_names": (1, "Check that the right files were given", "Some contig names repeat: the file may hold two assemblies joined together."),
    "confindr_contamination": (3, "Check for contamination", "ConFindr found reads from a second strain (or genus) mixed into this sample."),
    "amr_calls_incomplete": (5, "", "Some resistance genes were found only in part, so their calls need a check."),
    "low_breadth": (2, "Sequence more", "Part of the genome has too few reads, even if the average depth looks fine."),
    "breadth_unchecked": (5, "", "A depth report was found but how evenly the reads cover the genome could not be read."),
    "depth_below_agency": (5, "", "Enough reads for the EU guideline, but a database you may submit to asks for more."),
    "species_not_given": (5, "", "Without the species, Takivisor cannot tell whether the genome size and GC are normal."),
    "contamination_unchecked": (5, "", "Without CheckM2 or an identity report, low-level contamination can be missed."),
}
# Short human headings for the HTML report (the id stays in the technical details, JSON and markdown).
TITLES = {
    "species_unknown": "Species name not recognised", "species_mismatch": "Sample may be a different species",
    "taxonomy_mismatch": "Genome placed in the wrong group", "gc_range": "GC content unusual for this species",
    "report_mismatch": "A report belongs to another assembly", "report_unverified": "A report could not be matched",
    "tech_assembler_mismatch": "Reads do not fit the assembler", "read_gc_mismatch": "Reads may be from another sample",
    "long_reads_not_given": "Long reads not provided", "low_depth": "Not enough reads",
    "illumina_q30": "Low Illumina read quality", "fastp_low_pass": "Many reads failed quality filtering",
    "ont_short_reads": "Nanopore reads are short", "quality_simplified": "Read quality scores are missing",
    "dcs_reads": "Lambda control DNA in the reads", "checkm_contamination": "Contamination estimated by CheckM",
    "contig_contamination": "Contigs from another species", "contig_second_species": "Traces of a second species",
    "read_gc_outliers": "Possible reads from another organism", "identity_mixed": "Second species in the reads",
    "busco_duplicated": "Many genes appear twice", "low_cov_contigs": "Low-coverage contigs",
    "dcs_in_assembly": "Lambda control DNA in the assembly", "lambda_like_sequence": "Lambda-like sequence in the assembly",
    "indel_errors": "Sequence errors breaking genes", "n_content": "Many unknown bases (N)",
    "busco_low": "Core genes missing", "busco_fragmented": "Many core genes broken",
    "checkm_completeness": "Genome looks incomplete", "long_read_fragmented": "Genome in too many pieces",
    "short_read_fragmented": "Genome in very many pieces", "low_n50": "Contigs are short",
    "no_circular_chromosome": "Chromosome not closed", "quast_genome_fraction": "Part of the reference is missing",
    "size_range": "Genome size unusual for this species", "suspect_circular": "Suspicious circular contig",
    "canu_repeat_contigs": "Repeat contigs flagged by Canu", "assembler_unknown": "Assembler not recognised",
    "identity_unclear": "Species could not be confirmed", "contig_identity": "Contigs match the species",
    "contig_identity_unchecked": "Contigs not checked against a species database",
    "low_breadth": "Parts of the genome have few reads", "breadth_unchecked": "Depth breadth not read",
    "species_not_given": "Species not given", "depth_below_agency": "Below a submission database's depth",
    "duplicate_contig_names": "Repeated contig names",
    "confindr_contamination": "Mixed strains (ConFindr)", "amr_calls_incomplete": "Incomplete resistance gene calls",
"contamination_unchecked": "Contamination only partly checked",
}
GLOSSARY = {
    "Contigs": "Pieces of the assembled genome; a complete bacterial genome is usually 1 chromosome plus a few plasmids.",
    "N50 / N90": "Half (or 90%) of the assembly is in contigs at least this long; bigger is more complete.",
    "Largest": "Length of the longest contig.",
    "GC %": "Share of G and C bases; each species has a typical value.",
    "N %": "Share of unknown bases (gaps).",
    "Circular contigs": "Contigs closed into a circle, as bacterial chromosomes and plasmids are.",
    "Median contig cov (header)": "Typical number of reads covering each base, from the assembler.",
    "ORF N50 (frameshift check)": "Typical gene length in amino acids; short values mean errors are breaking genes.",
    "Read N50": "Half of the bases are in reads at least this long.",
    "Mean Q (error-averaged)": "Average read accuracy: Q10 = 90%, Q20 = 99%, Q30 = 99.9%.",
    "Depth breadth (mosdepth)": "Share of the genome read at least 10 (or 20) times; the EU guideline asks 95% at 10x.",
    "% bases >=Q30": "Bases with at least 99.9% accuracy.",
    "% bases >=Q20": "Bases with at least 99% accuracy.",
    "Est. depth (x)": "How many times, on average, each base of the genome was read.",
    "% reads with GC outside assembly range": "Reads whose GC does not fit the assembly; high values suggest another organism.",
    "% long reads that are lambda DCS": "ONT's lambda control DNA mixed into the reads.",
}


# ---------------- fitness for use (what the data can still be used for) ----------------
# Ratings follow the EU reference laboratories' dry-lab QC guideline (EURL 2025, short reads). Where it calls a
# metric critical but a verified study found the analysis still works, the rating is CAUTION and the note names
# both. Where it is silent (long reads, same-species mixes, genetic context) the studies decide. Full matrix and
# sources: ~/takivisor_validation/qc_standards_review_2026-10-05.md, parts 7-8.
ANALYSES = {"id": "Species identification", "mlst": "MLST / serotype", "genes": "AMR/virulence genes present",
            "absence": "Claims that a gene is absent", "relatedness": "Relatedness (cgMLST, SNP distances, outbreaks)",
            "context": "Genetic context (plasmid vs chromosome, complete genome)"}
# NOT_CHECKED: a check that would rate this use is known to be unreliable on this run (e.g. lambda DCS inflating
# a classifier's second species); any real CAUTION/NOT_RECOMMENDED from another check still wins.
RATING_ORDER = {"OK": 0, "NOT_CHECKED": 1, "CAUTION": 2, "NOT_RECOMMENDED": 3}
RATING_WORD = {"OK": "OK", "NOT_CHECKED": "Not checked", "CAUTION": "Use with caution", "NOT_RECOMMENDED": "Not recommended"}
SRC = {"EU": "EU EURL dry-lab QC guideline 2025", "EUCAST": "EUCAST WGS-AST report (Ellington 2017)",
       "PALMA": "Palma et al. 2022", "RAVEN": "Raven et al. 2020", "KOZ": "Kozyreva et al. 2017",
       "ROON": "Rooney et al. 2022", "RGI": "CARD RGI documentation", "PASTEUR": "Institut Pasteur BIGSdb accepted technologies", "GR25": "Genome Research 2025 (ONT R10.4.1 genotyping)",
       "SAND": "Sanderson et al. 2024", "NEU": "Neuenschwander et al. 2026", "LMJCM": "JCM 2024 (ONT Listeria)",
       "GOIG": "Goig et al. 2020", "PIGHT": "Pightling et al. 2019", "CONF": "Low et al. 2019 (ConFindr)",
       "ARR": "Arredondo-Alonso et al. 2017", "SYLPH": "Shaw & Yu 2025 (sylph)", "AQUAMIS": "Deneke et al. 2021 (AQUAMIS)",
       "ENTERO": "EnteroBase QC", "NCBI": "NCBI genome quality rules", "AMRF": "AMRFinderPlus documentation",
       "MLSTDOC": "mlst documentation", "WICK": "Wick et al. 2024", "TK23": "Takivisor real-data check (#23)",
       "ATB": "AllTheBacteria high-quality criteria", "WHOTB": "WHO TB mutation catalogue 2023",
       "RESF": "ResFinder 4.0 (Bortolaia 2020)", "RGI": "CARD RGI documentation",
       "TK27": "Takivisor real-data check (#27): 3 of 5 such calls were the isolate's own reads",
       "TKDCS": "Takivisor lambda DCS check; Kraken2 counts DCS as E. coli/A. baumannii (real run, 2026-10-06)"}
_FIT_CLOSE = ("Listeria", "Enterococcus")   # methylation-prone with PCR-free ONT libraries (LMJCM 8/78; NEU VRE 0-44)


def fitness(res):
    """Per analysis: worst rating among the reasons, each reason naming its basis. Checks that could not run
    are listed under not_checked and do not change the rating."""
    out = {k: {"rating": "OK", "reasons": [], "not_checked": []} for k in ANALYSES}
    ids = {f["id"]: f for f in res["findings"] if f["severity"] != "INFO" or f["id"] in ("no_circular_chromosome", "suspect_circular")}
    reads, rep, sp, a = res.get("reads") or {}, res.get("reports") or {}, res.get("species"), res["assembly"]

    def put(keys, rating, text, *src):
        for k in keys.split():
            o = out[k]
            o["reasons"].append({"rating": rating, "text": text, "basis": [SRC[x] for x in src]})
            if RATING_ORDER[rating] > RATING_ORDER[o["rating"]]:
                o["rating"] = rating

    def miss(keys, text):
        for k in keys.split():
            out[k]["not_checked"].append(text.rstrip("."))

    # --- depth, per sequencing type (hybrid: best of the two for typing, long reads for context)
    depth_rows = {}
    for tech, r in reads.items():
        d = r.get("est_depth")
        if d is None:
            continue
        rows = {}
        if tech == "illumina":
            below = f"Below the EU guideline (30x; depth {d:g}x)."
            rows["id"] = ("OK", f"Depth {d:g}x.", ["EU"]) if d >= 30 else \
                ("CAUTION", below + " Kozyreva 2017 validated species ID at 10x; sylph identifies species below 1x.", ["EU", "KOZ", "SYLPH"])
            rows["mlst"] = ("OK", f"Depth {d:g}x.", ["EU"]) if d >= 30 else \
                ("CAUTION", below + " Kozyreva 2017 validated MLST with a 10x minimum.", ["EU", "KOZ"]) if d >= 10 else \
                ("NOT_RECOMMENDED", f"Depth {d:g}x is below the 10x minimum used in validated MLST.", ["EU", "KOZ"])
            rows["genes"] = ("OK", f"Depth {d:g}x.", ["EU"]) if d >= 30 else \
                ("CAUTION", below + " Rooney 2022 found 15x enough to detect AMR genes in 948 E. coli isolates.", ["EU", "ROON"]) if d >= 15 else \
                ("CAUTION", below + " and below the 15x Rooney 2022 found enough for AMR genes.", ["EU", "ROON"]) if d >= 10 else \
                ("NOT_RECOMMENDED", f"Depth {d:g}x is below every validated minimum.", ["EU", "KOZ"])
            rows["absence"] = ("OK", f"Depth {d:g}x.", ["EU", "EUCAST"]) if d >= 30 else \
                ("CAUTION", below + " Rooney 2022: at 15x a few genes were missed in some isolates.", ["EU", "ROON"]) if d >= 15 else \
                ("NOT_RECOMMENDED", f"Depth {d:g}x: genes can be missed, so 'not detected' can't be trusted.", ["EU", "EUCAST", "ROON"])
            rows["relatedness"] = ("OK", f"Depth {d:g}x.", ["EU"]) if d >= 50 else \
                ("OK", f"Depth {d:g}x meets the EU routine level (30x); for outbreak SNP decisions the EU guideline asks 50x "
                       "(Palma 2022 and Raven 2020 found >= 40x equivalent to higher depth).", ["EU", "PALMA", "RAVEN"]) if d >= 30 else \
                ("CAUTION", below + " Palma 2022 found assembly-based cgMLST about 100% correct from 20x; confirm close links.", ["EU", "PALMA"]) if d >= 20 else \
                ("NOT_RECOMMENDED", f"Depth {d:g}x: allele and SNP errors appear below 20x.", ["EU", "PALMA", "RAVEN"])
            rows["context"] = ("CAUTION", "Short reads can't reliably place genes on plasmids or chromosomes "
                                          "(plasmid binning precision 0.75).", ["ARR"]) if d >= 10 else \
                ("NOT_RECOMMENDED", f"Depth {d:g}x is below every validated minimum.", ["EU", "KOZ"])
        elif tech == "ont":
            rows["id"] = ("OK", f"Depth {d:g}x.", ["SYLPH"]) if d >= 10 else ("CAUTION", f"Depth {d:g}x.", ["SYLPH"])
            rows["mlst"] = ("OK", f"Depth {d:g}x.", ["GR25", "NEU"]) if d >= 25 else \
                ("CAUTION", f"Depth {d:g}x: nanopore typing matched Illumina from 25x.", ["GR25"])
            rows["genes"] = ("OK", f"Depth {d:g}x.", ["SAND"]) if d >= 30 else \
                ("CAUTION", f"Depth {d:g}x: gene recovery above 99% was shown from 30x; lower depth untested.", ["SAND"])
            rows["absence"] = ("OK", f"Depth {d:g}x.", ["SAND", "EUCAST"]) if d >= 30 else \
                ("NOT_RECOMMENDED", f"Depth {d:g}x: below the 30x where complete gene recovery was shown.", ["SAND", "EUCAST"])
            rows["relatedness"] = ("OK", f"Depth {d:g}x.", ["GR25"]) if d >= 35 else \
                ("CAUTION", f"Depth {d:g}x: nanopore cgMLST was accurate from 25x for some genera, 35x for all tested.", ["GR25"]) if d >= 25 else \
                ("NOT_RECOMMENDED", f"Depth {d:g}x: nanopore cgMLST errors were not minimised below 25x.", ["GR25"])
            rows["context"] = ("OK", f"Depth {d:g}x.", ["SAND"]) if d >= 30 else \
                ("CAUTION", f"Depth {d:g}x: plasmid recovery was impaired below 30x.", ["SAND"])
        else:
            miss("relatedness genes absence", f"{TECH_NAME.get(tech, tech)} depth: no peer-reviewed bacterial benchmark, not rated.")
        if rows:
            depth_rows[tech] = rows
    if depth_rows:
        for k in ANALYSES:
            cands = [(t, rows[k]) for t, rows in depth_rows.items() if k in rows]
            if k == "context" and len(cands) > 1:
                cands = [c for c in cands if c[0] != "illumina"] or cands
            t, (rating, text, src) = min(cands, key=lambda c: RATING_ORDER[c[1][0]])
            put(k, rating, f"{TECH_NAME.get(t, t)}: {text}", *src)
    else:
        miss(" ".join(ANALYSES), "Sequencing depth (no reads given).")
    if "ont" in reads:
        m = re.search(r"model (\S+)", reads["ont"].get("tech_basis", ""))
        if m and re.search(r"hac|fast", m.group(1), re.I):
            put("relatedness mlst genes", "CAUTION", f"Basecalling model {m.group(1)}: high-accuracy/fast models gave "
                "outliers with hundreds of allele errors; the superior-accuracy (sup) model is recommended. Pasteur "
                "BIGSdb only accepts nanopore assemblies from R10.4.1 + Dorado sup at >= 50x.", "GR25", "NEU", "PASTEUR")
            put("absence", "NOT_RECOMMENDED", f"Basecalling model {m.group(1)}: indel errors break genes, so 'not detected' "
                "can't be trusted.", "EUCAST")
        elif not m:
            miss("relatedness", "Nanopore basecalling model (not in read names): accuracy figures assume Dorado sup.")
        if sp and sp["species"].split()[0] in _FIT_CLOSE:
            put("relatedness", "CAUTION", f"{sp['species'].split()[0]} with PCR-free nanopore libraries: methylation "
                "causes allele errors in some strains; confirm close links with Illumina or a PCR-based library.", "LMJCM", "NEU")
    # --- Illumina read quality (EU critical, by read length)
    r = reads.get("illumina")
    if r and r.get("pct_q30") is not None and not r.get("qual_simplified"):
        lim = next(m for L, m in T["q30_min_by_len"] if r.get("max_len", 150) <= L)
        q = r["pct_q30"]
        if q < 50:
            put(" ".join(ANALYSES), "NOT_RECOMMENDED", f"Only {q}% of bases >= Q30.", "EU", "KOZ")
        elif q < lim:
            put(" ".join(ANALYSES), "CAUTION", f"{q}% of bases >= Q30 is below the EU guideline ({lim}%); "
                "Kozyreva 2017 accepted runs with >= 50%.", "EU", "KOZ")
    # --- identity and other species
    idr = rep.get("identity")
    dcs = dcs_suspect(idr, reads)
    dcs_why = (f"Species shares in the {idr['tool']} report: {dcs:g}% of long reads are lambda DCS, which it counts as "
               f"{idr['second']}; re-run it after removing DCS reads." if dcs else "")
    if dcs and idr.get("matches_assembly") is not False:
        put("id", "NOT_CHECKED", dcs_why, "TKDCS")
    elif idr and idr.get("matches_assembly") is not False:
        top = idr["top_share"]
        if top >= 95:
            put("id", "OK", f"{idr['tool']}: {idr['top']} {top}% of species-level assignments.", "AQUAMIS")
        elif top >= 70:
            put("id", "CAUTION", f"{idr['tool']}: top species only {top}% (AQUAMIS expects > 95%; EnteroBase accepts >= 70%).", "AQUAMIS", "ENTERO")
        else:
            put("id", "NOT_RECOMMENDED", f"{idr['tool']}: top species only {top}% (< 70%).", "ENTERO")
    else:
        miss("id", "Species confirmation from reads (no Kraken2/sylph/GTDB-Tk report).")
    mixed = ids.get("identity_mixed")
    k = ((rep.get("checkm") or {}).get("contamination")) if "checkm" in rep else None
    share = idr["second_share"] if mixed and idr and not dcs else None
    if dcs and mixed:
        put("mlst genes absence relatedness context", "NOT_CHECKED", dcs_why, "TKDCS")
    lim = second_limit(sp["species"] if sp else (idr or {}).get("top"))
    if (share is not None and share >= lim) or (k is not None and k > 5):
        reads_hit = share is not None and share >= lim
        what = f"{idr['second']} at {share}% of reads" if reads_hit else f"CheckM2 contamination {k}%"
        eu = f"{lim:g}%" if reads_hit else "5%"
        put("mlst", "CAUTION", f"{what}: above the EU limit ({eu}); Pightling 2019 found contamination by another species "
            "caused few typing errors.", "EU", "PIGHT")
        depth = max((x.get("est_depth") or 0) for x in reads.values()) if reads else 0
        if share is not None and share * depth / 100 >= 10:
            put("genes", "NOT_RECOMMENDED", f"{what}: the contaminant has ~{share * depth / 100:.0f}x coverage, enough "
                "for its genes to be detected (from ~10-15x).", "EU", "ROON")
        else:
            put("genes", "NOT_RECOMMENDED", f"{what}: above the EU limit ({eu}); genes may come from the contaminant.", "EU")
        put("absence context", "NOT_RECOMMENDED", f"{what}: above the EU limit ({eu}).", "EU")
        put("relatedness", "CAUTION", f"{what}: above the EU limit ({eu}); Pightling 2019 found few SNP errors from "
            "another species, and Goig 2020 showed filtering contaminant reads removes them.", "EU", "PIGHT", "GOIG")
    elif (idr and not dcs and idr.get("matches_assembly") is not False and idr["second"]
          and idr["second_share"] >= T["id_second_caution"] and not close_relatives(idr["top"], idr["second"])
          and idr.get("n_assigned", T["id_min_reads"]) >= T["id_min_reads"]):
        # #27 real data: at this level the reads were often the isolate's own, misassigned (3 of 5 isolates).
        what = (f"{idr['second']} at {idr['second_share']}% of reads (may be the isolate's own reads misassigned; "
                "map them to the assembly to check)")
        put("mlst", "CAUTION", f"{what}: below the EU limit ({lim:g}%) but detected (AllTheBacteria requires >= 99% "
            "one species; Raven 2020 controls < 0.4% Kraken reads); look for mixed alleles.", "ATB", "RAVEN", "MLSTDOC", "TK27")
        put("genes", "CAUTION", f"{what}: genes may come from it; treat genes at < 1/10 of the isolate's depth "
            "as suspect (ResFinder 4.0).", "RESF", "ROON")
        put("absence context", "CAUTION", f"{what}: below the EU limit ({lim:g}%) but detected.", "EU", "ATB")
        put("relatedness", "CAUTION", f"{what}: filter contaminant reads before SNP calling (Goig 2020: errors even "
            "below 1%).", "GOIG")
    elif {"contig_contamination", "read_gc_outliers"} & ids.keys() or (
            sp and ids.get("size_range") and a["total_length"] > sp["max"]):
        put("mlst", "CAUTION", "Signs of another organism (level unknown): look for mixed alleles.", "MLSTDOC")
        put("genes absence context", "CAUTION", "Signs of another organism (level unknown): genes may come from it.", "EU")
        put("relatedness", "CAUTION", "Signs of another organism: filter contaminant reads before SNP calling.", "GOIG")
    if "busco_duplicated" in ids:
        put("mlst relatedness", "CAUTION", "BUSCO duplication suggests a mix of strains (level unknown); "
            "same-species mixes cause typing errors (Raven 2020 detects them above 10%).", "RAVEN", "PIGHT", "CONF")
        put("genes absence context", "CAUTION", "Possible mix of strains: genes may come from the second strain.", "ROON")
    cf = confindr_row(rep)
    if cf and cf["contaminated"]:
        what = f"ConFindr: {cf['snvs']} contaminating SNVs" + (f" (~{cf['pct']:g}%)" if cf["pct"] is not None else "")
        put("mlst relatedness", "CAUTION", f"{what}: mixed strains give mixed alleles and false SNPs (Pightling 2019: "
            "harm grows with the contaminant's share).", "CONF", "PIGHT")
        put("genes absence context", "CAUTION", f"{what}: genes may come from the second strain.", "CONF", "ROON")
    elif cf and cf["contaminated"] is False:
        put("relatedness", "OK", "ConFindr: no mixing of strains found.", "CONF")
    elif cf:
        miss("relatedness", f"Same-species mixing (ConFindr could not check this sample: genus {cf['genus']}, "
             f"{cf['bases'] if cf.get('bases') is not None else '?'} bases examined; other genera need its rMLST database)")
    else:
        miss("relatedness", "Same-species mixing directly (ConFindr report given but no row for this sample; pass "
             "--sample)." if rep.get("confindr") else
             "Same-species mixing directly (needs a ConFindr report; BUSCO duplication is indirect).")
    md = mosdepth_row(rep, a)
    if "low_breadth" in ids:
        put("absence relatedness", "CAUTION", f"Depth breadth below the EU guideline (95% at 10x): genes and SNPs in "
            "thin regions can be missed.", "EU", "WHOTB")
    elif md and md[0]["at"]:
        put("absence", "OK", f"mosdepth: {md[0]['at'][10]:g}% of {md[1]} at >= 10x.", "EU")
    else:
        miss("absence relatedness", "Depth breadth (needs a mosdepth report).")
    # --- AMR gene calls cut or broken by the assembly
    for key, src in (("amrfinder", "AMRF"), ("rgi", "RGI")):
        hits = amr_hits(rep[key], a) if rep.get(key) else None
        cut = sorted({h["gene"] for h in hits or [] if h["issue"] in ("edge", "partial", "stop")})
        if cut:
            tool = rep[key]["tool"]
            put("genes", "CAUTION", f"{tool}: {len(cut)} gene(s) cut at contig ends, partial or with internal stops "
                f"({', '.join(cut[:3])}{', ...' if len(cut) > 3 else ''}); confirm before reporting.", src)
            if any(h["issue"] in ("edge", "partial") for h in hits):
                put("absence", "CAUTION", f"{tool}: genes cut at contig ends or truncated were found here, so others "
                    "may be missed.", src)
    # --- species identity
    for i in ("species_mismatch", "taxonomy_mismatch", "gc_range", "species_unknown"):
        if i in ids:
            put("id mlst absence relatedness", "NOT_RECOMMENDED", f"Species identity unresolved ({i.replace('_', ' ')}): "
                "species-specific schemes and gene rules can't be chosen.", "EU", "MLSTDOC", "AMRF")
            put("genes context", "CAUTION", "Species identity unresolved; gene presence itself is still informative.", "AMRF")
            break
    if not sp:
        miss("mlst relatedness", "Species not given: species-specific schemes can't be checked.")
    # --- completeness
    ck = rep.get("checkm")
    small = sp and ids.get("size_range") and a["total_length"] < sp["min"]
    c = ck.get("completeness") if ck else None
    if (c is not None and c < 90) or small:
        why = f"CheckM2 completeness {c}%" if c is not None and c < 90 else "genome smaller than the species range"
        put("mlst", "CAUTION", f"{why}: some alleles may be missing.", "MLSTDOC")
        put("absence relatedness context", "NOT_RECOMMENDED", f"{why} (EU guideline 95%, NCBI 90%).", "EU", "NCBI")
    elif c is not None and c < 95:
        put("absence relatedness context", "CAUTION", f"CheckM2 completeness {c}%: below the EU guideline (95%); "
            "NCBI accepts >= 90%.", "EU", "NCBI")
    elif c is None:
        miss("absence relatedness", "Completeness (no CheckM2 report).")
    # --- assembly structure and errors
    if {"short_read_fragmented", "low_n50"} & ids.keys():
        put("absence", "CAUTION", "Fragmented assembly: genes split at contig ends can look partial or missing.", "AMRF", "EUCAST")
        put("relatedness", "CAUTION", "Fragmented assembly: check that >= 95% of cgMLST loci are called (EU guideline).", "EU")
        put("context", "NOT_RECOMMENDED", "Fragmented short-read assembly: plasmid reconstruction is limited by contig length.", "ARR")
    if "indel_errors" in ids:
        put("mlst genes", "CAUTION", "Frameshift (indel) errors: alleles and genes can look novel, partial or truncated.", "MLSTDOC", "AMRF")
        put("absence relatedness", "NOT_RECOMMENDED", "Frameshift (indel) errors: false allele differences and missed genes; "
            "polish first.", "AMRF", "NEU", "GR25", "LMJCM")
    ref = sp["expected"] if sp else a["total_length"]
    closed = any(c.get("circular") and c["length"] >= min(T["circ_min_len"], T["near_genome_frac"] * ref)
                 for c in a.get("contigs", []))
    if {"long_read_fragmented", "no_circular_chromosome", "suspect_circular"} & ids.keys() or not closed:
        put("context", "CAUTION", "Chromosome not closed (or a suspect circle): not a complete genome.", "NCBI")
    if "dcs_in_assembly" in ids:
        put("context", "CAUTION", "Lambda control (DCS) sequence in the assembly: remove it and reassemble.", "TK23")
    return out


def acts_of(f):
    """(lab, bioinformatics) actions: the finding's own pair if it has one, else ACTIONS by id."""
    return f.get("acts") or ACTIONS.get(f["id"])


# Lab and bioinformatics actions per finding (what to do at the bench vs at the computer).
ACTIONS = {
    "low_depth": ("Top up or re-sequence: the EU guideline asks 30x for routine work and 50x for outbreak SNP analysis.",
                  "Prefer read-mapping SNP calling with a per-site depth floor; call a gene absent only if every base "
                  "of it has >= 10x (WHO TB catalogue 2023)."),
    "low_breadth": ("Top up the reads; if the gaps stay, check for lost plasmids or a mixed culture.",
                    "Call a gene absent only if every base of it has >= 10x (WHO TB catalogue 2023); mask low-depth "
                    "regions in SNP calling."),
    "illumina_q30": ("Check the run (cluster density, PhiX error rate, % passing filter); re-run the library if Q30 stays low.",
                     "Quality-trim (e.g. fastp --cut_right) and use base-quality filters for SNP calling."),
    "fastp_low_pass": ("Check library and run quality; re-sequence if many reads fail.",
                       "Check what fastp removed (adapters, low quality, too short) and whether the settings were too strict."),
    "ont_short_reads": ("Use a gentler DNA extraction and size selection; rapid kits keep small plasmids.",
                        "Filter very short reads only if depth allows; repeats need reads longer than the repeat."),
    "identity_mixed": ("Re-streak a single colony, re-extract and re-sequence (Raven 2020: 22/25 failures passed after re-isolation).",
                       "Remove contaminant reads with a taxonomic classifier before SNP calling (Goig 2020); don't report genes "
                       "until resolved, and treat genes at < 1/10 of the isolate's depth as the contaminant's (ResFinder 4.0)."),
    "contig_contamination": ("Re-streak a single colony and re-sequence if the contamination is confirmed.",
                             "Remove the foreign contigs and screen with NCBI FCS-GX before submission; genes on those contigs "
                             "are the contaminant's."),
    "read_gc_outliers": ("If an identity report confirms another organism, re-streak a single colony and re-sequence.",
                         "Run Kraken2/sylph on the reads to name the outlying reads; filter them before SNP calling."),
    "confindr_contamination": ("Re-streak a single colony and re-sequence from it.",
                               "Don't use this genome for typing or SNP distances until a pure culture is sequenced; "
                               "ConFindr's contamination CSV lists the mixed sites."),
    "duplicate_contig_names": ("Check sample tracking: the assembly file may hold two samples.",
                               "Check how the FASTA was made (concatenated assemblies repeat names); rename "
                               "contigs uniquely before per-contig analysis."),
    "read_gc_mismatch": ("Check sample tracking: these reads may belong to another sample.",
                         "Confirm the read and assembly files belong together; run an identity report on the reads."),
    "species_mismatch": ("Check sample labels; confirm the species (MALDI-TOF or ANI).",
                         "Re-run species-specific steps (MLST/cgMLST scheme, AMRFinderPlus --organism) with the confirmed species."),
    "taxonomy_mismatch": ("Check sample labels; confirm the species (MALDI-TOF or ANI).",
                          "Re-run species-specific steps with the confirmed species."),
    "gc_range": ("Check sample labels and culture purity.", "Confirm the species with an identity report (Kraken2/sylph/ANI)."),
    "size_range": ("Too large: check the culture for a mix. Too small: check DNA extraction (plasmid loss).",
                   "Too large: look for contaminant or duplicated contigs. Too small: check assembly filters and depth."),
    "checkm_contamination": ("Re-streak a single colony and re-sequence.",
                             "Name the contaminant (identity report, low-coverage contigs) and remove it before analysis."),
    "checkm_completeness": ("Check depth and DNA quality; long reads can recover missing regions.",
                            "Check assembly filters (minimum contig length/coverage) and reassemble."),
    "busco_duplicated": ("Re-streak a single colony (possible mix of strains).",
                         "Check for duplicated contigs or untrimmed ends; run ConFindr on the reads."),
    "short_read_fragmented": ("Check DNA integrity and depth; add long reads to close the genome.",
                              "Reassemble and remove low-coverage contigs; for typing check that >= 95% of cgMLST loci are called; "
                              "for AMR run RGI with --low_quality (contigs < 20 kb) and count AMRFinderPlus PARTIAL_CONTIG_ENDX/ENDP hits."),
    "low_n50": ("Check DNA integrity and depth.", "Reassemble; prefer read-mapping methods for typing; for AMR run RGI "
                "with --low_quality (contigs < 20 kb) and count AMRFinderPlus PARTIAL_CONTIG_ENDX/ENDP hits."),
    "long_read_fragmented": ("Increase depth or read length.", "Try another assembler or Autocycler; check read N50 and depth."),
    "indel_errors": ("Basecall with Dorado sup (bacterial methylation model) at >= 40x, or add Illumina reads.",
                     "Polish with Medaka, then Polypolish + Pypolca if Illumina reads exist (most gain by 25x, Wick 2024)."),
    "no_circular_chromosome": ("More depth or longer reads.", "Try another assembler or Autocycler; rotate with dnaapler once closed."),
    "dcs_in_assembly": ("None needed: lambda DCS is ONT's control spike-in.", "Remove DCS reads (NanoLyse or chopper) and reassemble."),
    "dcs_reads": ("None needed: lambda DCS is ONT's control spike-in.", "Remove DCS reads (NanoLyse or chopper) before assembly."),
    "low_cov_contigs": ("If many, check culture purity.", "Remove contigs below ~20% of median coverage (EnteroBase) and re-check."),
}


# Bioinformatics-tab titles for "Start here" steps that are lab work; other steps keep their title.
BIO_STEPS = {
    "Confirm which organism this sample really is": "Confirm the species from the reads (identity report or ANI)",
    "Check the Illumina run quality": "Quality-trim the reads before analysis",
    "Get longer nanopore reads (DNA extraction, size selection)": "Work around the short reads",
    "Sequence more": "Work around the low depth",
    "Check for contamination": "Name the contaminant and remove its reads and contigs",
}


def start_here(findings):
    """Ordered to-do list from the WARN/FAIL findings: [(step text, [finding ids])]."""
    probs = [f for f in findings if f["severity"] != "INFO"]
    steps, later = {}, []
    for f in probs:
        stage, step, _ = GUIDE.get(f["id"], (4, "", ""))
        if stage == 5 or not step:
            later.append(f["id"])
        else:
            steps.setdefault((stage, step), []).append(f["id"])
    out = [(step, ids) for (_, step), ids in sorted(steps.items())]
    if out and later:
        out.append(("Then re-run Takivisor; these may improve after the steps above", later))
    elif later:   # only consequence-type problems: act on them directly
        out += [(f["action"], [f["id"]]) for f in probs if f["id"] in later]
    elif out:
        out.append(("Then re-run Takivisor to confirm", []))
    return out


# ---------------- HTML report (--html) ----------------
# One standalone file, stdlib only: verdict banner, finding cards, inline-SVG charts with hover tooltips
# and a data table each, then the markdown tables. Colours: dataviz reference palette, slots 1-2
# validated all-pairs in light and dark mode (2026-10-01); status colours always paired with icon + label.
import html as _html

_E = lambda x: _html.escape(str(x), quote=True)
_W, _H, _ML, _MR, _MT, _MB = 640, 250, 66, 16, 14, 40   # chart box and margins (px)
_STATUS = {"PASS": ("good", "✓"), "REVIEW": ("warning", "!"), "FAIL": ("critical", "✕"),
           "WARN": ("warning", "!"), "INFO": ("info", "i")}


def _ticks(lo, hi, n=5):
    """Clean linear ticks (1/2/5 x 10^k steps)."""
    span = (hi - lo) or 1
    raw = span / n
    step = min((m * 10 ** math.floor(math.log10(raw)) for m in (1, 2, 5, 10)), key=lambda x: abs(x - raw))
    t = math.ceil(lo / step) * step
    out = []
    while t <= hi + step * 1e-9:
        out.append(round(t, 10)); t += step
    return out


def _fmt_len(v):
    return f"{v/1e6:.3g} Mb" if v >= 1e6 else f"{v/1e3:.3g} kb" if v >= 1e3 else f"{v:.3g} bp"


def _axes(xs, ys, xlab, ylab, xfmt=str, yfmt=str):
    """Gridlines + tick labels. xs/ys: [(pixel, value)]."""
    o = []
    for py, v in ys:
        o.append(f'<line class="grid" x1="{_ML}" x2="{_W - _MR}" y1="{py:.1f}" y2="{py:.1f}"/>'
                 f'<text class="tick" x="{_ML - 6}" y="{py + 4:.1f}" text-anchor="end">{_E(yfmt(v))}</text>')
    for px, v in xs:
        o.append(f'<text class="tick" x="{px:.1f}" y="{_H - _MB + 16}" text-anchor="middle">{_E(xfmt(v))}</text>')
    o.append(f'<line class="axis" x1="{_ML}" x2="{_W - _MR}" y1="{_H - _MB}" y2="{_H - _MB}"/>'
             f'<text class="lab" x="{(_ML + _W - _MR) / 2}" y="{_H - 4}" text-anchor="middle">{_E(xlab)}</text>'
             f'<text class="lab" transform="translate(12,{(_MT + _H - _MB) / 2}) rotate(-90)" text-anchor="middle">{_E(ylab)}</text>')
    return "".join(o)


def _chart(title, sub, body, legend, head, rows):
    """Card: title, subtitle, legend (>= 2 series), SVG, collapsible data table."""
    leg = "".join(f'<span class="key"><i class="{k}"></i>{_E(t)}</span>' for k, t in legend) if len(legend) > 1 else ""
    tbl = "".join("<tr>" + "".join(f"<td>{_E(c)}</td>" for c in r) + "</tr>" for r in rows)
    return (f'<section class="card"><h3>{_E(title)}</h3><p class="sub">{_E(sub)}</p>'
            + (f'<div class="legend">{leg}</div>' if leg else "")
            + f'<svg viewBox="0 0 {_W} {_H}" role="img" aria-label="{_E(title)}">{body}</svg>'
            + f'<details><summary>Data table</summary><table><tr>{"".join(f"<th>{_E(h)}</th>" for h in head)}</tr>'
            + f"{tbl}</table></details></section>")


def _hit(x, y, w, h, tip):
    return f'<rect class="hit" x="{x:.1f}" y="{y:.1f}" width="{max(w, 1):.1f}" height="{h:.1f}" data-tip="{_E(tip)}"/>'


def _chart_cumulative(a, sp):
    cs = sorted((c["length"] for c in a["contigs"]), reverse=True)
    cum, t = [], 0
    for l in cs:
        t += l; cum.append(t)
    top = max(cum[-1], sp["max"] if sp else 0) * 1.05
    n = len(cs)
    X = lambda i: _ML + (i / max(1, n - 1) if n > 1 else 0.5) * (_W - _ML - _MR)
    Y = lambda v: _H - _MB - v / top * (_H - _MB - _MT)
    o = []
    if sp:   # NCBI species range as a wash band, expected size as a hairline
        o.append(f'<rect class="band" x="{_ML}" width="{_W - _ML - _MR}" y="{Y(sp["max"]):.1f}" '
                 f'height="{Y(sp["min"]) - Y(sp["max"]):.1f}"/><line class="ref" x1="{_ML}" x2="{_W - _MR}" '
                 f'y1="{Y(sp["expected"]):.1f}" y2="{Y(sp["expected"]):.1f}"/><text class="note" x="{_ML + 6}" '
                 f'y="{Y(sp["max"]) + 13:.1f}">NCBI size range for {_E(sp["species"])}</text>')
    o.append(_axes([(X(i), i + 1) for i in sorted({0, n - 1} | {round(k * (n - 1) / 4) for k in range(5)})],
                   [(Y(v), v) for v in _ticks(0, top)], "Contig rank (longest first)", "Cumulative length",
                   yfmt=lambda v: f"{v/1e6:g} Mb"))
    pts = " ".join(f"{X(i):.1f},{Y(v):.1f}" for i, v in enumerate(cum))
    o.append(f'<polyline class="s1 line" points="{pts}"/><circle class="s1 dot" cx="{X(n - 1):.1f}" cy="{Y(cum[-1]):.1f}" r="4"/>'
             f'<text class="note" x="{X(n - 1) - 6:.1f}" y="{Y(cum[-1]) - 8:.1f}" text-anchor="end">{cum[-1]/1e6:.2f} Mb</text>')
    bw = (_W - _ML - _MR) / max(1, n)
    o += [_hit(X(i) - bw / 2, _MT, bw, _H - _MB - _MT, f"Contig {i + 1}: {cs[i]:,} bp · cumulative {v/1e6:.3f} Mb")
          for i, v in enumerate(cum[:3000])]
    sub = (f"{n} contig{'s' * (n != 1)}, {cum[-1]/1e6:.2f} Mb" + (f"; expected {sp['expected']/1e6:.2f} Mb "
           f"(NCBI {sp['min']/1e6:.2f}-{sp['max']/1e6:.2f} Mb)" if sp else "; no species given, no expected size"))
    rows = [(i + 1, f"{l:,}", f"{c:,}") for i, (l, c) in enumerate(zip(cs, cum))][:500]
    return _chart("Assembly length", sub, "".join(o), [], ("Rank", "Length (bp)", "Cumulative (bp)"), rows)


def _chart_coverage(a):
    cs = [c for c in a["contigs"] if c["cov"] and c["length"]]
    med = a.get("median_contig_cov")
    if len(cs) < 2 or not med:
        return ""
    rel = [c["cov"] / med for c in cs]
    lx0, lx1 = math.floor(math.log10(min(c["length"] for c in cs))), math.ceil(math.log10(max(c["length"] for c in cs)))
    ly0, ly1 = math.floor(math.log10(min(min(rel), 0.1))), math.ceil(math.log10(max(max(rel), 10)))
    lx1, ly1 = max(lx1, lx0 + 1), max(ly1, ly0 + 1)
    X = lambda l: _ML + (math.log10(l) - lx0) / (lx1 - lx0) * (_W - _ML - _MR)
    Y = lambda r: _H - _MB - (math.log10(r) - ly0) / (ly1 - ly0) * (_H - _MB - _MT)
    o = [_axes([(X(10 ** k), 10 ** k) for k in range(lx0, lx1 + 1)], [(Y(10 ** k), 10 ** k) for k in range(ly0, ly1 + 1)],
               "Contig length (log)", "Coverage / median (log)", xfmt=_fmt_len, yfmt=lambda v: f"{v:g}x"),
         f'<line class="ref" x1="{_ML}" x2="{_W - _MR}" y1="{Y(1):.1f}" y2="{Y(1):.1f}"/>'
         f'<text class="note" x="{_W - _MR - 4}" y="{Y(1) - 5:.1f}" text-anchor="end">median {med:.3g}x</text>']
    circ_known = any(c["circular"] is not None for c in cs)
    for c, r in sorted(zip(cs, rel), key=lambda x: -x[0]["length"]):
        k = "s2" if c["circular"] else "s1"
        tip = f'{c["name"]}: {c["length"]:,} bp · {c["cov"]:g}x ({r:.2f}x median)' + (" · circular" if c["circular"] else "")
        o.append(f'<circle class="{k} dot" cx="{X(c["length"]):.1f}" cy="{Y(r):.1f}" r="4"/>'
                 f'<circle class="hit" cx="{X(c["length"]):.1f}" cy="{Y(r):.1f}" r="12" data-tip="{_E(tip)}"/>')
    legend = [("s1", "Linear"), ("s2", "Circular")] if circ_known else []
    rows = [(c["name"], f'{c["length"]:,}', f'{c["cov"]:g}', f"{r:.2f}", {True: "yes", False: "no"}.get(c["circular"], "?"))
            for c, r in sorted(zip(cs, rel), key=lambda x: -x[0]["length"])][:500]
    return _chart("Contig coverage vs length", "Contaminant contigs usually sit well below the median; "
                  "plasmids at or above it.", "".join(o), legend,
                  ("Contig", "Length (bp)", "Coverage", "x median", "Circular"), rows)


def _chart_bars(title, sub, bins, xlab, ylab, xfmt, tipfmt, xticks, band=None, marker=None, head=("Bin", "Value")):
    """Column histogram. bins: [(x0, x1, value)] on a numeric x; xticks: [(x, label)]; band: (x0, x1, label);
    marker: (x, label)."""
    x0, x1 = bins[0][0], bins[-1][1]
    top = max(v for *_, v in bins) * 1.1 or 1
    X = lambda v: _ML + (v - x0) / (x1 - x0) * (_W - _ML - _MR)
    Y = lambda v: _H - _MB - v / top * (_H - _MB - _MT)
    o = []
    if band:
        o.append(f'<rect class="band" x="{X(band[0]):.1f}" width="{X(band[1]) - X(band[0]):.1f}" y="{_MT}" '
                 f'height="{_H - _MB - _MT}"/><text class="note" x="{X(band[0]) + 4:.1f}" y="{_MT + 12}">{_E(band[2])}</text>')
    o.append(_axes([], [(Y(v), v) for v in _ticks(0, top)], xlab, ylab, yfmt=lambda v: f"{v:g}%"))
    o += [f'<text class="tick" x="{X(t):.1f}" y="{_H - _MB + 16}" text-anchor="middle">{_E(lab)}</text>'
          for t, lab in xticks if x0 <= t <= x1]
    for b0, b1, v in bins:
        w = X(b1) - X(b0)
        bw = min(24, max(1, w - 2))
        h = Y(0) - Y(v)
        x = X(b0) + (w - bw) / 2
        if v > 0:   # 4px rounded data end, square at the baseline
            r = min(4, bw / 2, h)
            o.append(f'<path class="s1 bar" d="M{x:.1f},{Y(0):.1f}V{Y(v) + r:.1f}q0,{-r:.1f} {r:.1f},{-r:.1f}'
                     f'H{x + bw - r:.1f}q{r:.1f},0 {r:.1f},{r:.1f}V{Y(0):.1f}Z"/>')
        o.append(_hit(X(b0), _MT, w, _H - _MB - _MT, tipfmt(b0, b1, v)))
    if marker:
        o.append(f'<line class="ref" x1="{X(marker[0]):.1f}" x2="{X(marker[0]):.1f}" y1="{_MT}" y2="{_H - _MB}"/>'
                 f'<text class="note" x="{X(marker[0]) + 4:.1f}" y="{_MT + 26}">{_E(marker[1])}</text>')
    return _chart(title, sub, "".join(o), [], head, [(xfmt(b0), f"{v:.2f}") for b0, b1, v in bins])


def _chart_read_length(tech, r):
    h = r.get("len_hist")
    if not h or len(h) < 3:   # Illumina: one or two bins, not worth a chart
        return ""
    ks = range(min(h), max(h) + 1)
    tot = sum(h.values())
    bins = [(10 ** (k / 10), 10 ** ((k + 1) / 10), 100 * h.get(k, 0) / tot) for k in ks]
    lb = [(math.log10(b0), math.log10(b1), v) for b0, b1, v in bins]
    return _chart_bars(f"Read length ({tech})", f"{r['n_reads']:,} reads, {r['bases']/1e6:,.1f} Mb; "
                       "share of bases per length bin", lb, "Read length (log)", "% of bases",
                       lambda v: _fmt_len(float(f"{10 ** v:.2g}")),
                       lambda b0, b1, v: f"{_fmt_len(round(10 ** b0))}-{_fmt_len(round(10 ** b1))}: {v:.1f}% of bases",
                       [(math.log10(m * 10 ** e), _fmt_len(m * 10 ** e)) for e in range(1, 8) for m in (1, 2, 5)],
                       marker=(math.log10(r["read_n50"]), f"N50 {_fmt_len(r['read_n50'])}"),
                       head=("Bin start", "% of bases"))


def _chart_read_gc(tech, r, a):
    rg, rng = r.get("read_gc"), a.get("gc_window_range") or {}
    if not rg:
        return ""
    h = Counter()
    for (c, g), n in rg.items():
        h[g] += n
    tot = sum(h.values())
    lo, hi = max(0, min(h) - 3), min(100, max(h) + 4)
    bins = [(g, g + 1, 100 * h.get(g, 0) / tot) for g in range(lo, hi)]
    cls = Counter()
    for (c, g), n in rg.items():
        cls[c] += n
    W = next((c for c, _ in cls.most_common() if c in rng), None)
    band = (rng[W][0], rng[W][1] + 1, f"assembly 99% range ({_fmt_len(W)} windows)") if W else None
    dcs = r.get("dcs_pct", 0) >= T["dcs_reads_pct"] and r["dcs_pct"]
    sub = (f"{tot:,} reads" + (f"; {r['gc_outlier_pct']}% outside the assembly's range for their length"
                               if "gc_outlier_pct" in r else "") + ". A second hump suggests another organism"
           + (f"; the line marks ONT's lambda DCS control ({DCS_GC:g}% GC, {dcs:g}% of reads)." if dcs else "."))
    return _chart_bars(f"Read GC ({tech})", sub, bins, "GC %", "% of reads", lambda v: f"{v:g}",
                       lambda b0, b1, v: f"GC {b0}%: {v:.2f}% of reads", [(g, str(g)) for g in range(0, 101, 5)],
                       band=band, head=("GC %", "% of reads"),
                       marker=(max(range(int(DCS_GC) - 2, int(DCS_GC) + 3), key=lambda g: h.get(g, 0)) + 0.5,
                               f"lambda DCS {dcs:g}%") if dcs and lo <= DCS_GC < hi else None)


def _inl(t):
    return re.sub(r"`([^`]+)`", r"<code>\1</code>", re.sub(r"\*\*([^*]+)\*\*", r"<b>\1</b>", _E(t)))


def _md_html(md):
    """Our own markdown subset (headings, tables, bullets, bold, code) -> HTML."""
    inl = _inl
    out, tbl = [], []
    for line in md.splitlines() + [""]:
        if line.startswith("|"):
            cells = [c.strip() for c in line.strip().strip("|").split("|")]
            if not set("".join(cells)) <= set("-: "):
                tbl.append(cells)
            continue
        if tbl:
            cell = lambda c, j: (f'<abbr title="{_E(GLOSSARY[c])}">{inl(c)}</abbr>' if j == 0 and c in GLOSSARY else inl(c))
            out.append("<table>" + "".join("<tr>" + "".join(f"<{'th' if i == 0 else 'td'}>{cell(c, j)}</{'th' if i == 0 else 'td'}>"
                                                          for j, c in enumerate(r)) + "</tr>" for i, r in enumerate(tbl)) + "</table>")
            tbl = []
        if line.startswith("## "):
            out.append(f"<h2>{inl(line[3:])}</h2>")
        elif line.startswith("- "):
            out.append(f"<p class='li'>{inl(line[2:])}</p>")
    return "".join(out)


_CSS = """
:root{color-scheme:light;--page:#f9f9f7;--surface:#fcfcfb;--ink:#0b0b0b;--ink2:#52514e;--muted:#898781;
--grid:#e1e0d9;--axis:#c3c2b7;--ring:rgba(11,11,11,.10);--s1:#2a78d6;--s2:#eb6834;--wash:rgba(42,120,214,.10);
--good:#0ca30c;--warning:#fab219;--critical:#d03b3b;--info:#898781}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){color-scheme:dark;--page:#0d0d0d;--surface:#1a1a19;
--ink:#fff;--ink2:#c3c2b7;--grid:#2c2c2a;--axis:#383835;--ring:rgba(255,255,255,.10);--s1:#3987e5;--s2:#d95926;
--wash:rgba(57,135,229,.14)}}
:root[data-theme="dark"]{color-scheme:dark;--page:#0d0d0d;--surface:#1a1a19;--ink:#fff;--ink2:#c3c2b7;--grid:#2c2c2a;
--axis:#383835;--ring:rgba(255,255,255,.10);--s1:#3987e5;--s2:#d95926;--wash:rgba(57,135,229,.14)}
*{box-sizing:border-box}body{margin:0;background:var(--page);color:var(--ink);font:15px/1.5 system-ui,-apple-system,"Segoe UI",sans-serif}
main{max-width:1100px;margin:0 auto;padding:24px 16px 48px}h1{font-size:22px;margin:0 0 4px}h2{font-size:17px;margin:28px 0 8px}
h3{font-size:15px;margin:0}.sub{color:var(--ink2);margin:2px 0 8px;font-size:13px}
.banner{display:flex;gap:14px;align-items:center;padding:14px 18px;border-radius:10px;background:var(--surface);
box-shadow:0 0 0 1px var(--ring);border-left:6px solid var(--st)}.banner .ic{flex:none;width:34px;height:34px;border-radius:50%;
background:var(--st);color:#fff;display:grid;place-items:center;font-weight:700;font-size:18px}
.banner b{font-size:20px}.meta{color:var(--ink2);font-size:13px;margin:10px 0 0}
.finds{display:grid;gap:8px}.find{background:var(--surface);box-shadow:0 0 0 1px var(--ring);border-radius:8px;
padding:10px 14px;border-left:4px solid var(--st)}.find .tag{font-weight:700;font-size:12px;letter-spacing:.03em}
.find .tag i{font-style:normal;display:inline-grid;place-items:center;width:18px;height:18px;border-radius:50%;
background:var(--st);color:#fff;margin-right:6px;font-size:11px}.find .act{color:var(--ink2);font-size:13px}
.charts{display:grid;grid-template-columns:repeat(auto-fit,minmax(min(100%,480px),1fr));gap:12px}
.card{background:var(--surface);box-shadow:0 0 0 1px var(--ring);border-radius:10px;padding:14px;min-width:0}
svg{width:100%;height:auto;display:block;overflow:visible}.grid{stroke:var(--grid);stroke-width:1}.axis{stroke:var(--axis);stroke-width:1}
.tick,.note{fill:var(--muted);font-size:11px;font-variant-numeric:tabular-nums}.lab{fill:var(--ink2);font-size:12px}
.ref{stroke:var(--muted);stroke-width:1}.band{fill:var(--wash)}.line{fill:none;stroke:var(--s1);stroke-width:2;stroke-linejoin:round;stroke-linecap:round}
.s1.dot,.s1.bar{fill:var(--s1)}.s2.dot{fill:var(--s2)}.dot{stroke:var(--surface);stroke-width:2}
.hit{fill:transparent;cursor:crosshair}.hit:hover,.hit:focus{fill:var(--ink);fill-opacity:.06;outline:none}
.legend{display:flex;gap:14px;font-size:12px;color:var(--ink2);margin-bottom:4px}.key i{display:inline-block;width:10px;height:10px;
border-radius:50%;margin-right:5px;vertical-align:-1px}.key i.s1{background:var(--s1)}.key i.s2{background:var(--s2)}
details{margin-top:6px;font-size:13px}summary{cursor:pointer;color:var(--ink2)}table{border-collapse:collapse;margin:6px 0;font-size:13px;
max-width:100%;display:block;overflow-x:auto}th,td{text-align:left;padding:4px 10px;border-bottom:1px solid var(--grid);
font-variant-numeric:tabular-nums}code{font-size:12px;background:var(--grid);padding:1px 4px;border-radius:4px}
.li{margin:4px 0}.find h3{margin:2px 0 4px;font-size:16px}.find p{margin:4px 0}.find details{margin-top:4px}
.tag{font-size:12px}.meta{display:flex;flex-wrap:wrap;gap:6px 22px;margin:12px 0 4px;color:var(--ink2);font-size:14px}
.meta b{color:var(--ink);font-weight:600;margin-right:4px}.pill{font-size:11px;padding:1px 7px;border-radius:9px;
box-shadow:0 0 0 1px var(--ring);color:var(--ink2);margin-left:4px}.pill.good{color:var(--ink);box-shadow:0 0 0 1px var(--good)}
.pill.critical{color:var(--ink);box-shadow:0 0 0 1px var(--critical)}.pill.warning{color:var(--ink);box-shadow:0 0 0 1px var(--warning)}.fit{background:var(--surface);box-shadow:0 0 0 1px var(--ring);border-radius:10px;padding:2px 18px 10px;margin-top:16px}.fit table{width:100%}.fit td,.fit th{vertical-align:top}.fit .pill{white-space:nowrap;display:inline-block;margin:0}@media(max-width:600px){.fit table,.fit tbody,.fit tr,.fit th,.fit td{display:block}.fit tr{padding:6px 0;border-bottom:1px solid var(--ring)}.fit th{border:0;padding:0}.fit td{border:0;padding:2px 0}}.fit .src{color:var(--ink2);font-size:12px}
.tiles{display:grid;grid-template-columns:repeat(auto-fit,minmax(min(100%,165px),1fr));gap:10px;margin:14px 0 4px}
.tile{background:var(--surface);box-shadow:0 0 0 1px var(--ring);border-radius:10px;padding:12px 14px;border-top:4px solid var(--st)}
.tl{font-size:13px;color:var(--ink2)}.tv{font-size:24px;font-weight:600;margin:2px 0}.ts{font-size:12px;color:var(--muted);min-height:2.4em}
.chip{font-size:12px;font-weight:600;margin-top:6px}.chip i{font-style:normal;display:inline-grid;place-items:center;width:16px;
height:16px;border-radius:50%;background:var(--st);color:#fff;margin-right:5px;font-size:10px}
.start a,.find a{color:var(--s1)}details.sec{background:var(--surface);box-shadow:0 0 0 1px var(--ring);border-radius:10px;
padding:10px 14px;margin:8px 0}details.sec>summary{font-weight:600;color:var(--ink)}.files li{margin:3px 0;word-break:break-all}
.foot{margin-top:24px}.plain{font-size:13px;margin:2px 0}abbr{text-decoration:underline dotted var(--muted);cursor:help}
.start{background:var(--surface);box-shadow:0 0 0 1px var(--ring);border-radius:10px;padding:2px 18px 10px;margin-top:16px}
.start h2{margin:14px 0 6px}.start ol{margin:0;padding-left:22px}.start li{margin:8px 0}.todo{font-size:14px;color:var(--ink2);margin:2px 0}
.tabs>input{position:absolute;opacity:0;pointer-events:none}.tabs>label{display:inline-block;padding:6px 14px;margin:0 4px 10px 0;
border-radius:999px;box-shadow:0 0 0 1px var(--ring);cursor:pointer;font-size:14px;color:var(--ink2)}
.tabs>input:checked+label{background:var(--s1);color:#fff;box-shadow:none;font-weight:600}
.tabs>input:focus-visible+label{outline:2px solid var(--s1);outline-offset:2px}.tab{display:none}
#st-lab:checked~#tab-lab,#st-bio:checked~#tab-bio{display:block}#tip{position:fixed;pointer-events:none;background:var(--surface);color:var(--ink);box-shadow:0 2px 10px rgba(0,0,0,.18),
0 0 0 1px var(--ring);border-radius:6px;padding:6px 9px;font-size:12px;max-width:320px;display:none;z-index:9}
"""

_JS = """
const tip=document.getElementById('tip');
function show(e,el){tip.textContent=el.dataset.tip;tip.style.display='block';
 const r=el.getBoundingClientRect(),x=e&&e.clientX!==undefined?e.clientX:r.left+r.width/2,y=e&&e.clientY!==undefined?e.clientY:r.top;
 const w=tip.offsetWidth;tip.style.left=Math.min(window.innerWidth-w-8,Math.max(8,x+12))+'px';tip.style.top=(y-tip.offsetHeight-12)+'px';}
document.querySelectorAll('.hit').forEach((el,i)=>{if(i<400)el.setAttribute('tabindex','0');
 el.addEventListener('pointermove',e=>show(e,el));el.addEventListener('focus',()=>show(null,el));
 ['pointerleave','blur'].forEach(t=>el.addEventListener(t,()=>tip.style.display='none'));});
"""


TECH_NAME = {"ont": "Oxford Nanopore", "illumina": "Illumina", "hifi": "PacBio HiFi", "clr": "PacBio CLR"}
ASM_NAME = {"spades": "SPAdes", "skesa": "SKESA", "hifiasm": "hifiasm", "unknown": "not recognised"}
_SEVWORD = {"FAIL": "Detected issue", "WARN": "Check", "INFO": "Note"}
_TILEWORD = {"good": ("✓", "OK"), "warning": ("!", "Check"), "critical": ("✕", "Detected issue"), "info": ("–", "Not checked")}


def _status(findings, ids):
    sev = {f["severity"] for f in findings if f["id"] in ids}
    return "critical" if "FAIL" in sev else "warning" if "WARN" in sev else "good"


def _tiles(res):
    """Overview: six key numbers, each with a status taken from the findings behind it."""
    a, sp, F = res["assembly"], res.get("species"), res["findings"]
    rep, reads = res.get("reports") or {}, res.get("reads") or {}
    ck = rep.get("checkm") if (rep.get("checkm") or {}).get("matches_assembly") is not False else None
    bu = rep.get("busco") if (rep.get("busco") or {}).get("matches_assembly") is not False else None
    depth = [r["est_depth"] for r in reads.values() if r.get("est_depth") is not None]
    stage3 = [k for k, v in GUIDE.items() if v[0] == 3]
    cont_ran = bool(ck or bu or rep.get("identity") or confindr_row(rep) or (rep.get("contig_kraken") or {}).get("contigs")
                    or any("gc_outlier_pct" in r for r in reads.values()))
    cont_st = _status(F, stage3) if cont_ran else "info"
    tiles = [
        ("Genome size", _fmt_len(a["total_length"]),
         f"expected {_fmt_len(sp['expected'])} (normal {_fmt_len(sp['min'])}–{_fmt_len(sp['max'])})" if sp
         else "species not given", _status(F, ["size_range", "gc_range"]) if sp else "info"),
        ("Pieces (contigs)", f"{a['n_contigs']:,}",
         f"longest {_fmt_len(a['largest'])}" + (f", {sum(1 for c in a['contigs'] if c['circular'])} circular"
                                                if any(c["circular"] is not None for c in a["contigs"]) else ""),
         _status(F, ["long_read_fragmented", "short_read_fragmented", "low_n50", "n_content", "no_circular_chromosome"])),
        ("Read depth", f"{max(depth):.0f}×" if depth else "–",
         ", ".join(TECH_NAME.get(t, t) for t in reads) if depth else "no reads given",
         _status(F, ["low_depth"]) if depth else "info"),
        ("Completeness", f"{ck['completeness']:g}%" if ck else f"{bu['complete']:g}%" if bu else "–",
         "CheckM2 estimate" if ck and ck.get("tool") == "CheckM2" else "CheckM estimate" if ck else
         "BUSCO core genes complete" if bu else "no BUSCO/CheckM2 report",
         _status(F, ["busco_low", "checkm_completeness", "busco_fragmented"]) if (ck or bu) else "info"),
        ("Contamination", f"{ck['contamination']:g}%" if ck else
         ("Possible" if cont_st != "good" else "No signs") if cont_ran else "–",
         "CheckM2 estimate" if ck else "from reads, identity and gene checks" if cont_ran else "no checks possible",
         cont_st),
        ("Sequence accuracy", "Errors" if any(f["id"] == "indel_errors" for f in F) else
         "No signs" if ("orf_n50_used" in a or bu or ck) else "–",
         "gene-length check for small errors", _status(F, ["indel_errors"])),
    ]
    return "".join(f'<div class="tile" style="--st:var(--{st})"><div class="tl">{_E(lab)}</div><div class="tv">{_E(v)}</div>'
                   f'<div class="ts">{_E(sub)}</div><div class="chip"><i>{_TILEWORD[st][0]}</i>{_TILEWORD[st][1]}</div></div>'
                   for lab, v, sub, st in tiles)


def _metrics_html(res):
    """All metrics with friendly units; label hover = glossary."""
    a, reads = res["assembly"], res.get("reads") or {}
    lab = lambda k: f'<abbr title="{_E(GLOSSARY[k])}">{_E(k)}</abbr>' if k in GLOSSARY else _E(k)
    table = lambda rows: "<table>" + "".join(f"<tr><th>{lab(k)}</th><td>{_E(v)}</td></tr>" for k, v in rows) + "</table>"
    rows = [("Contigs", f"{a['n_contigs']:,}"), ("Total length", _fmt_len(a["total_length"])),
            ("Largest", _fmt_len(a["largest"])), ("N50 / N90", f"{_fmt_len(a['n50'])} / {_fmt_len(a['n90'])}"),
            ("GC %", f"{a['gc_pct']}%"), ("N %", f"{a['n_pct']}%")]
    if any(c["circular"] is not None for c in a["contigs"]):
        rows.append(("Circular contigs", sum(1 for c in a["contigs"] if c["circular"])))
    if "median_contig_cov" in a:
        rows.append(("Median contig cov (header)", f"{a['median_contig_cov']:.3g}×"))
    if "orf_n50_used" in a:
        rows.append(("ORF N50 (frameshift check)", f"{a['orf_n50_used']} amino acids"))
    md = mosdepth_row(res.get("reports") or {}, a)
    if md and md[0]["at"]:
        rows.append(("Depth breadth (mosdepth)", f"{md[0]['at'][10]:g}% ≥10×, {md[0]['at'][20]:g}% ≥20×"))
    out = "<h3>Assembly</h3>" + table(rows)
    for tech, r in reads.items():
        rr = [(k, f(r[key])) for k, key, f in (
            ("Reads", "n_reads", lambda v: f"{v:,}"), ("Bases", "bases", _fmt_len),
            ("Median length", "median_len", _fmt_len), ("Read N50", "read_n50", _fmt_len),
            ("Mean Q (error-averaged)", "mean_q", lambda v: f"Q{v:g}"), ("% bases >=Q20", "pct_q20", lambda v: f"{v:g}%"),
            ("% bases >=Q30", "pct_q30", lambda v: f"{v:g}%"), ("% reads passing fastp", "pct_passed", lambda v: f"{v:g}%"),
            ("% reads with adapter", "adapter_pct", lambda v: f"{v:g}%"),
            ("% long reads that are lambda DCS", "dcs_pct", lambda v: f"{v:g}%"),
            ("% reads with GC outside assembly range", "gc_outlier_pct", lambda v: f"{v:g}%"),
            ("Est. depth (x)", "est_depth", lambda v: f"{v:g}×")) if key in r]
        out += f"<h3>Reads: {_E(TECH_NAME.get(tech, tech))}</h3>" + table(rr)
    return out


def _files_html(res):
    a, reads, rep = res["assembly"], res.get("reads") or {}, res.get("reports") or {}
    li = [f"Assembly: {a['path']}"] + [f"{TECH_NAME.get(t, t)} reads: {', '.join(r.get('files') or [r.get('source', '')])}"
                                     for t, r in reads.items()]
    li += [f"{k} report: {v['file']}" + (" (not for this assembly, ignored)" if v.get("matches_assembly") is False else "")
           for k, v in rep.items() if v.get("file")]
    li += [f"How it was chosen: {n}" for n in res.get("inputs") or []]
    return "<ul class='files'>" + "".join(f"<li>{_E(x)}</li>" for x in li) + "</ul>"


def _fitness_html(res):
    fit = res.get("fitness")
    if not fit:
        return ""
    cls = {"OK": "good", "NOT_CHECKED": "", "CAUTION": "warning", "NOT_RECOMMENDED": "critical"}
    rows = []
    for k, label in ANALYSES.items():
        o = fit[k]
        worst = [x for x in o["reasons"] if x["rating"] == o["rating"]]
        why = "".join(f'<div>{_E(x["text"])} <span class="src">({_E(", ".join(x["basis"]))})</span></div>'
                      for x in worst[:2]) or "No problem found by the checks that ran."
        rest = [x for x in o["reasons"] if x not in worst[:2]]
        more = (f'<details><summary>All reasons ({len(o["reasons"])})</summary>' + "".join(
            f'<div>{_E(RATING_WORD[x["rating"]])}: {_E(x["text"])} <span class="src">({_E(", ".join(x["basis"]))})'
            f'</span></div>' for x in o["reasons"]) + "</details>") if rest else ""
        rows.append(f'<tr><th>{_E(label)}</th><td><span class="pill {cls[o["rating"]]}">{_E(RATING_WORD[o["rating"]])}'
                    f'</span></td><td>{why}{more}</td></tr>')
    nc = sorted({n for o in fit.values() for n in o["not_checked"]})
    return ('<div class="fit"><h2>What this data is good for</h2><p class="sub">Ratings follow the EU reference '
            'laboratories\' dry-lab QC guideline (2025); where a study found data below it still usable, the rating '
            'is "use with caution" and both are named.</p><table>' + "".join(rows) + "</table>"
            + (f'<p class="sub">Not checked (ratings cover only the checks that ran): {_E("; ".join(nc))}</p>' if nc else "")
            + "</div>")


def html_report(res, md):
    a, sp, reads = res["assembly"], res.get("species"), res.get("reads") or {}
    F = res["findings"]
    st, ic = _STATUS[res["verdict"]]
    # problems in "Start here" order: severity first, then the guide stage (inputs, data, contamination, fixes)
    steps = start_here(F)
    order = [i for _, ids in steps for i in ids]   # cards follow the "Start here" plan
    probs = sorted((f for f in F if f["severity"] != "INFO"),
                   key=lambda f: order.index(f["id"]) if f["id"] in order else len(order))
    notes = [f for f in F if f["severity"] == "INFO"]
    title = lambda i: TITLES.get(i, i.replace("_", " ").capitalize())

    def card(f):
        stc, icon = _STATUS[f["severity"]]
        plain = GUIDE.get(f["id"], (0, "", ""))[2]
        return (f'<div class="find" id="f-{_E(f["id"])}" style="--st:var(--{stc})"><div class="tag"><i>{icon}</i>'
                f'{_SEVWORD[f["severity"]]}</div><h3>{_E(title(f["id"]))}</h3>'
                + (f'<p>{_E(plain)}</p>' if plain else "")
                + (f'<p class="act"><b>In the lab:</b> {_E(acts_of(f)[0])}</p>'
                   f'<p class="act"><b>In bioinformatics:</b> {_E(acts_of(f)[1])}</p>'
                   if acts_of(f) and f["severity"] != "INFO" else
                   f'<p class="act"><b>What to do:</b> {_E(f["action"])}</p>')
                + f'<details><summary>Technical details</summary><p>{_E(f["message"])}</p>'
                f'<p class="sub">Check id: <code>{_E(f["id"])}</code></p></details></div>')

    links = lambda ids: ", ".join(f'<a href="#f-{_E(i)}">{_E(title(i))}</a>' for i in ids)

    byid = {f["id"]: f for f in F}

    def todo(k):
        """Start-here steps with the lab (k=0) or bioinformatics (k=1) title and action of each finding; re-run steps stay plain."""
        items = []
        for t, ids in steps:
            acts = [] if t.startswith("Then re-run") else list(dict.fromkeys(acts_of(byid[i])[k] for i in ids if acts_of(byid[i])))
            items.append(f'<li><b>{_E(BIO_STEPS.get(t, t) if k else t)}</b>' + "".join(f'<div class="todo">{_E(x)}</div>' for x in acts)
                         + (f'<div class="sub">{links(ids)}</div>' if ids else "") + "</li>")
        return f'<ol>{"".join(items)}</ol>'
    charts = [_chart_cumulative(a, sp), _chart_coverage(a)]
    for tech, r in reads.items():
        charts += [_chart_read_length(TECH_NAME.get(tech, tech), r), _chart_read_gc(TECH_NAME.get(tech, tech), r, a)]
    sentence = ("No problems found by the checks." if not probs else
                f"{len(probs)} thing{'s' * (len(probs) > 1)} to check" +
                (f", starting with: {steps[0][0][0].lower() + steps[0][0][1:]}." if steps else "."))
    species = (f"<i>{_E(sp['species'])}</i> "
               + ('<span class="pill good">confirmed</span>' if sp.get("verified_by") else
                  '<span class="pill critical">conflict</span>' if sp.get("identity_conflict") else
                  '<span class="pill">from identity report</span>' if "report" in (res.get("species_match") or "") else
                  '<span class="pill">not verified</span>') if sp else "not given")
    meta = (f"<span><b>Species</b> {species}</span><span><b>Reads</b> "
            f"{_E(', '.join(TECH_NAME.get(t, t) for t in reads) or 'not given')}</span>"
            f"<span><b>Assembler</b> {_E(ASM_NAME.get(a['assembler'], a['assembler'].capitalize()))}</span>")
    tool_md = next((sec for sec in md.split("\n## ") if sec.startswith("Tool reports")), "")
    tool_md = "\n".join(l for l in tool_md.splitlines()[1:] if "_file_" not in l)
    sample = res.get("sample") or Path(a["path"]).name
    return (f'<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,'
            f'initial-scale=1"><title>{_E(sample)} · Takivisor</title><style>{_CSS}</style></head><body><main>'
            f'<p class="sub">Takivisor QC report · {datetime.date.today():%d %b %Y}</p><h1>{_E(sample)}</h1>'
            f'<div class="banner" style="--st:var(--{st})"><div class="ic">{ic}</div><div><b>{shown(res["verdict"])}</b>'
            f'<div>{_E(sentence)}</div></div></div><div class="meta">{meta}</div>'
            f'<div class="tiles">{_tiles(res)}</div>'
            + _fitness_html(res)
            + ('<div class="start"><h2>Start here</h2><div class="tabs">'
               '<input type="radio" name="st" id="st-lab" checked><label for="st-lab">In the lab</label>'
               '<input type="radio" name="st" id="st-bio"><label for="st-bio">In bioinformatics</label>'
               f'<div class="tab" id="tab-lab">{todo(0)}</div><div class="tab" id="tab-bio">{todo(1)}</div></div></div>'
               if steps else "")
            + (f'<h2>What was found</h2><div class="finds">{"".join(card(f) for f in probs)}</div>' if probs else "")
            + f'<h2>Charts</h2><div class="charts">{"".join(c for c in charts if c)}</div>'
            + "<h2>More details</h2>"
            + (f'<details class="sec"><summary>Notes ({len(notes)})</summary><div class="finds">'
               f'{"".join(card(f) for f in notes)}</div></details>' if notes else "")
            + f'<details class="sec"><summary>All metrics</summary>{_metrics_html(res)}</details>'
            + (f'<details class="sec"><summary>Tool reports used</summary>{_md_html(tool_md)}</details>'
               if tool_md.strip() else "")
            + f'<details class="sec"><summary>Files used</summary>{_files_html(res)}</details>'
            + '<details class="sec"><summary>Glossary</summary><table>'
            + "".join(f"<tr><th>{_E(k)}</th><td>{_E(v)}</td></tr>" for k, v in GLOSSARY.items()) + "</table></details>"
            + f'<p class="sub foot">Takivisor v{res["version"]} · research use only, not for clinical decisions.</p>'
            + f'</main><div id="tip" role="tooltip"></div><script>{_JS}</script></body></html>')


# ---------------- demo (--demo) ----------------
def make_demo(root):
    """Demo run folder: the real M. genitalium G37 genome as the 'assembly' plus SYNTHETIC ONT-like reads
    simulated here (15x, ~5% errors) with 3% fake contaminant reads at ~55% GC, so the report shows a
    depth warning, a contamination step and a second hump in the read-GC chart. Seeded: same every time."""
    import random
    rnd = random.Random(1)
    root = Path(root)
    (root / "assembly").mkdir(parents=True, exist_ok=True)
    (root / "reads").mkdir(exist_ok=True)
    src = _data("demo/mgenitalium_G37.fna.gz")
    with gzip.open(src, "rt") as f:
        text = f.read()
    (root / "assembly" / "assembly.fasta").write_text(text)
    genome = "".join(l for l in text.splitlines() if not l.startswith(">"))
    g = len(genome)
    (root / "assembly" / "assembly_info.txt").write_text(f"#seq_name\tlength\tcov.\tcirc.\trepeat\tmult.\talt_group\t"
                                                         f"graph_path\ncontig_1\t{g}\t15\tY\tN\t1\t*\t1\n")
    comp = str.maketrans("ACGT", "TGCA")

    def noisy(s):   # ~5% errors: substitutions, insertions, deletions in equal parts
        out = []
        for c in s:
            r = rnd.random()
            if r < 0.017:
                out.append(rnd.choice("ACGT"))
            elif r < 0.033:
                out += [c, rnd.choice("ACGT")]
            elif r >= 0.05:
                out.append(c)
        return "".join(out)

    with gzip.open(root / "reads" / "demo_ont.fastq.gz", "wt") as fq:
        bases, i = 0, 0
        while bases < 15 * g:
            n = min(int(rnd.lognormvariate(math.log(3000), 0.6)) + 200, 40_000)   # >= 2000 reads for read GC
            if rnd.random() < 0.03:   # fake contaminant: random sequence at ~55% GC
                s = "".join(rnd.choices("GCAT", weights=(27.5, 27.5, 22.5, 22.5), k=n))
            else:
                p = rnd.randrange(g)
                s = (genome + genome)[p:p + n]
                s = s if rnd.random() < 0.5 else s.translate(comp)[::-1]
            s = noisy(s)
            q = "".join(chr(33 + min(50, max(2, int(rnd.gauss(16, 6))))) for _ in s)
            fq.write(f"@demo_read_{i} synthetic\n{s}\n+\n{q}\n")
            bases += len(s); i += 1
    return root


# ---------------- run-folder discovery (takivisor RUN_FOLDER) ----------------
# Rules set from real run folders (TakiLine ONT/Illumina/polypolish runs, abby Flye+Racon): assembler
# work files and QC-tool FASTAs are skipped; polished/final names beat raw ones; filtered reads beat raw.
FASTA_EXT = (".fasta", ".fa", ".fna", ".fasta.gz", ".fa.gz", ".fna.gz")
FASTQ_EXT = (".fastq", ".fq", ".fastq.gz", ".fq.gz")
SKIP_DIR = re.compile(r"^(\d\d-[a-z]+|K\d+|misc|tmp|temp|busco.*|.*_busco|quast.*|checkm.*|prodigal_output|"
                      r"predicted_genes|run_.*|kraken.*|bakta.*|prokka.*|mob_?suite.*|plasmidfinder|abricate.*|"
                      r"amrfinder.*|mmseqs.*|\..*)$", re.I)
SKIP_FASTA = re.compile(r"before_rr|graph|draft_assembly|polished_\d|filtered_contigs|broken_scaffolds|"
                        r"repeat|predicted|unclassified|classified|resistance|plasmid", re.I)
FINAL = re.compile(r"final|consensus", re.I)
POLISHED = re.compile(r"polish|medaka|pilon|pypolca|racon|homopolish|dorado_polish|nextpolish", re.I)
PROCESSED = re.compile(r"filt|trim|clean|qc|chopper|fastp|porechop|nanofilt|filtlong|hq", re.I)
MAX_ASM_CANDIDATES = 8   # more than this: the folder is probably several runs


def _files(root, exts):
    out = []
    for p in sorted(Path(root).rglob("*")):
        rel = p.relative_to(root)
        if any(SKIP_DIR.match(part) for part in rel.parts[:-1]):
            continue
        if p.name.lower().endswith(exts) and p.is_file():
            out.append(p)
    return out


def _is_long(path, n=100):
    lens, it = [], _lines(path)
    for h in it:
        lens.append(len(next(it, b"").rstrip())); next(it, None); next(it, None)
        if len(lens) >= n:
            break
    return bool(lens) and sorted(lens)[len(lens) // 2] >= 500


def discover(root):
    """Pick the assembly, reads and Flye info inside a run folder. Returns (choices, notes)."""
    root = Path(root)
    fa = [p for p in _files(root, FASTA_EXT) if not SKIP_FASTA.search(p.name)
          and p.stat().st_size >= (20_000 if p.name.endswith(".gz") else 50_000)]
    if not fa:
        sys.exit(f"{root}: no assembly FASTA found (looked for {', '.join(FASTA_EXT)} >= 50 kb); pass --assembly")
    if len(fa) > MAX_ASM_CANDIDATES:
        sys.exit(f"{root}: {len(fa)} assembly-like FASTA files, so this looks like several runs. Point Takivisor "
                 f"at one run's folder, or pass --assembly. First ones: "
                 + ", ".join(str(p.relative_to(root)) for p in fa[:5]))
    rank = lambda p: (0 if FINAL.search(p.name) else 1 if POLISHED.search(str(p.relative_to(root))) else 2,
                      -p.stat().st_mtime)
    fa.sort(key=rank)
    asm, notes = fa[0], []
    why = ("final/consensus name" if rank(asm)[0] == 0 else "polished" if rank(asm)[0] == 1 else "newest")
    notes.append(f"assembly: {asm.relative_to(root)} ({why}"
                 + (f"; other candidates: {', '.join(str(p.relative_to(root)) for p in fa[1:])}" if len(fa) > 1 else "")
                 + ")")
    reads, skipped = [], []
    groups = {}
    for p in _files(root, FASTQ_EXT):
        groups.setdefault(_is_long(p), []).append(p)
    for long_, ps in groups.items():
        proc = [p for p in ps if PROCESSED.search(str(p.relative_to(root)))]
        keep = proc or ps   # filtered/trimmed beat raw, so the same reads are not counted twice
        reads += keep
        skipped += [p for p in ps if p not in keep]
    if reads:
        notes.append("reads: " + ", ".join(str(p.relative_to(root)) for p in reads)
                     + (f" (skipped raw: {', '.join(str(p.relative_to(root)) for p in skipped)})" if skipped else ""))
    infos = sorted((p for p in root.rglob("assembly_info.txt")),
                   key=lambda p: (p.parent != asm.parent, len(p.relative_to(root).parts)))
    flye = infos[0] if infos else None
    if flye:
        notes.append(f"Flye info: {flye.relative_to(root)}")
    return {"assembly": str(asm), "reads": [str(p) for p in reads], "flye_info": str(flye) if flye else None}, notes


# ---------------- report ----------------

# WARNs about the read data alone: they explain problems but are not problems with the genome, so
# on their own they don't turn a PASS into REVIEW (real MiSeq runs at Q30 51-73% can still give a
# good assembly). Low depth is NOT here: it can hide missing genome, especially without CheckM2.
DATA_ONLY_WARNS = {"illumina_q30", "fastp_low_pass", "ont_short_reads"}


def verdict(findings):
    sev = {f["severity"] for f in findings}
    if "FAIL" in sev:
        return "FAIL"
    return "REVIEW" if any(f["severity"] == "WARN" and f["id"] not in DATA_ONLY_WARNS for f in findings) else "PASS"


def _fitness_md(res):
    fit = res.get("fitness")
    if not fit:
        return []
    L = ["", "## What this data is good for", "", "| Analysis | Rating | Why (basis) |", "|---|---|---|"]
    for k, label in ANALYSES.items():
        o = fit[k]
        worst = [x for x in o["reasons"] if x["rating"] == o["rating"]][:2]
        why = "; ".join(f"{x['text']} ({', '.join(x['basis'])})" for x in worst) or "No problem found by the checks that ran."
        L.append(f"| {label} | {RATING_WORD[o['rating']]} | {why} |")
    nc = sorted({n for o in fit.values() for n in o["not_checked"]})
    if nc:
        L += ["", "Not checked (ratings cover only the checks that ran): " + "; ".join(nc)]
    return L


def markdown(res):
    a, reads, sp = res["assembly"], res.get("reads") or {}, res.get("species")
    techs = ", ".join(f"**{t}** ({r['tech_basis']})" for t, r in reads.items()) or "unknown (no reads/reports given)"
    L = [f"# Takivisor: {Path(a['path']).name}", "",
         f"**Verdict: {shown(res['verdict'])}**" + (
             f" ({n_data} read-data warning{'s' * (n_data > 1)}; assembly checks passed)"
             if res["verdict"] == "PASS" and (n_data := sum(f["severity"] == "WARN" for f in res["findings"])) else ""),
         "",
         *(f"- Auto-detected {n}" for n in res.get("inputs") or []),
         f"- Tech: {techs}" + (" → hybrid" if len(reads) > 1 else ""),
         f"- Assembler: **{a['assembler']}** ({a['assembler_confidence']:.0%} of headers match)"
         + (f", polishing: {', '.join(a['polishing_hints'])}" if a["polishing_hints"] else ""),
         f"- Species: " + (f"{sp['species']} (expected {sp['expected']/1e6:.2f} Mb; {res['species_match']})"
                           + (f"; identity confirmed by {sp['verified_by']}" if sp.get("verified_by") else
                              f"; identity CONFLICT: {sp['identity_conflict']}" if sp.get("identity_conflict") else
                              "" if "report" in (res.get("species_match") or "") else
                              "; identity not verified (identity report inconclusive)" if (res.get("reports") or {}).get("identity") else
                              "; identity not verified (no Kraken2/sylph/GTDB-Tk report)")
                           if sp else res.get("species_match") or "not given"),
         *(["", "## Start here", ""] + [f"{i}. {t}" + (f" ({', '.join(ids)})" if ids else "")
                                        for i, (t, ids) in enumerate(start_here(res["findings"]), 1)]
           if start_here(res["findings"]) else []),
         *_fitness_md(res),
         "", "## Assembly", "", "| Metric | Value |", "|---|---|",
         f"| Contigs | {a['n_contigs']} |", f"| Total length | {a['total_length']:,} |",
         f"| Largest | {a['largest']:,} |", f"| N50 / N90 | {a['n50']:,} / {a['n90']:,} |",
         f"| GC % | {a['gc_pct']} |", f"| N % | {a['n_pct']} |"]
    circ = sum(1 for c in a["contigs"] if c["circular"])
    if any(c["circular"] is not None for c in a["contigs"]):
        L.append(f"| Circular contigs | {circ} |")
    if "median_contig_cov" in a:
        L.append(f"| Median contig cov (header) | {a['median_contig_cov']} |")
    if "orf_n50_used" in a:
        L.append(f"| ORF N50 (frameshift check) | {a['orf_n50_used']} aa |")
    for tech, r in reads.items():
        L += ["", f"## Reads: {tech} ({r['source']})", "", "| Metric | Value |", "|---|---|"]
        L += [f"| {label} | {r[k]:,} |" for label, k in (
            ("Reads", "n_reads"), ("Bases", "bases"), ("Median length", "median_len"),
            ("Read N50", "read_n50"), ("Max length", "max_len"),
            ("Mean Q (error-averaged)", "mean_q"), ("Mean read Q (NanoPlot, per-read average)", "mean_read_q"),
            ("% bases >=Q20", "pct_q20"), ("% bases >=Q30", "pct_q30"),
            ("% reads passing fastp", "pct_passed"), ("Duplication % (fastp)", "dup_pct"),
            ("% reads with adapter", "adapter_pct"), ("% long reads that are lambda DCS", "dcs_pct"), ("% reads with GC outside assembly range", "gc_outlier_pct"),
            ("Est. depth (x)", "est_depth")) if k in r]
    rep = res.get("reports") or {}
    if rep:
        L += ["", "## Tool reports", "", "| Tool | Result |", "|---|---|"]
        if "busco" in rep:
            b = rep["busco"]
            L.append(f"| BUSCO ({b['lineage']}) | C:{b['complete']}% [S:{b['single']}%, D:{b['duplicated']}%], "
                     f"F:{b['fragmented']}%, M:{b['missing']}% |")
        if "checkm" in rep:
            c = rep["checkm"]
            extra = "".join(f", {lbl} {c[k]}" for lbl, k in (
                ("coding density", "coding_density"), ("mean gene aa", "avg_gene_length_aa"),
                ("strain heterogeneity", "strain_heterogeneity")) if c.get(k) is not None)
            L.append(f"| {c['tool']} | completeness {c['completeness']}%, contamination {c['contamination']}%{extra} |")
        if "identity" in rep:
            i = rep["identity"]
            L.append(f"| {i['tool']} (identity) | {i['top']} {i['top_share']}%"
                     + (f", then {i['second']} {i['second_share']}%" if i["second"] else "")
                     + " of bacterial species-level assignments |")
        if (rep.get("contig_kraken") or {}).get("contigs"):
            cs = rep["contig_kraken"]["contigs"]
            L.append(f"| Kraken2 (contigs) | {len(cs)} contigs; top species by length: "
                     + ", ".join(f"{c['name']} {c['top']} {c['share']}%" for c in sorted(cs, key=lambda c: -c["length"])[:3]) + " |")
        if "confindr" in rep and not confindr_row(rep):
            L.append(f"| ConFindr | {rep['confindr']['n_candidates']} samples, none named like this assembly (not used) |")
        elif "confindr" in rep:
            c = rep["confindr"]
            L.append(f"| ConFindr | {c['name']}: contamination "
                     + {True: "detected", False: "not detected", None: "not determined"}[c["contaminated"]]
                     + (f", {c['snvs']} SNVs" if c["snvs"] is not None else "")
                     + (f", ~{c['pct']:g}%" if c["pct"] is not None else "") + " |")
        if rep.get("mosdepth"):
            m = rep["mosdepth"]
            L.append(f"| mosdepth | {m['name']}: mean {m['mean']:g}x"
                     + (f", {m['at'][10]:g}% >= 10x, {m['at'][20]:g}% >= 20x" if m["at"] else ", no global.dist.txt") + " |")
        for key in ("amrfinder", "rgi"):
            if key in rep:
                r = rep[key]
                n_bad = sum(1 for h in r["hits"] if h["issue"])
                L.append(f"| {r['tool']} | {len(r['hits'])} gene calls, {n_bad} incomplete |")
        if "quast" in rep:
            q = rep["quast"]
            extra = "".join(f", {lbl} {q[k]:g}" for lbl, k in (
                ("rRNA genes", "rrna_genes"), ("genome fraction %", "genome_fraction"),
                ("misassemblies", "misassemblies")) if q.get(k) is not None)
            L.append(f"| QUAST | {q['n_contigs']} contigs, N50 {q['n50']:,}{extra} |")
        for v in rep.values():
            note = f" (best of {v['n_candidates']} candidates)" if v["n_candidates"] > 1 else ""
            note += " **not this assembly**" if v.get("matches_assembly") is False else ""
            L.append(f"| _file_ | {v['file']}{note} |")
    # Problems first; INFO notes go in their own compact section so they don't bury WARN/FAIL.
    problems = [f for f in res["findings"] if f["severity"] != "INFO"]
    notes = [f for f in res["findings"] if f["severity"] == "INFO"]
    L += ["", "## Findings", ""]
    if not problems:
        L.append("No problems found by the checks in this version.")
    plain = lambda f: f"  \n  In plain words: {GUIDE[f['id']][2]}" if f["id"] in GUIDE else ""
    for f in problems:
        acts = (res.get("actions") or {}).get(f["id"])
        L.append(f"- **{shown(f['severity'])}** `{f['id']}`: {f['message']}{plain(f)}  \n  → {f['action']}"
                 + (f"  \n  In the lab: {acts['lab']}  \n  In bioinformatics: {acts['bioinformatics']}" if acts else ""))
    if notes:
        L += ["", "## Notes", ""] + [f"- `{f['id']}`: {f['message']} {f['action']}" for f in notes]
    return "\n".join(L) + "\n"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("folder", nargs="?", help="a pipeline run folder: the assembly, reads, reports and Flye info "
                    "are found in it (what was picked is printed; any option below overrides it)")
    ap.add_argument("--assembly", help="assembly FASTA (gz ok); required unless a run folder is given")
    ap.add_argument("--reads", nargs="+", default=[],
                    help="FASTQ file(s), any mix of Illumina R1/R2 and ONT/HiFi; tech is detected per file")
    ap.add_argument("--species", help='e.g. "Klebsiella pneumoniae", "E. coli"')
    ap.add_argument("--tech", choices=["illumina", "ont", "hifi", "clr"], help="force tech for all --reads files")
    ap.add_argument("--assembler", help="override detection")
    ap.add_argument("--flye-info", help="Flye assembly_info.txt (circularity + coverage)")
    ap.add_argument("--sample", help="sample name used to match reports that carry no genome length "
                                     "(Kraken2, legacy CheckM); default: the assembly file name unless generic")
    ap.add_argument("--reports", nargs="+", default=[],
                    help="folder(s) searched (recursively) for BUSCO, CheckM2/CheckM, QUAST, fastp, NanoPlot, "
                         "Kraken2 / sylph-tax / GTDB-Tk reports, and Kraken2 per-contig output (+ its DB folder "
                         "for ktaxonomy.tsv)")
    ap.add_argument("--json", help="write full results JSON here")
    ap.add_argument("--html", help="where to write the HTML report with charts (default: takivisor_report.html "
                    "in the run folder, or <assembly name>.takivisor.html next to the assembly)")
    ap.add_argument("--no-html", action="store_true", help="don't write the HTML report")
    ap.add_argument("--demo", action="store_true", help="create a demo run folder (real M. genitalium genome + "
                    "synthetic reads) in FOLDER (default ./takivisor_demo) and review it")
    ap.add_argument("--exit-code", action="store_true",
                    help="exit with the verdict: 0 PASS, 3 REVIEW, 4 FAIL (1/2 = input/argument errors)")
    args = ap.parse_args(argv)
    found = []
    if args.demo:
        args.folder = str(make_demo(args.folder or "takivisor_demo"))
        args.species, args.tech = args.species or "Mycoplasmoides genitalium", args.tech or "ont"
        print(f"takivisor: demo data written to {args.folder} (real genome, SYNTHETIC reads with 3% fake "
              "contaminant reads)", file=sys.stderr)
    if args.folder:
        if not Path(args.folder).is_dir():
            ap.error(f"not a folder: {args.folder}")
        pick, found = discover(args.folder)
        args.assembly = args.assembly or pick["assembly"]
        args.reads = args.reads or pick["reads"]
        args.flye_info = args.flye_info or pick["flye_info"]
        args.reports = args.reports or [args.folder]
        found.append(f"reports: searched {', '.join(args.reports)}")
        for n in found:
            print(f"takivisor: using {n}", file=sys.stderr)
    elif not args.assembly:
        ap.error("give a run folder (takivisor RUN_FOLDER) or --assembly FILE")
    for f in [args.assembly, args.flye_info, *args.reads, *args.reports]:
        if f and not Path(f).exists():
            ap.error(f"not found: {f}")

    t0 = time.monotonic()
    print(f"takivisor: checking assembly {Path(args.assembly).name}...", file=sys.stderr, flush=True)
    asm = assembly_metrics(args.assembly, args.flye_info)
    print(f"takivisor: assembly done: {asm['n_contigs']:,} contigs, {_fmt_len(asm['total_length'])} "
          f"in {time.monotonic() - t0:.0f} s", file=sys.stderr, flush=True)
    if args.assembler:
        asm["assembler"], asm["assembler_confidence"] = args.assembler, 1.0
    asm_name = identifying(args.sample or Path(args.assembly).name.split(".")[0])
    rep = load_reports(args.reports, asm["total_length"], asm_name)
    if args.reports and not rep:
        print(f"warning: no QC reports recognised in {', '.join(args.reports)} (looked for BUSCO short_summary*.txt, "
              "CheckM2 quality_report.tsv, QUAST report.tsv, fastp *.json, NanoPlot NanoStats.txt, Kraken2/sylph/"
              "GTDB-Tk, ConFindr, AMRFinderPlus/RGI, mosdepth reports); continuing with the assembly and reads only.", file=sys.stderr)
    reads = read_sets(args.reads, args.tech) if args.reads else {}
    # read stats from reports (raw, before filtering) for techs not given as FASTQ
    if "fastp" in rep and "illumina" not in reads:
        reads["illumina"] = {**rep["fastp"], "source": "fastp report", "tech_basis": "fastp report"}
    elif "fastp" in rep:
        reads["illumina"].update(pct_passed=rep["fastp"]["pct_passed"], bases_after=rep["fastp"]["bases_after"])
    long_tech = args.tech if args.tech in LONG_TECHS else "ont"
    if "nanostats" in rep and not reads.keys() & LONG_TECHS:
        reads[long_tech] = {**rep["nanostats"], "source": "NanoPlot report",
                            "tech_basis": "NanoPlot report" + ("" if args.tech else "; pass --tech hifi if PacBio")}

    sp, how = species_range(args.species) if args.species else (None, None)
    if args.species and not sp:
        print(f"warning: species '{args.species}': {how}; size/GC checks skipped", file=sys.stderr)
    idr = rep.get("identity")
    if (not args.species and idr and idr["top_share"] >= T["id_min_share"]
            and (idr.get("n_assigned") is None or idr["n_assigned"] >= T["id_min_reads"])):
        sp, why = species_range(idr["top"])
        how = (f"from {idr['tool']} report, {idr['top_share']}% of species-level assignments"
               if sp else f"{idr['tool']} top species '{idr['top']}' {why}")

    if "contig_kraken" in rep:
        ckr, tax = rep["contig_kraken"], load_ktaxonomy(args.reports)
        ckr["contigs"] = contig_species(ckr.pop("rows"), tax) if tax else None
    findings = run_checks(asm, reads, sp, rep)
    if not args.species and not sp:   # beginners often skip --species; say what that costs
        findings.append({"severity": "INFO", "id": "species_not_given",
                         "message": "No species given and no usable identity report, so genome size and GC "
                                    "were not checked against the species.",
                         "action": 'Add --species "Genus species" (e.g. --species "Escherichia coli"), or give '
                                   "a Kraken2/sylph/GTDB-Tk report."})
    if args.species and not sp:  # #15: a stderr line was easy to miss and the run could still PASS
        sug = species_suggestions(args.species)
        findings.insert(0, {"severity": "WARN", "id": "species_unknown",
                            "message": f"--species '{args.species}' {how}; genome size and GC checks were skipped.",
                            "action": ("Did you mean: " + ", ".join(sug) + "? " if sug else "")
                                      + "Re-run with the exact NCBI species name."})
    sample = args.sample or (Path(args.folder).resolve().name if args.folder else Path(args.assembly).name.split(".")[0])
    # JSON contract (README "Output for pipelines"): schema, version, sample, verdict and findings[] with
    # severity/id/message/action are stable within a schema number; everything else is informational.
    res = {"schema": 1, "version": "1.0", "sample": sample, "inputs": found, "species": sp, "species_match": how, "assembly": asm,
           "reads": reads, "reports": rep, "findings": findings, "verdict": verdict(findings)}
    res["fitness"] = fitness(res)
    res["actions"] = {f["id"]: dict(zip(("lab", "bioinformatics"), acts_of(f)))
                      for f in findings if acts_of(f) and f["severity"] != "INFO"}
    md = markdown(res)
    print(md, end="")
    if not args.no_html:
        out = Path(args.html) if args.html else (Path(args.folder) / "takivisor_report.html" if args.folder else
                                                 Path(args.assembly).with_name(Path(args.assembly).name.split(".")[0]
                                                                               + ".takivisor.html"))
        try:
            out.write_text(html_report(res, md))
            print(f"takivisor: HTML report written to {out}", file=sys.stderr)
        except OSError as e:
            print(f"warning: could not write the HTML report ({e}); use --html PATH", file=sys.stderr)
    for r in reads.values():
        r.pop("read_gc", None)  # internal histogram with tuple keys; not for JSON
    if args.json:
        out_json = {**res, "findings": [{k: v for k, v in f.items() if k != "acts"} for f in findings]}   # acts -> "actions"
        Path(args.json).write_text(json.dumps(out_json, indent=1, default=sorted))  # sets (species lineage)
    return EXIT_CODES[res["verdict"]] if args.exit_code else 0


if __name__ == "__main__":
    sys.exit(main())
