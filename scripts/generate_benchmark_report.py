"""Generate the InstaNovo benchmark report markdown.

Reads results/<dataset>/peptide_precision_plot_data.csv and AA_precision_plot_data.csv,
ranks algorithms by AUC per dataset, and emits the report grouped into the same
sections used in the existing benchmark_report.md.
"""

import os
import csv
import glob
from collections import defaultdict
from datetime import date


RESULTS_DIR = "results"

# Section title -> emoji + ordered list of dataset names (must match dir names exactly)
SECTIONS = [
    ("🧬 Nine Species Benchmark", "Cross-species generalization on tryptic proteomics data from diverse organisms — the most widely used de novo sequencing benchmark.", [
        "9_species_apis_mellifera",
        "9_species_bacillus_subtilis",
        "9_species_candidatus_thiodiazotropha",
        "9_species_human",
        "9_species_methanosarcina_mazei",
        "9_species_mus_musculus",
        "9_species_saccharomyces_cerevisiae",
        "9_species_solanum_lycopersicum",
        "9_species_vigna_mungo",
    ]),
    ("📋 CAMPI Benchmark", "Community benchmarking standards from the Critical Assessment of Mass spectrometry-based Proteomics Informatics.", [
        "CAMPI_F07",
        "CAMPI_S11",
    ]),
    ("🔬 HeLa Single Cell", "Single-cell proteomics with limited input material.", [
        "hela_single_cell",
        "hela_single_cell_2",
    ]),
    ("🧑‍🔬 Human Standard Proteomics", "Standard human proteomics across instrument platforms and labeling strategies.", [
        "human_TMT",
        "human_agilent",
        "human_astral",
        "human_jurkat",
    ]),
    ("🛡️ Immunopeptidomics (HLA / MHC)", "Non-tryptic immunopeptides presented by HLA/MHC molecules. Critical for neoantigen discovery and immunotherapy.", [
        "human_HLAI",
        "human_HLAII",
        "human_MHCI",
        "mouse_MHCI",
        "mouse_MHCII",
    ]),
    ("🧪 Human Phosphoproteomics", "Phosphopeptide-enriched datasets with SILAC labeling.", [
        "human_phospho_Phos_SILAC",
        "human_phospho_pY_SILAC",
        "human_phospho_superSILAC",
    ]),
    ("⚠️ Human Multiprotease (Standard)", "Multi-enzyme digestion datasets using non-tryptic proteases (aspn, chymotrypsin, gluc, lysc, lysn) and trypsin.", [
        "human_multiprotease_aspn",
        "human_multiprotease_chymotrypsin",
        "human_multiprotease_gluc",
        "human_multiprotease_lysc",
        "human_multiprotease_lysn",
        "human_multiprotease_trypsin",
    ]),
    ("✨ Human Multiprotease with PTMs", "Multi-enzyme datasets combined with post-translational modifications.", [
        "human_multiprotease_ptm_argc",
        "human_multiprotease_ptm_aspn",
        "human_multiprotease_ptm_chymotrypsin",
        "human_multiprotease_ptm_gluc",
        "human_multiprotease_ptm_lysc",
        "human_multiprotease_ptm_trypsin",
    ]),
    ("💉 Human Herceptin (Antibody Digestion)", "Herceptin (trastuzumab) antibody digested with multiple enzymes.", [
        "human_herceptin_aspn",
        "human_herceptin_chymotrypsin",
        "human_herceptin_gluc",
        "human_herceptin_lysc",
        "human_herceptin_lysn",
        "human_herceptin_trypsin",
    ]),
    ("🧫 Monoclonal Antibody (mAb) Datasets", "Human and mouse monoclonal antibody samples with various proteases.", [
        "human_mAb_aspn",
        "human_mAb_chymotrypsin",
        "human_mAb_trypsin",
        "mouse_mAb_aspn",
        "mouse_mAb_chymotrypsin",
        "mouse_mAb_trypsin",
    ]),
    ("📊 Label-Free Quantification (LFQ)", "Datasets from Orbitrap, SCIEX, and timsTOF instruments.", [
        "LFQ_orbitrap",
        "LFQ_sciex",
        "LFQ_timstof",
    ]),
    ("⚗️ ProteomeTools — Orbitrap", "Synthetic peptide library datasets from the ProteomeTools project on Orbitrap instruments. Includes tryptic, non-tryptic, HLA, and TMT-labeled subsets.", [
        "PT_orbitrap_HLAI",
        "PT_orbitrap_HLAII",
        "PT_orbitrap_HLAI_TMT",
        "PT_orbitrap_aspn",
        "PT_orbitrap_aspn_TMT",
        "PT_orbitrap_lysn",
        "PT_orbitrap_lysn_TMT",
        "PT_orbitrap_tryptic_1",
        "PT_orbitrap_tryptic_2",
        "PT_orbitrap_tryptic_TMT",
    ]),
    ("🔧 ProteomeTools — timsTOF", "Synthetic peptide library datasets on timsTOF instruments.", [
        "PT_timstof_HLAI",
        "PT_timstof_HLAII",
        "PT_timstof_aspn",
        "PT_timstof_lysn",
        "PT_timstof_tryptic",
    ]),
    ("✏️ 21 Post-Translational Modifications", "Datasets testing recognition of 21 different PTMs, grouped by modified residue (K, P, R, Y).", [
        "21PTMs_Kmod",
        "21PTMs_Pmod",
        "21PTMs_Rmod",
        "21PTMs_Ymod",
    ]),
    ("🛠️ Multitool (Search Engine Comparison)", "Datasets for comparing de novo sequencing tools across Orbitrap and SCIEX platforms.", [
        "multitool_orbitrap",
        "multitool_sciex",
    ]),
    ("❓ iPRG Benchmark", "The iPRG (Proteome Informatics Research Group) benchmark datasets.", [
        "iPRG",
        "iPRG_2015",
    ]),
    ("🌍 Non-Human Organisms (Animal, Plant, Fungus, Bacteria)", "Cross-organism generalization tests.", [
        "animal_invertebrate",
        "animal_mammal_1",
        "animal_mammal_2",
        "fungus",
        "plant",
        "staphylococcus_aureus_timstof",
    ]),
    ("📦 Other Datasets", "", [
        "PT_nonnatural",
        "UPS",
    ]),
]


def load_auc_for_dataset(dataset, metric_file):
    """Return list of (algorithm, version, auc) tuples for a dataset/metric."""
    path = os.path.join(RESULTS_DIR, dataset, metric_file)
    if not os.path.exists(path):
        return []
    rows = []
    with open(path) as f:
        for r in csv.DictReader(f):
            algo = r["algorithm"]
            version = r["version"]
            try:
                auc = float(r["auc"])
            except (ValueError, KeyError):
                continue
            rows.append((algo, version, auc))
    return rows


def best_per_algorithm(entries):
    """Keep the best AUC per algorithm (across versions)."""
    best = {}
    for algo, version, auc in entries:
        if algo not in best or auc > best[algo][1]:
            best[algo] = (version, auc)
    return [(algo, ver, auc) for algo, (ver, auc) in best.items()]


def find_instanovo(entries, version):
    """Return (auc, ) or None for InstaNovo at a specific version."""
    for algo, ver, auc in entries:
        if algo == "instanovo" and ver == version:
            return auc
    return None


def medal(rank):
    return {1: "🥇", 2: "🥈", 3: "🥉"}.get(rank, "")


def rank_marker(rank):
    m = medal(rank)
    if m:
        return f"{m} **#{rank}**"
    return f"#{rank}"


def fmt_auc(x):
    if x is None:
        return "—"
    return f"{x:.4f}"


def fmt_delta(curr, prev):
    if curr is None or prev is None:
        return "—"
    return f"{(curr - prev) * 100:+.1f}%"


def section_table_rows(datasets, metric_file, v122_label, v112_label):
    """Compute table rows for a section/metric.

    Returns list of dicts with keys:
      dataset, v122_auc, v122_rank, n_algos, v112_auc, delta, best_other_algo, best_other_ver, best_other_auc
    """
    out = []
    for ds in datasets:
        entries = load_auc_for_dataset(ds, metric_file)
        if not entries:
            continue
        # All instanovo versions kept separately for v1.1.2 vs v1.2.2 comparison
        v122 = find_instanovo(entries, v122_label)
        v112 = find_instanovo(entries, v112_label)
        # For ranking: keep only best version per algorithm
        best = best_per_algorithm(entries)
        n_algos = len(best)
        # rank by AUC desc
        ranked = sorted(best, key=lambda r: -r[2])
        # Find v1.2.2 rank — look for instanovo entry in best (v1.2.2 should be the best instanovo if it improved)
        v122_rank = None
        for i, (algo, ver, auc) in enumerate(ranked, start=1):
            if algo == "instanovo":
                # Only count this rank if v1.2.2 is the best instanovo version
                if ver == v122_label:
                    v122_rank = i
                break
        # Best non-instanovo competitor
        best_other = None
        for algo, ver, auc in ranked:
            if algo != "instanovo":
                best_other = (algo, ver, auc)
                break
        out.append({
            "dataset": ds,
            "v122_auc": v122,
            "v122_rank": v122_rank,
            "n_algos": n_algos,
            "v112_auc": v112,
            "delta": (v122 - v112) if (v122 is not None and v112 is not None) else None,
            "best_other": best_other,
        })
    return out


def render_table(rows, label):
    lines = []
    lines.append(f"**{label}**\n")
    lines.append("| Dataset | v1.2.2 AUC | Rank | v1.1.2 AUC | Delta | Best Other Algorithm | Best Other AUC |")
    lines.append("|---|---|---|---|---|---|---|")
    for r in rows:
        rank_str = rank_marker(r["v122_rank"]) + f"/{r['n_algos']}" if r["v122_rank"] else "—"
        bo = r["best_other"]
        bo_str = f"{bo[0]} {bo[1]}" if bo else "—"
        bo_auc = fmt_auc(bo[2]) if bo else "—"
        lines.append(
            f"| {r['dataset']} | {fmt_auc(r['v122_auc'])} | {rank_str} | "
            f"{fmt_auc(r['v112_auc'])} | {fmt_delta(r['v122_auc'], r['v112_auc'])} | "
            f"{bo_str} | {bo_auc} |"
        )
    lines.append("")
    return "\n".join(lines)


def main():
    today = date.today().isoformat()

    # Discover overall stats
    all_datasets = sorted(
        d for d in os.listdir(RESULTS_DIR)
        if os.path.isdir(os.path.join(RESULTS_DIR, d))
    )
    # Count algorithms across the full dataset
    all_algos = set()
    max_per_dataset = 0
    for ds in all_datasets:
        entries = load_auc_for_dataset(ds, "peptide_precision_plot_data.csv")
        algos = {a for a, _, _ in entries}
        all_algos.update(algos)
        max_per_dataset = max(max_per_dataset, len(set(a for a, _, _ in best_per_algorithm(entries))))

    # Compute v1.2.2 vs v1.1.2 stats globally and regression list
    n_first_peptide = 0
    n_top3_peptide = 0
    deltas_peptide = []
    regressions = []
    for ds in all_datasets:
        for metric_file, metric_label in [
            ("peptide_precision_plot_data.csv", "Peptide"),
            ("AA_precision_plot_data.csv", "AA"),
        ]:
            entries = load_auc_for_dataset(ds, metric_file)
            if not entries:
                continue
            v122 = find_instanovo(entries, "1.2.2")
            v112 = find_instanovo(entries, "1.1.2")
            if v122 is not None and v112 is not None:
                delta = v122 - v112
                if metric_label == "Peptide":
                    deltas_peptide.append(delta)
                if delta < 0:
                    regressions.append((ds, metric_label, v112, v122, delta))
            if metric_label == "Peptide":
                # Rank check
                best = best_per_algorithm(entries)
                ranked = sorted(best, key=lambda r: -r[2])
                for i, (algo, ver, auc) in enumerate(ranked, start=1):
                    if algo == "instanovo" and ver == "1.2.2":
                        if i == 1:
                            n_first_peptide += 1
                        if i <= 3:
                            n_top3_peptide += 1
                        break

    n_datasets = len(all_datasets)
    mean_peptide_delta = sum(deltas_peptide) / len(deltas_peptide) if deltas_peptide else 0

    # Build output
    out = []
    out.append("# InstaNovo De Novo Benchmark Report\n")
    out.append("**Source:** [bittremieuxlab/denovo_benchmarks](https://github.com/bittremieuxlab/denovo_benchmarks/tree/main/results)  ")
    out.append(f"**Date:** {today}  ")
    out.append("**Metrics:** Area Under the Precision-Coverage Curve (AUC) for amino acid (AA) precision and peptide-level precision.  ")
    out.append(f"**Total datasets evaluated:** {n_datasets}  ")
    out.append(f"**Algorithms compared:** Up to {max_per_dataset} per dataset  ")
    out.append(f"**Algorithm pool:** {', '.join(sorted(all_algos))}\n")
    out.append("---\n")
    out.append("## Conclusion\n")
    out.append("### 🏆 Where InstaNovo v1.2.2 Leads\n")
    out.append(
        f"InstaNovo v1.2.2 ranks 🥇 on {n_first_peptide}/{n_datasets} datasets and top-3 on "
        f"{n_top3_peptide}/{n_datasets} for peptide-level precision. Compared to v1.1.2, it improves by a mean of "
        f"**{mean_peptide_delta * 100:+.1f}%** AUC across {len(deltas_peptide)} datasets where both versions ran.\n"
    )

    # Per-section
    for i, (title, blurb, datasets) in enumerate(SECTIONS, start=1):
        out.append(f"## {i}. {title}\n")
        if blurb:
            out.append(blurb + "\n")
        pep_rows = section_table_rows(datasets, "peptide_precision_plot_data.csv", "1.2.2", "1.1.2")
        aa_rows = section_table_rows(datasets, "AA_precision_plot_data.csv", "1.2.2", "1.1.2")
        if pep_rows:
            out.append(render_table(pep_rows, "Peptide Precision AUC"))
        if aa_rows:
            out.append(render_table(aa_rows, "Amino Acid Precision AUC"))
        out.append("---\n")

    # Regressions
    out.append(f"## {len(SECTIONS) + 1}. 📉 Regressions (v1.2.2 < v1.1.2)\n")
    if regressions:
        regressions.sort(key=lambda r: r[4])
        out.append(f"Out of {2 * sum(1 for ds in all_datasets if load_auc_for_dataset(ds, 'peptide_precision_plot_data.csv'))} dataset/metric combinations, "
                   f"{len(regressions)} show a decrease:\n")
        out.append("| Dataset | Metric | v1.1.2 AUC | v1.2.2 AUC | Delta |")
        out.append("|---|---|---|---|---|")
        for ds, metric, v112, v122, delta in regressions:
            out.append(f"| {ds} | {metric} | {fmt_auc(v112)} | {fmt_auc(v122)} | {delta * 100:+.1f}% |")
    else:
        out.append("No regressions: v1.2.2 ≥ v1.1.2 on every dataset/metric combination.")

    out.append("")

    # Datasets covered by sections vs all
    listed = set()
    for _, _, ds_list in SECTIONS:
        listed.update(ds_list)
    unlisted = sorted(set(all_datasets) - listed)
    if unlisted:
        out.append("\n## Appendix: Datasets not assigned to a section\n")
        for ds in unlisted:
            out.append(f"- {ds}")
        out.append("")

    return "\n".join(out)


if __name__ == "__main__":
    md = main()
    out_path = os.path.expanduser("~/Downloads/benchmark_report_updated.md")
    with open(out_path, "w") as f:
        f.write(md)
    print(f"Wrote {out_path} ({len(md)} bytes)")
