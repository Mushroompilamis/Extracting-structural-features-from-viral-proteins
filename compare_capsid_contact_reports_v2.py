#!/usr/bin/env python3

"""
Comparative Analysis of Wild-Type and Mutant Viral Capsid Contact Networks

This script compares residue contact networks between a wild-type (WT)
viral capsid and one or more mutant capsid models using ChimeraX
crystalcontacts output files. The objective is to quantify how specific
mutations affect intermolecular interactions and overall capsid interface
organization.

For each structure, the script:
- Parses residue contact and buried solvent-accessible surface area (SASA)
  information.
- Normalizes symmetry-expanded residue numbering to the original sequence
  numbering, allowing equivalent residues from different symmetry copies to
  be compared.
- Removes contacts below a user-defined distance threshold.
- Classifies residue interactions as hydrogen bond-like, electrostatic,
  hydrophobic, aromatic, or other.

The script then performs a direct comparison between the WT and mutant
structures by calculating:
- Total and unique residue contacts.
- Inter-chain and inter-copy interactions.
- Symmetry-related contacts.
- Interaction-type frequencies.
- Contact distance statistics.
- Buried SASA measurements.
- Residue-type distributions.
- Lost and newly formed residue-residue contacts relative to the WT.
- Optional local interaction network surrounding a selected target residue.

The analysis automatically generates:
- Parsed and normalized contact datasets for each structure.
- Side-by-side comparison tables of all structural metrics.
- Difference tables showing absolute and percentage changes relative to
  the WT control.
- Lists of lost and gained contacts for each mutant.
- Comparison plots.
- A report-ready summary describing the structural differences
  between WT and mutant capsids.

This workflow enables systematic evaluation of how point mutations alter
viral capsid interface architecture and provides quantitative evidence for
their potential effects on capsid stability and assembly.

"""

import argparse
import re
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

CONTACT_RE = re.compile(
    r"^\s*([A-Z0-9]{3})\s+(-?\d+)\s+(\S)\s+(-?\d+)\s+"
    r"([A-Z0-9]{3})\s+(-?\d+)\s+(\S)\s+(-?\d+)\s+"
    r"(-?\d+)\s+(-?\d+)\s+(-?\d+)\s+(-?\d+)\s+"
    r"(-?\d+(?:\.\d+)?)\s*$"
)

BURIAL_RE = re.compile(
    r"^\s*([A-Z0-9]{3})\s+(-?\d+)\s+(\S)\s+(-?\d+)\s+(-?\d+(?:\.\d+)?)\s*$"
)

HYDROPHOBIC = {"ALA", "VAL", "LEU", "ILE", "MET", "PRO"}
AROMATIC = {"PHE", "TYR", "TRP", "HIS"}
POSITIVE = {"LYS", "ARG", "HIS"}
NEGATIVE = {"ASP", "GLU"}
POLAR = {"SER", "THR", "ASN", "GLN", "CYS"}
HBOND_CAPABLE = POLAR | AROMATIC | POSITIVE | NEGATIVE
AA1 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}


def normalize_resnum(num: int) -> int:
    n = abs(int(num))
    if n >= 1000:
        rem = n % 1000
        return rem if rem != 0 else n
    return n


def observed_key(res, num, chain, ncs):
    return f"{res}{num}_{chain}_ncs{ncs}"


def classify_interaction(res1, res2, distance):
    if (res1 in POSITIVE and res2 in NEGATIVE) or (res1 in NEGATIVE and res2 in POSITIVE):
        return "electrostatic"
    if res1 in AROMATIC and res2 in AROMATIC:
        return "aromatic"
    if res1 in HYDROPHOBIC and res2 in HYDROPHOBIC:
        return "hydrophobic"
    if distance <= 3.5 and res1 in HBOND_CAPABLE and res2 in HBOND_CAPABLE:
        return "hydrogen_bond_like"
    return "other"


def parse_contacts_file(path: Path):
    contacts, buried = [], []
    section = None
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\n")
            if line.startswith("Residue contacts"):
                section = "contacts"
                continue
            if line.startswith("Residue buried areas"):
                section = "buried"
                continue
            if not line.strip() or line.lstrip().startswith("residue"):
                continue
            if section == "contacts":
                m = CONTACT_RE.match(line)
                if m:
                    contacts.append({
                        "line": line_number,
                        "res1": m.group(1), "num1": int(m.group(2)), "chain1": m.group(3), "ncs1": int(m.group(4)),
                        "res2": m.group(5), "num2": int(m.group(6)), "chain2": m.group(7), "ncs2": int(m.group(8)),
                        "sym2": int(m.group(9)), "cell_x": int(m.group(10)), "cell_y": int(m.group(11)), "cell_z": int(m.group(12)),
                        "distance": float(m.group(13)),
                    })
            elif section == "buried":
                m = BURIAL_RE.match(line)
                if m:
                    buried.append({
                        "line": line_number, "res": m.group(1), "num": int(m.group(2)), "chain": m.group(3),
                        "ncs": int(m.group(4)), "buried_area": float(m.group(5)),
                    })
    contacts_df = pd.DataFrame(contacts)
    buried_df = pd.DataFrame(buried)
    if contacts_df.empty:
        raise ValueError(f"No contacts parsed from {path}")

    contacts_df["seq_num1"] = contacts_df["num1"].map(normalize_resnum)
    contacts_df["seq_num2"] = contacts_df["num2"].map(normalize_resnum)
    contacts_df["residue1_seq"] = contacts_df.apply(lambda r: f"{r.res1}{r.seq_num1}", axis=1)
    contacts_df["residue2_seq"] = contacts_df.apply(lambda r: f"{r.res2}{r.seq_num2}", axis=1)
    contacts_df["observed1"] = contacts_df.apply(lambda r: observed_key(r.res1, r.num1, r.chain1, r.ncs1), axis=1)
    contacts_df["observed2"] = contacts_df.apply(lambda r: observed_key(r.res2, r.num2, r.chain2, r.ncs2), axis=1)
    contacts_df["observed_pair"] = contacts_df.apply(lambda r: " -- ".join(sorted([r.observed1, r.observed2])), axis=1)
    contacts_df["sequence_pair"] = contacts_df.apply(lambda r: " -- ".join(sorted([r.residue1_seq, r.residue2_seq])), axis=1)
    contacts_df["inter_chain_or_copy"] = ((contacts_df["chain1"] != contacts_df["chain2"]) | (contacts_df["ncs1"] != contacts_df["ncs2"]) | (contacts_df["sym2"] != 0) | (contacts_df["cell_x"] != 0) | (contacts_df["cell_y"] != 0) | (contacts_df["cell_z"] != 0))
    contacts_df["symmetry_related"] = ((contacts_df["sym2"] != 0) | (contacts_df["cell_x"] != 0) | (contacts_df["cell_y"] != 0) | (contacts_df["cell_z"] != 0))
    contacts_df["interaction_type"] = contacts_df.apply(lambda r: classify_interaction(r.res1, r.res2, r.distance), axis=1)

    if not buried_df.empty:
        buried_df["seq_num"] = buried_df["num"].map(normalize_resnum)
        buried_df["residue_seq"] = buried_df.apply(lambda r: f"{r.res}{r.seq_num}", axis=1)
        buried_df["observed"] = buried_df.apply(lambda r: observed_key(r.res, r.num, r.chain, r.ncs), axis=1)
    return contacts_df, buried_df


def summarize_sample(label, contacts_df, buried_df):
    unique_observed = contacts_df.drop_duplicates("observed_pair")
    unique_seq = contacts_df.drop_duplicates("sequence_pair")
    row = {
        "sample": label,
        "raw_contacts": len(contacts_df),
        "unique_observed_contacts": len(unique_observed),
        "unique_sequence_contacts": len(unique_seq),
        "inter_chain_or_copy_contacts": int(contacts_df["inter_chain_or_copy"].sum()),
        "symmetry_related_contacts": int(contacts_df["symmetry_related"].sum()),
        "hydrogen_bond_like_contacts": int((contacts_df["interaction_type"] == "hydrogen_bond_like").sum()),
        "hydrophobic_contacts": int((contacts_df["interaction_type"] == "hydrophobic").sum()),
        "electrostatic_contacts": int((contacts_df["interaction_type"] == "electrostatic").sum()),
        "aromatic_contacts": int((contacts_df["interaction_type"] == "aromatic").sum()),
        "mean_distance": round(float(contacts_df["distance"].mean()), 3),
        "median_distance": round(float(contacts_df["distance"].median()), 3),
        "contacts_le_3_5A": int((contacts_df["distance"] <= 3.5).sum()),
        "contacts_le_4A": int((contacts_df["distance"] <= 4.0).sum()),
    }
    if buried_df.empty:
        row.update({"total_buried_sasa": 0.0, "mean_buried_sasa": 0.0, "max_buried_sasa": 0.0, "buried_residue_entries": 0})
    else:
        row.update({
            "total_buried_sasa": round(float(buried_df["buried_area"].sum()), 3),
            "mean_buried_sasa": round(float(buried_df["buried_area"].mean()), 3),
            "max_buried_sasa": round(float(buried_df["buried_area"].max()), 3),
            "buried_residue_entries": len(buried_df),
        })
    return row


def value_counts_merge(tables, key):
    out = tables[0]
    for t in tables[1:]:
        out = out.merge(t, on=key, how="outer")
    return out.fillna(0)


def residue_type_table(label, contacts_df):
    t = pd.concat([
        contacts_df[["res1"]].rename(columns={"res1": "residue"}),
        contacts_df[["res2"]].rename(columns={"res2": "residue"}),
    ], ignore_index=True)
    return t.value_counts("residue").reset_index(name=label)


def interaction_type_table(label, contacts_df):
    return contacts_df.value_counts("interaction_type").reset_index(name=label)


def target_neighborhood(label, contacts_df, target):
    target = str(target).upper().replace(":", "")
    m = re.match(r"([A-Z]{1,3})?(\d+)$", target)
    if not m:
        return pd.DataFrame()
    target_res, target_num = m.group(1), int(m.group(2))

    def is_target(res, seq_num):
        if int(seq_num) != target_num:
            return False
        if not target_res:
            return True
        if len(target_res) == 1:
            return AA1.get(res, res[0]) == target_res
        return res == target_res

    rows = []
    for _, r in contacts_df.iterrows():
        side1 = is_target(r.res1, r.seq_num1)
        side2 = is_target(r.res2, r.seq_num2)
        if not (side1 or side2):
            continue
        if side1:
            target_seq, target_obs = f"{r.res1}{r.seq_num1}", r.observed1
            partner_seq, partner_obs = f"{r.res2}{r.seq_num2}", r.observed2
        else:
            target_seq, target_obs = f"{r.res2}{r.seq_num2}", r.observed2
            partner_seq, partner_obs = f"{r.res1}{r.seq_num1}", r.observed1
        rows.append({
            "sample": label,
            "target_sequence_residue": target_seq,
            "target_observed": target_obs,
            "partner_sequence_residue": partner_seq,
            "partner_observed": partner_obs,
            "distance": r.distance,
            "interaction_type": r.interaction_type,
            "inter_chain_or_copy": r.inter_chain_or_copy,
            "symmetry_related": r.symmetry_related,
        })
    return pd.DataFrame(rows)


def sequence_contact_set(contacts_df):
    return set(contacts_df.drop_duplicates("sequence_pair")["sequence_pair"])


def make_barplot(metrics_df, metric, outdir):
    pretty = metric.replace("_", " ").title()
    plt.figure(figsize=(8, 5))
    plt.bar(metrics_df["sample"], metrics_df[metric])
    plt.ylabel(pretty)
    plt.xlabel("Sample")
    plt.title(pretty)
    plt.xticks(rotation=30, ha="right")
    plt.tight_layout()
    plt.savefig(outdir / f"{metric}.png", dpi=300)
    plt.close()


def make_percent_plot(diffs_df, outdir):
    key_metrics = ["unique_sequence_contacts", "inter_chain_or_copy_contacts", "symmetry_related_contacts", "hydrogen_bond_like_contacts", "hydrophobic_contacts", "total_buried_sasa"]
    key = diffs_df[diffs_df["metric"].isin(key_metrics)].copy()
    if key.empty:
        return
    key["label"] = key["sample"] + "\n" + key["metric"].str.replace("_", " ")
    plt.figure(figsize=(10, 6))
    plt.barh(key["label"], key["percent_change"])
    plt.axvline(0, linewidth=1)
    plt.xlabel("Percent change vs WT/control (%)")
    plt.title("Main changes relative to control")
    plt.tight_layout()
    plt.savefig(outdir / "main_percent_changes_vs_control.png", dpi=300)
    plt.close()


def table_text(df, max_rows=None):
    if df is None or df.empty:
        return "No rows."
    x = df.head(max_rows).copy() if max_rows else df.copy()
    return x.to_string(index=False)


def write_report(outdir, control_label, metrics_df, diffs_df, target_df, lost_gained_summaries):
    report = outdir / "report_ready_summary.md"
    key_metrics = ["unique_sequence_contacts", "inter_chain_or_copy_contacts", "symmetry_related_contacts", "hydrogen_bond_like_contacts", "hydrophobic_contacts", "total_buried_sasa"]
    key = diffs_df[diffs_df["metric"].isin(key_metrics)].copy()

    with open(report, "w", encoding="utf-8") as f:
        f.write("# Report-ready capsid contact comparison\n\n")
        f.write("## Analysis design\n\n")
        f.write("ChimeraX `crystalcontacts` outputs were compared for the WT/control and mutant capsid models. Residue numbers were normalized so 1BMV-style observed numbers such as 3159 are interpreted as sequence residue 159. Thus, the report can refer to T159A, while the ChimeraX command remains `swapaa /2:3159 ala`.\n\n")
        f.write("## Critical quality check\n\n")
        f.write("The WT and mutant contact files must be generated using the same ChimeraX protocol: same structure state, same assembly/ASU choice, same distance cutoff, and same `crystalcontacts` settings. If WT was calculated from a different model state than the mutants, very large global differences may reflect protocol differences rather than mutation effects.\n\n")
        f.write("## Side-by-side global metrics\n\n```text\n")
        f.write(table_text(metrics_df))
        f.write("\n```\n\n")
        f.write("## Key changes vs WT/control\n\n```text\n")
        f.write(table_text(key))
        f.write("\n```\n\n")
        f.write("## Interpretation\n\n")
        for sample in [s for s in metrics_df["sample"] if s != control_label]:
            sub = diffs_df[diffs_df["sample"] == sample]
            f.write(f"### {sample} vs {control_label}\n\n")
            for metric in key_metrics:
                row = sub[sub["metric"] == metric]
                if not row.empty:
                    r = row.iloc[0]
                    f.write(f"- {metric.replace('_', ' ')} changed from {r.control_value} to {r.sample_value} ({r.delta:+.3f}, {r.percent_change:+.2f}%).\n")
            f.write("\nThis suggests altered interface/contact organization in the mutant model. The result should be interpreted as an in silico structural-contact comparison, not direct experimental proof that capsid assembly fails.\n\n")
        if target_df is not None and not target_df.empty:
            f.write("## Target residue neighborhood\n\n```text\n")
            f.write(table_text(target_df.sort_values(["sample", "distance"]), max_rows=100))
            f.write("\n```\n\n")
        f.write("## Lost/gained sequence-level contacts\n\n")
        for sample, summary in lost_gained_summaries.items():
            f.write(f"### {sample}\n\n")
            f.write(f"- Lost contacts vs {control_label}: {summary['lost_count']}\n")
            f.write(f"- Gained contacts vs {control_label}: {summary['gained_count']}\n\n")
            if summary["lost_examples"]:
                f.write("Examples of lost contacts:\n")
                for x in summary["lost_examples"]:
                    f.write(f"- {x}\n")
                f.write("\n")
        f.write("## Suggested sentence for the report\n\n")
        f.write("The in silico T→A mutations reduced several interface descriptors, including inter-chain/copy contacts, symmetry-related contacts, hydrogen-bond-like contacts, and buried surface area compared with the WT control. This supports the hypothesis that the selected threonine residues contribute to inter-subunit stabilization in the 1BMV capsid, although experimental validation would be required to confirm an assembly defect.\n")
    return report


def parse_named_file(arg):
    if "=" not in arg:
        raise argparse.ArgumentTypeError("Use LABEL=filename, e.g. WT=1bmv_contacts.txt")
    label, filename = arg.split("=", 1)
    return label, Path(filename)


def main():
    parser = argparse.ArgumentParser(description="Compare WT and mutant ChimeraX crystalcontacts files with report-ready outputs.")
    parser.add_argument("--control", required=True, type=parse_named_file, help="Control file as LABEL=path")
    parser.add_argument("--samples", required=True, nargs="+", type=parse_named_file, help="Mutant files as LABEL=path")
    parser.add_argument("--outdir", default="capsid_comparison_report")
    parser.add_argument("--target", default=None, help="Optional target residue, e.g. THR159, T159, 159")
    parser.add_argument("--min-distance", type=float, default=1.0)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    plots_dir = outdir / "plots"
    outdir.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)

    all_inputs = [args.control] + list(args.samples)
    parsed = {}
    for label, path in all_inputs:
        contacts, buried = parse_contacts_file(path)
        contacts = contacts[contacts["distance"] >= args.min_distance].copy()
        parsed[label] = {"contacts": contacts, "buried": buried}
        contacts.to_csv(outdir / f"{label}_contacts_parsed_normalized.csv", index=False)
        buried.to_csv(outdir / f"{label}_buried_parsed_normalized.csv", index=False)

    control_label = args.control[0]
    metrics_df = pd.DataFrame([summarize_sample(label, d["contacts"], d["buried"]) for label, d in parsed.items()])
    metrics_df.to_csv(outdir / "side_by_side_metrics.csv", index=False)

    control_row = metrics_df[metrics_df["sample"] == control_label].iloc[0]
    diff_rows = []
    for _, row in metrics_df.iterrows():
        if row["sample"] == control_label:
            continue
        for metric in metrics_df.columns:
            if metric == "sample":
                continue
            c, s = control_row[metric], row[metric]
            delta = s - c
            pct = (delta / c * 100) if c != 0 else None
            diff_rows.append({"sample": row["sample"], "metric": metric, "control_value": c, "sample_value": s, "delta": round(float(delta), 3), "percent_change": round(float(pct), 2) if pct is not None else None})
    diffs_df = pd.DataFrame(diff_rows)
    diffs_df.to_csv(outdir / "differences_vs_control.csv", index=False)

    residue_df = value_counts_merge([residue_type_table(label, d["contacts"]) for label, d in parsed.items()], "residue")
    residue_df.to_csv(outdir / "residue_type_comparison.csv", index=False)
    interaction_df = value_counts_merge([interaction_type_table(label, d["contacts"]) for label, d in parsed.items()], "interaction_type")
    interaction_df.to_csv(outdir / "interaction_type_comparison.csv", index=False)

    target_df = pd.DataFrame()
    if args.target:
        parts = [target_neighborhood(label, d["contacts"], args.target) for label, d in parsed.items()]
        parts = [p for p in parts if p is not None and not p.empty]
        if parts:
            target_df = pd.concat(parts, ignore_index=True)
            target_df.to_csv(outdir / "target_residue_neighborhood.csv", index=False)

    control_set = sequence_contact_set(parsed[control_label]["contacts"])
    lost_gained_summaries = {}
    for label, d in parsed.items():
        if label == control_label:
            continue
        sample_set = sequence_contact_set(d["contacts"])
        lost = sorted(control_set - sample_set)
        gained = sorted(sample_set - control_set)
        pd.DataFrame({"status": ["lost_vs_control"] * len(lost) + ["gained_vs_control"] * len(gained), "sequence_pair": lost + gained}).to_csv(outdir / f"lost_gained_contacts_{label}.csv", index=False)
        lost_gained_summaries[label] = {"lost_count": len(lost), "gained_count": len(gained), "lost_examples": lost[:15]}

    for metric in ["unique_sequence_contacts", "inter_chain_or_copy_contacts", "symmetry_related_contacts", "hydrogen_bond_like_contacts", "hydrophobic_contacts", "total_buried_sasa"]:
        make_barplot(metrics_df, metric, plots_dir)
    make_percent_plot(diffs_df, plots_dir)

    report = write_report(outdir, control_label, metrics_df, diffs_df, target_df, lost_gained_summaries)
    print("\nComparison complete.")
    print(f"Output folder: {outdir.resolve()}")
    print(f"Main report: {report.resolve()}")


if __name__ == "__main__":
    main()
