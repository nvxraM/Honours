import os
import re
import pandas as pd
from pathlib import Path

# --- Optional GC3 dependencies (Biopython) ---
try:
    from Bio import SeqIO
    from Bio.SeqUtils import gc_fraction
except Exception:
    SeqIO = None
    gc_fraction = None

# ====================== CONFIG ======================
ROOT_DIR = Path("sequences/Species_POP_Mega_Haploid_split")   # DnaSP .out files
MEGA_BASE_DIR = Path("sequences/Interactive_group_editor")    # MEGA metrics tree

OUTPUT_DIR = Path("sequences/Species_POP_Mega_CSVOutput")
OUTPUT_XLSX = OUTPUT_DIR / "All_species_population_summary.xlsx"

CHECK_CONFLICTS = True  # warn if multiple files disagree on Tajima/Fu values

# ====================================================
COLUMNS = [
    "Species","Common name","Family","Population","Combined Common (Pop)",
    "Number of Sequences used","Variable sites (VS)","Haplotype number (HN)",
    "Haplotype diversity (Hd)","Haplotype diversity Standard Deviation (Hd.SD)",
    "Nucleotide diversity (Pi)","Nucleotide diversity standard deviation (Pi.SD)",
    "Tajima's D","Tajima's D significance","Fu and Li's D*","Fu and Li's D* significance",
    "d","d S.E.","dN","dN S.E.","dS","dS S.E.","dN/dS Ratio","Avg GC3 (%)","Sequence No.",
]

INT_COLS = ["Number of Sequences used","Variable sites (VS)","Haplotype number (HN)","Sequence No."]

TEXT_DECIMAL_COLS = [
    "Haplotype diversity (Hd)","Haplotype diversity Standard Deviation (Hd.SD)",
    "Nucleotide diversity (Pi)","Nucleotide diversity standard deviation (Pi.SD)",
    "Tajima's D","Fu and Li's D*","d","d S.E.","dN","dN S.E.","dS","dS S.E.","dN/dS Ratio","Avg GC3 (%)",
]

# Number pattern (handles scientific notation)
DEC = r"([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)"

# =================== LOOKUP TABLES ===================
COMMON_NAME = {
    "Balaena_mysticetus":"Bowhead whale","Balaenoptera_musculus":"Blue whale",
    "Balaenoptera_physalus":"Fin whale","Delphinapterus_leucas":"Beluga whale",
    "Eubalaena_japonica":"North Pacific right whale","Globicephala_macrorhynchus":"Short-finned pilot whale",
    "Hyperoodon_ampullatus":"Northern bottlenose whale","Mesoplodon_grayi":"Gray's beaked whale",
    "Monodon_monoceros":"Narwhal","Orcaella_brevirostris":"Irrawaddy dolphin",
    "Orcinus_orca":"Killer whale (Orca)","Peponocephala_electra":"Melon-headed whale",
    "Phocoena_phocoena":"Harbour porpoise","Phocoena_sinus":"Vaquita",
    "Physeter_macrocephalus":"Sperm whale","Pseudorca_crassidens":"False killer whale",
    "Stenella_longirostris":"Spinner dolphin","Tursiops_aduncus":"Indo-Pacific bottlenose dolphin",
    "Tursiops_truncatus":"Common bottlenose dolphin","Ziphius_cavirostris":"Cuvier's beaked whale",
}
FAMILY_BY_SPECIES = {
    "Balaena_mysticetus":"Balaenidae","Eubalaena_japonica":"Balaenidae",
    "Balaenoptera_musculus":"Balaenopteridae","Balaenoptera_physalus":"Balaenopteridae",
    "Delphinapterus_leucas":"Monodontidae","Monodon_monoceros":"Monodontidae",
    "Globicephala_macrorhynchus":"Delphinidae","Orcaella_brevirostris":"Delphinidae",
    "Orcinus_orca":"Delphinidae","Peponocephala_electra":"Delphinidae",
    "Pseudorca_crassidens":"Delphinidae","Stenella_longirostris":"Delphinidae",
    "Tursiops_aduncus":"Delphinidae","Tursiops_truncatus":"Delphinidae",
    "Phocoena_phocoena":"Phocoenidae","Phocoena_sinus":"Phocoenidae",
    "Physeter_macrocephalus":"Physeteridae","Hyperoodon_ampullatus":"Ziphiidae",
    "Mesoplodon_grayi":"Ziphiidae","Ziphius_cavirostris":"Ziphiidae",
}
FAMILY_BY_GENUS = {
    "Balaena":"Balaenidae","Eubalaena":"Balaenidae","Balaenoptera":"Balaenopteridae",
    "Delphinapterus":"Monodontidae","Monodon":"Monodontidae",
    "Globicephala":"Delphinidae","Orcaella":"Delphinidae","Orcinus":"Delphinidae",
    "Peponocephala":"Delphinidae","Pseudorca":"Delphinidae","Stenella":"Delphinidae","Tursiops":"Delphinidae",
    "Phocoena":"Phocoenidae","Physeter":"Physeteridae",
    "Hyperoodon":"Ziphiidae","Mesoplodon":"Ziphiidae","Ziphius":"Ziphiidae",
}
WORD_TO_NUM = {"one":"1","two":"2","three":"3","four":"4","five":"5","six":"6","seven":"7","eight":"8","nine":"9","ten":"10"}

# =================== HELPERS ===================
def normalize_text(s: str) -> str:
    """Normalize quotes/spaces/newlines to make regex robust."""
    return (s.replace("\r\n","\n").replace("\r","\n")
             .replace("\u2019","'").replace("\u2018","'")
             .replace("\u2013","-").replace("\u2014","-")
             .replace("\xa0"," "))

def map_population_to_group(pop_text: str) -> str:
    if not pop_text:
        return ""
    m = re.search(r"\bPopulation_(\d+|one|two|three|four|five|six|seven|eight|nine|ten)\b",
                  pop_text, flags=re.IGNORECASE)
    if not m:
        return pop_text
    token = m.group(1)
    num_str = WORD_TO_NUM.get(token.lower(), token)
    return f"Gp{num_str}"

def common_name_for(scientific: str) -> str:
    return COMMON_NAME.get(scientific, scientific.replace("_"," ").title())

def family_for(scientific: str) -> str:
    if scientific in FAMILY_BY_SPECIES:
        return FAMILY_BY_SPECIES[scientific]
    genus = scientific.split("_",1)[0]
    return FAMILY_BY_GENUS.get(genus,"Unknown")

def shorten_scientific(scientific: str) -> str:
    if not scientific:
        return ""
    parts = scientific.split("_")
    if len(parts) >= 2:
        genus, species = parts[0], parts[1]
        return f"{genus[0]}. {species.lower()}"
    name = scientific.replace("_"," ").strip()
    return name and f"{name.split()[0][0]}. {' '.join(name.split()[1:])}".strip() or ""

def build_combined(common: str, species: str, pop: str) -> str:
    common = ("" if pd.isna(common) else str(common).strip())
    species = ("" if pd.isna(species) else str(species).strip())
    pop = ("" if pd.isna(pop) else str(pop).strip())
    short = shorten_scientific(species) if species else ""
    if common and short and pop:  return f"{common} ({short}) {pop}"
    if common and short:          return f"{common} ({short})"
    if common and pop:            return f"{common} {pop}"
    return common or short or pop or ""

# =================== DNASP PARSERS ===================
def parse_tests(text: str) -> dict:
    """
    Extract Tajima's D (+sig) and Fu & Li's D* (+sig) from normalized text.
    Handles variations like: "Fu and Li’s D* test statistic, FLD*: 0.22909"
    """
    t = normalize_text(text)
    out = {"Tajima's D":"", "Tajima's D significance":"", "Fu and Li's D*":"", "Fu and Li's D* significance":""}

    # Accept straight or curly apostrophes (already normalized), flexible spaces
    taj = re.search(rf"Tajima's\s*D\s*:\s*{DEC}", t, flags=re.IGNORECASE)
    if taj:
        out["Tajima's D"] = taj.group(1)
        # inline or next-line significance
        # first try after the Tajima match; fallback to anywhere
        sig_after = re.search(r"Statistical\s+significance:\s*([^\r\n]+)", t[taj.end():], flags=re.IGNORECASE)
        if not sig_after:
            sig_after = re.search(r"Statistical\s+significance:\s*([^\r\n]+)", t, flags=re.IGNORECASE)
        if sig_after:
            out["Tajima's D significance"] = sig_after.group(1).strip()

    # Fu & Li's D*:
    #  - "Fu and Li's D* test statistic: 0.22909"
    #  - "Fu and Li's D* test statistic, FLD*: 0.22909"
    #  - (fallback) "FLD*: 0.22909"
    fu1 = re.search(
        rf"Fu\s+and\s+Li's\s*D\*\s*(?:test\s*statistic(?:,?\s*FLD\*)?)?\s*:\s*{DEC}",
        t, flags=re.IGNORECASE
    )
    fu_val = None
    if fu1:
        fu_val = fu1.group(1)
    else:
        fu2 = re.search(rf"\bFLD\*\s*[:=]\s*{DEC}", t, flags=re.IGNORECASE)
        if fu2:
            fu_val = fu2.group(1)

    if fu_val is not None:
        out["Fu and Li's D*"] = fu_val
        sig_after_fu = re.search(r"Statistical\s+significance:\s*([^\r\n]+)", t[fu1.end() if fu1 else 0:], flags=re.IGNORECASE)
        if not sig_after_fu:
            sig_after_fu = re.search(r"Statistical\s+significance:\s*([^\r\n]+)", t, flags=re.IGNORECASE)
        if sig_after_fu:
            out["Fu and Li's D* significance"] = sig_after_fu.group(1).strip()

    return out

def extract_fields_from_dnasp_out(text: str, filename: str, species: str) -> dict:
    """Pull counts + diversities + Tajima/Fu&Li from a DnaSP .out file."""
    t = normalize_text(text)

    # Population from content or filename
    pop_match = re.search(r"Population\s+used:\s*(\S+)", t, flags=re.IGNORECASE)
    if pop_match:
        population_raw = pop_match.group(1)
    else:
        mfile = re.search(r"(Population_(\d+|one|two|three|four|five|six|seven|eight|nine|ten))",
                          filename, re.IGNORECASE)
        population_raw = mfile.group(1) if mfile else ""
    population = map_population_to_group(population_raw)

    patterns = {
        "Number of Sequences used": r"Number of sequences used:\s*(\d+)",
        "Variable sites (VS)": r"Number of polymorphic \(segregating\) sites,\s*S:\s*(\d+)",
        "Haplotype number (HN)": r"Number of Haplotypes,\s*h:\s*(\d+)",
        "Haplotype diversity (Hd)": rf"Haplotype\s*\(gene\)\s*diversity,\s*Hd:\s*{DEC}",
        "Haplotype diversity Standard Deviation (Hd.SD)": rf"Standard\s*Deviation\s*of\s*Haplotype\s*diversity:\s*{DEC}",
        "Nucleotide diversity (Pi)": rf"Nucleotide\s*diversity,\s*Pi:\s*{DEC}",
        "Nucleotide diversity standard deviation (Pi.SD)": rf"Standard\s*deviation\s*of\s*Pi:\s*{DEC}",
    }

    row = {
        "Species": species,
        "Common name": common_name_for(species),
        "Family": family_for(species),
        "Population": population,
        "Tajima's D": "",
        "Tajima's D significance": "",
        "Fu and Li's D*": "",
        "Fu and Li's D* significance": "",
    }

    for field, pat in patterns.items():
        m2 = re.search(pat, t, flags=re.IGNORECASE)
        row[field] = m2.group(1) if m2 else row.get(field, "")

    row.update(parse_tests(t))
    return row

# --------- Multiple .out files: priority + conflict checking ----------
def out_priority(path: Path) -> int:
    n = path.name.lower()
    if "tajima" in n: return 0
    if "fu" in n:     return 1
    return 2

def _merge_rows(base: dict, new: dict, src: str = "") -> dict:
    for k, v in new.items():
        if v is None or str(v).strip() == "":
            continue
        if k not in base or str(base[k]).strip() == "":
            base[k] = v
        elif CHECK_CONFLICTS and k in ("Tajima's D","Fu and Li's D*","Tajima's D significance","Fu and Li's D* significance"):
            if str(base[k]).strip() != str(v).strip():
                print(f"⚠️ Conflict for {k} in {src}: kept '{base[k]}', ignored '{v}'")
    return base

# =================== MEGA PARSERS ===================
_MEGA_LINE = re.compile(r'^(\S+)\s+(\S+)\s+Population\s+(\d+)\s+(\S+)\s+(\S+)$')

def _safe_str(x):
    return "" if x is None else str(x)

def _fmt_ratio(dn: str, ds: str) -> str:
    try:
        dn_f = float(dn); ds_f = float(ds)
        if ds_f == 0: return ""
        return f"{dn_f/ds_f:.12g}"
    except Exception:
        return ""

def _gc3_percent_from_seqs(seq_list):
    if not seq_list or gc_fraction is None:
        return ""
    thirds = [s[2::3] for s in seq_list]
    all3 = "".join(thirds)
    if not all3:
        return ""
    return f"{gc_fraction(all3)*100:.12g}"

def collect_mega_metrics(base_dir: Path) -> pd.DataFrame:
    rows = []
    if not base_dir.exists():
        return pd.DataFrame(columns=["Species","Population","d","d S.E.","dN","dN S.E.","dS","dS S.E.","dN/dS Ratio","Avg GC3 (%)","Sequence No."])

    for genus in os.listdir(base_dir):
        genus_path = base_dir / genus
        if not genus_path.is_dir(): continue

        for species_dir in os.listdir(genus_path):
            species_path = genus_path / species_dir
            if not species_path.is_dir(): continue

            pop_dir = species_path / "MEGA_Populations"
            if not pop_dir.is_dir(): continue

            d_file      = pop_dir / f"{species_dir}_d_value.txt"
            nonsyn_file = pop_dir / f"{species_dir}_nonsynonymous_only.txt"
            syn_file    = pop_dir / f"{species_dir}_Synonymous_only.txt"
            grp_file    = pop_dir / f"{species_dir}_population.grp"
            fasta_path  = base_dir / genus / species_dir / f"{species_dir}_concatenated_gapped_sequences.fasta"

            pop_metrics = {}

            def parse_value_file(path, val_key, se_key):
                if not path.exists(): return
                with open(path,"r",encoding="utf-8",errors="ignore") as f:
                    for line in f:
                        line = line.strip()
                        if not line or line.startswith("#"): continue
                        m = _MEGA_LINE.match(line)
                        if not m: continue
                        _, _, pop_str, value_str, se_str = m.groups()
                        pnum = int(pop_str)
                        pop_metrics.setdefault(pnum, {})
                        pop_metrics[pnum][val_key] = _safe_str(value_str)
                        pop_metrics[pnum][se_key]  = _safe_str(se_str)

            parse_value_file(d_file, "d", "d S.E.")
            parse_value_file(nonsyn_file, "dN", "dN S.E.")
            parse_value_file(syn_file, "dS", "dS S.E.")

            # grp mapping
            pop_to_acc = {}
            if grp_file.exists():
                with open(grp_file,"r",encoding="utf-8",errors="ignore") as gf:
                    for line in gf:
                        line = line.strip()
                        if not line or line.startswith("#") or "=" not in line: continue
                        left, right = line.split("=",1)
                        m = re.search(r"_Population_(\d+)$", right)
                        if not m: continue
                        pnum = int(m.group(1))
                        pop_to_acc.setdefault(pnum, []).append(left)

            # GC3
            if SeqIO is not None and fasta_path.exists() and fasta_path.stat().st_size > 0:
                pop_to_seqs = {p: [] for p in pop_to_acc}
                for rec in SeqIO.parse(str(fasta_path), "fasta"):
                    for pnum, accs in pop_to_acc.items():
                        if rec.id in accs:
                            pop_to_seqs[pnum].append(str(rec.seq))
                for pnum, seqs in pop_to_seqs.items():
                    gc3 = _gc3_percent_from_seqs(seqs)
                    pop_metrics.setdefault(pnum, {})
                    pop_metrics[pnum]["Avg GC3 (%)"] = gc3
                    pop_metrics[pnum]["Sequence No."] = str(len(seqs))
            else:
                for pnum in pop_to_acc:
                    pop_metrics.setdefault(pnum, {})
                    pop_metrics[pnum]["Avg GC3 (%)"] = ""
                    pop_metrics[pnum]["Sequence No."] = ""

            # Ratio
            for pnum, met in pop_metrics.items():
                dn = met.get("dN","")
                ds = met.get("dS","")
                met["dN/dS Ratio"] = _fmt_ratio(dn, ds)

            for pnum, met in pop_metrics.items():
                rows.append({
                    "Species": species_dir,
                    "Population": f"Gp{pnum}",
                    "d": met.get("d",""),
                    "d S.E.": met.get("d S.E.",""),
                    "dN": met.get("dN",""),
                    "dN S.E.": met.get("dN S.E.",""),
                    "dS": met.get("dS",""),
                    "dS S.E.": met.get("dS S.E.",""),
                    "dN/dS Ratio": met.get("dN/dS Ratio",""),
                    "Avg GC3 (%)": met.get("Avg GC3 (%)",""),
                    "Sequence No.": met.get("Sequence No.",""),
                })
    return pd.DataFrame(rows)

# =================== MAIN ===================
def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # 1) Parse & merge all DnaSP .out per (Species, Population)
    by_key = {}
    for species_folder in ROOT_DIR.iterdir():
        if not species_folder.is_dir():
            continue
        species = species_folder.name

        for out_file in sorted(species_folder.glob("*.out"), key=out_priority):
            try:
                raw = out_file.read_text(encoding="utf-8", errors="ignore")
            except Exception as e:
                print(f"Skip {out_file}: {e}")
                continue

            parsed = extract_fields_from_dnasp_out(raw, out_file.name, species)
            key = (parsed.get("Species", species), parsed.get("Population",""))

            if key not in by_key:
                by_key[key] = {
                    "Species": parsed.get("Species", species),
                    "Common name": common_name_for(species),
                    "Family": family_for(species),
                    "Population": parsed.get("Population",""),
                }
            by_key[key] = _merge_rows(by_key[key], parsed, src=out_file.name)

    df = pd.DataFrame(list(by_key.values()))

    # 2) Combined label
    df["Combined Common (Pop)"] = df.apply(
        lambda r: build_combined(r.get("Common name",""), r.get("Species",""), r.get("Population","")),
        axis=1
    )

    # 3) MEGA metrics + merge
    mega_df = collect_mega_metrics(MEGA_BASE_DIR)
    if not mega_df.empty:
        df = df.merge(mega_df, on=["Species","Population"], how="left")

    # 4) Ensure columns & order
    for col in COLUMNS:
        if col not in df.columns:
            df[col] = ""
    df = df[COLUMNS]

    # 5) Cast integer columns; keep the rest as text
    for col in INT_COLS:
        df[col] = pd.to_numeric(df[col], errors="coerce").astype("Int64")

    # 6) Write Excel with explicit text formatting for decimal/test/ratio columns
    try:
        with pd.ExcelWriter(
            OUTPUT_XLSX,
            engine="xlsxwriter",
            engine_kwargs={"options": {"strings_to_numbers": False}},
        ) as writer:
            df.to_excel(writer, index=False, sheet_name="Summary")
            wb = writer.book
            ws = writer.sheets["Summary"]

            # Autosize columns
            for i, col in enumerate(df.columns):
                ws.set_column(i, i, max(14, min(48, len(col) + 2)))

            # Integers as numbers
            int_fmt = wb.add_format({"num_format": "0"})
            for col_name in INT_COLS:
                col_idx = df.columns.get_loc(col_name)
                ws.set_column(col_idx, col_idx, None, int_fmt)

            # Force text for decimal/test/ratio + combined label
            text_fmt = wb.add_format({"num_format": "@"})
            for col_name in TEXT_DECIMAL_COLS + ["Combined Common (Pop)"]:
                col_idx = df.columns.get_loc(col_name)
                ws.set_column(col_idx, col_idx, None, text_fmt)
    except Exception as e:
        df.to_excel(OUTPUT_XLSX, index=False)
        print(f"Note: wrote without column formats ({e})")

    print(f"Done! Output written to {OUTPUT_XLSX}")

if __name__ == "__main__":
    main()
