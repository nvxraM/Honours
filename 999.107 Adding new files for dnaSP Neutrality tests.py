from pathlib import Path

ROOT = Path("sequences/Species_POP_Mega_Haploid_split")
SUFFIXES = [" Fu's Fs.out", " Tajima's D.out"]

created = skipped = 0

for pop_file in ROOT.rglob("*_Population_*.out"):
    name = pop_file.name
    if name.endswith(" Fu's Fs.out") or name.endswith(" Tajima's D.out"):
        continue  # don't reprocess the files we create
    stem = pop_file.with_suffix("").name
    species_dir = pop_file.parent
    for suf in SUFFIXES:
        target = species_dir / f"{stem}{suf}"
        if target.exists():
            skipped += 1
        else:
            target.touch()
            created += 1

print(f"Done. Created {created} files, skipped {skipped} existing.")
