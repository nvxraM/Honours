#!/usr/bin/env python3
"""
find_nearest.py
"""
from Bio import Phylo  # import Biopython's tree parsing module

def main():
    # Load the Newick tree from file
    tree_file = "sequences/CDS_Reference/Phylogeny/ML_tree_with_brlen.nwk"
    tree = Phylo.read(tree_file, "newick")  # parse tree structure into a Phylo object

    # List of species (without accession prefixes) for which to find neighbours
    targets = [
        "Delphinus_delphis", "Tursiops_aduncus", "Tursiops_truncatus",
        "Stenella_longirostris", "Stenella_attenuata", "Tursiops_australis",
        "Orcaella_brevirostris", "Orcinus_orca", "Pseudorca_crassidens",
        "Globicephala_macrorhynchus", "Peponocephala_electra", "Delphinapterus_leucas",
        "Monodon_monoceros", "Neophocaena_asiaeorientalis", "Neophocaena_phocaenoides",
        "Phocoena_sinus", "Phocoena_spinipinnis", "Phocoena_phocoena",
        "Phocoenoides_dalli", "Ziphius_cavirostris", "Hyperoodon_ampullatus",
        "Mesoplodon_mirus", "Mesoplodon_europaeus", "Mesoplodon_grayi",
        "Mesoplodon_densirostris", "Platanista_gangetica", "Balaenoptera_musculus",
        "Balaenoptera_physalus", "Balaena_mysticetus", "Eubalaena_japonica",
        "Eubalaena_australis", "Eubalaena_glacialis", "Physeter_macrocephalus"
    ]

    # Extract all tip (leaf) nodes from the tree
    terminals = tree.get_terminals()

    # Build mapping: species name -> corresponding tip node
    species_map = {}
    for tip in terminals:
        # tip.name is like 'NC_012058.1_Tursiops_aduncus'; split into 3 parts
        parts = tip.name.split('_', 2)
        if len(parts) == 3:
            species_name = parts[2]  # keep only the species portion
            species_map[species_name] = tip

    # Dictionary to hold nearest neighbour info: species -> (neighbour, distance)
    nearest_info = {}
    for sp in targets:
        focal = species_map.get(sp)
        if focal is None:
            print(f"⚠️ WARNING: target '{sp}' not found in tree tips.")
            continue

        # Initialize minimum distance tracking
        min_dist = float('inf')
        min_neighbour = None

        # Compare focal species to every other target
        for other in targets:
            if other == sp:
                continue  # skip self
            other_tip = species_map.get(other)
            if other_tip is None:
                continue  # skip missing tips

            # Compute path length between the two tips
            dist = tree.distance(focal, other_tip)
            if dist < min_dist:
                min_dist = dist
                min_neighbour = other

        if min_neighbour is None:
            print(f"⚠️ No neighbours computed for '{sp}'")
            continue

        nearest_info[sp] = (min_neighbour, min_dist)

    # Print a human-readable list of nearest neighbours with distances
    print("# Closest neighbours with branch-length distances:")
    for sp, (nbr, dist) in sorted(nearest_info.items()):
        print(f"{sp} -> {nbr} (distance = {dist:.4f})")

    # Print a Python dict literal mapping species to its closest relative
    print("\n# Mapping from species to its closest relative:")
    print("closest_relatives = {")
    for sp, (nbr, _) in sorted(nearest_info.items()):
        print(f'    "{sp}": "{nbr}",')
    print("}")

if __name__ == "__main__":
    main()