import os
import math
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

# 1) Your existing species-to-family mapping
species_family_mapping = {
    "Balaena_mysticetus": "Balaenidae (Right Whales)",
    "Eubalaena_australis": "Balaenidae (Right Whales)",
    "Eubalaena_glacialis": "Balaenidae (Right Whales)",
    "Eubalaena_japonica": "Balaenidae (Right Whales)",
    "Balaenoptera_musculus": "Balaenopteridae (Rorquals)",
    "Balaenoptera_physalus": "Balaenopteridae (Rorquals)",
    "Delphinus_delphis": "Delphinidae (Oceanic Dolphins)",
    "Globicephala_macrorhynchus": "Delphinidae (Oceanic Dolphins)",
    "Orcaella_brevirostris": "Delphinidae (Oceanic Dolphins)",
    "Orcinus_orca": "Delphinidae (Oceanic Dolphins)",
    "Peponocephala_electra": "Delphinidae (Oceanic Dolphins)",
    "Pseudorca_crassidens": "Delphinidae (Oceanic Dolphins)",
    "Stenella_attenuata": "Delphinidae (Oceanic Dolphins)",
    "Stenella_longirostris": "Delphinidae (Oceanic Dolphins)",
    "Tursiops_aduncus": "Delphinidae (Oceanic Dolphins)",
    "Tursiops_australis": "Delphinidae (Oceanic Dolphins)",
    "Tursiops_truncatus": "Delphinidae (Oceanic Dolphins)",
    "Phocoena_phocoena": "Phocoenidae (Porpoises)",
    "Phocoena_sinus": "Phocoenidae (Porpoises)",
    "Phocoena_spinipinnis": "Phocoenidae (Porpoises)",
    "Phocoenoides_dalli": "Phocoenidae (Porpoises)",
    "Neophocaena_asiaeorientalis": "Phocoenidae (Porpoises)",
    "Neophocaena_phocaenoides": "Phocoenidae (Porpoises)",
    "Delphinapterus_leucas": "Monodontidae (Beluga whale & Narwhal)",
    "Monodon_monoceros": "Monodontidae (Beluga whale & Narwhal)",
    "Hyperoodon_ampullatus": "Ziphiidae (Beaked Whales)",
    "Ziphius_cavirostris": "Ziphiidae (Beaked Whales)",
    "Mesoplodon_densirostris": "Ziphiidae (Beaked Whales)",
    "Mesoplodon_europaeus": "Ziphiidae (Beaked Whales)",
    "Mesoplodon_grayi": "Ziphiidae (Beaked Whales)",
    "Mesoplodon_mirus": "Ziphiidae (Beaked Whales)",
    "Physeter_catodon": "Physeteridae (Sperm Whales)",
    "Platanista_gangetica": "Platanistidae (River Dolphins)"
}

# 2) Paths to the source images and the combined outputs
SOURCE_ROOT = "sequences/Interactive_group_editor"
DEST_ROOT = "sequences/PNG_populations"


def main():
    # Invert the mapping to get family -> list of species
    family_dict = {}
    for species, family in species_family_mapping.items():
        family_dict.setdefault(family, []).append(species)

    # Ensure the output folder exists
    os.makedirs(DEST_ROOT, exist_ok=True)

    # Process each family
    for family_name, species_list in family_dict.items():
        images_for_family = []

        # Load each species' PNG
        for sp in species_list:
            genus = sp.split('_')[0]
            png_path = os.path.join(
                SOURCE_ROOT,
                genus,
                sp,
                f"{sp}_2d_scatterplot.png"
            )
            if not os.path.isfile(png_path):
                print(f"[Missing] {png_path}")
                continue

            try:
                img_data = mpimg.imread(png_path)
                images_for_family.append(img_data)
            except Exception as e:
                print(f"Could not read {png_path}: {e}")

        # If no images, skip
        if not images_for_family:
            print(f"No images found for family '{family_name}', skipping.")
            continue

        # Create exactly one combined figure for *all* images in this family
        combined_png_path = os.path.join(DEST_ROOT, f"{family_name}.png")
        create_combined_figure(images_for_family, family_name, combined_png_path)


def create_combined_figure(images, family_name, out_path):
    """
    Given a list of image arrays, arrange them in a grid.
    - If there are exactly 4 images, force a 2x2 layout.
    - Otherwise, use up to 3 columns and however many rows are needed.
    - Increase resolution by saving with dpi=300.
    """
    num_images = len(images)

    # If exactly 4 images, do a 2x2 grid
    if num_images == 4:
        rows, cols = 2, 2
    else:
        # Up to 3 columns otherwise
        cols = 3
        rows = math.ceil(num_images / cols)

    fig, axs = plt.subplots(
        nrows=rows,
        ncols=cols,
        figsize=(5 * cols, 4 * rows)  # You can adjust figsize as you like
    )

    # Title for entire figure
    fig.suptitle(f"Family: {family_name}", fontsize=16, y=0.98)

    # Handle shape issues if there's only one row or one column
    if rows == 1 and cols == 1:
        axs = [[axs]]
    elif rows == 1 or cols == 1:
        if rows == 1:
            axs = [axs]
        else:
            axs = [[ax] for ax in axs]

    # Flatten the 2D list of axes
    ax_list = [ax for row in axs for ax in row]

    # Place each image in a subplot
    for i, img_data in enumerate(images):
        ax_list[i].imshow(img_data)
        ax_list[i].axis("off")

    # Any leftover subplots get turned off
    for j in range(i + 1, len(ax_list)):
        ax_list[j].axis("off")

    # Adjust the layout so the suptitle doesn’t overlap subplots
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    # Increase the DPI to get higher resolution
    plt.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"Saved combined PNG for family '{family_name}' to: {out_path}")


if __name__ == "__main__":
    main()
