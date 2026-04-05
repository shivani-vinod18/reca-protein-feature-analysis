# Can Stress-Survival Proteins Carry an Extremophile Signature?
## Sequence-Feature Analysis of Bacterial RecA Proteins

## Objective
This project explores whether RecA protein sequences from extremophile bacteria carry detectable sequence-derived compositional patterns when compared with non-extremophile controls. The goal was to move beyond standard phylogenetic workflow and test whether a stress-relevant protein could show an environmental signature in feature space.

## Biological Motivation
RecA is a core bacterial protein involved in DNA repair and recombination. Because DNA repair is especially important under radiation, thermal, and other stress conditions, RecA was selected as a biologically meaningful candidate for exploratory analysis in the context of extremophile adaptation and astrobiology.

## Dataset
The dataset contains 20 bacterial RecA proteins:
- 10 extremophiles
- 10 non-extremophile controls

The extremophiles include organisms associated with:
- radiation resistance
- thermophily
- halophily
- acidophily

## Methods
1. Retrieved RecA protein sequences from NCBI Protein in FASTA format.
2. Curated a dataset of 20 full-length bacterial RecA proteins.
3. Wrote a Python script to extract sequence-derived features from each protein, including:
   - protein length
   - acidic residue fraction
   - basic residue fraction
   - hydrophobic residue fraction
   - aromatic residue fraction
   - glycine fraction
   - proline fraction
   - full amino acid composition
4. Compared mean core feature values between extremophiles and controls.
5. Applied PCA to the feature matrix to visualize whether the proteins showed structure or separation in reduced-dimensional space.

## Outputs
- data/ contains the individual protein FASTA files and the combined FASTA file.
- scripts/feature_extraction.py computes sequence-derived protein features.
- scripts/group_feature_summary.py compares feature means between groups and generates a summary plot.
- scripts/pca_analysis.py performs PCA on the feature matrix.
- results/protein_metadata.csv contains organism and category labels.
- results/protein_feature_matrix.csv contains the extracted protein feature table.
- results/group_feature_summary.csv contains group-level mean feature comparisons.
- results/core_feature_comparison.png shows group-level comparison of selected core features.
- results/pca_coordinates.csv contains PCA coordinates for each protein.
- results/pca_plot.png visualizes the proteins in PCA space.
- results/pca_summary.txt contains explained variance information for the PCA.

## Key Findings
- Extremophile RecA proteins showed modest but directional differences in sequence-derived composition relative to controls.
- Hydrophobic fraction showed the strongest mean shift between groups in this dataset.
- PCA of the RecA feature space showed partial structure, with some extremophiles occupying more distinct regions while others overlapped with controls.
- The results suggest that RecA may carry weak-to-moderate environmental signal, but not a clean universal extremophile signature.
- Different extremophile classes did not collapse into a single cluster, suggesting that molecular adaptation patterns may be heterogeneous.

## Interpretation
This was an exploratory feature-based protein analysis rather than a claim of strong classification or mechanistic proof. The main result is that a stress-relevant protein can show detectable compositional trends across environmental groups, while still preserving substantial overlap and biological complexity.

## Relevance
Extremophiles are useful analogs for life under extraterrestrial stress conditions. A protein-level approach is relevant because adaptation to radiation, heat, salinity, and acidity may be reflected in sequence-derived molecular properties. This project helped connect computational sequence analysis with questions about biological survival in extreme environments.

## Limitations
This is a small exploratory study using one protein family. Stronger future work would include:
- additional stress-related proteins
- larger and more balanced datasets
- stronger statistical testing
- broader environmental classes
- integration with structural or functional annotation

## Future Directions
A natural next step would be to expand beyond RecA and test whether a multi-protein panel provides stronger separation between environmental groups than a single stress-relevant protein alone.