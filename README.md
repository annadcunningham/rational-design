# Rational Design
For DMR lab.

## Installation
Coming soon... planning to get everything set up in a Vagrant VM.

## General Workflow
1. Get a list of all possible peptide alignments with homology scores
    - TODO: Try a few different amino acid lengths (5-10), specified by user
    - Keep only peptide pairs with homology score greater than or equal to
      2 stdev above the mean
2. For each homologous peptide pair:
    - Calculate the conservation of each peptide across the list of organisms
    - Filter out the ones that are poorly conserved?
        This is currently not done, but can be done by changing
        filter_by_cons=True in design.py
3. For each homologous peptide pair:
    - Blast to get other proteins that contain that peptide
        (exclude putative and fragment proteins)
        TODO: Allow peptides that are not exactly the same, for long peptides
    - Calculate the conservation of that peptide in that protein across
        the list of organisms
    - Make a heatmap showing the 10 proteins with the most conserved
        homologous peptide
4. TODO: Find where the peptides are on the protein (linear map)
    - TODO: Show this summary on the first page of the report
    - TODO: Show it on each individual peptide page
    - TODO: Also calculate and show secondary structure?
5. The report includes:
    - Summary Page
        - Table of peptides
            - Homology score
            - Peptide sequence
            - Peptide position
            - Peptide conservation
    - Peptide pages
            - Alignment of each peptide
            - Heatmap
            - Structural location?
            - Table of protein descriptions
