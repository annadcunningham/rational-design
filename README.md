# Rational Design
For DMR lab.

This program (`design.py`) generates a list of candidate peptide
sequences for inhibition of a desired protein-protein interaction, using the
rational design method developed in the Mochly-Rosen lab.

Inputs: Two .fasta protein sequence files from UniProt.
Output: A .pdf report containing homology, alignment, and heatmap information
for the specified number of candidate peptides.

####Usage:
```
$ python3 design.py -p1 name_of_protein_1 -p2 name_of_protein_2 -l peptide_length -o name_of_output_report.pdf protein1.fasta protein2.fasta
```
For example:
```
$ python3 design.py -p1 PDK -p2 DPKC -l 5 -o pdk_dpkc_report.pdf test/pdk.fasta test/dpkc.fasta
```

## General Workflow
1. Get a list of all possible peptide alignments with homology scores
    - Try a few different amino acid lengths (ex. 5-10), specified by user
    - Keep only peptide pairs with homology score greater than or equal to
      2 stdev above the mean
    - Combine overlapping peptides into one entry
2. For each homologous peptide pair:
    - Calculate the conservation of each peptide across the list of organisms
    - Filter out the ones that are poorly conserved?
        This is currently not done, but can be done by changing
        filter_by_cons=True in design.py
3. For each homologous peptide pair:
    - Blast to get other proteins that contain that peptide
        (exclude putative and fragment proteins)
    - Calculate the conservation of that peptide in that protein across
        the list of organisms
    - Make a heatmap showing the 10 proteins with the most conserved
        homologous peptide
4. TODO: Find where the peptides are on the protein (linear map)
    - TODO: Show this summary on the first page of the report
    - TODO: Show it on each individual peptide page
    - TODO: Also calculate and show secondary structure?
    - TODO: Show a .pdb structure if available?
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
        - TODO: Structural location
        - Table of protein descriptions

## Installation
This program is meant to be run in a virtual machine (VM) with all the
prerequisites already installed. This allows anyone on any operating system to
run the program without complicated and OS-specific installation issues.
To start the VM, you will need Vagrant and VirtualBox.

1. Download and install Vagrant here: https://www.vagrantup.com/downloads.html

2. Download and install VirtualBox here: https://www.virtualbox.org/wiki/Downloads

3. Clone this repository and initialize the VM:
```
$ cd rational-design/vagrant
$ vagrant up
```
Initializing the VM may take a while (10-20 minutes).

4. Enter the VM and navigate to the `design` directory:
```
$ vagrant ssh
$ cd /design/
```
You should save your .fasta input files somewhere in this `design` directory so
you can access them within the VM.
