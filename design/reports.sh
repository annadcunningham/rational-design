#!/bin/bash
python3 design.py -p1 Drp1 -p2 Mff -l 5-10 -o drp1_mff.pdf -pdb1 4BEJ test/drp1.fasta test/mff.fasta
python3 design.py -p1 PDK -p2 DPKC -l 5-10 -o pdk_dpkc.pdf -pdb1 4MP2 -pdb2 1BDY test/pdk.fasta test/dpkc.fasta
python3 design.py -p1 Drp1 -p2 Fis1 -l 5-10 -o drp1_fis1.pdf -pdb1 4BEJ -pdb2 1NZN test/drp1.fasta test/fis1.fasta
