#!/bin/bash
# commands for installing databases for ncbi-blast+
makeblastdb -in alligator8496.fasta -parse_seqids -dbtype prot -title alligator -out alligator -taxid 8496
makeblastdb -in ant443821.fasta -parse_seqids -dbtype prot -title ant -out ant -taxid 443821
makeblastdb -in candida237561.fasta -parse_seqids -dbtype prot -title candida -out candida -taxid 237561
makeblastdb -in cat9685.fasta -parse_seqids -dbtype prot -title cat -out cat -taxid 9685
makeblastdb -in cavefish7994.fasta -parse_seqids -dbtype prot -title cavefish -out cavefish -taxid 7994
makeblastdb -in chicken9031.fasta -parse_seqids -dbtype prot -title chicken -out chicken -taxid 9031
makeblastdb -in cow9913.fasta -parse_seqids -dbtype prot -title cow -out cow -taxid 9913
makeblastdb -in dog9615.fasta -parse_seqids -dbtype prot -title dog -out dog -taxid 9615
makeblastdb -in fly7227.fasta -parse_seqids -dbtype prot -title fly -out fly -taxid 7227
makeblastdb -in frog8364.fasta -parse_seqids -dbtype prot -title frog -out frog -taxid 8364
makeblastdb -in hamster10029.fasta -parse_seqids -dbtype prot -title hamster -out hamster -taxid 10029
makeblastdb -in human9606.fasta -parse_seqids -dbtype prot -title human -out human -taxid 9606
makeblastdb -in killifish8090.fasta -parse_seqids -dbtype prot -title killifish -out killifish -taxid 8090
makeblastdb -in mouse10090.fasta -parse_seqids -dbtype prot -title mouse -out mouse -taxid 10090
makeblastdb -in pig9823.fasta -parse_seqids -dbtype prot -title pig -out pig -taxid 9823
makeblastdb -in rat10116.fasta -parse_seqids -dbtype prot -title rat -out rat -taxid 10116
makeblastdb -in turtle13735.fasta -parse_seqids -dbtype prot -title turtle -out turtle -taxid 13735
makeblastdb -in worm135651.fasta -parse_seqids -dbtype prot -title worm -out worm -taxid 135651
makeblastdb -in yeast559292.fasta -parse_seqids -dbtype prot -title yeast -out yeast -taxid 559292
makeblastdb -in zebrafish7955.fasta -parse_seqids -dbtype prot -title zebrafish -out zebrafish -taxid 7955
