""" Functions and classes for peptides and proteins. """

import os.path
import re
from Bio import AlignIO, SeqRecord, Seq, Align

from blast import blast_protein_local, get_top_blast_subject
from alignments import MultipleSequenceAlignment, clustal_MSA_local, get_human_seq_from_alignment
from basebio import organisms, PeptideBase, ProteinBase, join_subdir


class Protein():
    def __init__(self, ProteinBase):
        self.name = ProteinBase.name
        self.proteinseq = ProteinBase.proteinseq
        self.fasta = ProteinBase.fasta
        self.proteinlength = len(self.proteinseq)

    def write_blast_fasta(self, subdir='TEMP'):
        """ Blasts the protein sequence against organisms in the database and
            writes a fasta file of the results. Returns the fasta filename.
        """
        fastafile = '{}_blast.fasta'.format(self.name)
        fastafile = join_subdir(fastafile, subdir)
        if not os.path.isfile(fastafile):
            with open(fastafile, 'w') as F:
                F.write('>HUMAN\n{}\n'.format(self.proteinseq))
                for organism in organisms:
                    hit_name, hit_seq = get_top_blast_subject(
                                            blast_protein_local(
                                                self.fasta, dbname=organism))
                    if (hit_name is not None) and (hit_seq is not None):
                        F.write('>{}\n{}\n'.format(hit_name, hit_seq))
        return fastafile

    def get_protein_MSA(self, organism_list=None):
        """ Uses clustalw2 to make a MSA for this protein.
            Generates the MSA and returns the MultipleSequenceAlignment object.
        """
        fastafile = self.write_blast_fasta()
        alnfile = fastafile.replace('fasta', 'aln')
        clustal_MSA_local(fastafile)
        try:
            aln = AlignIO.read(alnfile, 'clustal')
        except: # if alnfile is empty (only the human sequence)
            record = SeqRecord.SeqRecord(Seq.Seq(self.proteinseq), description='HUMAN')
            aln = Align.MultipleSeqAlignment([record])
        aln = MultipleSequenceAlignment(aln)
        if organism_list:
            aln = aln.trim_by_organism(organism_list)
        return aln


class Peptide():
    def __init__(self, PeptideBase):
        self.protein = PeptideBase.protein
        self.name = PeptideBase.protein.name
        self.proteinseq = PeptideBase.protein.proteinseq
        self.proteinfasta = PeptideBase.protein.fasta
        self.peptideseq = PeptideBase.peptideseq
        self.peptidelength = len(self.peptideseq)
        self.position = self.proteinseq.find(self.peptideseq) + 1

    def get_peptide_MSA(self, organism_list=None, trim_by_gaps=True):
        """ Gets the protein MSA and crops it to just the relevant peptide region.
            Returns the peptide MultipleSequenceAlignment object.
        """
        aln = self.protein.get_protein_MSA(organism_list=organism_list)
        if trim_by_gaps:
            aln_trimmed = aln .trim_by_gaps()
        else:
            aln_trimmed = aln
        # find coordinates of peptide in human sequence
        humanseq = get_human_seq_from_alignment(aln_trimmed.alignment)
        re_pattern = self.peptideseq.replace('', '-*')[2:-2]
        try:
            (start, end) = re.search(re_pattern, str(humanseq.seq)).span()
        except:
            print('Could not find seq {} in {}'.format(self.peptideseq, str(humanseq.seq)))
        # crop the alignment
        cropped_MSA = aln_trimmed.alignment[:, start:end]
        return MultipleSequenceAlignment(cropped_MSA)

    def get_conservation_score(self):
        """ Gets the MSA score from the cropped peptide MSA. """
        return self.get_peptide_MSA().calc_alignment_score()

    def get_relative_conservation_score(self, organism_list=None):
        """ Peptide MSA score divided by the protein MSA score. """
        protein_conservation = self.protein.get_protein_MSA(organism_list=organism_list).calc_alignment_score()
        return self.get_conservation_score() / float(protein_conservation)


class HomologousPeptidePair():
    def __init__(self, PairwiseAlignment):
        self.homology_alignment = PairwiseAlignment
        self.peptide1 = PairwiseAlignment.pep1
        self.peptide2 = PairwiseAlignment.pep2
        self.protein1 = PairwiseAlignment.pep1.protein
        self.protein2 = PairwiseAlignment.pep2.protein
        self.seq1 = PairwiseAlignment.pep1.peptideseq
        self.seq2 = PairwiseAlignment.pep2.peptideseq
        self.homology_score = PairwiseAlignment.score

    def get_organisms_with_both_proteins(self):
        """ Returns a list of organisms that are present in both protein MSAs.
        """
        p1_alignment = self.protein1.get_protein_MSA()
        p2_alignment = self.protein2.get_protein_MSA()
        p1_orglist = p1_alignment.get_organism_list()
        p2_orglist = p2_alignment.get_organism_list()
        shared_organisms = list(set(p1_orglist) & set(p2_orglist))
        return shared_organisms

    @property
    def shared_organisms(self):
        return self.get_organisms_with_both_proteins()

    @property
    def peptide1_MSA(self):
        """ Gets the peptide1 MSA with no gaps in the human sequence. """
        peptide1_MSA = self.peptide1.get_peptide_MSA(organism_list=self.shared_organisms)
        return peptide1_MSA.remove_gaps_from_human()

    @property
    def peptide2_MSA(self):
        """ Gets the peptide2 MSA with no gaps in the human sequence. """
        peptide2_MSA = self.peptide2.get_peptide_MSA(organism_list=self.shared_organisms)
        return peptide2_MSA.remove_gaps_from_human()

    @property
    def conservationscore1(self):
        """ Calculates the alignment score of the peptide1 MSA divided by the
            alignment score for the protein1 MSA.
        """
        prot1_MSA = self.protein1.get_protein_MSA(organism_list=self.peptide1_MSA.get_organism_list())
        prot1_cons = prot1_MSA.calc_alignment_score()
        pep1_cons = self.peptide1_MSA.calc_alignment_score()
        return float(pep1_cons) / prot1_cons

    @property
    def conservationscore2(self):
        """ Calculates the alignment score of the peptide2 MSA divided by the
            alignment score for the protein2 MSA.
        """
        prot2_MSA = self.protein2.get_protein_MSA(organism_list=self.peptide2_MSA.get_organism_list())
        prot2_cons = prot2_MSA.calc_alignment_score()
        pep2_cons = self.peptide2_MSA.calc_alignment_score()
        return float(pep2_cons) / prot2_cons
