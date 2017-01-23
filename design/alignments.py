""" Functions and classes for sequence alignments. """

import os.path
import re
import numpy as np
from itertools import combinations
from Bio.Align.Applications import ClustalwCommandline
from Bio.SubsMat import MatrixInfo

from basebio import organism_dict, PairwiseAlignment, PeptideBase, ProteinBase


class MultipleSequenceAlignment():
    def __init__(self, alignment):
        self.alignment = alignment
        self.length = alignment.get_alignment_length()
        self.num_seq = len(alignment)

    def calc_alignment_score(self, matrix=MatrixInfo.blosum80):
        """ Calculates the alignment score, normalized to alignment length.
            The alignment score is calculated using the matrix score for each
            possible pair of amino acids in each column.
        """
        score = 0
        for n in range(self.length):
            column = self.alignment[:, n]
            score += calc_pairwise_column_score(column, matrix=matrix)
        return float(score) / self.length

    def get_organism_list(self):
        """ Returns a list of organisms that are in the alignment. """
        organism_list = []
        for record in self.alignment:
            organism = record.id.split('|')[-1].split('_')[-1]
            organism_list.append(organism)
        return organism_list

    def trim_by_organism(self, organism_list, usedict=False):
        """ Trims the alignment to only include the organisms in the specified
            list. Returns trimmed MultipleSequenceAlignment object.
            If specified, will use the dictionary to convert organism names to
            organism abbreviations.
        """
        def sortkey(list_, x):
            for item in list_:
                if x.endswith(item):
                    return True
            return False
        num_organisms = len(organism_list)
        if usedict:
            organism_ids = [organism_dict[organism] for organism in organism_list]
        else:
            organism_ids = organism_list
        self.alignment.sort(lambda x: sortkey(organism_ids, x.id), reverse=True)
        new_alignment = self.alignment[:num_organisms]
        return MultipleSequenceAlignment(new_alignment)

    def trim_by_gaps(self, gap_percent=10):
        """ Removes sequences for which the gap percent is higher than 10%.
            Returns a MultipleSequenceAlignment object.
        """
        humanseq = str(get_human_seq_from_alignment(self.alignment).seq)
        humangaps = humanseq.count('-')
        gapcutoff = (len(humanseq) - humangaps) * gap_percent/100. + humangaps
        self.alignment.sort(lambda x: str(x.seq).count('-') < gapcutoff, reverse=True)
        for n in range(len(self.alignment)):
            if str(self.alignment[n].seq).count('-') > gapcutoff:
                break
        return MultipleSequenceAlignment(self.alignment[:n])

    def remove_gaps_from_human(self):
        """ Removes columns that are gaps in the human sequence.
            Returns a MultipleSequenceAlignment object.
        """
        aln = self.alignment
        humanseq = str(get_human_seq_from_alignment(aln).seq)
        while re.search('-', humanseq):
            start, stop = re.search('-', humanseq).span()
            aln = aln[:, 0:start] + aln[:, stop:]
            humanseq = str(get_human_seq_from_alignment(aln).seq)
        return MultipleSequenceAlignment(aln)

### Functions relating to MSAs ###
def calc_matrix_score(aa1, aa2, matrix=MatrixInfo.blosum80):
    """ Gets the score of two amino acids in the given matrix.
        If there is a gap, 'X' is substituted.
    """
    aa1 = 'X' if aa1 == '-' else aa1
    aa2 = 'X' if aa2 == '-' else aa1
    try:
        return matrix[(aa1, aa2)]
    except:
        return matrix[(aa2, aa1)]

def calc_pairwise_column_score(column, matrix=MatrixInfo.blosum80):
    """ Averages the matrix score for all the possible pairs of amino acids
        in a column.
    """
    pairs = list(combinations(list(column), 2))
    scores = [calc_matrix_score(aa1, aa2, matrix=matrix) for aa1, aa2 in pairs]
    return float(sum(scores)) / len(column) # should this be len(pairs)? or not divide by anything?

def get_human_seq_from_alignment(alignment):
    """ Returns the human SeqRecord object from an MultipleSeqAlignment object.
    """
    for record in alignment:
        if record.description == 'HUMAN':
            return record
    raise ValueError('HUMAN sequence not found: {}'.format(alignment))

def clustal_MSA_local(fastafile):
    """ Uses locally installed clustalw2 to generate a MSA from the
        sequences in the given fastafile.
    """
    alignfile = fastafile.replace('fasta', 'aln')
    if not os.path.isfile(alignfile):
        clustalw2_cline = ClustalwCommandline("clustalw", infile=fastafile, pwmatrix='blosum')
        stdout, stderr = clustalw2_cline()
    return None

def keep_top_n_alignments(alignmentlist, num2keep):
    """ Returns the top num2keep alignments from a list of PairwiseAlignments.
    >>> alist = [Alignment('', '', 1), Alignment('', '', 10), Alignment('', '', 3)]
    >>> keep_top_n_alignments(alist, 1)
    [Alignment(pep1='', pep2='', score=10)]
    >>> alist.append(Alignment('', '', 10))
    >>> keep_top_n_alignments(alist, 1)
    [Alignment(pep1='', pep2='', score=10), Alignment(pep1='', pep2='', score=10)]
    """
    alignmentlist.sort(key=lambda x: x.score, reverse=True)
    min_score = alignmentlist[num2keep-1].score
    top_alignments = []
    for alignment in alignmentlist:
        if alignment.score >= min_score:
            top_alignments.append(alignment)
        else:
            break
    return top_alignments

def keep_top_std_alignments(alignmentlist, stdev=2):
    """ Returns the top alignments from a list of PairwiseAlignments, with
        scores higher than 2 stdevs above the mean.
    """
    alignmentlist.sort(key=lambda x: x.score, reverse=True)
    scoreslist = [a.score for a in alignmentlist]
    min_score = np.mean(scoreslist) + stdev*np.std(scoreslist)
    top_alignments = []
    for alignment in alignmentlist:
        if alignment.score >= min_score:
            top_alignments.append(alignment)
        else:
            break
    return top_alignments
