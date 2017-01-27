import re
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from Bio import SeqIO, pairwise2
from Bio.SubsMat import MatrixInfo

from basebio import organisms, organism_dict2, PeptideBase, ProteinBase, PairwiseAlignment, join_subdir
from peptides import Peptide, Protein, HomologousPeptidePair
from alignments import keep_top_n_alignments, keep_top_std_alignments
from blast import blast_peptide_local, get_homologous_blast_subjects, write_protein_fastas_from_accession_numbers

def get_protein_name_and_seq(fastafile):
    """ Parses a fasta file (containing a single entry) into a Protein namedtuple """
    fasta_parsed = list(SeqIO.parse(fastafile, 'fasta'))
    assert len(fasta_parsed) == 1, 'Fasta file contains multiple entries: {}'.format(fastafile)
    protein_name = fasta_parsed[0].id.replace('|', '-')
    protein_seq = str(fasta_parsed[0].seq)
    return Protein(ProteinBase(protein_name, protein_seq, fastafile))

def get_full_protein_name(fastafile):
    """ Get the full name of the protein from a fasta file. """
    with open(fastafile) as F:
        fasta = F.read()
    nameline = fasta.split('\n')[0]
    assert nameline.startswith('>'), 'Not a valid fasta file (line must start with \'>\'): {}'.format(nameline)
    names = re.split('Name:|Full=', nameline)[2::2]
    return [name.split(';')[0] for name in names]

def filter_homologous_pairs_by_conservation(list_of_homologous_pairs, use_mean=False):
    """ From a list of HomologousPeptidePairs, keeps only the pairs for which
        both conservation scores are > 1.
        If use_mean=True, keeps the pairs for which both conservation scores
        are above the mean score.
    """
    filtered_pairs = []
    if use_mean:
        mean_cons_1 = np.mean([pair.conservationscore1 for pair in list_of_homologous_pairs])
        mean_cons_2 = np.mean([pair.conservationscore2 for pair in list_of_homologous_pairs])
        for pair in list_of_homologous_pairs:
            if pair.conservationscore1 > mean_cons_1 and pair.conservationscore2 > mean_cons_2:
                filtered_pairs.append(pair)
    else:
        filtered_pairs = list_of_homologous_pairs

    filtered_pairs2 = []
    for pair in filtered_pairs:
        if pair.conservationscore1 > 1 and pair.conservationscore2 > 1:
            filtered_pairs2.append(pair)
    return filtered_pairs2

def align_peptide_along_protein(peptidebase, protein, matrix=MatrixInfo.blosum80, gap_penalty=-6):
    """ Aligns the peptide to every peptide-length peptide in sequence.
        Returns the top-scoring alignment.
    #>>> align_peptide_along_protein(testpeptide, testprotein)
    #Alignment(pep1=PeptideBase(name='pdk', seq='ALSTD'), pep2=PeptideBase(name='dpkc', seq='ALSTE'), score=20)
    """
    num_iter = protein.proteinlength - len(peptidebase.peptideseq)
    alignments = [PairwiseAlignment('', '', 0)]
    for n in range(num_iter):
        temp_peptide = PeptideBase(protein, protein.proteinseq[n : n + len(peptidebase.peptideseq)])
        temp_alignment = pairwise2.align.globalds(peptidebase.peptideseq, temp_peptide.peptideseq, matrix,
                                                  gap_penalty, gap_penalty)
        score = temp_alignment[0][2] / float(len(peptidebase.peptideseq))
        alignment = PairwiseAlignment(Peptide(peptidebase), Peptide(temp_peptide), score)
        alignments.append(alignment)
    topalignment = keep_top_n_alignments(alignments, 1)
    return topalignment[0]

def find_homologous_peptides(fasta1, fasta2, peptide_length, use_mean=False,
            matrix=MatrixInfo.blosum80, filter_by_cons=False):
    """ Returns a list of the top 10 homologous peptide Alignments.
    """
    protein1 = get_protein_name_and_seq(fasta1)
    protein2 = get_protein_name_and_seq(fasta2)
    alignments = []
    for n in range(protein1.proteinlength - peptide_length):
        peptide1 = PeptideBase(protein1, protein1.proteinseq[n : n + peptide_length])
        alignment1 = align_peptide_along_protein(peptide1, protein2, matrix=matrix)
        alignments.append(alignment1)
    top_pairwise_alignments = keep_top_std_alignments(alignments)
    top_homologous_pairs = [HomologousPeptidePair(aln) for aln in top_pairwise_alignments]
    if filter_by_cons:
        top_homologous_pairs = filter_homologous_pairs_by_conservation(top_homologous_pairs, use_mean=use_mean)
    return top_homologous_pairs

def concatenate_homologous_peptides(list_of_top_homologous_pairs_lists):
    """ Reduce multiple lists of homologous peptide pairs down to a single list
        with little repetition.
        When making a list of top homologous peptide pairs, longer peptides tend
        to have lower homology scores. This is why I opted to get a list of top
        homologous pairs of each length, then join them (rather than filtering
        all lengths by the top stdev together).
    """
    homologous_pairs_concat = []
    for top_homologous_pairs in list_of_top_homologous_pairs_lists:
        for pair1 in top_homologous_pairs:
            add_peptide = True
            for pair2 in homologous_pairs_concat:
                # check if the sequence is already included
                if pair1.seq1 in pair2.seq1:
                    pair2.similar_pairs.append((pair1.seq1, pair1.seq2))
                    add_peptide = False
                # if there is major overlap between this and a previous one
                elif pair1.seq1[1:] in pair2.seq1 or pair1.seq1[:-1] in pair2.seq1:
                    pair2.similar_pairs.append((pair1.seq1, pair1.seq2))
                    add_peptide = False
            if add_peptide:
                homologous_pairs_concat.append(pair1)
    return homologous_pairs_concat

def calculate_conservation_for_heatmap(fastafile, peptideseq, gap_penalty=-6, matrix=MatrixInfo.blosum80):
    """ Makes a MSA of the given protein in the fastafile, prunes it to the
        peptide seq, then calculates conservation of the peptide in each
        organism compared to human. Returns a dictionary of values.
    """
    protein = get_protein_name_and_seq(fastafile)
    peptide = Peptide(PeptideBase(protein, peptideseq))
    # get the peptide alignment object
    aln = peptide.get_peptide_MSA(trim_by_gaps=False).alignment
    # go through each entry and add the alignment score to the dictionary
    scores_dict = {}
    for record in aln:
        temp_alignment = pairwise2.align.globalds(peptideseq,
                                                  str(record.seq).replace('-', 'X'),
                                                  matrix, gap_penalty, gap_penalty)
        score = temp_alignment[0][2]
        try:
            organism = organism_dict2[record.id.split('_')[-1]]
        except:
            organism = record.id
        scores_dict[organism] = score
    return scores_dict

def calculate_heatmap(peptide, accessions_to_skip=None, axis=None, num2show=10):
    """ Calculates a heatmap of conservation of the peptide in the given list of
        proteins. Peptide is a peptide object.
        Returns a dataframe of conservation scores and a dictionary of protein
        accessions and their descriptions.
    """
    fname = blast_peptide_local(peptide.peptideseq)
    accession_list, peptide_list = get_homologous_blast_subjects(fname, peptide.peptideseq)
    print('Fetching {} proteins with homologous peptides'.format(len(accession_list)))
    fastas = write_protein_fastas_from_accession_numbers(accession_list, peptide_list, accessions_to_skip)
    dictdict = {}
    protein_name_dict = {}
    # calculate for protein of interest
    accession_label = peptide.protein.fasta.split('/')[-1].split('.fasta')[0]
    scores_dict = calculate_conservation_for_heatmap(peptide.protein.fasta, peptide.peptideseq)
    dictdict[accession_label] = scores_dict
    # calculate for other homologous proteins
    print('Assembling heatmap for homologous peptides in {} proteins'.format(len(fastas)))
    for fasta in fastas:
        print('Calculating conservation for {}'.format(fasta))
        fname = fasta.split('/')[-1].split('.fasta')[0]
        accession_label, peptide_seq = fname.split('_')
        scores_dict = calculate_conservation_for_heatmap(fasta, peptide_seq)
        dictdict[accession_label] = scores_dict
        protein_name_dict[accession_label] = [peptide_seq, get_full_protein_name(fasta)]
    heatmap_df = pd.DataFrame(dictdict)
    return (heatmap_df, protein_name_dict)

def plot_heatmap(df1, df2, figname, tempdir, num2show=30, pname1=None, pname2=None):
    df3 = pd.concat([df1, df2], axis=1)
    # weed out repeats
    accessions = df3.columns.values.tolist()
    double_accessions = []
    for accession in accessions:
        if (accessions.count(accession) > 1) and (accession not in double_accessions):
            df3['{}_'.format(accession)] = df3[accession].max(axis=1)
            double_accessions.append(accession)
    for accession in double_accessions:
        df3.pop(accession)
    # rename the columns
    df3.columns = [col.replace('_', '') for col in df3.columns.values.tolist()]
    # sort
    df3_sorted = df3.ix[['human'] + organisms, :]
    df3_sorted = df3_sorted.reindex_axis(df3_sorted.mean().sort_values(ascending=False).index, axis=1)
    df3_sorted[df3_sorted < 0] = 0
    if pname1 is not None and pname2 is not None:
        df3_columns = df3.columns.values.tolist()
        df3_new_columns = [pname1, pname2]
        df3_new_columns.extend([acc for acc in df3_columns if acc not in df3_new_columns])
        df3_sorted = df3_sorted.ix[:, df3_new_columns]
    df3_sorted = df3_sorted.ix[:, :num2show]
    # make the figure
    fig = plt.figure(figsize=(11,5))
    ax1 = fig.add_subplot(1,1,1)
    sns.heatmap(df3_sorted, vmin=0, vmax=max(df3.max()), ax=ax1)
    plt.yticks(rotation=0)
    plt.xticks(rotation=45, ha='right')
    print('Saving heatmap image {}.png'.format(figname))
    plt.savefig(join_subdir('{}.png'.format(figname), tempdir))
    return df3_sorted, df3_sorted.columns.values.tolist()[:num2show]


#if __name__ == '__main__':
    """
    NOTES:
    - make a map showing the position? possibly with a predicted secondary structure?
    - user input for how many peptides to show?
    - sort the heatmap first by human sequence?
    - show the two peptides on top of each other in the summary table
    - try reversing peptide sequence?
    - something to calculate the heatmap for fewer peptides when the peptides are short...
    """
