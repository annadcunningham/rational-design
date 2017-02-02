""" Command line tool for rational peptide design. """

from argparse import ArgumentParser, FileType
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

from basebio import join_subdir, organisms
from blast import try_Entrez_efetch, write_efetch_fasta
from dmr import find_homologous_peptides, concatenate_homologous_peptides, calculate_heatmap, plot_heatmap
from report import build_report


def argument_parser():
    """ Parses input arguments. """
    parser = ArgumentParser()

    # PROTEIN INPUTS
    # User can either input:
    # 1. Protein names and accession numbers from Uniprot
    parser.add_argument('-p1', type=str, nargs=2, required=False, default=None,
                        help='Protein name and Uniprot accession number for protein 1.')
    parser.add_argument('-p2', type=str, nargs=2, required=False, default=None,
                        help='Protein name and Uniprot accession number for protein 2.')
    # 2. Fasta files downloaded from uniprot
    parser.add_argument('--fasta1', type=FileType('r'), required=False, default=None,
                        help='Fasta file for protein 1 (please download from Uniprot).')
    parser.add_argument('--fasta2', type=FileType('r'), required=False, default=None,
                        help='Fasta file for protein 2 (please download from Uniprot).')

    # PEPTIDE INPUTS
    parser.add_argument('-l', '--peptidelength', type=str,
                        required=True, help=('Length of desired peptide. Input '
                        'a single integer (i.e. 6) or a range (i.e. 5-10).'))
    # STRUCTURE INPUTS
    parser.add_argument('-pdb1', type=str, required=False, default=None,
                        help='PDB ID of protein 1.')
    parser.add_argument('-pdb2', type=str, required=False, default=None,
                        help='PDB ID of protein 2.')
    # ENTREZ OPTION
    parser.add_argument('--no-entrez', action='store_false', required=False,
                        help='Use this flag to skip Entrez fetching if all files are already fetched.')

    # OUTPUT SPECIFICATIONS
    parser.add_argument('-o', '--out', type=FileType('w'), required=False,
                        default='output_report.pdf',
                        help='Name of output file, as a .pdf (default is output_report.pdf)')
    parser.add_argument('-temp', type=str, required=False, default='TEMP',
                        help='Name of directory for writing temporary files.')

    return parser.parse_args()

if __name__ == '__main__':
    args = argument_parser()

    ## 1. Get the protein names and fasta files
    fasta1, fasta2 = None, None
    protein1_name, protein2_name = None, None
    accession1, accession2 = None, None
    #  If protein names and accessions are specified:
    if args.p1 is not None:
        protein1_name = args.p1[0]
        accession1 = args.p1[1]
        fasta1 = write_efetch_fasta(
                    try_Entrez_efetch(accession1),
                    '{}.fasta'.format(protein1_name)
                    )
    if args.p2 is not None:
        protein2_name = args.p2[0]
        accession2 = args.p2[1]
        fasta2 = write_efetch_fasta(
                    try_Entrez_efetch(accession2),
                    '{}.fasta'.format(protein2_name)
                    )
    #  If fasta files are specified:
    if args.fasta1 is not None:
        fasta1 = args.fasta1.name
        protein1_name = fasta1.split('/')[-1].split('.fasta')[0]
    if args.fasta2 is not None:
        fasta2 = args.fasta2.name
        protein2_name = fasta2.split('/')[-1].split('.fasta')[0]
    #  Check that proteins and fastas are defined
    assert fasta1 != None, 'Error: No fasta file for protein 1'
    assert fasta2 != None, 'Error: No fasta file for protein 2'
    assert protein1_name != None, 'Error: No name for protein 1'
    assert protein2_name != None, 'Error: No name for protein 1'

    ## 2. Parse the peptide_length argument
    if '-' in args.peptidelength:
        length1, length2 = [int(n) for n in args.peptidelength.split('-')]
        peptide_lengths = range(length2, length1-1, -1)
    else:
        peptide_lengths = [int(args.peptidelength)]

    ## 3. Generate a list of homologous peptide pair candidates
    list_of_top_homologous_pairs_lists = []
    for peptide_length in peptide_lengths:
        top_homologous_pairs = find_homologous_peptides(fasta1,
                                                        fasta2,
                                                        peptide_length,
                                                        filter_by_cons=False)
        list_of_top_homologous_pairs_lists.append(top_homologous_pairs)
        print('{} peptides of length {}'.format(len(top_homologous_pairs), peptide_length))

    top_homologous_pairs = concatenate_homologous_peptides(list_of_top_homologous_pairs_lists)

    print('{} candidate peptides identified'.format(len(top_homologous_pairs)))

    for pair in top_homologous_pairs:
        print('{}\t{}\t{}'.format(round(pair.homology_score, 2), pair.seq1, pair.seq2))

    #top_homologous_pairs = top_homologous_pairs[:2]

    ## 3. Generate a heatmap for each candidate
    #  Try to get the accession numbers for the proteins so they don't repeat in the heatmap
    if accession1 is None:
        try:
            accession1 = top_homologous_pairs[0].peptide1.name.split('-')[1]
        except:
            pass
    if accession2 is None:
        try:
            accession2 = top_homologous_pairs[0].peptide2.name.split('-')[1]
        except:
            pass
    accessions_to_skip = [A for A in [accession1, accession2] if A is not None]

    #  Make the heatmap
    protein_names_for_report = {}
    for n in range(len(top_homologous_pairs)):
        pair = top_homologous_pairs[n]
        print('Generating heatmap for peptide {}'.format(n+1))
        df1, protein_names_1 = calculate_heatmap(pair.peptide1, entrez=args.no_entrez)
        df2, protein_names_2 = calculate_heatmap(pair.peptide2, entrez=args.no_entrez)
        heatmap_df, ordered_keys = plot_heatmap(
                df1, df2, n, args.temp,
                pname1=protein1_name,
                pname2=protein2_name
                )
        protein_names_1.update(protein_names_2)
        protein_names_for_report[n] = (protein_names_1, ordered_keys)

    ## 4. Generate the report
    build_report(protein1_name, protein2_name,
                 top_homologous_pairs,
                 protein_names_for_report,
                 ['human'] + organisms,
                 PDB1=args.pdb1, PDB2=args.pdb2,
                 out=args.out.name,
                 subdir=args.temp)

    """
    NOTES:
    - make a linear map showing the position of each peptide?
    - show the two peptides on top of each other in the summary table
    - try reversing peptide sequence during initial alignment
    """
