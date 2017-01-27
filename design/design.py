""" Command line tool for rational peptide design. """

from argparse import ArgumentParser, FileType
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

from basebio import join_subdir, organisms
from dmr import find_homologous_peptides, concatenate_homologous_peptides, calculate_heatmap, plot_heatmap
from report import build_report


def argument_parser():
    """ Parses input arguments. """
    parser = ArgumentParser()

    # PROTEIN INPUTS
    parser.add_argument('fasta1', type=FileType('r'),
                        help='Fasta file for protein 1.')
    parser.add_argument('fasta2', type=FileType('r'),
                        help='Fasta file for protein 2.')

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

    # Parse the protein names from the fasta files
    protein1_name = args.fasta1.name.split('/')[-1].split('.fasta')[0]
    protein2_name = args.fasta2.name.split('/')[-1].split('.fasta')[0]

    # Parse the peptide_length argument
    if '-' in args.peptidelength:
        length1, length2 = [int(n) for n in args.peptidelength.split('-')]
        peptide_lengths = range(length2, length1-1, -1)
    else:
        peptide_lengths = [int(args.peptidelength)]

    # Generate a list of homologous peptide pair candidates
    list_of_top_homologous_pairs_lists = []
    for peptide_length in peptide_lengths:
        top_homologous_pairs = find_homologous_peptides(args.fasta1.name,
                                                        args.fasta2.name,
                                                        peptide_length,
                                                        filter_by_cons=False)
        list_of_top_homologous_pairs_lists.append(top_homologous_pairs)
        print('{} peptides of length {}'.format(len(top_homologous_pairs), peptide_length))

    top_homologous_pairs = concatenate_homologous_peptides(list_of_top_homologous_pairs_lists)

    print('{} candidate peptides identified'.format(len(top_homologous_pairs)))

    for pair in top_homologous_pairs:
        print('{}\t{}\t{}'.format(round(pair.homology_score, 2), pair.seq1, pair.seq2))

    # top_homologous_pairs = top_homologous_pairs[:2]

    # Try to get the accession numbers for the proteins so they don't repeat
    # in the heatmap
    try:
        accession1 = top_homologous_pairs[0].peptide1.name.split('-')[1]
        accession2 = top_homologous_pairs[0].peptide2.name.split('-')[1]
        accessions_to_skip = [accession1, accession2]
    except:
        accessions_to_skip = None

    # Generate a heatmap for each candidate
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

    # Generate the report
    build_report(protein1_name, protein2_name,
                 top_homologous_pairs,
                 protein_names_for_report,
                 ['human'] + organisms,
                 PDB1=args.pdb1, PDB2=args.pdb2,
                 out=args.out.name,
                 subdir=args.temp)
