from subprocess import call
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

from Bio.PDB import PDBParser, PDBList
from Bio.PDB.DSSP import make_dssp_dict, DSSP
from Bio.SubsMat import MatrixInfo
from Bio import pairwise2


ASA = {
        'A': 106.0, 'R': 248.0, 'N': 157.0, 'D': 163.0, 'C': 135.0, 'Q': 198.0,
        'E': 194.0, 'G': 84.0,  'H': 184.0, 'I': 169.0, 'L': 164.0, 'K': 205.0,
        'M': 188.0, 'F': 197.0, 'P': 136.0, 'S': 130.0, 'T': 142.0, 'W': 227.0,
        'Y': 222.0, 'V': 142.0
      }


def fetch_PDB_file(PDB_ID):
    """ Fetches the given PDB file.
    """
    pdb1 = PDBList()
    outfile = pdb1.retrieve_pdb_file(PDB_ID, pdir='PDB')
    return outfile

def get_PDB_structure(PDB_ID, pdbfile):
    """ Gets the PDB structure object from a PDB file.
    """
    p = PDBParser()
    structure = p.get_structure(PDB_ID, pdbfile)
    return structure

def get_PDB_DSSP_dict(structure, pdbfile):
    """ Uses Biopython to get dssp information from a PDB structure object and
        the PDB file.
    """
    dssp_dict = {}
    model = structure[0]
    dssp = DSSP(model, pdbfile)
    for key in list(dssp.keys()):
        pos, res, ss, sa, rsa, phi, psi, *ignore = dssp[key]
        dssp_dict[pos] = (res, ss, rsa)
    return dssp_dict

def generate_DSSP_file(pdbfile):
    """ Uses subprocess.call to run DSSP on a PDB file.
        Returns the .dssp file name.
    """
    pdbname = pdbfile.split('.')[0]
    dssp_output = pdbname + '.dssp'
    call(['dssp', '-i', pdbfile, '-o', dssp_output])
    return dssp_output

def parse_DSSP_dict(dsspfile):
    """ Parses a .dssp file into a dictionary containing relevant information:
            aa  amino acid (1-letter abbreviation)
            ss  secondary structure
            rsa relative accessible surface area
    """
    dssp_dict = make_dssp_dict(dsspfile)[0]
    parsed_dict = {}
    n = 0
    for key, value in sorted(dssp_dict.items()):
        chain, (_, pos, _) = key
        (aa, ss, sa, phi, psi, dssp_index, *ignore) = value
        rsa = sa / ASA[aa]
        n += 1
        parsed_dict[n] = (aa, ss, rsa, pos, phi, psi, dssp_index)
    return parsed_dict

def trydict(dict_, key, ind):
    try:
        return dict_[key][ind]
    except:
        return '-'

def make_structure_map(peptideseq, peptidepos, proteinseq, PDB_ID, filename):
    """ Makes a .png plotting the relative accessible surface area and
        secondary structure of the given peptide.
        Secondary structure abbreviations are as follows:
            H  Alpha helix (4-12)
            B  Isolated beta-bridge residue
            E  Strand
            G  3-10 helix
            I  pi helix
            T  Turn
            S  Bend
            -  None
        I simplified to:
            A  Alpha helix
            B  Beta strand
            L  Loop/turn
    """
    pdbfile = fetch_PDB_file(PDB_ID)
    dssp = parse_DSSP_dict(generate_DSSP_file(pdbfile))
    peptiderange = [r+peptidepos for r in range(len(peptideseq))]

    prot_to_PDB = align_PDB_to_seq(proteinseq, dssp)
    pdb_peptiderange = [prot_to_PDB[r] for r in peptiderange]

    pdb_sequence = [trydict(dssp, r, 0) for r in pdb_peptiderange]
    secondary_structure = [trydict(dssp, r, 1) for r in pdb_peptiderange]
    ss_ = ''.join(secondary_structure)
    ss_ = ss_.replace('H', 'A').replace('G', 'A').replace('I', 'A') # Alpha-helices
    ss_ = ss_.replace('E', 'B') # Beta-sheets
    ss_ = ss_.replace('T', 'L').replace('S', 'L') # Loops
    secondary_structure = list(ss_)
    accessible_surface_area = [trydict(dssp, r, 2) for r in pdb_peptiderange]
    accessible_surface_area = [-0.1 if a == '-' else a for a in accessible_surface_area]

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    # Plot RASA
    nrange = range(len(peptideseq))
    plt.bar(nrange, accessible_surface_area)
    # Plot secondary structure
    markerdict = {'A': 'o', 'B': 's', 'L': '_', '-': ''}
    for n in nrange:
        plt.scatter(n, -0.2, c='black', marker=markerdict[secondary_structure[n]], s=200)
    # Add labels
    ticks = ['{}\n{}'.format(peptideseq[n], pdb_sequence[n]) for n in nrange]
    ax.set_xticklabels(ticks)
    ax.set_xticks(nrange)
    plt.ylim([-0.25, 1])
    plt.xlim([-0.5, len(peptideseq)-0.5])
    plt.ylabel('Relative accessible surface area')
    ax.spines['bottom'].set_position('zero')
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    # Add legend
    plt.scatter(20, 20, c='black', marker='o', s=100, label='alpha-helix')
    plt.scatter(20, 20, c='black', marker='s', s=100, label='beta-sheet')
    plt.scatter(20, 20, c='black', marker='_', s=100, label='loop')
    ax.legend()

    plt.savefig(filename)
    return ''.join(pdb_sequence)

def align_PDB_to_seq(proteinseq, PDB_dssp_dict):
    """ The PDB structures are sometimes indexed weird.
        This will return a dictionary of protein_pos: PDB_pos
    """
    prot_to_PDB = {}

    pdb_seq = ''
    for pos in PDB_dssp_dict:
        pdb_seq += PDB_dssp_dict[pos][0]

    aln = pairwise2.align.globalds(proteinseq, pdb_seq, MatrixInfo.blosum80, -5, -1)[0]
    protein_aln = aln[0]
    pdb_aln = aln[1]

    # cycle through the alignment
    proteinseq_pos = 0
    pdbseq_pos = 0
    for n in range(len(protein_aln)):
        pdb_pdb = '-'
        if pdb_aln[n] != '-':
            pdbseq_pos += 1
            pdb_pdb = pdbseq_pos
        if protein_aln[n] != '-':
            proteinseq_pos += 1
            prot_to_PDB[proteinseq_pos] = pdb_pdb

    return prot_to_PDB
