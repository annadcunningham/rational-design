""" Base definitions and functions. """

import os.path
from collections import namedtuple


mammals = ['mouse', 'rat', 'hamster', 'dog', 'cat', 'cow', 'pig']
chicken = ['chicken']
reptiles = ['alligator', 'frog', 'turtle']
fish = ['cavefish', 'killifish', 'zebrafish']
insect = ['worm', 'fly', 'ant']
yeast = ['yeast', 'candida']
organisms = mammals + chicken + reptiles + fish + insect + yeast
# for now - just use a few organisms
# organisms = ['mouse', 'rat', 'chicken', 'frog', 'fly', 'yeast']
proteome_names = ['MOUSE', 'RAT', 'CRIGR', 'CANLF', 'FELCA', 'BOVIN', 'PIG',
                 'CHICK', 'ALLMI', 'XENTR', 'PELSI', 'ASTMX', 'ORYLA', 'DANRE',
                 'CAEBE', 'DROVI', 'CERBI', 'YEAST', 'CANAL']
organism_dict = dict(zip(organisms, proteome_names))
organism_dict2 = dict(zip(proteome_names, organisms))
organism_dict['human'] = 'HUMAN'
organism_dict2['HUMAN'] = 'human'

ProteinBase = namedtuple('ProteinBase', ['name', 'proteinseq', 'fasta'])
PeptideBase = namedtuple('PeptideBase', ['protein', 'peptideseq'])
PairwiseAlignment = namedtuple('PairwiseAlignment', ['pep1', 'pep2', 'score']) # pep1 and pep2 are PeptideBases

def join_subdir(filename, subdir):
    if filename.startswith(subdir):
        return filename
    else:
        basename = os.path.basename(filename)
        newfilename = os.path.join(subdir, basename)
        return newfilename

def ndiff(str1, str2):
    """ Calculates the number of differences between the two strings.
    """
    diffs = 0
    for m, n in zip(list(str1), list(str2)):
        if m != n:
            diffs += 1
    return diffs
