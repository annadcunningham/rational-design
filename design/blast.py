""" Functions for NCBI blast. """

import os.path
import time
from Bio import Entrez
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastxCommandline, NcbiblastpCommandline

from basebio import join_subdir, ndiff

Entrez.email = 'annadcunningham@gmail.com'
ncbi_dbpath = 'ncbi-proteomes/'


def blast_protein_local(fastafile, dbname='human', evalue=0.001, subdir='TEMP'):
    """ BLASTs the protein sequence in the fastafile using local database.
        Writes a .xml file containing the BLAST results.
    """
    dbpath = ncbi_dbpath + dbname
    outfile = '{}_{}.xml'.format(fastafile.split('.fasta')[0], dbname)
    outfile = join_subdir(outfile, subdir)
    if not os.path.isfile(outfile):
        blastp_cline = NcbiblastxCommandline(cmd='blastp', query=fastafile, db=dbpath,
                                             evalue=evalue, outfmt=5, out=outfile)
        stdout, stderr = blastp_cline()
    return outfile

def blast_peptide_local(peptideseq, dbname='human', evalue=1000, subdir='TEMP'):
    """ BLASTs the peptide sequence (input as string) using local database.
        Writes a .xml file containing the BLAST results.
    """
    dbpath = ncbi_dbpath + dbname
    infile = '{}.fasta'.format(peptideseq)
    infile = join_subdir(infile, subdir)
    with open(infile, 'w') as F:
        F.write('>peptide\n')
        F.write(peptideseq)
    outfile = '{}_{}.xml'.format(peptideseq, dbname)
    outfile = join_subdir(outfile, subdir)
    if not os.path.isfile(outfile):
        blastp_cline = NcbiblastpCommandline(task='blastp-short', query=infile, db=dbpath,
                                             evalue=evalue, outfmt=5, out=outfile)
        stdout, stderr = blastp_cline()
    return outfile

def get_top_blast_subject(filename):
    """ From a .xml blast record, return the name and sequence of the top alignment.
    """
    with open(filename) as F:
        blast_record = NCBIXML.read(F)
        try:
            alignment = blast_record.alignments[0]
            subject_title = alignment.title
            subject_seq = alignment.hsps[0].sbjct
            return subject_title, subject_seq
        except:
            return None, None

def get_homologous_blast_subjects(filename, peptideseq):
    """ From a .xml blast record of short input sequence, return a list of
        accession numbers with only the proteins that are a
        perfect match to the query sequence.
    """
    with open(filename) as F:
        blast_record = NCBIXML.read(F)
    alignment_list = []
    peptide_list = []
    for alignment in blast_record.alignments:
        # first, check whether the protein is putative or a fragment
        if 'putative' in alignment.hit_def.lower() or \
           'fragment' in alignment.hit_def.lower():
            continue
        # this is a hacky way to get the match sequence (n.match) :/
        for n in alignment.hsps:
            pass
        if (len(n.match) == len(peptideseq)) and ('-' not in n.sbjct):
            alignment_list.append(alignment)
            peptide_list.append(n.sbjct)
    return [a.accession for a in alignment_list], peptide_list

def try_Entrez_efetch(accessionID):
    """ Tries Entrez.efetch until it works.
    """
    print('Entrez fetching {}'.format(accessionID), end='')
    try:
        handle = Entrez.efetch(db='protein', rettype='fasta', retmode='text', id=accessionID)
        print('\t\t.')
    except:
        handle = None
        print('\t\tx')
    return handle

def write_efetch_fasta(handle, fasta_filename):
    try:
        fasta = handle.read()
    except:
        fasta = ''
    if fasta != '':
        with open(fasta_filename, 'w') as F:
            F.write(fasta)
    return fasta_filename

def write_protein_fastas_from_accession_numbers(accession_list, peptide_list, accessions_to_skip=[], subdir='TEMP', entrez=True):
    """ Uses Entrez to search for the protein sequences of a given list of
        accession numbers. Writes a fasta if a sequence is found.
        Returns a list of fasta files.
    """
    accessions_with_proteins = []
    for accession, peptide in zip(accession_list, peptide_list):
        if accessions_to_skip != []:
            if accession in accessions_to_skip:
                continue
        fasta_filename = join_subdir('{}_{}.fasta'.format(accession, peptide), subdir)
        if entrez:
            if not os.path.isfile(fasta_filename):
                handle = try_Entrez_efetch(accession)
                write_efetch_fasta(handle, fasta_filename)
        if os.path.isfile(fasta_filename):
            print('Entrez {} fetched'.format(accession))
            accessions_with_proteins.append(fasta_filename)
    return accessions_with_proteins
