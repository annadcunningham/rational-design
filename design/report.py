import pandas as pd
from reportlab.lib.enums import TA_JUSTIFY
from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Frame, Image, Table, PageBreak, FrameBreak, PageTemplate, NextPageTemplate
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch

from peptides import Peptide, Protein, HomologousPeptidePair
from pdb import make_structure_map
from basebio import join_subdir

styles = getSampleStyleSheet()
styles.add(ParagraphStyle(name='Align', fontName='Courier', fontSize=8))
marginsize = 0.6*inch

def format_table(table, rows=[0], bold=True, style='Normal'):
    for row in rows:
        header = table[row]
        if bold:
            bold_header = [Paragraph('<b>{}</b>'.format(n), styles[style]) for n in header]
        else:
            bold_header = [Paragraph('{}'.format(n), styles[style]) for n in header]
        table[row] = bold_header
    return table

def parse_top_homologous_pairs_for_report(protein1_name, protein2_name, list_of_HomologousPeptidePairs):
    table = []
    table_header_0 = ['', protein1_name, '', '', protein2_name, '', '']
    table_header_1 = ['Homology Score', 'Sequence', 'Position', 'Conservation',
                      'Sequence', 'Position', 'Conservation']
    table.append(table_header_0)
    table.append(table_header_1)
    for pair in list_of_HomologousPeptidePairs:
        p1, p2 = pair.peptide1, pair.peptide2
        row = [str(round(pair.homology_score, 2)),
               p1.peptideseq, p1.position, str(round(pair.conservationscore1, 2)),
               p2.peptideseq, p2.position, str(round(pair.conservationscore2, 2))]
        table.append(row)
    return table

def align_similar_peptides(peptide, similar_peptides):
    aligned_peptide = '-{}-'.format(peptide)
    aligned_similar_peptides = []
    for pep in similar_peptides:
        length_diff = len(aligned_peptide) - len(pep)
        aligned_pep = ''
        ref_pep = aligned_peptide[:]
        if aligned_peptide.find(pep) == -1:
            if aligned_peptide.find(pep[1:]) != -1:
                aligned_pep = pep + '-'*length_diff
            else:
                aligned_pep = '-'*length_diff + pep
        else:
            aligned_pep = '-'*aligned_peptide.find(pep) + pep + '-'*(length_diff-1)
        aligned_similar_peptides.append(aligned_pep)
    return [aligned_peptide] + aligned_similar_peptides

def position_strings(position, peptide_length):
    pos1 = str(position)
    pos2 = str(position+peptide_length-1)
    spaces = peptide_length + 2 - len(pos1) - len(pos2)
    p = ' '*(peptide_length+2)
    pp = ' |' + ' '*(peptide_length-2) + '| '
    if len(pos1) == 1:
        p = p[0] + pos1 + p[2:]
    elif len(pos1) == 2:
        p = p[0] + pos1 + p[3:]
    else:
        p = pos1 + p[len(pos1):]
    if len(pos2) == 1:
        p = p[:-1] + pos2
    else:
        p = p[:-len(pos2)] + pos2
    return p, pp

def build_report(protein1_name, protein2_name,
                 list_of_HomologousPeptidePairs,
                 dict_of_heatmap_protein_names,
                 list_of_organisms,
                 PDB1=None,
                 PDB2=None,
                 out='test_report.pdf',
                 subdir='TEMP'):

    doc = SimpleDocTemplate(out, pagesize=letter,
                        rightMargin=marginsize, leftMargin=marginsize,
                        topMargin=marginsize, bottomMargin=marginsize)
                        #showBoundary=1)

    Story = []
    Story.append(NextPageTemplate('SummaryPage'))
    doc.addPageTemplates(
             PageTemplate(
                id='SummaryPage',
                frames=Frame(doc.leftMargin, doc.bottomMargin, doc.width, doc.height)
                )

        )
    title = 'Candidate peptides for {} and {}'.format(protein1_name, protein2_name)
    Story.append(Paragraph(title, styles['Title']))
    Story.append(Spacer(1, 12))

    Story.append(Paragraph('Summary of candidate peptides', styles['Heading3']))
    Story.append(Spacer(1, 12))

    table_ = parse_top_homologous_pairs_for_report(protein1_name, protein2_name,
                                                   list_of_HomologousPeptidePairs)
    table_ = format_table(table_, rows=[0,1], style='Align')
    table_ = format_table(table_, rows=range(2, len(table_)), bold=False, style='Align')
    Story.append(Table(table_, hAlign='LEFT'))

    for n in range(len(list_of_HomologousPeptidePairs)):
        Story.append(NextPageTemplate('PeptidePage'))
        Story.append(PageBreak())

        pair = list_of_HomologousPeptidePairs[n]

        # Put a title on the page
        title = 'Peptide {}: {} / {}'.format(n+1, pair.seq1, pair.seq2)
        Story.append(Paragraph(title, styles['Title']))
        Story.append(Spacer(1, 10))
        Story.append(Paragraph('Similar peptides:', styles['Heading3']))
        Story.append(FrameBreak())

        # Show the similar peptides
        peptide_list_1 = align_similar_peptides(pair.seq1, [p[0] for p in pair.similar_pairs])
        peptide_list_2 = align_similar_peptides(pair.seq2, [p[1] for p in pair.similar_pairs])
        similar_pairs_len = len(peptide_list_1) + 4

        for peptide_list, protein_name, peptide in zip([peptide_list_1, peptide_list_2],
                                                       [protein1_name, protein2_name],
                                                       [pair.peptide1, pair.peptide2]):
            Story.append(Paragraph('<u>{}</u>'.format(protein_name), styles['Align']))
            Story.append(Spacer(1, 5))
            positions, markers = position_strings(peptide.position, len(peptide.peptideseq))
            Story.append(Paragraph(positions.replace(' ', '&nbsp;'), styles['Align']))
            Story.append(Paragraph(markers.replace(' ', '&nbsp;'), styles['Align']))
            Story.append(Paragraph('<b>{}</b>'.format(peptide_list[0]), styles['Align']))
            for aligned_pep in peptide_list[1:]:
                Story.append(Paragraph(aligned_pep, styles['Align']))
            Story.append(FrameBreak())

        # Add structure information
        Story.append(Paragraph('Secondary structure and rASA', styles['Heading3']))
        if PDB1 is not None:
            Story.append(Paragraph('<u>{}</u>'.format(protein1_name), styles['Align']))
            Story.append(Spacer(1, 5))
            filename = join_subdir('{}_{}.png'.format(PDB1, n), subdir)
            pdb_seq1 = make_structure_map(pair.seq1, pair.peptide1.position, PDB1, filename)
            Story.append(Paragraph('Peptide seq: {}'.format(pair.seq1), styles['Align']))
            Story.append(Paragraph('PDB str seq: {}'.format(pdb_seq1), styles['Align']))
            image_ = Image(filename)
            image_._restrictSize(3 * inch, 3 * inch)
            Story.append(image_)

        if PDB2 is not None:
            Story.append(Paragraph('<u>{}</u>'.format(protein2_name), styles['Align']))
            Story.append(Spacer(1, 5))
            filename = join_subdir('{}_{}.png'.format(PDB2, n), subdir)
            pdb_seq2 = make_structure_map(pair.seq2, pair.peptide2.position, PDB2, filename)
            Story.append(Paragraph('Peptide seq: {}'.format(pair.seq2), styles['Align']))
            Story.append(Paragraph('PDB str seq: {}'.format(pdb_seq2), styles['Align']))
            image_ = Image(filename)
            image_._restrictSize(3 * inch, 3 * inch)
            Story.append(image_)

        Story.append(NextPageTemplate('HeatmapPage'))
        Story.append(PageBreak())

        # Add the two alignments
        Story.append(Spacer(1, 17))
        for organism in list_of_organisms:
            Story.append(Paragraph('<b>{}</b>'.format(organism), styles['Align']))
        Story.append(FrameBreak())

        alignment_len = max(pair.peptide1_MSA.num_seq, pair.peptide2_MSA.num_seq)
        alignment_names_dict_list = []
        for protein_name, alignment_, protein_homologs in zip(
                [protein1_name, protein2_name],
                [pair.peptide1_MSA, pair.peptide2_MSA],
                [pair.protein1.get_protein_homologs(), pair.protein2.get_protein_homologs()]
                ):
            Story.append(Paragraph('<u>{}</u>'.format(protein_name), styles['Align']))
            Story.append(Spacer(1, 5))
            sorted_indices = alignment_.get_sorted_indices(list_of_organisms)
            for si in sorted_indices:
                if si == -1:
                    Story.append(Spacer(1, 10))
                else:
                    line = alignment_.alignment[si]
                    line_id = '|'.join(line.id.split('|')[1:])
                    if line_id != '':
                        for alignment_name, description in protein_homologs.items():
                            if line_id in alignment_name:
                                alignment_names_dict_list.append((line_id.split('|')[0], description))
                        if len(line_id) > 32 - len(line.seq):
                            line_id = line_id[:32-len(line.seq)]
                    Story.append(Paragraph('{} {}'.format(line.seq, line_id), styles['Align']))
            Story.append(FrameBreak())
        # Add the heatmap
        image_ = Image(join_subdir('{}.png'.format(n), subdir))
        _, image_height = image_._restrictSize(7.5 * inch, 7.5 * inch)
        Story.append(image_)

        Story.append(NextPageTemplate('TablePage'))
        Story.append(PageBreak())

        # Add the list of proteins / protein names
        table_ = [['Accession ID', 'Homologous sequence', 'Protein description']]
        for accession, description in alignment_names_dict_list:
            table_.append([accession, '   ', Paragraph(description, styles['Normal'])])
        (protein_name_dict, ordered_keys) = dict_of_heatmap_protein_names[n]
        for accession in ordered_keys:
            if accession.upper() not in [protein1_name.upper(), protein2_name.upper()]:
                [peptide_seq, descriptions] = protein_name_dict[accession]
                table_.append([accession,
                               Paragraph(peptide_seq, styles['Align']),
                               Paragraph(', '.join(descriptions), styles['Normal'])
                               ])
        Story.append(Table(
                        format_table(table_),
                        hAlign='LEFT',
                        colWidths=[1.2*inch, 1.8*inch, 4.3*inch]
                        ))

        doc.addPageTemplates(
                [PageTemplate(
                    id='PeptidePage',
                    frames=build_frames_peptidepage(doc, similar_pairs_len)
                    ),
                 PageTemplate(
                    id='HeatmapPage',
                    frames=build_frames_heatmappage(doc, alignment_len)
                    ),
                 PageTemplate(
                    id='TablePage',
                    frames=Frame(doc.leftMargin, doc.bottomMargin, doc.width, doc.height)
                    )
                 ]
            )

    doc.build(Story)


def build_frames_peptidepage(doc, peptides_height):
    frame_height = 16*(peptides_height+6) + 12
    title_frame = Frame(doc.leftMargin, doc.height, doc.width, 1.1*inch)
    peptide_frame_1 = Frame(
                        doc.leftMargin,
                        doc.height - frame_height,
                        doc.width/2,
                        frame_height
                        )
    peptide_frame_2 = Frame(
                        doc.leftMargin + doc.width/2,
                        doc.height - frame_height,
                        doc.width/2,
                        frame_height
                        )
    structure_frame_1 = Frame(
                        doc.leftMargin,
                        doc.bottomMargin + 2.3*inch,
                        doc.width/2,
                        2.3*inch
                        )
    structure_frame_2 = Frame(
                        doc.leftMargin + doc.width/2,
                        doc.bottomMargin + 2.3*inch,
                        doc.width/2,
                        2.3*inch
                        )
    structure_frame_3 = Frame(
                        doc.leftMargin,
                        doc.bottomMargin,
                        doc.width/2,
                        2.3*inch
                        )
    structure_frame_4 = Frame(
                        doc.leftMargin + doc.width/2,
                        doc.bottomMargin,
                        doc.width/2,
                        2.3*inch
                        )
    return [title_frame, peptide_frame_1, peptide_frame_2,
            structure_frame_1, structure_frame_2, structure_frame_3, structure_frame_4]

def build_frames_heatmappage(doc, alignment_height):
    frame_height = 12*(alignment_height+4) + 14
    align_frame_0 = Frame(
                        doc.leftMargin,
                        doc.height - frame_height + marginsize,
                        1*inch,
                        frame_height
                        )
    align_frame_1 = Frame(
                        doc.leftMargin + 1*inch,
                        doc.height - frame_height + marginsize,
                        (doc.width-1*inch)/2,
                        frame_height
                        )
    align_frame_2 = Frame(
                        doc.leftMargin + 1*inch + (doc.width-1*inch)/2,
                        doc.height - frame_height + marginsize,
                        (doc.width-1*inch)/2,
                        frame_height
                        )
    heatmap_frame = Frame(
                        doc.leftMargin,
                        doc.bottomMargin,
                        doc.width,
                        doc.height - frame_height
                        )
    return [align_frame_0, align_frame_1, align_frame_2, heatmap_frame]
