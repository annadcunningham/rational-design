import pandas as pd
from reportlab.lib.enums import TA_JUSTIFY
from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Frame, Image, Table, PageBreak, FrameBreak, PageTemplate, NextPageTemplate
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch

from peptides import Peptide, Protein, HomologousPeptidePair
from basebio import join_subdir

styles = getSampleStyleSheet()
styles.add(ParagraphStyle(name='Align', fontName='Courier'))

def make_table_header_bold(table, rows=[0]):
    for row in rows:
        header = table[row]
        bold_header = [Paragraph('<b>{}</b>'.format(n), styles['Normal']) for n in header]
        table[row] = bold_header
    return table

def make_table_wrap(table):
    newtable = [table[0]]
    for line in table[1:]:
        line = [Paragraph(n, styles['Normal']) for n in line]
        newtable.append(line)
    return newtable

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

def build_report(protein1_name, protein2_name,
                 list_of_HomologousPeptidePairs,
                 dict_of_heatmap_protein_names,
                 out='test_report.pdf',
                 subdir='TEMP'):

    marginsize = 0.6*inch
    doc = SimpleDocTemplate(out, pagesize=letter,
                        rightMargin=marginsize, leftMargin=marginsize,
                        topMargin=marginsize, bottomMargin=marginsize)

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

    heading1 = 'Table of candidate peptides'
    Story.append(Paragraph(heading1, styles['Heading3']))
    Story.append(Spacer(1, 12))

    table_ = parse_top_homologous_pairs_for_report(protein1_name, protein2_name,
                                                   list_of_HomologousPeptidePairs)
    Story.append(Table(make_table_header_bold(table_, rows=[0,1]), hAlign='LEFT'))

    for n in range(len(list_of_HomologousPeptidePairs)):
        Story.append(NextPageTemplate('PeptidePage'))
        Story.append(PageBreak())

        pair = list_of_HomologousPeptidePairs[n]
        # Put a title on the page
        title = 'Peptide {}: {} / {}'.format(n+1, pair.seq1, pair.seq2)
        Story.append(Paragraph(title, styles['Title']))
        Story.append(FrameBreak())
        # Add the two alignments
        alignment_len = max(pair.peptide1_MSA.num_seq, pair.peptide2_MSA.num_seq)
        for alignment_ in [str(pair.peptide1_MSA.alignment), str(pair.peptide2_MSA.alignment)]:
            for line in alignment_.split('\n')[1:]:
                Story.append(Paragraph(line, styles['Align']))
            Story.append(FrameBreak())
        # Add the heatmap
        image_ = Image(join_subdir('{}.png'.format(n), subdir))
        _, image_height = image_._restrictSize(6 * inch, 6 * inch)
        Story.append(image_)

        Story.append(NextPageTemplate('TablePage'))
        Story.append(PageBreak())

        # Add the list of proteins / protein names
        table_ = [['Accession ID', 'Homologous sequence', 'Protein description']]
        for (protein_name_dict, ordered_keys) in dict_of_heatmap_protein_names[n]:
            for accession in ordered_keys:
                if accession.upper() not in [protein1_name.upper(), protein2_name.upper()]:
                    [peptide_seq, descriptions] = protein_name_dict[accession]
                    table_.append([accession, peptide_seq, ', '.join(descriptions)])
        Story.append(Table(
                        make_table_header_bold(make_table_wrap(table_)),
                        hAlign='LEFT',
                        colWidths=[1.2*inch, 1.8*inch, 4.3*inch]
                        ))

        doc.addPageTemplates(
                [PageTemplate(
                    id='PeptidePage',
                    frames=build_frames(doc, alignment_len, image_height)
                    ),
                 PageTemplate(
                    id='TablePage',
                    frames=Frame(doc.leftMargin, doc.bottomMargin, doc.width, doc.height)
                    )
                 ]
            )

    doc.build(Story)


def build_frames(doc, alignment_len, image_height):
    alignframeheight = 12*alignment_len + 12
    title_frame = Frame(doc.leftMargin, doc.height, doc.width, 0.5*inch)
    align_frame_1 = Frame(
                        doc.leftMargin,
                        doc.height - alignframeheight,
                        doc.width/2,
                        alignframeheight
                        )
    align_frame_2 = Frame(
                        doc.leftMargin + doc.width/2,
                        doc.height - alignframeheight,
                        doc.width/2,
                        alignframeheight
                        )
    heatmap_frame = Frame(
                        doc.leftMargin,
                        doc.bottomMargin,
                        doc.width,
                        6*inch
                        )
    return [title_frame, align_frame_1, align_frame_2, heatmap_frame]
