"""PDF report generation for Brimer-PLAST.

Produces the primary design report with:
1. Genome-view diagram (Overview + Zoom panes)
2. Primer pair table
3. Technical record (full sequences and run parameters)

Uses ReportLab.
"""

from __future__ import annotations

import os
import sys
from datetime import datetime

from reportlab.lib import colors
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.platypus import (
    HRFlowable,
    PageBreak,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)

from brimer_plast.diagram import (
    TRANSCRIPT_CAP,
    _GeneDiagram,
    _SequenceRow,
)
from brimer_plast.genome import exons_in_template_order, reverse_complement

# ── Layout constants ────────────────────────────────────────────────────────
MARGIN = 0.75 * inch
PAGE_WIDTH, PAGE_HEIGHT = landscape(A4)
CONTENT_WIDTH = PAGE_WIDTH - 2 * MARGIN

# ReportLab's Frame adds 6pt internal padding on each side beyond the
# user-set margins.  This is the actual space available for flowables.
_FRAME_PAD = 12  # 6 left + 6 right, likewise top + bottom
FRAME_W = PAGE_WIDTH - 2 * MARGIN - _FRAME_PAD
FRAME_H = PAGE_HEIGHT - 2 * MARGIN - _FRAME_PAD

# Pagination: upper bound on primer pairs shown per genome-view page.
# The actual count per page is computed dynamically from available frame
# height so the diagram never overflows the frame.
DIAGRAM_PAIR_CAP = 1000


# ── Report Builder ──────────────────────────────────────────────────────────

def build_pdf_report(output_path, chains, locus, filtered_pairs, target_gene, target_transcript, genome_path, annotations_path, genome_md5="", annotations_md5="", version_str="0.1.0", cli_args=None):
    doc = SimpleDocTemplate(output_path, pagesize=landscape(A4), leftMargin=MARGIN, rightMargin=MARGIN, topMargin=MARGIN, bottomMargin=MARGIN)
    styles = getSampleStyleSheet()

    story = [
        Paragraph("Brimer-PLAST Primer Design Report", styles["Title"]),
        HRFlowable(width="100%", color=colors.grey),
        Spacer(1, 8),
        Paragraph(f"Date: {datetime.now():%Y-%m-%d %H:%M:%S}", styles["Normal"]),
        Paragraph(f"Target: <i>{target_gene or target_transcript}</i>", styles["Normal"]),
        Paragraph(f"Genome database: {os.path.basename(genome_path)}", styles["Normal"]),
        Spacer(1, 12),
        Paragraph("1. Genome-view Diagram", styles["Heading2"]),
        Paragraph(
            "<b>Color legend</b> - Blue: primer3-designed positions on the "
            "conserved exon chain (Panel A). Red: tnBLAST-verified genomic "
            "amplicon from the specificity screen, split by exon boundaries "
            "(Panel B).",
            styles["Normal"],
        ),
        Spacer(1, 4),
    ]

    if locus:
        # Gene diagram draws primer pairs across multiple pages.  The
        # number of pairs per page is computed dynamically from the
        # available frame height so the diagram never overflows.
        #
        # Drawing geometry (from draw_gene_diagram trace):
        #   base overhead = 82 + n_trans*11 pt (with "+N others" line)
        #                  = 74 + n_trans*11 pt (without "+N others")
        #   each pair     = 20 pt
        n_trans = min(len(locus.transcripts), TRANSCRIPT_CAP)
        has_others = len(locus.transcripts) > TRANSCRIPT_CAP
        diagram_base_h = (82 if has_others else 74) + n_trans * 11

        # Frame geometry: actual space inside ReportLab's frame
        # (which adds 6pt internal padding, hence FRAME_H vs PAGE_HEIGHT).
        # Overhead measured empirically — later pages have none.
        PAGE1_OVERHEAD = 138  # title, date, heading, legend, spacers

        def _max_pairs(overhead: float) -> int:
            return max(1, int((FRAME_H - overhead - diagram_base_h) // 20))

        first_page_pairs = _max_pairs(PAGE1_OVERHEAD)
        later_page_pairs = _max_pairs(0)  # full frame for page 2+

        total = len(filtered_pairs)
        page_start = 0
        page_num = 0
        while page_start < total:
            cap = first_page_pairs if page_num == 0 else later_page_pairs
            page_end = min(page_start + cap, total)
            page_pairs = filtered_pairs[page_start:page_end]
            n_pairs_this_page = len(page_pairs)
            h = diagram_base_h + n_pairs_this_page * 20

            if page_num > 0:
                story.append(PageBreak())

            story.append(
                _GeneDiagram(locus, page_pairs, chains, FRAME_W, h, target_transcript=target_transcript)
            )
            page_start = page_end
            page_num += 1

    story.append(PageBreak())
    story.append(Paragraph("2. Filtered Primer Pairs", styles["Heading2"]))

    if filtered_pairs:
        table_cell_style = ParagraphStyle(
            "TableCell",
            parent=styles["Normal"],
            fontSize=7,
            leading=8,
            fontName="Courier",
        )

        data: list[list[str | Paragraph]] = [["Pair", "Forward", "Tm", "GC%", "Reverse", "Tm", "GC%", "Size"]]
        for p in filtered_pairs:
            rc_rev = reverse_complement(p.reverse_seq or "")
            rev_p = Paragraph(f"&nbsp;{p.reverse_seq}<br/>({rc_rev})", table_cell_style)
            data.append([
                str(p.pair_number),
                Paragraph(p.forward_seq or "", table_cell_style),
                f"{p.forward_tm:.1f}",
                f"{p.forward_gc:.0f}",
                rev_p,
                f"{p.reverse_tm:.1f}",
                f"{p.reverse_gc:.0f}",
                str(p.product_size)
            ])

        cw = [CONTENT_WIDTH * x for x in [0.06, 0.22, 0.08, 0.08, 0.22, 0.08, 0.08, 0.18]]
        t = Table(data, colWidths=cw)
        t.setStyle(TableStyle([
            ('BACKGROUND', (0,0), (-1,0), colors.lightgrey),
            ('GRID', (0,0), (-1,-1), 0.5, colors.grey),
            ('FONTSIZE', (0,0), (-1,-1), 7),
            ('VALIGN', (0,0), (-1,-1), 'MIDDLE'),
            ('LEFTPADDING', (0,0), (-1,-1), 2),
            ('RIGHTPADDING', (0,0), (-1,-1), 2),
        ]))
        story.append(t)
    else:
         story.append(Paragraph("No primer pairs passed specificity filtering.", styles["Normal"]))

    # ── Page 3: Run Information ──────────────────────────────────────────────
    story.append(PageBreak())
    story.append(Paragraph("3. Run Information", styles["Heading2"]))
    story.append(Spacer(1, 10))

    file_data = [
        [Paragraph("<b>File type</b>", styles["Normal"]), Paragraph("<b>Location and Checksum</b>", styles["Normal"])],
        ["Genome", Paragraph(f"{os.path.basename(genome_path)}<br/>(md5: {genome_md5})", styles["Normal"])],
        ["Annotations", Paragraph(f"{os.path.basename(annotations_path)}<br/>(md5: {annotations_md5})", styles["Normal"])]
    ]
    ft = Table(file_data, colWidths=[CONTENT_WIDTH*0.25, CONTENT_WIDTH*0.75])
    ft.setStyle(TableStyle([
        ('BACKGROUND', (0,0), (-1,0), colors.lightgrey),
        ('GRID', (0,0), (-1,-1), 0.5, colors.grey),
        ('VALIGN', (0,0), (-1,-1), 'TOP'),
    ]))
    story.append(ft)
    story.append(Spacer(1, 15))

    if cli_args:
        story.append(Paragraph("CLI Parameters:", styles["Heading3"]))
        param_data = [["Parameter", "Value"]]
        for k, v in sorted(cli_args.items()):
            param_data.append([f"--{k}", str(v)])
        pt = Table(param_data, colWidths=[CONTENT_WIDTH*0.4, CONTENT_WIDTH*0.6])
        pt.setStyle(TableStyle([
            ('BACKGROUND', (0,0), (-1,0), colors.lightgrey),
            ('GRID', (0,0), (-1,-1), 0.5, colors.grey),
            ('FONTSIZE', (0,0), (-1,-1), 8),
        ]))
        story.append(pt)

    story.append(Spacer(1, 15))
    story.append(Paragraph("Software Versions:", styles["Heading3"]))
    ver_lines = [
        f"Brimer-PLAST: {version_str}",
        f"Python: {sys.version.split()[0]}",
    ]
    for vl in ver_lines:
        story.append(Paragraph(vl, styles["Normal"]))

    # ── Page 4+: Conserved Exon Chains ───────────────────────────────────────
    if chains:
        story.append(PageBreak())
        story.append(Paragraph("4. Conserved Exon Chains", styles["Heading2"]))
        story.append(Spacer(1, 10))

        for chain in chains:
            story.append(Paragraph(f"Chain: {chain.id}", styles["Heading3"]))

            ordered_exons = exons_in_template_order(chain.exons)

            cumulative_pos = 0
            for i, ex in enumerate(ordered_exons, start=1):
                exon_len = ex.end - ex.start + 1

                if ex.strand == "-":
                    coord_str = f"{ex.seqid}:{ex.end}-{ex.start} (-)"
                else:
                    coord_str = f"{ex.seqid}:{ex.start}-{ex.end} (+)"

                story.append(Paragraph(f"Exon {i} ({coord_str})", styles["Heading4"]))

                exon_seq = chain.template[cumulative_pos : cumulative_pos + exon_len]

                for r in range(0, len(exon_seq), 100):
                    row_seq = exon_seq[r : r + 100]
                    chain_pairs = [p for p in filtered_pairs if p.chain_id == chain.id]
                    story.append(_SequenceRow(row_seq, cumulative_pos + r + 1, chain_pairs))

                story.append(Spacer(1, 10))

                cumulative_pos += exon_len

    doc.build(story)