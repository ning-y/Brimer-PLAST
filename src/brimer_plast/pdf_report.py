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
from typing import Any

from reportlab.lib import colors
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.platypus import (
    Flowable,
    HRFlowable,
    PageBreak,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
    Preformatted,
)

from brimer_plast.genome import _exons_in_template_order, _reverse_complement
from brimer_plast.models import ConservedExonChain, ExonInfo, GeneLocus, GenomicFragment, PrimerPair

# ── Colour palette ──────────────────────────────────────────────────────────
PANEL_A_COLOR = colors.HexColor("#2166ac")  # blue — primer3
PANEL_B_COLOR = colors.HexColor("#b2182b")  # reddish — tnBLAST
PANEL_A_FILL = colors.HexColor("#d1e5f0")
PANEL_B_FILL = colors.HexColor("#fddbc7")
GENE_FILL = colors.Color(0.8, 0.8, 0.8)
VIEWBOX_FILL = colors.Color(1, 0, 0, alpha=0.1)

# ── Layout constants ────────────────────────────────────────────────────────
MARGIN = 0.75 * inch
PAGE_WIDTH, PAGE_HEIGHT = landscape(A4)
CONTENT_WIDTH = PAGE_WIDTH - 2 * MARGIN

LEFT_MARGIN = 60  
RIGHT_MARGIN = 60 # increased for labels on right
DRAW_W = CONTENT_WIDTH - LEFT_MARGIN - RIGHT_MARGIN

EXON_HEIGHT = 8   
TRACK_SPACING = 3 
PRIMER_ROW_H = 14 
PAIR_GAP = 6
LABEL_FONT_SIZE = 6
TICK_FONT_SIZE = 6
TRANSCRIPT_CAP = 10

OVERVIEW_H = 15 

# ── Diagram drawing helpers ─────────────────────────────────────────────────

def _draw_intron(canvas, x1, x2, y, h, strand):
    canvas.setStrokeColor(colors.grey)
    canvas.setLineWidth(0.5)
    canvas.line(x1, y + h/2, x2, y + h/2)
    span = abs(x2-x1)
    if span > 10:
        n = max(1, int(span/20))
        dist = span / (n+1)
        for i in range(1, n+1):
            cx = min(x1,x2) + i*dist
            cy = y + h/2
            s = 2
            if strand == "+":
                canvas.line(cx-s, cy-s, cx, cy)
                canvas.line(cx-s, cy+s, cx, cy)
            else:
                canvas.line(cx+s, cy-s, cx, cy)
                canvas.line(cx+s, cy+s, cx, cy)

def _draw_fragments(canvas, fragments, y, h, view_min, view_span, origin_x, draw_w, color, fill, label, alignment):
    """Draw one or more fragments for a primer, connecting with dashed lines if segmented."""
    def to_x(g):
        f = (g - view_min) / view_span
        return origin_x + f * draw_w

    drawn_xs = []
    
    for i, frag in enumerate(fragments):
        gx1, gx2 = frag.start, frag.end
        x1, x2 = to_x(gx1), to_x(gx2)
        rx1, rx2 = max(origin_x, min(x1, x2)), min(origin_x + draw_w, max(x1, x2))
        
        if rx1 < rx2:
            canvas.setFillColor(fill)
            canvas.setStrokeColor(color)
            canvas.setLineWidth(0.5)
            canvas.rect(rx1, y, rx2-rx1, h, fill=1, stroke=1)
            drawn_xs.append((rx1, rx2))
            
            # If segmented, connect to previous with dashed bridge
            if i > 0 and drawn_xs:
                prev_rx1, prev_rx2 = drawn_xs[-2]
                canvas.setLineWidth(0.5)
                canvas.setDash(1, 1)
                canvas.line(prev_rx2, y + h/2, rx1, y + h/2)
                canvas.setDash()

    # Draw Label (once per multi-segment set)
    if drawn_xs:
        min_x = min(x for x, _ in drawn_xs)
        max_x = max(x for _, x in drawn_xs)
        canvas.setFont("Helvetica", 5)
        canvas.setFillColor(color)
        if alignment == "left":
            canvas.drawRightString(min_x - 3, y + h/2 - 2, label)
        else:
            canvas.drawString(max_x + 3, y + h/2 - 2, label)


def draw_gene_diagram(canvas, x, y, width, locus, filtered_pairs, chains):
    # Overall genomic range
    g_min, g_max = locus.min_start, locus.max_end
    g_span = max(1, g_max - g_min)
    
    origin_x = x + LEFT_MARGIN
    draw_w = width - LEFT_MARGIN - RIGHT_MARGIN
    
    def to_x_over(g): return origin_x + ((g-g_min)/g_span)*draw_w

    # Determine Zoom bounds from all fragment lists
    p_coords = []
    for p in filtered_pairs:
        for frag in list(p.primer3_forward_fragments) + list(p.primer3_reverse_fragments) \
                     + list(p.tnblast_forward_fragments) + list(p.tnblast_reverse_fragments):
            p_coords.extend([frag.start, frag.end])
    
    if p_coords:
        pmin, pmax = min(p_coords), max(p_coords)
        pad = max(500, int((pmax-pmin)*0.4))
        v_min, v_max = max(g_min, pmin-pad), min(g_max, pmax+pad)
    else:
        v_min, v_max = g_min, g_max
    
    v_span = max(1, v_max - v_min)
    def to_x_zoom(g): return origin_x + ((g-v_min)/v_span)*draw_w

    # ── 1. Header ────────────────────────────────────────────────────────────
    canvas.setFont("Helvetica-Bold", 8)
    canvas.setFillColor(colors.black)
    canvas.drawString(x, y + 10, f"Gene: {locus.gene_name} | Region: {locus.seqid}:{g_min:,}-{g_max:,}")
    
    # ── 2. Overview Pane ─────────────────────────────────────────────────────
    over_y = y - OVERVIEW_H
    
    # Consensus Intron Line
    _draw_intron(canvas, origin_x, origin_x + draw_w, over_y + 5, 8, locus.strand)
    
    # Collapsed Gene Model
    all_exons = []
    for exons in locus.transcripts.values():
        all_exons.extend(exons)
    # Deduplicate exons by coords for cleaner union view
    unique_exons = {(e.start, e.end) for e in all_exons}
    for start, end in unique_exons:
        lx, rx = to_x_over(start), to_x_over(end)
        if rx > lx:
            canvas.setFillColor(GENE_FILL)
            canvas.setStrokeColor(colors.black)
            canvas.setLineWidth(0.1)
            canvas.rect(lx, over_y + 5, rx-lx, 8, fill=1, stroke=1)
            
    # Viewbox
    vx1, vx2 = to_x_over(v_min), to_x_over(v_max)
    canvas.setFillColor(VIEWBOX_FILL)
    canvas.setStrokeColor(colors.red)
    canvas.setLineWidth(0.5)
    canvas.rect(vx1, over_y + 3, max(2, vx2-vx1), 12, fill=1, stroke=1)
    
    # Axis labels for overview
    canvas.setFont("Helvetica", 5)
    canvas.setFillColor(colors.grey)
    canvas.drawString(origin_x, over_y + OVERVIEW_H + 2, f"{g_min:,}")
    canvas.drawRightString(origin_x+draw_w, over_y + OVERVIEW_H + 2, f"{g_max:,}")

    # Visual connectors
    zoom_top_y = over_y - 25
    canvas.setStrokeColor(colors.red)
    canvas.setLineWidth(0.15)
    canvas.setDash(1, 2)
    canvas.line(vx1, over_y + 3, origin_x, zoom_top_y)
    canvas.line(vx2, over_y + 3, origin_x + draw_w, zoom_top_y)
    canvas.setDash()

    # ── 3. Zoom Pane Axis ────────────────────────────────────────────────────
    canvas.setStrokeColor(colors.grey)
    canvas.setLineWidth(0.4)
    canvas.line(origin_x, zoom_top_y, origin_x+draw_w, zoom_top_y)
    canvas.setFont("Helvetica", TICK_FONT_SIZE)
    canvas.drawString(origin_x, zoom_top_y+2, f"{v_min:,}")
    canvas.drawRightString(origin_x+draw_w, zoom_top_y+2, f"{v_max:,}")

    # ── 4. Transcripts (Zoom) ────────────────────────────────────────────────
    curr_y = zoom_top_y - EXON_HEIGHT - 5
    sorted_tids = sorted(locus.transcripts.keys())
    for tid in sorted_tids[:TRANSCRIPT_CAP]:
        exons = locus.transcripts[tid]
        canvas.setFont("Helvetica", LABEL_FONT_SIZE)
        canvas.setFillColor(colors.black)
        canvas.drawString(x, curr_y + 2, tid[:15])
        
        for i in range(len(exons)-1):
            x1, x2 = to_x_zoom(exons[i].end), to_x_zoom(exons[i+1].start)
            if x1 < origin_x+draw_w and x2 > origin_x:
                _draw_intron(canvas, max(origin_x, x1), min(origin_x+draw_w, x2), curr_y, EXON_HEIGHT, locus.strand)
        
        for ex in exons:
            lx, rx = to_x_zoom(ex.start), to_x_zoom(ex.end)
            if lx < origin_x+draw_w and rx > origin_x:
                canvas.setFillColor(GENE_FILL)
                canvas.setStrokeColor(colors.black)
                canvas.setLineWidth(0.2)
                canvas.rect(max(origin_x, lx), curr_y, min(origin_x+draw_w, rx)-max(origin_x, lx), EXON_HEIGHT, fill=1, stroke=1)
        curr_y -= (EXON_HEIGHT + TRACK_SPACING)

    if len(sorted_tids) > TRANSCRIPT_CAP:
        canvas.setFont("Helvetica-Oblique", 5)
        canvas.drawString(x, curr_y+2, f"+ {len(sorted_tids)-TRANSCRIPT_CAP} others")
        curr_y -= 8

    # ── 5. Primer Plots ──────────────────────────────────────────────────────
    curr_y -= 10
    
    for pair in filtered_pairs:
        pnum = pair.pair_number or "?"
        canvas.setFont("Helvetica-Bold", LABEL_FONT_SIZE)
        canvas.setFillColor(colors.black)
        canvas.drawString(x, curr_y + 4, f"Pair {pnum}")

        # Panel A (Blue) — pre-computed primer3 fragments
        if pair.primer3_forward_fragments or pair.primer3_reverse_fragments:
            all_pts = []
            if pair.primer3_forward_fragments:
                all_pts.append((pair.primer3_forward_fragments, "F"))
            if pair.primer3_reverse_fragments:
                all_pts.append((pair.primer3_reverse_fragments, "R"))
            all_pts.sort(key=lambda item: item[0][0].start)

            for index, (frags, label) in enumerate(all_pts):
                alignment = "left" if index == 0 else "right"
                _draw_fragments(canvas, frags, curr_y+4, 4, v_min, v_span, origin_x, draw_w,
                                PANEL_A_COLOR, PANEL_A_FILL, f"{pnum}{label}", alignment)

        # Panel B (Red) — pre-computed tnBLAST fragments
        if pair.tnblast_forward_fragments or pair.tnblast_reverse_fragments:
            b_items = []
            if pair.tnblast_forward_fragments:
                b_items.append((pair.tnblast_forward_fragments, "F'"))
            if pair.tnblast_reverse_fragments:
                b_items.append((pair.tnblast_reverse_fragments, "R'"))
            b_items.sort(key=lambda item: item[0][0].start)

            for index, (frags, label) in enumerate(b_items):
                alignment = "left" if index == 0 else "right"
                _draw_fragments(canvas, frags, curr_y-1, 4, v_min, v_span, origin_x, draw_w,
                                PANEL_B_COLOR, PANEL_B_FILL, f"{pnum}{label}", alignment)

        curr_y -= (PRIMER_ROW_H + PAIR_GAP)

    return curr_y - 5


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
            "<b>Color legend</b> — Blue: primer3-designed positions on the "
            "conserved exon chain (Panel A). Red: tnBLAST-verified genomic "
            "amplicon from the specificity screen, split by exon boundaries "
            "(Panel B).",
            styles["Normal"],
        ),
        Spacer(1, 4),
    ]

    if locus:
        n_pairs = len(filtered_pairs)
        n_trans = min(len(locus.transcripts), TRANSCRIPT_CAP)
        h = 130 + n_trans*11 + n_pairs*22
        story.append(_GeneDiagram(locus, filtered_pairs, chains, CONTENT_WIDTH, h))
    
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
        
        data = [["Pair", "Forward", "Tm", "GC%", "Reverse", "Tm", "GC%", "Size"]]
        for p in filtered_pairs:
            rc_rev = _reverse_complement(p.reverse_seq or "")
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
            
            ordered_exons = _exons_in_template_order(chain.exons)
            
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


class _SequenceRow(Flowable):
    """Draws a single 100-base row of DNA sequence with indices and highlights."""
    def __init__(self, sequence: str, start_index: int, pairs: list[PrimerPair]):
        super().__init__()
        self.sequence = sequence
        self.start_index = start_index
        self.pairs = pairs
        self.line_h = 10
        self.height = 3.5 * self.line_h
        self.width = 600
        
    def wrap(self, w, h): return (self.width, self.height)
    
    def draw(self):
        canv = self.canv
        canv.setFont("Courier", 8)
        char_w = canv.stringWidth("A", "Courier", 8)
        
        seq = self.sequence
        idx = self.start_index
        
        f_hits: dict[int, set[int]] = {}
        r_hits: dict[int, set[int]] = {}
        
        for p in self.pairs:
            pnum = p.pair_number or 0
            if p.forward_start is not None:
                for pos in range(p.forward_start + 1, p.forward_start + p.forward_len + 1):
                    f_hits.setdefault(pos, set()).add(pnum)
            if p.reverse_start is not None:
                for pos in range(p.reverse_start - p.reverse_len + 2, p.reverse_start + 2):
                    r_hits.setdefault(pos, set()).add(pnum)

        for i in range(len(seq)):
            pos = idx + i
            f_ids = f_hits.get(pos, set())
            r_ids = r_hits.get(pos, set())
            
            if not f_ids and not r_ids:
                continue
            
            tx = 40 + i * char_w + (i // 20) * char_w
            
            if f_ids:
                canv.setFillColor(colors.Color(1, 1, 0, alpha=0.5))
                canv.rect(tx, 0.5 * self.line_h - 2, char_w, self.line_h, fill=1, stroke=0)
            
            if r_ids:
                canv.setFillColor(colors.Color(0, 1, 0, alpha=0.5))
                canv.rect(tx, 0.5 * self.line_h - 2, char_w, self.line_h, fill=1, stroke=0)

        canv.setFont("Helvetica", 4)
        for i in range(len(seq)):
            pos = idx + i
            ids = sorted(f_hits.get(pos, set()) | r_hits.get(pos, set()))
            if not ids:
                continue
                
            label = ",".join(map(str, ids))
            
            next_pos = pos + 1
            next_ids = sorted(f_hits.get(next_pos, set()) | r_hits.get(next_pos, set()))
            next_label = ",".join(map(str, next_ids))
            
            if label != next_label:
                tx = 40 + i * char_w + (i // 20) * char_w
                canv.setFillColor(colors.black)
                canv.drawString(tx + char_w - 1, 1.5 * self.line_h - 1, label)

        canv.setFont("Courier", 8)
        canv.setFillColor(colors.grey)
        for b in range(0, len(seq), 20):
            if b + 19 < len(seq):
                label_val = idx + b + 19
                tx = 40 + (b + 19) * char_w + (b // 20) * char_w
                canv.drawRightString(tx + char_w, 2.5 * self.line_h, str(label_val))

        canv.setFillColor(colors.black)
        canv.drawRightString(35, 0.5 * self.line_h, str(idx))
        
        blocks = []
        for b in range(0, len(seq), 20):
            blocks.append(seq[b:b+20])
        canv.drawString(40, 0.5 * self.line_h, " ".join(blocks))
        
        end_idx = idx + len(seq) - 1
        canv.drawString(40 + 100 * char_w + 7 * char_w, 0.5 * self.line_h, str(end_idx))


class _GeneDiagram(Flowable):
    def __init__(self, locus, filtered_pairs, chains, width, height):
        super().__init__()
        self.locus, self.filtered_pairs, self.chains, self.width, self.height = locus, filtered_pairs, chains, width, height
    def wrap(self, w, h): return (self.width, self.height)
    def draw(self): draw_gene_diagram(self.canv, 0, self.height-30, self.width, self.locus, self.filtered_pairs, self.chains)
