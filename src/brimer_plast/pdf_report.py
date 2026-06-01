"""PDF report generation for Brimer-PLAST.

Produces an archival PDF report with:
1. Genome-view diagram (Overview + Zoom panes)
2. Primer pair table
3. Run information (archival record)

Uses ReportLab.
"""

from __future__ import annotations

import os
import sys
from datetime import datetime
from typing import Any

from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
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
)

from brimer_plast.genome import template_to_genomic
from brimer_plast.models import ConservedExonChain, ExonInfo, GeneLocus, PrimerPair

# ── Colour palette ──────────────────────────────────────────────────────────
PANEL_A_COLOR = colors.HexColor("#2166ac")  # blue — primer3
PANEL_B_COLOR = colors.HexColor("#b2182b")  # reddish — tnBLAST
PANEL_A_FILL = colors.HexColor("#d1e5f0")
PANEL_B_FILL = colors.HexColor("#fddbc7")
GENE_FILL = colors.Color(0.8, 0.8, 0.8)
VIEWBOX_FILL = colors.Color(1, 0, 0, alpha=0.1)

# ── Layout constants ────────────────────────────────────────────────────────
MARGIN = 0.75 * inch
PAGE_WIDTH, PAGE_HEIGHT = letter
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
        gx1, gx2 = frag[1], frag[2]
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

    # Determine Zoom bounds
    p_coords = []
    for p in filtered_pairs:
        if p.tntblast_amplicon_start:
            p_coords.extend([p.tntblast_amplicon_start, p.tntblast_amplicon_end])
    
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
    _draw_intron(canvas, origin_x, origin_x + draw_w, over_y + 5, 2, locus.strand)
    
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
    chain_map = {c.id: c for c in chains}
    
    for pair in filtered_pairs:
        pnum = pair.pair_number or "?"
        canvas.setFont("Helvetica-Bold", LABEL_FONT_SIZE)
        canvas.setFillColor(colors.black)
        canvas.drawString(x, curr_y + 4, f"Pair {pnum}")

        # Panel A (Blue)
        if pair.chain_id in chain_map:
            exons = chain_map[pair.chain_id].exons
            
            # Gather fragments for F and R
            f_frags = template_to_genomic(pair.forward_start, pair.forward_len, exons) if pair.forward_start is not None else []
            r_frags = template_to_genomic(pair.reverse_start, pair.reverse_len, exons) if pair.reverse_start is not None else []
            
            # Sort by genomic start to decide alignment
            all_pts = []
            if f_frags: all_pts.append((f_frags, "F", "blue"))
            if r_frags: all_pts.append((r_frags, "R", "blue"))
            
            all_pts.sort(key=lambda item: item[0][0][1]) # sort by start of first fragment
            
            for index, (frags, label, _) in enumerate(all_pts):
                alignment = "left" if index == 0 else "right"
                _draw_fragments(canvas, frags, curr_y+4, 4, v_min, v_span, origin_x, draw_w, PANEL_A_COLOR, PANEL_A_FILL, f"{pnum}{label}", alignment)

        # Panel B (Red)
        if pair.tntblast_amplicon_start:
            # We treatPanel B hits as single-segment points for the "Sanity Check"
            # as they are re-derived from mRNA hits which already indicate validity.
            
            # Determine ordering for labels
            b_items = [
                ([("", pair.tntblast_amplicon_start, pair.tntblast_amplicon_start + (pair.forward_len or 20) - 1, "")], "F'"),
                ([("", pair.tntblast_amplicon_end - (pair.reverse_len or 20) + 1, pair.tntblast_amplicon_end, "")], "R'")
            ]
            b_items.sort(key=lambda item: item[0][0][1])
            
            for index, (frags, label) in enumerate(b_items):
                alignment = "left" if index == 0 else "right"
                _draw_fragments(canvas, frags, curr_y-1, 4, v_min, v_span, origin_x, draw_w, PANEL_B_COLOR, PANEL_B_FILL, f"{pnum}{label}", alignment)

        curr_y -= (PRIMER_ROW_H + PAIR_GAP)

    return curr_y - 5

# ── Report Builder ──────────────────────────────────────────────────────────

def build_pdf_report(output_path, chains, locus, filtered_pairs, target_gene, target_transcript, genome_path, annotations_path, genome_md5="", annotations_md5="", version_str="0.1.0", cli_args=None):
    doc = SimpleDocTemplate(output_path, pagesize=letter, leftMargin=MARGIN, rightMargin=MARGIN, topMargin=MARGIN, bottomMargin=MARGIN)
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
    ]

    if locus:
        n_pairs = len(filtered_pairs)
        n_trans = min(len(locus.transcripts), TRANSCRIPT_CAP)
        h = 130 + n_trans*11 + n_pairs*22
        story.append(_GeneDiagram(locus, filtered_pairs, chains, CONTENT_WIDTH, h))
    
    story.append(PageBreak())
    story.append(Paragraph("2. Filtered Primer Pairs", styles["Heading2"]))
    
    if filtered_pairs:
        data = [["Pair", "Forward (5->3)", "Tm", "GC%", "Reverse (5->3)", "Tm", "GC%", "Size"]]
        for p in filtered_pairs:
            data.append([str(p.pair_number), p.forward_seq, f"{p.forward_tm:.1f}", f"{p.forward_gc:.0f}", p.reverse_seq, f"{p.reverse_tm:.1f}", f"{p.reverse_gc:.0f}", str(p.product_size)])
        
        cw = [CONTENT_WIDTH * x for x in [0.06, 0.22, 0.08, 0.08, 0.22, 0.08, 0.08, 0.18]]
        t = Table(data, colWidths=cw)
        t.setStyle(TableStyle([
            ('BACKGROUND', (0,0), (-1,0), colors.lightgrey),
            ('GRID', (0,0), (-1,-1), 0.5, colors.grey),
            ('FONTSIZE', (0,0), (-1,-1), 7),
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
    
    # Files table (2 columns: File, Name + MD5)
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
    
    # CLI Parameters
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

    # Software versions
    story.append(Spacer(1, 15))
    story.append(Paragraph("Software Versions:", styles["Heading3"]))
    ver_lines = [
        f"Brimer-PLAST: {version_str}",
        f"Python: {sys.version.split()[0]}",
    ]
    for vl in ver_lines:
        story.append(Paragraph(vl, styles["Normal"]))

    doc.build(story)

class _GeneDiagram(Flowable):
    def __init__(self, locus, filtered_pairs, chains, width, height):
        super().__init__()
        self.locus, self.filtered_pairs, self.chains, self.width, self.height = locus, filtered_pairs, chains, width, height
    def wrap(self, w, h): return (self.width, self.height)
    def draw(self): draw_gene_diagram(self.canv, 0, self.height-30, self.width, self.locus, self.filtered_pairs, self.chains)
