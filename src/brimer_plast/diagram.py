"""Genome-view diagram and sequence row flowables for PDF reports.

Each genome view is drawn per conserved exon chain (rather than one view
for the whole target).  Contributing transcripts — those whose exons
contain every exon in the chain — are drawn with normal fill and solid
stroke.  Non-contributing transcripts are drawn with no fill, dashed exon
borders, and dashed intron arrows.
"""

from __future__ import annotations

from reportlab.lib import colors
from reportlab.platypus import Flowable

from brimer_plast.models import ConservedExonChain, ExonInfo, GeneLocus, PrimerPair

# ── Colour palette ──────────────────────────────────────────────────────────
PANEL_A_COLOR = colors.HexColor("#2166ac")  # blue - primer3
PANEL_B_COLOR = colors.HexColor("#b2182b")  # reddish - tnBLAST
PANEL_A_FILL = colors.HexColor("#d1e5f0")
PANEL_B_FILL = colors.HexColor("#fddbc7")
GENE_FILL = colors.Color(0.8, 0.8, 0.8)
VIEWBOX_FILL = colors.Color(1, 0, 0, alpha=0.1)

# ── Layout constants ────────────────────────────────────────────────────────
LEFT_MARGIN = 76
RIGHT_MARGIN = 60  # increased for labels on right
EXON_HEIGHT = 8
TRACK_SPACING = 3
PRIMER_ROW_H = 14
PAIR_GAP = 6
# Upper bound on primer pairs shown per genome-view page.  The actual
# count per page is computed dynamically from available frame height so
# the diagram never overflows the frame.
DIAGRAM_PAIR_CAP = 1000
LABEL_FONT_SIZE = 6
TICK_FONT_SIZE = 6
TRANSCRIPT_CAP = 10

OVERVIEW_H = 15


# ── Helpers ─────────────────────────────────────────────────────────────────


def _transcript_contains_chain(
    transcript_exons: list[ExonInfo],
    chain_exons: list[ExonInfo],
) -> bool:
    """True if every exon in *chain_exons* (by start/end coords) appears in
    *transcript_exons*."""
    keys = {(e.start, e.end) for e in transcript_exons}
    return all((e.start, e.end) in keys for e in chain_exons)


def compute_contributing_tids(
    locus: GeneLocus,
    chain: ConservedExonChain,
) -> set[str]:
    """Return the set of transcript IDs whose exons contain every exon
    in *chain*."""
    return {
        tid
        for tid, exons in locus.transcripts.items()
        if _transcript_contains_chain(exons, chain.exons)
    }


def compute_zoom_bounds(
    locus: GeneLocus,
    chain_pairs: list[PrimerPair],
    contributing_tids: set[str],
) -> tuple[int, int]:
    """Return ``(v_min, v_max)`` for the zoom pane.

    If *chain_pairs* is non-empty, zoom to the primer fragment
    coordinates.  Otherwise fall back to the exon range of the
    contributing transcripts (with padding).
    """
    p_coords: list[int] = []
    for p in chain_pairs:
        for frag in (
            list(p.primer3_forward_fragments)
            + list(p.primer3_reverse_fragments)
            + list(p.tnblast_forward_fragments)
            + list(p.tnblast_reverse_fragments)
        ):
            p_coords.extend([frag.start, frag.end])

    if p_coords:
        pmin, pmax = min(p_coords), max(p_coords)
        pad = max(500, int((pmax - pmin) * 0.4))
        return max(locus.min_start, pmin - pad), min(locus.max_end, pmax + pad)

    # Fallback: zoom to contributing transcripts' exon range
    all_starts: list[int] = []
    all_ends: list[int] = []
    for tid in contributing_tids:
        for ex in locus.transcripts[tid]:
            all_starts.append(ex.start)
            all_ends.append(ex.end)
    if not all_starts:
        return locus.min_start, locus.max_end
    cmin, cmax = min(all_starts), max(all_ends)
    pad = max(500, int((cmax - cmin) * 0.4))
    return max(locus.min_start, cmin - pad), min(locus.max_end, cmax + pad)


# ── Diagram drawing helpers ─────────────────────────────────────────────────


def _draw_intron(canvas, x1, x2, y, h, strand, dashed=False):
    canvas.setStrokeColor(colors.grey)
    canvas.setLineWidth(0.5)
    if dashed:
        canvas.setDash(2, 2)
    canvas.line(x1, y + h / 2, x2, y + h / 2)
    span = abs(x2 - x1)
    if span > 10:
        n = max(1, int(span / 20))
        dist = span / (n + 1)
        for i in range(1, n + 1):
            cx = min(x1, x2) + i * dist
            cy = y + h / 2
            s = 2
            if strand == "+":
                canvas.line(cx - s, cy - s, cx, cy)
                canvas.line(cx - s, cy + s, cx, cy)
            else:
                canvas.line(cx + s, cy - s, cx, cy)
                canvas.line(cx + s, cy + s, cx, cy)
    canvas.setDash()


def _draw_fragments(
    canvas,
    fragments,
    y,
    h,
    view_min,
    view_span,
    origin_x,
    draw_w,
    color,
    fill,
    label,
    alignment,
):
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
            canvas.rect(rx1, y, rx2 - rx1, h, fill=1, stroke=1)
            drawn_xs.append((rx1, rx2))

            # If segmented, connect to previous with dashed bridge
            if i > 0 and len(drawn_xs) >= 2:
                prev_rx1, prev_rx2 = drawn_xs[-2]
                canvas.setLineWidth(0.5)
                canvas.setDash(1, 1)
                canvas.line(prev_rx2, y + h / 2, rx1, y + h / 2)
                canvas.setDash()

    # Draw Label (once per multi-segment set)
    if drawn_xs:
        min_x = min(x for x, _ in drawn_xs)
        max_x = max(x for _, x in drawn_xs)
        canvas.setFont("Helvetica", 5)
        canvas.setFillColor(color)
        if alignment == "left":
            canvas.drawRightString(min_x - 3, y + h / 2 - 2, label)
        else:
            canvas.drawString(max_x + 3, y + h / 2 - 2, label)


# ── Main drawing function ──────────────────────────────────────────────────


def draw_gene_diagram(
    canvas,
    x,
    y,
    width,
    locus: GeneLocus,
    chain: ConservedExonChain,
    chain_pairs: list[PrimerPair],
    contributing_tids: set[str],
    v_min: int,
    v_max: int,
):
    """Draw one genome view for a single conserved exon chain.

    Contributing transcripts (those whose exons contain the chain's
    exons) are drawn with normal fill and solid stroke.  Non-contributing
    transcripts are drawn with no fill, dashed exon borders, and dashed
    intron arrows.

    When *chain_pairs* is empty, a "no pairs" message is shown instead of
    primer tracks.
    """
    # Overall genomic range
    g_min, g_max = locus.min_start, locus.max_end
    g_span = max(1, g_max - g_min)

    origin_x = x + LEFT_MARGIN
    draw_w = width - LEFT_MARGIN - RIGHT_MARGIN

    def to_x_over(g):
        return origin_x + ((g - g_min) / g_span) * draw_w

    v_span = max(1, v_max - v_min)

    def to_x_zoom(g):
        return origin_x + ((g - v_min) / v_span) * draw_w

    # ── 1. Header ────────────────────────────────────────────────────────────
    canvas.setFont("Helvetica-Bold", 8)
    canvas.setFillColor(colors.black)
    canvas.drawString(
        x,
        y + 10,
        f"Gene: {locus.gene_name}  |  Chain: {chain.id}"
        f"  |  Region: {locus.seqid}:{g_min:,}-{g_max:,}",
    )

    # ── 2. Overview Pane (contributing transcripts' exons only) ──────────────
    over_y = y - OVERVIEW_H

    # Consensus Intron Line
    _draw_intron(canvas, origin_x, origin_x + draw_w, over_y + 5, 8, locus.strand)

    # Collapsed Gene Model — only exons from contributing transcripts
    contrib_exons: list[ExonInfo] = []
    for tid in sorted(contributing_tids):
        contrib_exons.extend(locus.transcripts[tid])
    unique_contrib = {(e.start, e.end) for e in contrib_exons}
    for start, end in unique_contrib:
        lx, rx = to_x_over(start), to_x_over(end)
        if rx > lx:
            canvas.setFillColor(GENE_FILL)
            canvas.setStrokeColor(colors.black)
            canvas.setLineWidth(0.1)
            canvas.rect(lx, over_y + 5, rx - lx, 8, fill=1, stroke=1)

    # Viewbox
    vx1, vx2 = to_x_over(v_min), to_x_over(v_max)
    canvas.setFillColor(VIEWBOX_FILL)
    canvas.setStrokeColor(colors.red)
    canvas.setLineWidth(0.5)
    canvas.rect(vx1, over_y + 3, max(2, vx2 - vx1), 12, fill=1, stroke=1)

    # Axis labels for overview
    canvas.setFont("Helvetica", 5)
    canvas.setFillColor(colors.grey)
    canvas.drawString(origin_x, over_y + OVERVIEW_H + 2, f"{g_min:,}")
    canvas.drawRightString(origin_x + draw_w, over_y + OVERVIEW_H + 2, f"{g_max:,}")

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
    canvas.line(origin_x, zoom_top_y, origin_x + draw_w, zoom_top_y)
    canvas.setFont("Helvetica", TICK_FONT_SIZE)
    canvas.drawString(origin_x, zoom_top_y + 2, f"{v_min:,}")
    canvas.drawRightString(origin_x + draw_w, zoom_top_y + 2, f"{v_max:,}")

    # ── 4. Transcripts (Zoom) ────────────────────────────────────────────────
    curr_y = zoom_top_y - EXON_HEIGHT - 5
    sorted_tids = sorted(locus.transcripts.keys())
    for tid in sorted_tids[:TRANSCRIPT_CAP]:
        exons = locus.transcripts[tid]
        is_contributing = tid in contributing_tids

        # Label
        canvas.setFont("Helvetica", LABEL_FONT_SIZE)
        canvas.setFillColor(colors.black if is_contributing else colors.Color(0, 0, 0, alpha=0.5))
        max_label_chars = max(10, int((LEFT_MARGIN - 4) / 3.8))
        canvas.drawString(x, curr_y + 2, tid[:max_label_chars])

        # Introns (dashed for non-contributing)
        for i in range(len(exons) - 1):
            x1, x2 = to_x_zoom(exons[i].end), to_x_zoom(exons[i + 1].start)
            if x1 < origin_x + draw_w and x2 > origin_x:
                _draw_intron(
                    canvas,
                    max(origin_x, x1),
                    min(origin_x + draw_w, x2),
                    curr_y,
                    EXON_HEIGHT,
                    locus.strand,
                    dashed=not is_contributing,
                )

        # Exons
        for ex in exons:
            lx, rx = to_x_zoom(ex.start), to_x_zoom(ex.end)
            if lx < origin_x + draw_w and rx > origin_x:
                if is_contributing:
                    canvas.setFillColor(GENE_FILL)
                    canvas.setStrokeColor(colors.black)
                    canvas.setLineWidth(0.2)
                    canvas.rect(
                        max(origin_x, lx),
                        curr_y,
                        min(origin_x + draw_w, rx) - max(origin_x, lx),
                        EXON_HEIGHT,
                        fill=1,
                        stroke=1,
                    )
                else:
                    # Non-contributing: no fill, dashed border
                    canvas.setStrokeColor(colors.Color(0.5, 0.5, 0.5))
                    canvas.setLineWidth(0.2)
                    canvas.setDash(2, 2)
                    canvas.rect(
                        max(origin_x, lx),
                        curr_y,
                        min(origin_x + draw_w, rx) - max(origin_x, lx),
                        EXON_HEIGHT,
                        fill=0,
                        stroke=1,
                    )
                    canvas.setDash()
        curr_y -= EXON_HEIGHT + TRACK_SPACING

    if len(sorted_tids) > TRANSCRIPT_CAP:
        canvas.setFont("Helvetica-Oblique", 5)
        canvas.drawString(x, curr_y + 2, f"+ {len(sorted_tids)-TRANSCRIPT_CAP} others")
        curr_y -= 8

    # ── 5. Primer Plots ──────────────────────────────────────────────────────
    curr_y -= 10

    if not chain_pairs:
        # No pairs designed for this chain
        canvas.setFont("Helvetica-Oblique", 7)
        canvas.setFillColor(colors.Color(0.5, 0.5, 0.5))
        canvas.drawString(
            x,
            curr_y + 4,
            "No specificity-filtered primer pairs designed for this chain.",
        )
        curr_y -= PRIMER_ROW_H + PAIR_GAP
    else:
        for pair in chain_pairs:
            pname = pair.pair_name or f"Pair {pair.pair_number or '?'}"
            canvas.setFont("Helvetica-Bold", LABEL_FONT_SIZE)
            canvas.setFillColor(colors.black)
            canvas.drawString(x, curr_y + 4, pname)

            # Panel A (Blue) - pre-computed primer3 fragments
            if pair.primer3_forward_fragments or pair.primer3_reverse_fragments:
                all_pts = []
                if pair.primer3_forward_fragments:
                    all_pts.append((pair.primer3_forward_fragments, "F"))
                if pair.primer3_reverse_fragments:
                    all_pts.append((pair.primer3_reverse_fragments, "R"))
                all_pts.sort(key=lambda item: item[0][0].start)

                for index, (frags, label) in enumerate(all_pts):
                    alignment = "left" if index == 0 else "right"
                    _draw_fragments(
                        canvas,
                        frags,
                        curr_y + 4,
                        4,
                        v_min,
                        v_span,
                        origin_x,
                        draw_w,
                        PANEL_A_COLOR,
                        PANEL_A_FILL,
                        label,
                        alignment,
                    )

            # Panel B (Red) - pre-computed tnBLAST fragments
            if pair.tnblast_forward_fragments or pair.tnblast_reverse_fragments:
                b_items = []
                if pair.tnblast_forward_fragments:
                    b_items.append((pair.tnblast_forward_fragments, "F'"))
                if pair.tnblast_reverse_fragments:
                    b_items.append((pair.tnblast_reverse_fragments, "R'"))
                b_items.sort(key=lambda item: item[0][0].start)

                for index, (frags, label) in enumerate(b_items):
                    alignment = "left" if index == 0 else "right"
                    _draw_fragments(
                        canvas,
                        frags,
                        curr_y - 1,
                        4,
                        v_min,
                        v_span,
                        origin_x,
                        draw_w,
                        PANEL_B_COLOR,
                        PANEL_B_FILL,
                        label,
                        alignment,
                    )

            curr_y -= PRIMER_ROW_H + PAIR_GAP

    return curr_y - 5


# ── Flowables ────────────────────────────────────────────────────────────────


class _SequenceRow(Flowable):
    """Draws a single 100-base row of DNA sequence with indices."""

    def __init__(self, sequence: str, start_index: int, pairs: list[PrimerPair]):
        super().__init__()
        self.sequence = sequence
        self.start_index = start_index
        self.pairs = pairs
        self.line_h = 10
        self.height = 3.5 * self.line_h
        self.width = 600

    def wrap(self, aW, aH):  # noqa: N803
        return (self.width, self.height)

    def draw(self):
        canv = self.canv
        canv.setFont("Courier", 8)
        char_w = canv.stringWidth("A", "Courier", 8)

        seq = self.sequence
        idx = self.start_index

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
            blocks.append(seq[b : b + 20])
        canv.drawString(40, 0.5 * self.line_h, " ".join(blocks))

        end_idx = idx + len(seq) - 1
        canv.drawString(40 + 100 * char_w + 7 * char_w, 0.5 * self.line_h, str(end_idx))


class _GeneDiagram(Flowable):
    """A flowable that wraps :func:`draw_gene_diagram` for one chain."""

    def __init__(
        self,
        locus: GeneLocus,
        chain: ConservedExonChain,
        chain_pairs: list[PrimerPair],
        contributing_tids: set[str],
        v_min: int,
        v_max: int,
        width: float,
        height: float,
    ):
        super().__init__()
        self.locus = locus
        self.chain = chain
        self.chain_pairs = chain_pairs
        self.contributing_tids = contributing_tids
        self.v_min = v_min
        self.v_max = v_max
        self.width = width
        self.height = height

    def wrap(self, aW, aH):  # noqa: N803
        # Clamp to the actual available width so we never exceed the
        # frame (which has 6pt internal padding on each side).
        if aW < self.width:
            self.width = aW
        return (self.width, self.height)

    def draw(self):
        draw_gene_diagram(
            self.canv,
            0,
            self.height - 30,
            self.width,
            self.locus,
            self.chain,
            self.chain_pairs,
            self.contributing_tids,
            self.v_min,
            self.v_max,
        )
