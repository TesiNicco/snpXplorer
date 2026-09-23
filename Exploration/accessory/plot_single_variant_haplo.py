# Tentative script to generate a figure recapitulating the single-variant locus
"""Plot CADD scores across a locus for a query variant and its LD partners.

Offline companion to ``scatterplot_plotly`` in ``snpxplorer_functions.py``. It
reuses the same visual grammar (stacked Plotly panels, shared x-axis in Mb,
gene track underneath) but answers a different question: instead of showing
GWAS significance per dataset, it shows the *deleteriousness* of every variant
in LD with the query SNP.

  * panel 1 - CADD PHRED per variant, coloured by LD r2 with the query SNP and
              shaped by whether the variant is coding and/or a QTL
  * panel 2 - gene track for the same window

Input is the ``info`` payload produced by ``SNPbot.routes.run_variant_query``
(as stored in a pickled workspace); genes are read straight from the SNPXplorer
tabix-indexed gene file, so this script does not import the Flask app.

Usage
-----
    python3 plot_single_variant_haplo.py --pkl ~/Downloads/temp_workspace.pkl \
        --out rs7412_cadd.html
"""

from __future__ import annotations
import argparse
import html
import math
import os
import pickle
import types
import sys
from pathlib import Path
from typing import Any, Optional, Sequence
import plotly.graph_objs as go
from plotly.subplots import make_subplots

# ---------------------------------------------------------------------------
# One-time helper that produced the example workspace pickle. Preserved for
# reference; it is inert during a normal run so executing this script can never
# overwrite an existing pickle. Flip the flag (or paste the body into the
# interactive session that built `info`) to regenerate the snapshot.
# ---------------------------------------------------------------------------

SAVE_WORKSPACE_SNAPSHOT = False
WORKSPACE_SNAPSHOT_PATH = "/Users/nicco/Downloads/temp_workspace.pkl"

if SAVE_WORKSPACE_SNAPSHOT:
    # Save manually made example workspace -- ONE TIME OPERATION
    workspace = {}
    for k, v in globals().copy().items():
        if k.startswith("_"):
            continue
        if isinstance(v, (types.FunctionType, types.ModuleType, type)):
            continue
        if k == "workspace":
            continue
        try:
            pickle.dumps(v)
            workspace[k] = v
        except Exception:
            pass
    with open(WORKSPACE_SNAPSHOT_PATH, "wb") as f:
        pickle.dump(workspace, f)

# ---------------------------------------------------------------------------
# Styling constants
# ---------------------------------------------------------------------------

# LD colouring, identical to assign_color() in snpxplorer_functions.py so the
# two figures stay legible side by side.
LD_BINS = [
    (0.8, 'red', 'r² ≥ 0.8'),
    (0.6, 'orange', '0.6 ≤ r² < 0.8'),
    (0.4, 'green', '0.4 ≤ r² < 0.6'),
    (0.2, 'lightblue', '0.2 ≤ r² < 0.4'),
]
LD_COLOR_BELOW = 'grey'
LD_LABEL_BELOW = 'r² < 0.2'

# Marker shape encodes the functional class of the variant.
SHAPE_CODING_QTL = 'star'
SHAPE_CODING = 'diamond'
SHAPE_QTL = 'triangle-up'
SHAPE_OTHER = 'circle'
SHAPE_LEGEND = [
    (SHAPE_CODING_QTL, 'Coding + QTL'),
    (SHAPE_CODING, 'Coding'),
    (SHAPE_QTL, 'QTL'),
    (SHAPE_OTHER, 'Other'),
]

TARGET_COLOR = 'black'

# Column layout of Data/databases/Genes/genes_hg38.txt.gz, one transcript per
# line: name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount,
# exonStarts, exonEnds, geneName.
GENE_COL_TRANSCRIPT = 0
GENE_COL_STRAND = 2
GENE_COL_TXSTART = 3
GENE_COL_TXEND = 4
GENE_COL_EXON_N = 7
GENE_COL_EXON_STARTS = 8
GENE_COL_EXON_ENDS = 9
GENE_COL_NAME = 10

# How many QTL lines to spell out in a tooltip before collapsing the rest.
MAX_QTL_HOVER_ROWS = 6

# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

def r2_color(r2):
    """Map an LD r2 value onto the SNPXplorer colour bins."""
    if r2 is None or (isinstance(r2, float) and math.isnan(r2)):
        return LD_COLOR_BELOW
    for threshold, color, _label in LD_BINS:
        if r2 >= threshold:
            return color
    return LD_COLOR_BELOW


def variant_shape(is_coding, is_qtl):
    """Map the coding/QTL flags onto a Plotly marker symbol."""
    if is_coding and is_qtl:
        return SHAPE_CODING_QTL
    if is_coding:
        return SHAPE_CODING
    if is_qtl:
        return SHAPE_QTL
    return SHAPE_OTHER


def is_coding_annotation(annotypes, consequences=None):
    """CADD calls a variant coding via AnnoType; mirror routes.py's test."""
    if 'coding' in str(annotypes or '').lower():
        return True
    # Belt and braces: some rows only carry the consequence.
    coding_terms = ('NON_SYNONYMOUS', 'SYNONYMOUS', 'STOP_GAINED', 'STOP_LOST', 'FRAME_SHIFT', 'INFRAME')
    cons = str(consequences or '').upper()
    return any(term in cons for term in coding_terms)


def _to_float(value):
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(out) else out


def _fmt(value, digits=3):
    num = _to_float(value)
    if num is None:
        return 'NA'
    if num != 0 and (abs(num) < 1e-3 or abs(num) >= 1e5):
        return '%.*e' % (digits - 1, num)
    return ('%.*f' % (digits, num)).rstrip('0').rstrip('.')


def _esc(value):
    return html.escape(str(value)) if value not in (None, '') else 'NA'


def _alleles_from_uniq(partner_uniq):
    """'44909976:G:T' -> ('G', 'T'); anything unexpected -> (None, None)."""
    parts = str(partner_uniq or '').split(':')
    if len(parts) >= 3:
        return parts[-2], parts[-1]
    return None, None


# ---------------------------------------------------------------------------
# Reshaping the `info` payload
# ---------------------------------------------------------------------------

def index_cadd_by_position(cadd_rows):
    """Group CADD rows by position, keyed by (ref, alt) where available.

    CADD reports one row per alternate allele, so a single position can carry
    several PHRED scores. We keep them all and let the caller pick the row that
    matches the observed allele, falling back to the most damaging one.
    """
    by_pos = {}
    for row in cadd_rows or []:
        pos = _to_float(row.get('pos'))
        if pos is None:
            continue
        by_pos.setdefault(int(pos), []).append(row)
    return by_pos


def pick_cadd_row(rows, ref=None, alt=None):
    """Prefer the row matching the observed alleles, else the highest PHRED."""
    if not rows:
        return None
    if ref and alt:
        for row in rows:
            if str(row.get('ref')) == str(ref) and str(row.get('alt')) == str(alt):
                return row
    return max(rows, key=lambda r: _to_float(r.get('phred_max')) or float('-inf'))


def index_qtl_by_position(qtl_rows):
    """Group eQTL/sQTL rows by position, strongest association first."""
    by_pos = {}
    for row in qtl_rows or []:
        pos = _to_float(row.get('pos'))
        if pos is None:
            continue
        by_pos.setdefault(int(pos), []).append(row)
    for rows in by_pos.values():
        rows.sort(key=lambda r: _to_float(r.get('pval_nominal')) if _to_float(r.get('pval_nominal')) is not None else 1.0)
    return by_pos


def _qtl_hover_block(title, rows):
    """Render the eQTL/sQTL section of a tooltip."""
    if not rows:
        return ''
    lines = ['<br><b>%s</b> (%d)' % (title, len(rows))]
    for row in rows[:MAX_QTL_HOVER_ROWS]:
        target = row.get('gene') or row.get('ensemble') or row.get('phenotype_id') or 'NA'
        lines.append(
            '<br>  <i>%s</i> · %s · p=%s · slope=%s'
            % (_esc(target), _esc(row.get('tissue')), _fmt(row.get('pval_nominal')), _fmt(row.get('slope')))
        )
    if len(rows) > MAX_QTL_HOVER_ROWS:
        lines.append('<br>  … and %d more' % (len(rows) - MAX_QTL_HOVER_ROWS))
    return ''.join(lines)


def build_hover(rsid, chrom, pos, ref, alt, phred, r2, dist_bp, cadd_row, eqtl_rows, sqtl_rows, is_target=False):
    """Assemble the full tooltip for one variant."""
    header = '<b>%s</b>' % _esc(rsid or '%s:%s' % (chrom, pos))
    if is_target:
        header += '  <i>(query variant)</i>'
    parts = [
        header,
        '<br>chr%s:%s %s→%s' % (_esc(chrom), '{:,}'.format(int(pos)), _esc(ref), _esc(alt)),
        '<br>CADD PHRED: <b>%s</b>' % _fmt(phred, 2),
    ]
    if is_target:
        parts.append('<br>LD: query variant (r² = 1)')
    else:
        parts.append('<br>r² with query: <b>%s</b>' % _fmt(r2, 3))
        if dist_bp is not None:
            parts.append(' · %s bp away' % '{:,}'.format(int(dist_bp)))
    if cadd_row:
        parts.append('<br>Gene(s): %s' % _esc(cadd_row.get('genes')))
        parts.append('<br>Annotation: %s' % _esc(cadd_row.get('annotypes')))
        parts.append('<br>Consequence: %s' % _esc(cadd_row.get('consequences')))
    parts.append(_qtl_hover_block('eQTL', eqtl_rows))
    parts.append(_qtl_hover_block('sQTL', sqtl_rows))
    return ''.join(parts)


def build_variant_table(info, include_target=True):
    """Flatten `info` into one plottable record per variant.

    Returns (records, n_missing_cadd). Variants without a CADD score cannot be
    placed on the y-axis and are reported back to the caller rather than
    silently dropped.
    """
    chrom = info.get('Chromosome')
    target_pos = info.get('Position_hg38')

    cadd_by_pos = index_cadd_by_position(info.get('ld_cadd'))
    cadd_by_pos.update(index_cadd_by_position(info.get('cadd')))
    eqtl_by_pos = index_qtl_by_position(info.get('ld_eqtl'))
    for pos, rows in index_qtl_by_position(info.get('eqtl')).items():
        eqtl_by_pos.setdefault(pos, []).extend(rows)
    sqtl_by_pos = index_qtl_by_position(info.get('ld_sqtl'))
    for pos, rows in index_qtl_by_position(info.get('sqtl')).items():
        sqtl_by_pos.setdefault(pos, []).extend(rows)

    records = []
    missing_cadd = []
    seen_positions = set()

    def add_record(pos, rsid, ref, alt, r2, dist_bp, is_target):
        pos = int(pos)
        if pos in seen_positions:
            return
        cadd_row = pick_cadd_row(cadd_by_pos.get(pos), ref, alt)
        phred = _to_float(cadd_row.get('phred_max')) if cadd_row else None
        if phred is None:
            missing_cadd.append(rsid or pos)
            return
        seen_positions.add(pos)
        eqtl_rows = eqtl_by_pos.get(pos, [])
        sqtl_rows = sqtl_by_pos.get(pos, [])
        is_qtl = bool(eqtl_rows or sqtl_rows)
        is_coding = is_coding_annotation(cadd_row.get('annotypes'), cadd_row.get('consequences'))
        records.append({
            'pos': pos,
            'pos_mb': pos / 1_000_000,
            'rsid': rsid or cadd_row.get('ld_rsid') or '%s:%s' % (chrom, pos),
            'phred': phred,
            'r2': 1.0 if is_target else _to_float(r2),
            'color': TARGET_COLOR if is_target else r2_color(_to_float(r2)),
            'symbol': variant_shape(is_coding, is_qtl),
            'is_target': is_target,
            'hover': build_hover(
                rsid, chrom, pos, cadd_row.get('ref') or ref, cadd_row.get('alt') or alt,
                phred, r2, dist_bp, cadd_row, eqtl_rows, sqtl_rows, is_target=is_target,
            ),
        })

    # The query variant first, so it always wins the position de-duplication.
    if include_target and target_pos is not None:
        add_record(target_pos, info.get('rsID'), info.get('REF'), info.get('ALT'), 1.0, 0, True)

    for partner in info.get('ld') or []:
        pos = partner.get('partner_pos')
        if pos is None:
            continue
        ref, alt = _alleles_from_uniq(partner.get('partner_uniq'))
        add_record(pos, partner.get('rsid'), ref, alt, partner.get('r2'), partner.get('dist_bp'), False)

    records.sort(key=lambda r: r['pos'])
    return records, missing_cadd


# ---------------------------------------------------------------------------
# Gene track
# ---------------------------------------------------------------------------

def extract_genes(data_path, chrom, start_pos, end_pos):
    """hg38 gene records overlapping the window, one transcript per gene.

    Mirrors extract_genes() in snpxplorer_functions.py (GRCh38 branch) but is
    inlined so this script stays importable without the Flask stack.
    """
    data_path = str(data_path)
    gene_file = os.path.join(data_path, 'databases/Genes/genes_hg38.txt.gz')
    if not os.path.exists(gene_file):
        print('warning: gene file not found at %s - gene track will be empty' % gene_file, file=sys.stderr)
        return []
    cmd = 'tabix %s chr%s:%s-%s' % (
        gene_file.replace(' ', r'\ '),
        str(chrom).replace('chr', ''),
        str(max(int(start_pos) - 1000, 1)),
        str(int(end_pos) + 1000),
    )
    genes = [x.rstrip().split('\t') for x in os.popen(cmd)]
    genes = [g for g in genes if len(g) > GENE_COL_NAME]

    # Several transcripts can share a gene name; keep the longest one.
    unique_genes = {}
    for g in genes:
        name = g[GENE_COL_NAME]
        span = int(g[GENE_COL_TXEND]) - int(g[GENE_COL_TXSTART])
        if name not in unique_genes or span > unique_genes[name][0]:
            unique_genes[name] = (span, g)
    genes = [g for _span, g in unique_genes.values()]

    genes = [g for g in genes if 'LOC' not in g[GENE_COL_NAME]]
    genes.sort(key=lambda g: int(g[GENE_COL_TXSTART]))
    # One lane per gene, as the main app does, so labels never collide.
    return [[*g, i] for i, g in enumerate(genes)]


def add_gene_track(fig, genes, start_pos, end_pos, row, show_exons=True):
    """Draw the gene panel: transcript body, strand arrow, exons, label."""
    if not genes:
        fig.add_annotation(
            x=(start_pos + (end_pos - start_pos) / 2) / 1_000_000, y=0.5,
            text="<b><span style='font-size:18px'>No genes to plot</span></b>",
            showarrow=False, font=dict(size=18, color='black', family='Arial'),
            row=row, col=1,
        )
        return

    max_y = max(g[-1] for g in genes)
    if len(genes) < 8:
        fsize, line_w = 15, 6
    elif len(genes) < 20:
        fsize, line_w = 11, 4
    else:
        fsize, line_w = 7, 2
    # Keep labels visibly above the transcript/exon strokes. A fixed offset is
    # more reliable than scaling by the lane index, especially for one-gene
    # loci where max_y is zero.
    yoff = 0.32 if len(genes) < 8 else (0.25 if len(genes) < 20 else 0.20)

    for g in genes:
        x0 = int(g[GENE_COL_TXSTART]) / 1_000_000
        x1 = int(g[GENE_COL_TXEND]) / 1_000_000
        lane = int(g[-1])
        strand = g[GENE_COL_STRAND]
        color = 'grey' if strand == '+' else 'black'

        fig.add_trace(
            go.Scatter(
                x=[x0, x1], y=[lane, lane], mode='lines',
                line=dict(color=color, width=line_w),
                hoverinfo='skip', showlegend=False,
            ), row=row, col=1)

        forward = strand == '+'
        fig.add_trace(
            go.Scatter(
                x=[x1 if forward else x0], y=[lane], mode='markers',
                marker=dict(symbol='triangle-right' if forward else 'triangle-left',
                            size=10, color=color, line=dict(width=0)),
                hoverinfo='skip', showlegend=False,
            ), row=row, col=1)

        if show_exons:
            starts = [int(x) / 1_000_000 for x in str(g[GENE_COL_EXON_STARTS]).split(',') if x]
            ends = [int(x) / 1_000_000 for x in str(g[GENE_COL_EXON_ENDS]).split(',') if x]
            for xs, xe in zip(starts, ends):
                fig.add_trace(
                    go.Scatter(
                        x=[xs, xe], y=[lane, lane], mode='lines',
                        line=dict(color=color, width=line_w * 2),
                        hoverinfo='skip', showlegend=False,
                    ), row=row, col=1)

        fig.add_annotation(
            x=x0 + (x1 - x0) / 2, y=lane + yoff,
            text='<b><i>%s</i></b>' % _esc(g[GENE_COL_NAME]),
            showarrow=False, yshift=12, font=dict(size=fsize),
            hovertext=('<b>%s</b><br>Transcript: %s<br>Strand: %s<br>TxStart: %s<br>TxEnd: %s<br>Exons: %s'
                       % (_esc(g[GENE_COL_NAME]), _esc(g[GENE_COL_TRANSCRIPT]), _esc(strand),
                          '{:,}'.format(int(g[GENE_COL_TXSTART])), '{:,}'.format(int(g[GENE_COL_TXEND])),
                          _esc(g[GENE_COL_EXON_N]))),
            row=row, col=1)

    fig.update_yaxes(
        range=[-0.35, max_y + 0.70],
        row=row,
        col=1,
    )


# ---------------------------------------------------------------------------
# Legends
# ---------------------------------------------------------------------------

def add_legend_proxies(fig, records, include_target):
    """Plotly has one legend, so colour and shape get their own legend groups.

    Data traces carry per-point colours and symbols and stay out of the legend;
    these invisible proxy traces carry the keys instead.
    """
    used_colors = {r['color'] for r in records}
    used_symbols = {r['symbol'] for r in records}

    for threshold, color, label in LD_BINS:
        if color not in used_colors:
            continue
        fig.add_trace(go.Scatter(
            x=[None], y=[None], mode='markers',
            marker=dict(size=12, color=color, symbol='circle', line=dict(width=1, color='#444')),
            name=label, legendgroup='ld', legendgrouptitle_text='LD with query',
            hoverinfo='skip', showlegend=True,
        ), row=1, col=1)
    if LD_COLOR_BELOW in used_colors:
        fig.add_trace(go.Scatter(
            x=[None], y=[None], mode='markers',
            marker=dict(size=12, color=LD_COLOR_BELOW, symbol='circle', line=dict(width=1, color='#444')),
            name=LD_LABEL_BELOW, legendgroup='ld', legendgrouptitle_text='LD with query',
            hoverinfo='skip', showlegend=True,
        ), row=1, col=1)

    for symbol, label in SHAPE_LEGEND:
        if symbol not in used_symbols:
            continue
        fig.add_trace(go.Scatter(
            x=[None], y=[None], mode='markers',
            marker=dict(size=12, color='#777', symbol=symbol, line=dict(width=1, color='#444')),
            name=label, legendgroup='class', legendgrouptitle_text='Variant class',
            hoverinfo='skip', showlegend=True,
        ), row=1, col=1)

    target = next((r for r in records if r['is_target']), None)
    if include_target and target is not None:
        # Mirror the query variant's own shape so the key matches the marker.
        fig.add_trace(go.Scatter(
            x=[None], y=[None], mode='markers',
            marker=dict(size=15, color=TARGET_COLOR, symbol=target['symbol'], line=dict(width=1.5, color='#444')),
            name='Query variant', legendgroup='query', legendgrouptitle_text=' ',
            hoverinfo='skip', showlegend=True,
        ), row=1, col=1)


# ---------------------------------------------------------------------------
# Main figure
# ---------------------------------------------------------------------------

def plot_single_variant_cadd(info, data_path, include_target=True, show_exons=True,
                             window=None, label_top_n=5, genes=None):
    """Two-panel CADD locus plot for a query variant and its LD partners.

    Parameters
    ----------
    info : dict
        The `info` payload from SNPbot.routes.run_variant_query().
    data_path : str or Path
        SNPXplorer Data/ directory, used to read the hg38 gene track.
    include_target : bool
        Also plot the query variant itself (in black).
    window : int, optional
        Half-width in bp around the query variant. Defaults to the span of the
        plotted variants plus 10 kb of padding.
    label_top_n : int
        Annotate this many variants with the highest CADD PHRED.
    genes : list, optional
        Pre-extracted gene records; extracted from `data_path` when omitted.
    """
    chrom = info.get('Chromosome')
    target_pos = info.get('Position_hg38')

    records, missing_cadd = build_variant_table(info, include_target=include_target)
    if not records:
        raise ValueError('no variant in this locus carries a CADD score - nothing to plot')
    if missing_cadd:
        print('note: %d variant(s) skipped, no CADD score: %s'
              % (len(missing_cadd), ', '.join(str(x) for x in missing_cadd[:10])), file=sys.stderr)

    positions = [r['pos'] for r in records]
    if window:
        start_pos = int(target_pos) - int(window)
        end_pos = int(target_pos) + int(window)
    else:
        pad = 10_000
        start_pos = min(positions) - pad
        end_pos = max(positions) + pad

    if genes is None:
        genes = extract_genes(data_path, chrom, start_pos, end_pos)

    fig = make_subplots(rows=2, cols=1, shared_xaxes=True, vertical_spacing=0.04,
                        row_heights=[0.64, 0.36])
    for r in (1, 2):
        fig.update_xaxes(showgrid=True, gridcolor='#e0e0e0', zeroline=False, row=r, col=1)
        fig.update_yaxes(showgrid=True, gridcolor='#e0e0e0', zeroline=False, row=r, col=1)
    fig.update_layout(plot_bgcolor='white', paper_bgcolor='white')

    # ---------- panel 1: CADD ----------
    partners = [r for r in records if not r['is_target']]
    if partners:
        # Stems make it easy to read a score off the axis in a dense locus.
        stem_x, stem_y = [], []
        for r in partners:
            stem_x.extend([r['pos_mb'], r['pos_mb'], None])
            stem_y.extend([0, r['phred'], None])
        fig.add_trace(go.Scatter(
            x=stem_x, y=stem_y, mode='lines',
            line=dict(color='#cccccc', width=1),
            hoverinfo='skip', showlegend=False,
        ), row=1, col=1)

        fig.add_trace(go.Scatter(
            x=[r['pos_mb'] for r in partners],
            y=[r['phred'] for r in partners],
            mode='markers',
            marker=dict(
                size=13,
                color=[r['color'] for r in partners],
                symbol=[r['symbol'] for r in partners],
                line=dict(width=1, color='#333333'),
            ),
            customdata=[[r['hover']] for r in partners],
            hovertemplate='%{customdata[0]}<extra></extra>',
            showlegend=False, name='LD partners',
        ), row=1, col=1)

    target = next((r for r in records if r['is_target']), None)
    if target is not None:
        fig.add_trace(go.Scatter(
            x=[target['pos_mb'], target['pos_mb']], y=[0, target['phred']],
            mode='lines', line=dict(color=TARGET_COLOR, width=1.5, dash='dot'),
            hoverinfo='skip', showlegend=False,
        ), row=1, col=1)
        fig.add_trace(go.Scatter(
            x=[target['pos_mb']], y=[target['phred']], mode='markers',
            marker=dict(size=20, color=TARGET_COLOR, symbol=target['symbol'],
                        line=dict(width=1.5, color='#333333')),
            customdata=[[target['hover']]],
            hovertemplate='%{customdata[0]}<extra></extra>',
            showlegend=False, name='Query variant',
        ), row=1, col=1)

    # Label the most deleterious variants; the query is always labelled.
    to_label = sorted(records, key=lambda r: r['phred'], reverse=True)[:max(label_top_n, 0)]
    if target is not None and target not in to_label:
        to_label.append(target)
    for r in to_label:
        fig.add_annotation(
            x=r['pos_mb'], y=r['phred'], xanchor='left', yanchor='middle',
            ax=25, ay=-18, showarrow=True, arrowwidth=1, arrowhead=2,
            text='<b>%s</b>' % _esc(r['rsid']),
            font=dict(size=13, family='Arial', color='black'),
            row=1, col=1)

    max_phred = max(r['phred'] for r in records)
    fig.update_yaxes(
        range=[0, max(max_phred * 1.20, 5)], row=1, col=1,
        title_text='CADD PHRED',
        title_font=dict(size=18, family='Arial', color='black', weight='bold'),
        tickfont=dict(size=14, family='Arial', color='black'))

    # CADD's conventional "top 1% most deleterious" cut-off.
    if max_phred >= 15:
        fig.add_hline(y=20, line=dict(color='#999999', width=1, dash='dash'),
                      annotation_text='PHRED 20', annotation_position='top left',
                      annotation_font=dict(size=11, color='#666666'),
                      row=1, col=1)

    add_legend_proxies(fig, records, include_target)

    # ---------- panel 2: genes ----------
    add_gene_track(fig, genes, start_pos, end_pos, row=2, show_exons=show_exons)
    fig.update_yaxes(title_text='Gene track', row=2, col=1, showticklabels=False,
                     title_font=dict(size=18, family='Arial', color='black', weight='bold'))
    fig.update_xaxes(
        title_text='Genomic position on chr%s (Mb, hg38)' % _esc(chrom), row=2, col=1,
        title_font=dict(size=18, family='Arial', color='black', weight='bold'),
        tickfont=dict(size=15, family='Arial', color='black'),
        showticklabels=True, ticks='outside', showline=True, automargin=True)

    xr = [start_pos / 1_000_000, end_pos / 1_000_000]
    fig.update_xaxes(range=xr, row=1, col=1)
    fig.update_xaxes(range=xr, row=2, col=1)

    fig.update_layout(
        hovermode='closest',
        hoverlabel=dict(align='left', bgcolor='white', font=dict(size=13, family='Arial')),
        legend=dict(font=dict(size=14), itemsizing='constant', groupclick='toggleitem',
                    yanchor='top', y=1.0, xanchor='left', x=1.01),
        margin=dict(l=60, r=40, t=50, b=60),
        height=700,
        autosize=True,
    )
    return fig


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def load_info_from_pickle(pkl_path):
    """Read a pickled run_variant_query() workspace and return (info, data_path)."""
    with open(pkl_path, 'rb') as handle:
        workspace = pickle.load(handle)
    if isinstance(workspace, dict) and 'info' in workspace:
        return workspace['info'], workspace.get('DATA_PATH')
    # Also accept a bare `info` dict.
    if isinstance(workspace, dict) and 'Position_hg38' in workspace:
        return workspace, None
    raise ValueError('%s does not look like a variant-query workspace' % pkl_path)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--pkl', required=True, help='pickled run_variant_query() workspace')
    parser.add_argument('--out', default='single_variant_cadd.html', help='output HTML file')
    parser.add_argument('--data-path', default=None, help='SNPXplorer Data/ directory (default: from the pickle)')
    parser.add_argument('--window', type=int, default=None, help='half-window in bp around the query variant')
    parser.add_argument('--no-target', action='store_true', help='plot only the LD partners')
    parser.add_argument('--no-exons', action='store_true', help='draw gene bodies without exons')
    parser.add_argument('--label-top', type=int, default=5, help='label this many top-CADD variants')
    parser.add_argument('--png', default=None, help='also write a static PNG (needs kaleido)')
    args = parser.parse_args(argv)

    info, pkl_data_path = load_info_from_pickle(args.pkl)
    data_path = args.data_path or pkl_data_path
    if data_path is None:
        parser.error('no Data/ directory in the pickle; pass --data-path')

    fig = plot_single_variant_cadd(
        info, data_path,
        include_target=not args.no_target,
        show_exons=not args.no_exons,
        window=args.window,
        label_top_n=args.label_top,
    )

    out = Path(args.out)
    fig.write_html(str(out), include_plotlyjs='cdn', full_html=True)
    print('wrote %s' % out.resolve())
    if args.png:
        fig.write_image(args.png, width=1400, height=900, scale=2)
        print('wrote %s' % Path(args.png).resolve())
    return 0


if __name__ == '__main__':
    sys.exit(main())
