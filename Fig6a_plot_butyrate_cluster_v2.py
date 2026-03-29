#!/usr/bin/env python3
"""
Butyrate生合成遺伝子クラスター ゲノム配置図 v2
5遺伝子の連続配置を示す

Howto:
    python3 plot_butyrate_cluster_v2.py [width] [height] [fontsize]

example:
    python3 plot_butyrate_cluster_v2.py
    python3 plot_butyrate_cluster_v2.py 16 5 9
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.font_manager as fm
import sys
import os

# ============================================================
# font setting
# ============================================================
available = [f.name for f in fm.fontManager.ttflist]
if 'Helvetica' in available:
    plt.rcParams['font.family'] = 'Helvetica'
elif 'Arial' in available:
    plt.rcParams['font.family'] = 'Arial'
else:
    plt.rcParams['font.family'] = 'DejaVu Sans'

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype']  = 42

OUTPUT_DIR = 'butyrate_cluster'
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================
# Argument handling
# ============================================================
fig_w    = float(sys.argv[1]) if len(sys.argv) > 1 else 16.0
fig_h    = float(sys.argv[2]) if len(sys.argv) > 2 else 5.0
fontsize = float(sys.argv[3]) if len(sys.argv) > 3 else 9.0

# ============================================================
# gene cluster
# ============================================================
GENES = [
    {
        'id'      : 'BES00121.1',
        'product' : 'Acetyl-CoA\nC-acetyltransferase',
        'start'   : 2843472,
        'end'     : 2844656,
        'strand'  : -1,
        'color'   : '#8E44AD',   # violet：proteome confirmation
        'identity': '100%',
        'ec'      : 'EC:2.3.1.9',
        'proteome': True,
    },
    {
        'id'      : 'BES00122.1',
        'product' : 'acyl-CoA\ndehydrogenase',
        'start'   : 2844665,
        'end'     : 2845819,
        'strand'  : -1,
        'color'   : '#BDC3C7',   # grey：genome only
        'identity': '100%',
        'ec'      : '',
        'proteome': False,
    },
    {
        'id'      : 'BES00123.1',
        'product' : 'enoyl-CoA\nhydratase',
        'start'   : 2845850,
        'end'     : 2846644,
        'strand'  : -1,
        'color'   : '#E74C3C',   # red：proteome confirmation
        'identity': '100%',
        'ec'      : 'EC:4.2.1.17',
        'proteome': True,
    },
    {
        'id'      : 'BES00124.1',
        'product' : '3-hydroxybutyryl-CoA\ndehydrogenase',
        'start'   : 2846631,
        'end'     : 2847476,
        'strand'  : -1,
        'color'   : '#BDC3C7',   # grey：genome only
        'identity': '100%',
        'ec'      : 'EC:1.1.1.157',
        'proteome': False,
    },
    {
        'id'      : 'BES00125.1',
        'product' : 'CoA transferase\nsubunit B',
        'start'   : 2847590,
        'end'     : 2848237,
        'strand'  : -1,
        'color'   : '#BDC3C7',   # grey：genome only
        'identity': '99.5%',
        'ec'      : '',
        'proteome': False,
    },
]

# Range of figure
PADDING = 400
X_MIN   = min(g['start'] for g in GENES) - PADDING
X_MAX   = max(g['end']   for g in GENES) + PADDING
X_RANGE = X_MAX - X_MIN

# ============================================================
# Draw Block Arrow
# ============================================================
def draw_gene_arrow(ax, gene, y_center, height, fs):
    start  = gene['start']
    end    = gene['end']
    strand = gene['strand']
    color  = gene['color']

    gene_len = end - start
    head_len = min(gene_len * 0.15, X_RANGE * 0.02)

    if strand == 1:
        arrow = mpatches.FancyArrow(
            start, y_center, gene_len, 0,
            width=height, head_width=height * 1.1,
            head_length=head_len,
            fc=color, ec='white', linewidth=1.0,
            length_includes_head=True, zorder=3)
    else:
        arrow = mpatches.FancyArrow(
            end, y_center, -gene_len, 0,
            width=height, head_width=height * 1.1,
            head_length=head_len,
            fc=color, ec='white', linewidth=1.0,
            length_includes_head=True, zorder=3)
    ax.add_patch(arrow)

    x_center = (start + end) / 2

    # Top Label (Product Name + Proteome Mark/Logo)
    label = gene['product']
    if gene['proteome']:
        label += '\n★'
    ax.text(x_center, y_center + height * 0.65,
            label,
            ha='center', va='bottom',
            fontsize=fs, fontweight='bold', zorder=5)

    # Lower Label (ID/Match Rate)
    ax.text(x_center, y_center - height * 0.65,
            f"{gene['id']}\n({gene['identity']})",
            ha='center', va='top',
            fontsize=fs - 1.5,
            color='#555555', zorder=5)

    # ECnumber（arrows）
    if gene['ec']:
        ax.text(x_center, y_center,
                gene['ec'],
                ha='center', va='center',
                fontsize=fs - 2,
                color='white', fontweight='bold',
                zorder=6)

# ============================================================
# main visualization
# ============================================================
fig, ax = plt.subplots(figsize=(fig_w, fig_h))

Y_CENTER = 0.5
HEIGHT   = 0.28

# genom backbone line
ax.plot([X_MIN, X_MAX], [Y_CENTER, Y_CENTER],
        color='#7F8C8D', linewidth=2.5,
        zorder=1, solid_capstyle='round')

# visualization of each gene
for gene in GENES:
    draw_gene_arrow(ax, gene, Y_CENTER,
                    HEIGHT, fontsize)

# gap between genes
gaps = [
    (GENES[0]['end'], GENES[1]['start'], '9 bp'),
    (GENES[1]['end'], GENES[2]['start'], '31 bp'),
    (GENES[2]['end'], GENES[3]['start'], '13 bp'),  # ← 修正
    (GENES[3]['end'], GENES[4]['start'], '114 bp'),
]
for g_s, g_e, label in gaps:
    x_mid = (g_s + g_e) / 2
    ax.text(x_mid, Y_CENTER - HEIGHT * 1.8,
            label, ha='center', va='top',
            fontsize=fontsize - 2,
            color='#888888', style='italic')
    for x in [g_s, g_e]:
        ax.plot([x, x],
                [Y_CENTER - HEIGHT * 0.5,
                 Y_CENTER + HEIGHT * 0.5],
                color='#AAAAAA', linewidth=0.8,
                linestyle='--', zorder=2)

# scale bar
SCALE_LEN = 500
scale_x   = X_MAX - SCALE_LEN - 100
scale_y   = Y_CENTER - HEIGHT * 2.2
ax.plot([scale_x, scale_x + SCALE_LEN],
        [scale_y, scale_y],
        'k-', linewidth=2,
        solid_capstyle='round', zorder=5)
for x in [scale_x, scale_x + SCALE_LEN]:
    ax.plot([x, x], [scale_y - 0.02, scale_y + 0.02],
            'k-', linewidth=1.5, zorder=5)
ax.text(scale_x + SCALE_LEN / 2, scale_y - 0.05,
        '500 bp', ha='center', va='top',
        fontsize=fontsize - 1)

# – strand option
ax.annotate('',
            xy=(X_MIN + 150, Y_CENTER + 0.01),
            xytext=(X_MIN + 550, Y_CENTER + 0.01),
            arrowprops=dict(arrowstyle='->',
                           color='#555555', lw=1.2),
            zorder=4)
ax.text(X_MIN + 350, Y_CENTER + 0.06,
        '– strand', ha='center', va='bottom',
        fontsize=fontsize - 2,
        color='#555555', style='italic')

# title
ax.set_title(
    'Butyrate biosynthetic gene cluster (AP028807.1)\n'
    'Caldifermentibacillus hisashii homologs '
    '(99.5-100% identity)',
    fontsize=fontsize + 1,
    fontweight='bold', pad=12)

# legend
legend_elements = [
    mpatches.Patch(
        facecolor='#8E44AD', edgecolor='white',
        label='Acetyl-CoA C-acetyltransferase '
              '(BES00121.1, 100%) ★Proteome'),
    mpatches.Patch(
        facecolor='#BDC3C7', edgecolor='white',
        label='acyl-CoA dehydrogenase '
              '(BES00122.1, 100%) Genome only'),
    mpatches.Patch(
        facecolor='#E74C3C', edgecolor='white',
        label='enoyl-CoA hydratase '
              '(BES00123.1, 100%) ★Proteome'),
    mpatches.Patch(
        facecolor='#BDC3C7', edgecolor='white',
        label='3-hydroxybutyryl-CoA dehydrogenase '
              '(BES00124.1, 100%) Genome only'),
    mpatches.Patch(
        facecolor='#BDC3C7', edgecolor='white',
        label='CoA transferase subunit B '
              '(BES00125.1, 99.5%) Genome only'),
]
ax.legend(
    handles=legend_elements,
    loc='upper center',
    bbox_to_anchor=(0.5, -0.15),
    ncol=2,
    fontsize=fontsize - 1,
    frameon=True, framealpha=0.9,
    edgecolor='lightgray')

# axis setting
ax.set_xlim(X_MIN, X_MAX)
ax.set_ylim(Y_CENTER - HEIGHT * 3.5,
            Y_CENTER + HEIGHT * 3.5)
ax.axis('off')

plt.tight_layout()

# ============================================================
# save
# ============================================================
pdf_path = os.path.join(OUTPUT_DIR,
                        'butyrate_cluster_v2.pdf')
png_path = os.path.join(OUTPUT_DIR,
                        'butyrate_cluster_v2.png')

plt.savefig(pdf_path, dpi=300, bbox_inches='tight')
plt.savefig(png_path, dpi=300, bbox_inches='tight')
plt.close()

print(f"フォント      : {plt.rcParams['font.family']}")
print(f"図サイズ      : {fig_w} x {fig_h} inch")
print(f"フォントサイズ: {fontsize} pt")
print(f"\n✅ 保存完了:")
print(f"   PDF: {pdf_path}")
print(f"   PNG: {png_path}")
print(f"\nクラスター情報（5遺伝子）:")
for g in GENES:
    mark = '★Proteome' if g['proteome'] else 'Genome only'
    print(f"  {g['id']}: "
          f"{g['start']:,}-{g['end']:,} bp "
          f"({g['identity']}) {mark}")
