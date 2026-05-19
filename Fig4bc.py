#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fig. 4b and 4c: Feature Metabolite Selection (Vertical layout)

Fig. 4b: Individual-level fold change heatmap (metabolites on y-axis)
Fig. 4c: Biosynthetic specificity evaluation table (metabolites on y-axis)

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec

plt.rcParams['font.family'] = 'Helvetica'
plt.rcParams['font.size'] = 11

# ============================================================
# Data
# ============================================================

individuals = ['No. 1293', 'No. 1295', 'No. 1296', 'No. 1298']
metabolites = ['Butyrate', 'AIB', 'Glyceraldehyde', 'Octadecanol',
               'Hydroquinone', 'Pyrogallol', 'Hydroxyisovalerate']

before = {
    'Water':    [84.389, 82.612, 85.015, 87.680],
    'Butyrate': [214.964, 166.606, 222.954, 213.756],
    'AIB':      [0.1561, 0.1485, 0.1398, 0.1254],
    'Glyceraldehyde': [0.0789, 0.0704, 0.1170, 0.0583],
    'Octadecanol':    [28.292, 29.789, 40.025, 38.777],
    'Hydroquinone':   [0.01132, 0.02253, 0.02867, 0.01155],
    'Pyrogallol':     [0.0, 0.01273, 0.02378, 0.01302],
    'Hydroxyisovalerate': [12.714, 10.233, 9.316, 5.379],
}

after = {
    'Water':    [82.833, 80.388, 84.910, 84.300],
    'Butyrate': [283.308, 426.861, 370.561, 713.751],
    'AIB':      [0.1614, 0.1705, 0.2027, 0.1898],
    'Glyceraldehyde': [0.1062, 0.1305, 0.1180, 0.1263],
    'Octadecanol':    [36.413, 47.374, 72.542, 61.656],
    'Hydroquinone':   [0.02514, 0.04593, 0.03509, 0.02384],
    'Pyrogallol':     [0.01419, 0.02934, 0.02808, 0.01820],
    'Hydroxyisovalerate': [13.749, 11.932, 10.506, 9.641],
}

# Calculate fold changes: rows=metabolites, cols=individuals
keys = ['Butyrate', 'AIB', 'Glyceraldehyde', 'Octadecanol',
        'Hydroquinone', 'Pyrogallol', 'Hydroxyisovalerate']

fc_data = np.zeros((7, 4))  # 7 metabolites x 4 individuals
for j, key in enumerate(keys):
    for i in range(4):
        if before[key][i] > 0:
            fc_data[j, i] = after[key][i] / before[key][i]
        else:
            fc_data[j, i] = np.nan

water_fc = [after['Water'][i] / before['Water'][i] for i in range(4)]


# ============================================================
# Fig. 4b: Vertical Heatmap (metabolites on y-axis)
# ============================================================

def plot_fig4b_vertical(ax):
    cmap = plt.cm.RdBu_r
    valid_fc = fc_data[~np.isnan(fc_data)]
    vmin = min(0.8, np.min(valid_fc))
    vmax = max(3.5, np.max(valid_fc))
    norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=1.0, vmax=vmax)

    im = ax.imshow(fc_data, cmap=cmap, norm=norm, aspect='auto')

    # FC text
    for j in range(7):
        for i in range(4):
            val = fc_data[j, i]
            if np.isnan(val):
                text = 'N/A'
                color = 'gray'
            else:
                text = f'{val:.2f}'
                color = 'white' if (val > 2.5 or val < 0.85) else 'black'
            ax.text(i, j, text, ha='center', va='center',
                    fontsize=11, fontweight='bold', color=color)

    # Highlight inconsistent cells
    # Glyceraldehyde (row 2) No.1296 (col 2)
    rect = plt.Rectangle((2 - 0.5, 2 - 0.5), 1, 1,
                          linewidth=2.5, edgecolor='black',
                          facecolor='none', linestyle='--')
    ax.add_patch(rect)

    # Pyrogallol (row 5) No.1293 (col 0)
    rect2 = plt.Rectangle((0 - 0.5, 5 - 0.5), 1, 1,
                           linewidth=2.5, edgecolor='black',
                           facecolor='none', linestyle='--')
    ax.add_patch(rect2)

    # Axis labels
    ax.set_xticks(range(4))
    ax.set_xticklabels(individuals, fontsize=11)
    ax.set_yticks(range(7))
    ax.set_yticklabels(metabolites, fontsize=11)

    # Water FC on top
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(range(4))
    ax2.set_xticklabels([f'Water FC\n{wfc:.2f}' for wfc in water_fc],
                         fontsize=9, color='steelblue', fontstyle='italic')
    ax2.tick_params(length=0)

    ax.set_title('Individual-level fold change\n(Final / Before)',
                 fontsize=14, fontweight='bold', pad=40)

    cbar = plt.colorbar(im, ax=ax, shrink=0.7, pad=0.03)
    cbar.set_label('Fold Change', fontsize=11)

    ax.text(0.5, -0.08,
            'Dashed boxes: inconsistent or incalculable',
            transform=ax.transAxes, ha='center',
            fontsize=9, fontstyle='italic', color='gray')


# ============================================================
# Fig. 4c: Vertical Evaluation Table (metabolites on y-axis)
# ============================================================

def plot_fig4c_vertical(ax):
    ax.axis('off')

    col_labels = ['Cliff\'s\nδ=1.0', 'Consistent\nFC', 'Biosyn.\nspecificity',
                  'Genome\nBGC', 'Proteome', 'Primary\norigin', 'Feature']

    row_labels = ['Butyrate', 'AIB', 'Glyceraldehyde', 'Octadecanol',
                  'Hydroquinone', 'Pyrogallol', 'Hydroxy-\nisovalerate']

    cells = [
        # Butyrate
        ['✓', '✓', 'High', '✓', '✓', 'Specific BGC', '★'],
        # AIB
        ['✓', '✓', 'High', '✓', '✓', 'Specific BGC', '★'],
        # Glyceraldehyde
        ['✓', '✗', 'Low', '–', '–', 'Glycolysis', ''],
        # Octadecanol
        ['✓', '✓', 'Low', '–', '–', 'General FAS', ''],
        # Hydroquinone
        ['✓', '✓', 'Low', '–', '–', 'Dietary/quinone', ''],
        # Pyrogallol
        ['✓', '✗', 'Low', '–', '–', 'Dietary polyphenol', ''],
        # Hydroxyisovalerate
        ['✓', '✓', 'Low', '–', '–', 'BCAA catabolism', ''],
    ]

    # Cell colors
    cell_colors = []
    for r in range(7):
        row_colors = []
        for c in range(7):
            if c == 6:  # Feature column
                if cells[r][c] == '★':
                    row_colors.append('#EF9A9A')
                else:
                    row_colors.append('#F5F5F5')
            elif cells[r][c] in ['✓', 'High'] and r <= 1:
                row_colors.append('#C8E6C9')
            elif cells[r][c] in ['✓']:
                row_colors.append('#E8F5E9')
            elif cells[r][c] in ['✗', '–', 'Low']:
                row_colors.append('#FFF3E0')
            else:
                row_colors.append('#FAFAFA')
        cell_colors.append(row_colors)

    table = ax.table(
        cellText=cells,
        rowLabels=row_labels,
        colLabels=col_labels,
        cellColours=cell_colors,
        rowColours=['#E3F2FD'] * 7,
        colColours=['#BBDEFB'] * 7,
        cellLoc='center',
        loc='center',
        bbox=[0.0, 0.0, 1.0, 1.0]
    )

    table.auto_set_font_size(False)
    table.set_fontsize(10)

    for key, cell in table.get_celld().items():
        row, col = key
        cell.set_edgecolor('#CCCCCC')
        cell.set_linewidth(0.5)
        cell.set_height(0.13)

        if row == 0:
            cell.set_text_props(fontweight='bold', fontsize=10)
            cell.set_facecolor('#BBDEFB')
        if col == -1:
            cell.set_text_props(fontweight='bold', fontsize=10)
            cell.set_facecolor('#BBDEFB')
            cell.set_width(0.18)

        if row > 0 and col >= 0:
            text = cells[row - 1][col]
            if text == '✓':
                cell.set_text_props(color='#2E7D32', fontweight='bold',
                                    fontsize=13)
            elif text == '✗':
                cell.set_text_props(color='#C62828', fontweight='bold',
                                    fontsize=13)
            elif text == '–':
                cell.set_text_props(color='#757575', fontsize=12)
            elif text == '★':
                cell.set_text_props(color='#C62828', fontweight='bold',
                                    fontsize=16)
            elif text == 'High':
                cell.set_text_props(color='#2E7D32', fontweight='bold')
            elif text == 'Low':
                cell.set_text_props(color='#E65100', fontweight='bold')

    ax.set_title('Evaluation of candidate metabolites (Cliff\'s δ = 1.0)',
                 fontsize=14, fontweight='bold', pad=15)


# ============================================================
# Combined vertical figure
# ============================================================

fig = plt.figure(figsize=(10, 16))
gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1], hspace=0.3)

ax_b = fig.add_subplot(gs[0])
plot_fig4b_vertical(ax_b)
ax_b.text(-0.12, 1.08, 'b', transform=ax_b.transAxes,
          fontsize=22, fontweight='bold', va='top')

ax_c = fig.add_subplot(gs[1])
plot_fig4c_vertical(ax_c)
ax_c.text(-0.12, 1.06, 'c', transform=ax_c.transAxes,
          fontsize=22, fontweight='bold', va='top')

plt.savefig('Fig4bc_vertical.pdf', bbox_inches='tight', dpi=300)
plt.savefig('Fig4bc_vertical.png', bbox_inches='tight', dpi=300)
plt.close()
print("Saved: Fig4bc_vertical.pdf / .png")

# ============================================================
# Separate panels
# ============================================================

fig_b, ax_b2 = plt.subplots(figsize=(7, 8))
plot_fig4b_vertical(ax_b2)
plt.savefig('Fig4b_vertical.pdf', bbox_inches='tight', dpi=300)
plt.savefig('Fig4b_vertical.png', bbox_inches='tight', dpi=300)
plt.close()
print("Saved: Fig4b_vertical.pdf / .png")

fig_c, ax_c2 = plt.subplots(figsize=(9, 7))
plot_fig4c_vertical(ax_c2)
plt.savefig('Fig4c_vertical.pdf', bbox_inches='tight', dpi=300)
plt.savefig('Fig4c_vertical.png', bbox_inches='tight', dpi=300)
plt.close()
print("Saved: Fig4c_vertical.pdf / .png")
