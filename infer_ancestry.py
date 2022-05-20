#!/usr/bin/env python

import argparse
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
import math

argparser = argparse.ArgumentParser(description = 'Generates plots from the summary file.')
argparser.add_argument('-s', '--summary', metavar = 'file', dest = 'summary_filename', required = True, help = 'Name of the summary file.')
argparser.add_argument('-pca1', '--pca-1000g', metavar = 'file', dest = 'pca_1000g_filename', required = True, help = 'File with the PCs from 1000G.')
argparser.add_argument('-pop1', '--populations-1000g', metavar = 'file', dest = 'pop_1000g_filename', required = True, help = 'File with the populations from 1000G.')
argparser.add_argument('-pca2', '--pca-hgdp', metavar = 'file', dest = 'pca_hgdp_filename', required = False, help = 'File with the PCs from HGDP.')
argparser.add_argument('-pop2', '--populations-hgdp', metavar = 'file', dest = 'pop_hgdp_filename', required = False, help = 'File with the populations from HGDP.')

def infer_genetic_ancestry(summary, reference, reference_name):
    counts = {
        6: Counter(),
        7: Counter(),
        8: Counter(),
        9: Counter(),
        10: Counter()
    }

    def euclidean_distance(v1, v2):
        return math.sqrt(sum((x - y)**2 for x, y in zip(v1, v2)))

    for index, study_sample in summary.iterrows():
        study_coords = [study_sample[f'{reference_name}_PC1_INTENDED'], study_sample[f'{reference_name}_PC2_INTENDED'], study_sample[f'{reference_name}_PC3_INTENDED'], study_sample[f'{reference_name}_PC4_INTENDED']]
        distances = reference[['PC1', 'PC2', 'PC3', 'PC4', 'POP']].copy()
        distances['DIST'] = distances.apply(lambda row: euclidean_distance(study_coords, [row['PC1'], row['PC2'], row['PC3'], row['PC4']]), axis = 1)
        pop, n = Counter(distances.sort_values(by = 'DIST')['POP'][:10]).most_common(1)[0]
        if n >= 6:
            for i in range(6, n + 1):
                counts[i][pop] += 1
            for i in range(n + 1, 11):
                counts[i]['Undefined'] += 1
        else:
            for i in range(6, 11):
                counts[i]['Undefined'] += 1

    return counts

  def tabulate_inferred_genetic_ancestry(inferred_ancestries, reference_colors, reference_name, output_filename):
    inferred_ancestries_table = []
    for i in range(6, 11):
        entry = { 'k-NN': i }
        for pop in list(reference_colors.keys()) + ['Undefined']:
            entry[pop] = inferred_ancestries[i][pop]
        inferred_ancestries_table.append(entry)
    inferred_ancestries_table = pd.DataFrame(inferred_ancestries_table)

    fig = plt.figure(figsize=(13, 3), dpi = 300)
    ax = plt.subplot(111, frame_on = False)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    table = ax.table(cellText = inferred_ancestries_table.values, colLabels = inferred_ancestries_table.columns, loc = "center", colColours = ['lightgrey'] * (len(reference_colors) + 2), fontsize = 14)
    ax.set_title(f'Inferred genetic ancestry based on k-nearest neighbours (k-NN) in {reference_name}', y = 0.75)
    table.auto_set_font_size(False)
    plt.tight_layout()
    plt.savefig(output_filename, bbox_inches='tight')


if __name__ == '__main__':
    args = argparser.parse_args()

    summary = pd.read_csv(args.summary_filename, header = 0, sep = ' ')

    pca_1000g = pd.read_csv(args.pca_1000g_filename, sep = '\t', header = None, names = ['SAMPLE', 'PC1', 'PC2', 'PC3', 'PC4'], usecols = [0, 1, 2, 3, 4])
    pop_1000g = pd.read_csv(args.pop_1000g_filename, sep = '\t')
    assert len(pca_1000g) == len(pop_1000g), f'{len(pca_hgdp)} != {len(samples_hgdp)}'
    samples_1000g = pca_1000g.merge(pop_1000g)
    assert len(samples_1000g) == len(pca_1000g), f'{len(pca_hgdp)} != {len(samples_hgdp)}'

    inferred_ancestry = infer_genetic_ancestry(summary, samples_1000g, '1000G')
    tabulate_inferred_genetic_ancestry(inferred_ancestry, colors_1000g, '1000G', '1000G_inferred_ancestry_table.jpeg')
 
    if args.pop_hgdp_filename is not None and args.pca_hgdp_filename is not None:
        pca_hgdp = pd.read_csv(args.pca_hgdp_filename, sep = '\t', header = None, names = ['SAMPLE', 'PC1', 'PC2', 'PC3', 'PC4'], usecols = [0, 1, 2, 3, 4])
        pop_hgdp = pd.read_csv(args.pop_hgdp_filename, sep = '\t')
        samples_hgdp = pca_hgdp.merge(pop_hgdp)
        assert len(samples_hgdp) == len(pca_hgdp), f'{len(pca_hgdp)} != {len(samples_hgdp)}'
        inferred_ancestry = infer_genetic_ancestry(summary, samples_hgdp, 'HGDP')
        tabulate_inferred_genetic_ancestry(inferred_ancestry, colors_hgdp, 'HGDP', 'HGDP_inferred_ancestry_table.jpeg')
