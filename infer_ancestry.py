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
    output = []

    def euclidean_distance(v1, v2):
        return math.sqrt(sum((x - y)**2 for x, y in zip(v1, v2)))

    for index, study_sample in summary.iterrows():
        sample_id = study_sample[0]
        study_coords = [study_sample[f'{reference_name}_PC1_INTENDED'], study_sample[f'{reference_name}_PC2_INTENDED'], study_sample[f'{reference_name}_PC3_INTENDED'], study_sample[f'{reference_name}_PC4_INTENDED']]
        distances = reference[['PC1', 'PC2', 'PC3', 'PC4', 'POP']].copy()
        distances['DIST'] = distances.apply(lambda row: euclidean_distance(study_coords, [row['PC1'], row['PC2'], row['PC3'], row['PC4']]), axis = 1)
        pop, n = Counter(distances.sort_values(by = 'DIST')['POP'][:10]).most_common(1)[0]
        output.append([sample_id, pop, n])

    data = pd.DataFrame(output)
    data.columns = ['Sample', 'Inferred ancestry', '# nearest neighbours']
    return data


if __name__ == '__main__':
    args = argparser.parse_args()

    summary = pd.read_csv(args.summary_filename, header = 0, sep = ' ')

    pca_1000g = pd.read_csv(args.pca_1000g_filename, sep = '\t', header = None, names = ['SAMPLE', 'PC1', 'PC2', 'PC3', 'PC4'], usecols = [0, 1, 2, 3, 4])
    pop_1000g = pd.read_csv(args.pop_1000g_filename, sep = '\t')
    assert len(pca_1000g) == len(pop_1000g), f'{len(pca_hgdp)} != {len(samples_hgdp)}'
    samples_1000g = pca_1000g.merge(pop_1000g)
    assert len(samples_1000g) == len(pca_1000g), f'{len(pca_hgdp)} != {len(samples_hgdp)}'

    inferred_ancestry = infer_genetic_ancestry(summary, samples_1000g, '1000G')
    inferred_ancestry.to_csv('inferred_1000g.csv', index=False)
 
    if args.pop_hgdp_filename is not None and args.pca_hgdp_filename is not None:
        pca_hgdp = pd.read_csv(args.pca_hgdp_filename, sep = '\t', header = None, names = ['SAMPLE', 'PC1', 'PC2', 'PC3', 'PC4'], usecols = [0, 1, 2, 3, 4])
        pop_hgdp = pd.read_csv(args.pop_hgdp_filename, sep = '\t')
        samples_hgdp = pca_hgdp.merge(pop_hgdp)
        assert len(samples_hgdp) == len(pca_hgdp), f'{len(pca_hgdp)} != {len(samples_hgdp)}'
        inferred_ancestry = infer_genetic_ancestry(summary, samples_hgdp, 'HGDP')
        inferred_ancestry.to_csv('inferred_HGDP.csv', index=False)
