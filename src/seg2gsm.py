import pandas as pd
import matlab_functions as mf
import numpy as np
import argparse
from datetime import datetime
import time
# record start time
start = time.time()

TODAY = datetime.now().strftime("%d%b%Y")


parser = argparse.ArgumentParser()
parser.add_argument('-i', '--id',
                    help='Cohort set name.',
                    required=True, type=str)
parser.add_argument('-s', '--sample_set',
                    help='Sample file listing all samples\' in cohort, even w/o SV event.',
                    required=True, type=str)
parser.add_argument('-v', '--seg_file',
                    help='The seg file (CBS/IGV format) to find your samples\' CNV events.',
                    required=True, type=str)
parser.add_argument('-x', '--cnv_blacklist_file',
                    help='CNV regions to exclude.',
                    required=False, type=str, default='../../data_tables/additional_gsm_inputs/CNV.hg19.bypos.111213.CR1_event_added.bed')
parser.add_argument('-a', '--arm_significance_file',
                    help='CNV arm-level regions to include.',
                    required=False, type=str, default='../../data_tables/additional_gsm_inputs/DLBCL_broad_significance.19Aug2024.tsv')
parser.add_argument('-f', '--focal_file',
                    help='CNV foca; regions to include.',
                    required=False, type=str,default='../../data_tables/additional_gsm_inputs/DLBCL_focal_peaks.18Aug2024.tsv')
parser.add_argument('-o', '--output_dir',
                    help='Output directory.',
                    required=False, type=str,default='./')
parser.add_argument('-g','--genome_build',
                    help='Genome build: hg19, hg38.',
                    required=False, type=str, default='hg19')

args = parser.parse_args()

segs = pd.read_csv(args.seg_file, sep='\t')
log2CR_field = segs.columns[-1]
segs.rename(columns={"Start.bp":"Start", "End.bp":"End"},inplace=True)

S = pd.read_csv(args.sample_set, sep='\t', index_col=0)
sample_set = sorted(list(S.index))

single_amp_threshold = 0.1
single_del_threshold = -0.1
double_amp_threshold = 0.9
double_del_threshold = -0.9
arm_length_fraction_threshold = 2

AL = pd.read_csv(args.focal_file, sep='\t')
arm_level_significance = pd.read_csv(args.arm_significance_file, sep='\t')

if args.genome_build == 'hg19':
    segs.loc[:, 'gstart'] = mf.xhg19(segs['Chromosome'], segs['Start'])
    segs.loc[:, 'gend'] = mf.xhg19(segs['Chromosome'], segs['End'])
else:
    segs.loc[:, 'gstart'] = mf.xhg38(segs['Chromosome'], segs['Start'])
    segs.loc[:, 'gend'] = mf.xhg38(segs['Chromosome'], segs['End'])

segs.loc[:, 'length'] = segs['End'] - segs['Start']

cnv_blacklist = pd.read_csv(args.cnv_blacklist_file, sep='\t')
segs = mf.apply_cnv_blacklist(segs, cnv_blacklist, AL, arm_level_significance)

old_seg_means = segs[log2CR_field].copy(deep=True)
# samples in seg file:
samples = sorted(segs['Sample'].unique())
# check for gaps in sample_set
if not (set(samples) == set(sample_set)):
    print('Warning: sample in seg file not identical to input sample_set: GSM will have gaps')
    sseg = set(samples)
    ssegx = [x for x in sample_set if x not in sseg]
    print(ssegx)

segs['log_segment_mean'] = segs[log2CR_field].copy(deep=True)
segs['Segment_Mean'] = np.power(2, (segs[log2CR_field] + 1)) - 2

for samp in sorted(segs['Sample'].unique()):
    sampseg = segs.loc[segs['Sample'] == samp].copy(deep=True)
    sample_median = mf.calc_region_median(sampseg, min(arm_level_significance['x1']), max(arm_level_significance['x2']), 2)
    segs.loc[segs['Sample'] == samp, log2CR_field] = segs.loc[segs['Sample'] == samp, log2CR_field] - sample_median

#segs['log_segment_mean'] = segs[log2CR_field].copy(deep=True)
#segs['Segment_Mean'] = np.power(2, (segs[log2CR_field] + 1)) - 2

BA = arm_level_significance.loc[(arm_level_significance['significant_amplification'] == 1)] #&
#                                (arm_level_significance['amplification_cohort'].str.contains('DLBCL'))]

BD = arm_level_significance.loc[(arm_level_significance['significant_deletion'] == 1)]  #&
#                                (arm_level_significance['deletion_cohort'].str.contains('DLBCL'))]

# Arm dels
arm_del_df = pd.DataFrame(0, index=BD['arm'] + '.DEL', columns=sorted(sample_set)) #segs['Sample'].unique()))

print('Arm DELs')
for sample in arm_del_df.columns:
    for arm in BD['arm']:
        if ((sample == 'DLBCL10925') | (sample == 'DLBCL11206') ) & (arm  == '6q'):
            print(sample, arm)

        boundix = (BD['arm'] == arm)
        boundix = boundix[boundix].index

        bd_x1 = BD.loc[boundix, 'x1'].values[0]
        bd_x2 = BD.loc[boundix, 'x2'].values[0]

        c1 = ((segs['gstart'] <= bd_x1) & (segs['gend'] >= bd_x1))
        c2 = ((segs['gstart'] < bd_x2) & (segs['gend'] > bd_x2))
        c3 = ((segs['gstart'] >= bd_x1) & (segs['gend'] <= bd_x2))
        c4 = (segs['Sample'] == sample)

        segix = (c1 | c2 | c3) & c4
        if sum(segix) == 0:
            continue

        seg1 = segs.loc[segix]

        # sum(seg1.length) < (BD.x2(boundix) - BD.x1(boundix))/arm_length_fraction_threshold

        if (seg1['length'].sum() < ((bd_x2 - bd_x1) / arm_length_fraction_threshold)):
            continue

        region_median = mf.calc_region_median(seg1, bd_x1, bd_x2, 2)

        # region_median < single_del_threshold & region_median >= double_del_threshold

        if (region_median < single_del_threshold) & (region_median >= double_del_threshold):
            arm_del_df.loc[arm + '.DEL', sample] = 1
        elif region_median < double_del_threshold:
            arm_del_df.loc[arm + '.DEL', sample] = 2


# Arm amps
arm_amp_df = pd.DataFrame(0, index=BA['arm'] + '.AMP', columns=sorted(sample_set)) #columns=sorted(segs['Sample'].unique()))

# Use loops for now just to ensure code exact replication.
# This is VERY slow, and should be vectorized.
print('Arm AMPs')
for sample in arm_amp_df.columns:
    for arm in BA['arm']:

        if (sample == 'DLBCL10904') & (arm  == '6p'):
            print(sample, arm)
            
        boundix = (BA['arm'] == arm)
        boundix = boundix[boundix].index

        ba_x1 = BA.loc[boundix, 'x1'].values[0]
        ba_x2 = BA.loc[boundix, 'x2'].values[0]

        c1 = ((segs['gstart'] <= ba_x1) & (segs['gend'] >= ba_x1))
        c2 = ((segs['gstart'] < ba_x2) & (segs['gend'] > ba_x2))
        c3 = ((segs['gstart'] >= ba_x1) & (segs['gend'] <= ba_x2))
        c4 = (segs['Sample'] == sample)

        segix = (c1 | c2 | c3) & c4
        if sum(segix) == 0:
            continue

        seg1 = segs.loc[segix]

        # sum(seg1.length) < (BD.x2(boundix) - BD.x1(boundix))/arm_length_fraction_threshold
        if seg1['length'].sum() < ((ba_x2 - ba_x1) / arm_length_fraction_threshold):
            continue

        region_median = mf.calc_region_median(seg1, ba_x1, ba_x2, 2)

        if (region_median > single_amp_threshold) & (region_median <= double_amp_threshold):
            arm_amp_df.loc[arm + '.AMP', sample] = 1
        elif region_median > double_amp_threshold:
            arm_amp_df.loc[arm + '.AMP', sample] = 2

# Focals
focal_df = pd.DataFrame(0, index=AL['Descriptor'],columns=sorted(sample_set)) # columns=sorted(segs['Sample'].unique()))
#focal_df.insert(0, 'cohort', AL['cohort'].values)

print('Focals')
#for sample in focal_df.columns[1::]:
i=1
for sample in focal_df.columns:
    if (i % 10) == 0:
        print('\n',i,'.',end='')
    print(sample, end=' ')
    i += 1

    seg1 = segs.loc[segs['Sample'] == sample]

    if seg1.shape[0] == 0:
        continue

    for peak in AL['Descriptor']:
        if peak[1] in {'p', 'q'}:
            arm = peak[0:2]
        else:
            arm = peak[0:3]

        gstart_alix = AL.loc[AL['Descriptor'] == peak, 'gstart'].values[0]
        gend_alix = AL.loc[AL['Descriptor'] == peak, 'gend'].values[0]

        c1 = ((seg1['gstart'] < gstart_alix) & (seg1['gend'] > gstart_alix))
        c2 = ((seg1['gstart'] < gend_alix) & (seg1['gend'] > gend_alix))
        c3 = ((seg1['gstart'] > gstart_alix) & (seg1['gend'] < gend_alix))

        segix = c1 | c2 | c3
        segix = segix[segix]

        seg2 = seg1.loc[segix.index]

        if seg2.shape[0] == 0:
            continue

        if ((arm in BA['arm'].values) & ('AMP' in peak)) | ((arm in BD['arm'].values) & ('DEL' in peak)):
            bound1 = arm_level_significance.loc[arm_level_significance['arm'] == arm, 'x1'].values[0]
            bound2 = arm_level_significance.loc[arm_level_significance['arm'] == arm, 'x2'].values[0]
            arm_median = mf.calc_region_median(seg1, bound1, bound2, 2)
            seg2['Segment_Mean'] = seg2['Segment_Mean'] - arm_median

        if 'AMP' in peak:
            seg2['Segment_Mean'] = seg2['Segment_Mean'] * -1

        region_median = mf.calc_region_median(seg2, gstart_alix, gend_alix, 5)

        if 'AMP' in peak:
            region_median = region_median * -1


        if (('AMP' in peak) and (region_median > double_amp_threshold)) or (('DEL' in peak) and (region_median < double_del_threshold)):
            focal_df.loc[peak, sample] = 2
        elif (('AMP' in peak) and (region_median <= double_amp_threshold) and (region_median > single_amp_threshold)) or \
             (('DEL' in peak) and (region_median >= double_del_threshold) and (region_median) < single_del_threshold):
            focal_df.loc[peak, sample] = 1

focal_df.index = focal_df.index.str.replace(':', '.')
#focal_df = focal_df.loc[focal_df['cohort'].str.contains('DLBCL')]
# end time
end = time.time()
print("Execution time :", (end-start), "sec")

scna_df = pd.concat([arm_del_df, arm_amp_df, focal_df])
scna_df.index = scna_df.index.str.upper()
gsm = scna_df.reset_index().copy()
gsm = gsm.rename(columns={"index":"classifier_name"}).copy()
# GSM.index = GSM.index.str.upper()
outfile = args.output_dir + args.id + '.' + TODAY + '.CNV.GSM.tsv'
print('output :', outfile)
gsm.to_csv(outfile, sep='\t', index=False)
