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
                    help='CNV focal regions to include.',
                    required=False, type=str,default='../../data_tables/additional_gsm_inputs/DLBCL_focal_peaks.18Aug2024.tsv')
parser.add_argument('-o', '--output_dir',
                    help='Output directory.',
                    required=False, type=str,default='./')
parser.add_argument('-g','--genome_build',
                    help='Genome build: hg19, hg38.',
                    required=False, type=str, default='hg19')
parser.add_argument('-t', '--thresholds',
                    help='Threshold for CNV events: lenient (default) or stringent. lenient: single amp=0.1, single del=-0.1, double amp=0.9, double del=-0.9; stringent: single amp=0.321, single del=-0.415, double amp=0.81, double del=-2.0',
                    required=False, type=str, default='lenient')

args = parser.parse_args()

segs = pd.read_csv(args.seg_file, sep='\t')
# standardize seg column names 
c = segs.columns
segs.rename(columns={c[0]: 'Sample',c[1]:'Chromosome',c[2]:'Start',c[3]:'End',c[-1]:'Segment_Mean'}, inplace=True)        
log2CR_field = segs.columns[-1]
sample_field = segs.columns[0]
sample_set = sorted(list(set(segs[sample_field].unique())))
if len(args.sample_set) > 0:
    S = pd.read_csv(args.sample_set, sep="\t", header=None, index_col=0)
    sample_set = list(S.index)
# remove X and Y - seg2gsm file is for autosomes
segs = segs.loc[~segs['Chromosome'].isin(['X', 'Y'])].copy(deep=True)

# subset to sample_set 
segs = segs.loc[segs['Sample'].isin(sample_set)].copy(deep=True) # reset_index()

# Whether a segment is called gained or deleted depends on whether CN-2 exceeds one of these thresholds
if args.thresholds == 'lenient': 
    single_amp_threshold = 0.1
    single_del_threshold = -0.1
    double_amp_threshold = 0.9
    double_del_threshold = -0.9
elif args.thresholds == 'stringent':
    single_amp_threshold = 0.5
    single_del_threshold = -0.5
    double_amp_threshold = 1.5
    double_del_threshold = -1.5
else: 
    raise ValueError("Invalid value for --thresholds. Please specify 'lenient' or 'stringent'.")

# arm length fraction threshold of 2 corresponds to the 50% median 
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

# Segment_Mean in log2CR units
old_seg_means = segs[log2CR_field].copy(deep=True)
segs['Segment_Mean'] = segs[log2CR_field].copy(deep=True)

for samp in sorted(segs['Sample'].unique()):
    sampseg = segs.loc[segs['Sample'] == samp].copy(deep=True)
    sample_median = mf.calc_region_median(sampseg, min(arm_level_significance['x1']), max(arm_level_significance['x2']), 2)
    segs.loc[segs['Sample'] == samp, 'Segment_Mean'] = segs.loc[segs['Sample'] == samp, 'Segment_Mean'] - sample_median

segs['log_segment_mean'] = segs['Segment_Mean'].copy(deep=True)
segs['Segment_Mean'] = np.power(2, (segs['log_segment_mean'] + 1)) - 2


BA = arm_level_significance.loc[(arm_level_significance['significant_amplification'] == 1)] 

BD = arm_level_significance.loc[(arm_level_significance['significant_deletion'] == 1)] 

# Arm dels
arm_del_df = pd.DataFrame(0, index=BD['arm'] + '.DEL', columns=sorted(sample_set)) #segs['Sample'].unique()))

print('Arm DELs')
for sample in arm_del_df.columns:
    for arm in BD['arm']:

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

print('Arm AMPs')
for sample in arm_amp_df.columns:
    for arm in BA['arm']:
            
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
focal_df = pd.DataFrame(0, index=AL['Descriptor'],columns=sorted(sample_set)) 

print('Focals')
#i=1
for sample in focal_df.columns:
    #print(sample, end=' ')
    #i += 1

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

end = time.time()
print("Execution time :", (end-start), "sec")

scna_df = pd.concat([arm_del_df, arm_amp_df, focal_df])
scna_df = scna_df.reset_index()
scna_df = scna_df.rename(columns={"index":"classifier_name"}).copy()
scna_df.classifier_name = scna_df.classifier_name.str.upper()

# fix wonky focal region name 
scna_df["classifier_name"] = scna_df["classifier_name"].str.replace( "19Q13.32_1.DEL","19Q13.32.DEL")   

# DLBclass CNV region lables
CNVS = """10Q23.31.DEL 11P.AMP 11Q.AMP 11Q23.3.AMP 12P.AMP 12P13.2.DEL 
12Q.AMP 13Q.AMP 13Q14.2.DEL 13Q34.DEL 14Q32.31.DEL 15Q15.3.DEL 16Q12.1.DEL 
17P.DEL 17Q24.3.AMP 17Q25.1.DEL 18P.AMP 18Q.AMP 18Q21.32.AMP 18Q22.2.AMP 
18Q23.DEL 19P13.2.DEL 19P13.3.DEL 19Q.AMP 19Q13.32.DEL 19Q13.42.AMP 
1P13.1.DEL 1P31.1.DEL 1P36.11.DEL 1P36.32.DEL 1Q.AMP 1Q23.3.AMP 1Q32.1.AMP 
1Q42.12.DEL 21Q.AMP 2P16.1.AMP 2Q22.2.DEL 3P.AMP 3P21.31.DEL 3Q.AMP 3Q28.AMP 
3Q28.DEL 4Q21.22.DEL 4Q35.1.DEL 5P.AMP 5Q.AMP 6P.AMP 6P21.1.AMP 6P21.33.DEL 
6Q.DEL 6Q14.1.DEL 6Q21.DEL 7P.AMP 7Q.AMP 7Q22.1.AMP 8Q12.1.DEL 8Q24.22.AMP 
9P21.3.DEL 9P24.1.AMP 9Q.AMP 9Q21.13.DEL"""
CNVS = CNVS.replace('\n', '')
CNVS = CNVS.split(' ')
print(*CNVS)
idxKEEP = []
for idx1, row in scna_df.iterrows():
   if (row["classifier_name"] in CNVS):
   		idxKEEP = idxKEEP + [idx1]
scna_df = scna_df.loc[idxKEEP,:]   
scna_df = scna_df.sort_values(by='classifier_name')
scna_df = scna_df.reset_index(drop=True)

outfile = args.output_dir + "/" + args.id + '.' + TODAY + '.CNV.GSM.tsv'
scna_df.to_csv(outfile, sep='\t',index=False)
print('complete')

print('output :', outfile)