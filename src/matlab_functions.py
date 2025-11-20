import numpy as np
import pandas as pd


def xhg19(chromosome_v, start_positions, opt='hg19'):
    L = [249250621, 243199373, 198022430, 191154276, 180915260,
         171115067, 159138663, 146364022, 141213431, 135534747,
         135006516, 133851895, 115169878, 107349540, 102531392,
         90354753, 81195210, 78077248, 59128983, 63025520, 48129895,
         51304566, 155270560, 59373566, 16569]

    L18 = [247249719, 242951149, 199501827, 191273063, 180857866, 170899992,
           158821424, 146274826, 140273252, 135374737, 134452384, 132349534,
           114142980, 106368585, 100338915, 88827254, 78774742, 76117153, 63811651,
           62435964, 46944323, 49691432, 154913754, 57772954, 16571]

    Lmm9 = [197195432, 197195432, 181748087, 181748087, 159599783,
            159599783, 155630120, 155630120, 152537259, 152537259, 149517037,
            149517037, 152524553, 152524553, 131738871, 131738871, 124076172,
            124076172, 129993255, 129993255, 121843856, 121843856, 121257530,
            121257530, 120284312, 120284312, 125194864, 125194864, 103494974,
            103494974, 98319150, 98319150, 95272651, 95272651, 90772031, 90772031,
            61342430, 61342430, 0, 0, 0, 166650296, 166650296, 15902555, 16299]

    if opt == 'hg18':
        L = L18
    elif opt == 'mm9':
        L = Lmm9

    L = pd.Series(L)

    if opt == 'pad':
        PAD = 50e6
        L = L + PAD

    C = pd.concat([pd.Series([1]), L.cumsum()])
    C.index = np.arange(1, len(C) + 1)

    chromosome_v = chromosome_v.replace({'M': 25, 'MT': 25, 'X': 23, 'Y': 24}).astype(int)
    global_position = C.loc[chromosome_v] + start_positions.values

    return global_position.values


def apply_cnv_blacklist(segs, cnv_blacklist, AL, arm_level_significance):
    cnv_blacklist['gstart'] = xhg19(cnv_blacklist['Chromosome'], cnv_blacklist['Start'])
    cnv_blacklist['gend'] = xhg19(cnv_blacklist['Chromosome'], cnv_blacklist['End'])
    cnv_blacklist['keep_event'] = 0

    arm_level_significance = arm_level_significance.loc[(arm_level_significance['significant_amplification'].astype(bool)) |
                                                        (arm_level_significance['significant_deletion'].astype(bool))]

    AL['length'] = AL['gend'] - AL['gstart']
    cnv_blacklist['length'] = cnv_blacklist['gend'] - cnv_blacklist['gstart']
    arm_level_significance['length'] = arm_level_significance['x2'] - arm_level_significance['x1']

    # ((cnv_blacklist.gstart <= AL.gstart(i) & cnv_blacklist.gend >= AL.gstart(i)) |
    # (cnv_blacklist.gstart <= AL.gend(i) & cnv_blacklist.gend >= AL.gend(i)) |
    # (cnv_blacklist.gstart >= AL.gstart(i) & cnv_blacklist.gend <= AL.gend(i))) &
    # (cnv_blacklist.length > AL.length(i)*0.01)
    for idx, row in AL.iterrows():
        gstart_i = row['gstart']
        gend_i = row['gend']
        length_i = row['length']

        c1 = ((cnv_blacklist['gstart'] <= gstart_i) & (cnv_blacklist['gend'] >= gstart_i))
        c2 = ((cnv_blacklist['gstart'] <= gend_i) & (cnv_blacklist['gend'] >= gend_i))
        c3 = ((cnv_blacklist['gstart'] >= gstart_i) & (cnv_blacklist['gend'] <= gend_i))
        c4 = (cnv_blacklist['length'] > (length_i * 0.01))

        overlapix = (c1 | c2 | c3) & c4
        if sum(overlapix) != 0:
            overlapix = overlapix[overlapix].index
            cnv_blacklist.loc[overlapix, 'keep_event'] = 1

    # ((cnv_blacklist.gstart <= arm_level_significance.x1(i) & cnv_blacklist.gend >= arm_level_significance.x2(i)) |
    # (cnv_blacklist.gstart <= arm_level_significance.x1(i) & cnv_blacklist.gend >= arm_level_significance.x2(i)) |
    # (cnv_blacklist.gstart >= arm_level_significance.x1(i) & cnv_blacklist.gend <= arm_level_significance.x2(i))) &
    # ((cnv_blacklist.length > arm_level_significance.length(i)*0.01))
    for idx, row in arm_level_significance.iterrows():
        gstart_i = row['x1']
        gend_i = row['x2']
        length_i = row['length']

        c1 = ((cnv_blacklist['gstart'] <= gstart_i) & (cnv_blacklist['gend'] >= gend_i))
        c2 = ((cnv_blacklist['gstart'] <= gstart_i) & (cnv_blacklist['gend'] >= gend_i))
        c3 = ((cnv_blacklist['gstart'] >= gstart_i) & (cnv_blacklist['gend'] <= gend_i))
        c4 = (cnv_blacklist['length'] > (length_i * 0.01))

        overlapix = (c1 | c2 | c3) & c4
        if sum(overlapix) != 0:
            overlapix = overlapix[overlapix].index
            cnv_blacklist.loc[overlapix, 'keep_event'] = 1

    cnv_blacklist = cnv_blacklist.loc[cnv_blacklist['keep_event'] == 1]

    segs['remove_seg'] = 0

    for idx, bl_region in cnv_blacklist.iterrows():
        # (segs.gstart <= cnv_blacklist.gstart(i) & segs.gend >= cnv_blacklist.gstart(i)) |
        # (segs.gstart <= cnv_blacklist.gend(i) & segs.gend >= cnv_blacklist.gend(i)) |
        # (segs.gstart >= cnv_blacklist.gstart(i) & segs.gend <= cnv_blacklist.gend(i)) & ~segs.remove_seg)

        gstart_i = bl_region['gstart']
        gend_i = bl_region['gend']

        c1 = ((segs['gstart'] <= gstart_i) & (segs['gend'] >= gstart_i))
        c2 = ((segs['gstart'] <= gend_i) & (segs['gend'] >= gend_i))
        c3 = ((segs['gstart'] >= gstart_i) & (segs['gend'] <= gend_i))
        c4 = ~segs['remove_seg'].astype(bool)

        segix = c1 | c2 | c3 & c4
        if sum(segix) != 0:
            segix = segix[segix].index

        for ix in segix:
            if ix is False: 
                continue
            gstart_seg = segs.loc[ix, 'gstart']
            gend_seg = segs.loc[ix, 'gend']
            # (segs.gstart(segix(j)) < cnv_blacklist.gstart(i) & segs.gend(segix(j)) > cnv_blacklist.gend(i))
            if (gstart_seg < gstart_i) & (gend_seg > gend_i):
                seg_length = segs.loc[ix, 'length']
                cnv_length = gend_i - gstart_i
                if (cnv_length / seg_length) > 0.80:
                    segs.loc[ix, 'remove_seg'] = 1
            # (segs.gstart(segix(j)) < cnv_blacklist.gstart(i) & segs.gend(segix(j)) > cnv_blacklist.gstart(i))
            elif (gstart_seg < gstart_i) & (gend_seg > gstart_i):
                segs.loc[ix, 'gend'] = gstart_i
                segs.loc[ix, 'End'] = bl_region['Start']
            # (segs.gstart(segix(j))   < cnv_blacklist.gend(i) & segs.gend(segix(j)) > cnv_blacklist.gend(i))
            elif (gstart_seg < gend_i) & (gend_seg > gend_i):
                segs.loc[ix, 'gstart'] = gend_i
                segs.loc[ix, 'Start'] = bl_region['End']
            # (segs.gstart(segix(j)) >= cnv_blacklist.gstart(i) & segs.gend(segix(j)) <= cnv_blacklist.gend(i))
            elif (gstart_seg >= gstart_i) & (gend_seg <= gend_i):
                segs.loc[ix, 'remove_seg'] = 1

    segs = segs.loc[segs['remove_seg'] == 0]
    return segs


def calc_region_median(segs, bound1, bound2, frac):
    sorted_segs = segs.sort_values(by='Segment_Mean')
    sorted_segs['length'] = np.maximum(np.minimum(sorted_segs['gend'], bound2) - np.maximum(sorted_segs['gstart'], bound1),
                                       0)
    sorted_segs = sorted_segs.loc[sorted_segs['length'] != 0]
    sorted_segs = sorted_segs.reset_index()
    medianix = sum(sorted_segs['length']) / frac

    pos = 0
    for idx, row in sorted_segs.iterrows():
        pos = pos + row['length']
        if pos >= medianix:
            region_median = row['Segment_Mean']
            return region_median

    print('done')

