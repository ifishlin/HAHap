#!/usr/bin/env python3
"""
Search SNPs blocks and then output to stdout or to main prog.
"""
import logging
import sys
import pysam
from bisect import bisect_left

logger = logging.getLogger()


def b_search_while(lst, target):
    left = 0
    right = len(lst)-1
    while left <= right:
        avg = (left + right)//2
        mid = lst[avg]

        if mid < target:
            left = avg + 1
        elif mid > target:
            right = avg - 1
        else:
            return avg

    return left


def make_pair_connected(cc, locus_listed_dict):
    l = list(cc)
    for i in l:
        for j in l[l.index(i)+1:]:
            locus_listed_dict[i].add(j)
            locus_listed_dict[j].add(i)

def trans_loc_to_alignidx(loc, r_st, cigartuples):  
    # only work with get_aligned_pairs(matches_only=True)
    #print("+++", loc, r_st)
    idx = 0
    d = loc - r_st + 1
    for c, l in cigartuples: 
        #print(c, l, d)
        if c == 0:
            if l >= d:
                return idx + d - 1
            else:
               d = d - l
               idx += l
        elif c == 1:
            idx += l
        elif c == 2:
            if l >= d:
                return None
            d = d - l 
        elif c == 4:
            idx += l
        #print("idx", idx, "d", d)
    return None

def trans_loc_to_cigaridx(loc, region):
    # only work with get_aligned_pairs(matches_only=True)
    idx = 0
    rs = region[0][0]
    for r in region:
        if loc in range(r[0], r[1]+1):
            idx += loc - r[0]
            return idx
        idx += (r[1] - r[0] + 1)
    return None

def trans_cigartuples_to_region(start, cigartuples):
    """
    """
    #BAM_CMATCH = 0
    #BAM_CINS = 1
    #BAM_CDEL = 2
    #BAM_CREF_SKIP = 3
    #BAM_CSOFT_CLIP = 4
    #BAM_CHARD_CLIP = 5

    region = []
    end = 0
    for c, l in cigartuples:

        if c == 0:
            end = start + l - 1
            region.append([start, end])
            start = end + 1
        #elif c == BAM_CINS:
        #    pass
        elif c == 2:
            start = start + l  # next
        #elif c == BAM_CREF_SKIP:
        #    pass

    return region


def build_cc(locus_listed_dict, i, cc, walked):
    if i not in walked:
        cc.append(i)
        walked.add(i)
        for j in locus_listed_dict[i]:
            build_cc(locus_listed_dict, j, cc, walked)


def add_arguments(parser):
    arg = parser.add_argument
    arg('variant_file', metavar='VCF', help='VCF file with variants needed to be phased')
    arg('input_file', metavar='BAM', help='BAM file')


def main(args, chrom, vars_loc, timer):
    timer.start('000.head') 
    vars_loc.append(sys.maxsize)
    len_locus = len(vars_loc)
    locus_listed_dict = [set() for i in range(len_locus - 1)]

    timer.stop('000.head')
    timer.start('000.open')
    samfile = pysam.AlignmentFile(args.bam_file, "rb")
    timer.stop('000.open')
    timer.start('000.head')
    header = samfile.header

    PG = header['PG']
    SQ = header['SQ']
    RG = header['RG']
    HD = header['HD']

    reference_name = [q['SN'] for q in SQ]

    read_map = dict()

    # 0-based, end value not included
    iter1 = samfile.fetch(reference=chrom, start=vars_loc[0]-1, end=vars_loc[-2])

    pre_locus = 0
    pre_c = 0
    c = 0
    timer.stop('000.head')
    timer.start('000.loop')
    for read in iter1:
        timer.start('000.if')
        if not read.is_proper_pair or read.mapping_quality < args.mms:
            continue
        timer.stop('000.if')

        timer.start('000.connect')
        connected = set()
        timer.stop('000.connect')
        timer.start('000.pre')
        s = read.reference_start + 1  # transfer 0-base to 1-base
        timer.stop('000.pre')

        timer.start('000.pre2')
        key = read.query_name

        timer.stop('000.pre2')
        timer.start('000.bisect')
        if s > pre_locus:
            c = bisect_left(vars_loc, s)
            # c = b_search_while(vars_loc, s)
            pre_c = c
            pre_locus = vars_loc[c]
        else:
            c = pre_c
        timer.stop('000.bisect')

        timer.start('003.trans')
        region = trans_cigartuples_to_region(s, read.cigartuples)
        timer.stop('003.trans')
        timer.start('001.search')
        for s, e in region:
            while c < len_locus and vars_loc[c] <= e:
                connected.add(c)
                c += 1
        timer.stop('001.search')

        '''
        timer.start('000.pre2')
        if len(connected) > 0:
            key = read.query_name
        timer.stop('000.pre2')
        '''

        timer.start('002.putin')
        if key in read_map:
            uset = connected.union(read_map[key])
            make_pair_connected(uset, locus_listed_dict)
            del read_map[key]
        else:
            read_map[key] = connected
        timer.stop('002.putin')
    timer.stop('000.loop')

    timer.start("000.end")
    for k, v in read_map.items():
        make_pair_connected(v, locus_listed_dict)

    timer.stop("000.end")

    timer.start("000.close")
    samfile.close()
    timer.stop("000.close")

    connected_component = []
  
    timer.start('003.walk') 
    walked = set()
    for i in range(len(vars_loc) - 1):
        cc = []
        build_cc(locus_listed_dict, i, cc, walked)
        if len(cc) > 1:
            connected_component.append(sorted(cc))
    timer.stop('003.walk')
    return connected_component
