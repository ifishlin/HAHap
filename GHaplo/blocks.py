#!/usr/bin/env python3
"""
Search SNPs blocks and then output to stdout or to main prog.
"""
import logging
import sys
import pysam
#from concurrent.futures import ProcessPoolExecutor

#def func(i):
#    return i[1]


BAM_CMATCH = 0
BAM_CINS   = 1
BAM_CDEL   = 2
BAM_CREF_SKIP  = 3
BAM_CSOFT_CLIP = 4
BAM_CHARD_CLIP = 5

logger = logging.getLogger()

def b_search_while(lst, target):
    left = 0
    right = len(lst)-1
    while(left <= right):
        avg = (left + right)//2
        mid = lst[avg]

        if (mid < target):
            left = avg + 1
        elif (mid > target):
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


def trans_CigarTuples_To_Region(idx, cigarstuples):
    region = []
    start = idx
    end   = 0
    for c, l in cigarstuples:

        if c == BAM_CMATCH:
            end = start + l - 1
            region.append([start, end])
            start = end + 1
        elif c == BAM_CINS:
            pass
        elif c == BAM_CDEL:
            start = start + l #next
        elif c == BAM_CREF_SKIP:
            pass
        else:
            pass

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


def main(args, chrom, locus_list, locus_Mm_list, locus_idx_list, timer):

    #ee = ProcessPoolExecutor()
 
    locus_list.append(sys.maxsize)
    len_locus = len(locus_list)
    locus_listed_dict = [set() for i in range(len_locus - 1)]

    samfile = pysam.AlignmentFile(args.input_file, "rb" )

    header = samfile.header

    PG = header['PG']
    SQ = header['SQ']
    RG = header['RG']
    HD = header['HD']

    reference_name = [q['SN'] for q in SQ]

    read_map = dict()

    #print('fetch from ',chrom, locus_list[0]-1, locus_list[-2])
    #0-based, end value not included
    iter1 = samfile.fetch(reference=chrom, start=locus_list[0]-1, end=locus_list[-2])

    logger.info("hash implement")
    pre_locus = 0
    pre_c = 0
    c = 0
    for read in iter1:
        if read.is_proper_pair:
        #if True:
            timer.start('pre')
            connected = set()
            s = read.reference_start + 1 # transfer 1-base to 0-base
            if read.is_read1 == True:
                key = read.query_name + "_" + str(read.next_reference_start)
            else:
                key = read.query_name + "_" + str(read.reference_start)

            #if read.cigartuples == None:
            #    continue



            target = 'test1_2_105435_106338_0:0:0_0:0:0_13f7_16176089'

            flag = False
            if key == target:
                flag = True
                print(read.query_name, read.mapping_quality)

            if read.mapping_quality < args.mms:
                continue

            if s > pre_locus:
                c = b_search_while(locus_list, s)
                pre_c = c
                pre_locus = locus_list[c]
            else:
                c = pre_c
            timer.stop('pre')

            #print(aligned_locations)
            region = trans_CigarTuples_To_Region(s, read.cigartuples)
            timer.start('search')
            for s, e in region:
                while c < len_locus and  locus_list[c] <= e:
                    connected.add(c)
                    c += 1
            if flag:
                print(connected) 
            timer.stop('search')
            timer.start('putin')
            if key in read_map:
                uset = connected.union(read_map[key])
                '''
                if len(uset) > 1:
                    print("pair-end", key, uset, connected, read_map[key])
                '''
                make_pair_connected(uset, locus_listed_dict)
                del read_map[key]
            else:
                read_map[key] = connected
            timer.stop('putin')
        else: 
            pass
            #print(read.query_name, read.cigartuples)
            #print(read.get_aligned_pairs())
 
    # print("single")
    for k, v in read_map.items():
        make_pair_connected(v, locus_listed_dict)

    logger.info("end hash implement")

    samfile.close()

    connected_component = []
  
    timer.start('walk') 
    walked = set()
    for i in range(len(locus_list) - 1):
        cc = []
        build_cc(locus_listed_dict, i, cc, walked)
        if len(cc) > 1:
            connected_component.append(sorted(cc))
    timer.stop('walk')

    return connected_component, locus_Mm_list, locus_idx_list
 
