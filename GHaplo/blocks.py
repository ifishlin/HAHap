#!/usr/bin/env python3
"""
Search SNPs blocks and then output to stdout or to main prog.


"""
import logging
import sys
import pysam

BAM_CMATCH = 0
BAM_CINS   = 1
BAM_CDEL   = 2
BAM_CREF_SKIP  = 3
BAM_CSOFT_CLIP = 4
BAM_CHARD_CLIP = 5

def b_search_while(lst, target):
    left = 0
    right = len(lst)-1
    while(left <= right):
        avg = (left + right)//2
        mid = lst[avg]
        #print(left, right, avg, mid, target)
        if (mid < target):
            left = avg + 1
        elif (mid > target):
            right = avg - 1
        else:
            return avg
        #print("== ", left, right, ' ==')
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
        #print(c, l , start, end)
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


def main(args, chrom, locus_list, locus_Mm_list, locus_idx_list):
  
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
    for r in reference_name:

        iter1 = samfile.fetch(reference=chrom, start=locus_list[0], end=locus_list[-2]+1)

        logging.info("hash implement")
        count = 0
        proper_count = 0
        pre_locus = 0
        pre_c = 0
        c = 0
        singletons_count = 0
        for read in iter1:
            count += 1
            '''
            if read.is_unmapped == False and read.mate_is_unmapped == True:
                singletons_count += 1
                continue
            '''
            if read.is_proper_pair:
                proper_count += 1
                connected = set()
                s = read.reference_start
                if read.is_read1 == True:
                    key = read.query_name + "_" + str(read.next_reference_start)
                else:
                    key = read.query_name + "_" + str(read.reference_start)

                if s > pre_locus:
                    c = b_search_while(locus_list, s)
                    pre_c = c
                    pre_locus = locus_list[c]
                    '''
                    cp1 = c + 1
                    if cp1 < len_locus and s < locus_list[cp1]:
                        pre_c = cp1
                        pre_locus = locus_list[cp1]
                        c = cp1
                    else:
                        c = b_search_while(locus_list, s)
                        pre_c = c
                        pre_locus = locus_list[c]
                    '''
                else:
                    c = pre_c

                region = trans_CigarTuples_To_Region(s, read.cigartuples)
                for s, e in region:
                    while c < len_locus and  locus_list[c] <= e:
                        connected.add(c)
                        c += 1

                if key in read_map:
                    uset = connected.union(read_map[key])
                    make_pair_connected(uset, locus_listed_dict)
                    del read_map[key]
                else:
                    read_map[key] = connected

        logging.info("end hash implement")
        break

    samfile.close()

    single_count = 0
    for k, v in read_map.items():
        single_count += 1


    connected_component = []
    
    walked = set()
    for i in range(len(locus_list) - 1):
        cc = []
        build_cc(locus_listed_dict, i, cc, walked)
        if len(cc) > 5:
            connected_component.append(cc)
            print ("=== " ,len(cc) , " ===")
            for c in cc:
                print(c, locus_Mm_list[c])

    print (count, proper_count, single_count, singletons_count)

    print('in blocks')
    
    return connected_component, locus_Mm_list, locus_idx_list
 
## backup

def main2(args, indel=False):
   
    locus_list = []
    locus_Mm_list = []
    locus_idx_list = []
    variants_vcf = open(args.variant_file,'r')
    for line in variants_vcf:
        if line[0] == '#':
            continue
        e = line.split("\t")
        c, i, _, m, n = e[0:5]

        if not indel and (len(m) > 1 or len(n) > 1):
            continue
            
        locus_list.append(int(i))
        locus_Mm_list.append([m, n])
        locus_idx_list.append(i)

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
    for r in reference_name:

        iter1 = samfile.fetch(reference='1', start=locus_list[0], end=locus_list[-2]+1)

        logging.info("hash implement")
        count = 0
        proper_count = 0
        pre_locus = 0
        pre_c = 0
        c = 0
        singletons_count = 0
        for read in iter1:
            count += 1
            '''
            if read.is_unmapped == False and read.mate_is_unmapped == True:
                singletons_count += 1
                continue
            '''
            if read.is_proper_pair:
                proper_count += 1
                connected = set()
                s = read.reference_start
                if read.is_read1 == True:
                    key = read.query_name + "_" + str(read.next_reference_start)
                else:
                    key = read.query_name + "_" + str(read.reference_start)

                if s > pre_locus:
                    c = b_search_while(locus_list, s)
                    pre_c = c
                    pre_locus = locus_list[c]
                    '''
                    cp1 = c + 1
                    if cp1 < len_locus and s < locus_list[cp1]:
                        pre_c = cp1
                        pre_locus = locus_list[cp1]
                        c = cp1
                    else:
                        c = b_search_while(locus_list, s)
                        pre_c = c
                        pre_locus = locus_list[c]
                    '''
                else:
                    c = pre_c

                region = trans_CigarTuples_To_Region(s, read.cigartuples)
                for s, e in region:
                    while c < len_locus and  locus_list[c] <= e:
                        connected.add(c)
                        c += 1

                if key in read_map:
                    uset = connected.union(read_map[key])
                    make_pair_connected(uset, locus_listed_dict)
                    del read_map[key]
                else:
                    read_map[key] = connected

        logging.info("end hash implement")
        break

    samfile.close()

    sys.exit()

    single_count = 0
    for k, v in read_map.items():
        single_count += 1


    connected_component = []
    
    walked = set()
    for i in range(len(locus_list) - 1):
        cc = []
        build_cc(locus_listed_dict, i, cc, walked)
        if len(cc) > 5:
            connected_component.append(cc)
            print ("=== " ,len(cc) , " ===")
            for c in cc:
                print(c, locus_Mm_list[c])

    print (count, proper_count, single_count, singletons_count)

    print('in blocks')
    
    return connected_component, locus_Mm_list, locus_idx_list
