#!/usr/bin/env python3
"""
Function phase


"""
import logging
import sys
from GHaplo.blocks import main as blocks_main
from GHaplo.vcf import output_CCs2VCF, output_CC2VCF
from GHaplo.variants import InputVcfReader
from GHaplo.entropy import create_sv_freq_dict, create_pv_freq_dict, calc_score_matrix
from GHaplo.assembly import HiMergeNode, hc_merge, rebuild_haplo

def print_var_matrix(var_matrix, v_matrix_ids, encoding_table, prange):
    """
    """

    max_list = []

    #max len in each variance
    if(var_matrix.shape[0] == 0):
        sys.exit()
    for j in range(var_matrix.shape[1]):
        max_var = 0
        for i in range(var_matrix.shape[0]):
            max_var = len(encoding_table[var_matrix[i,j]]) if len(encoding_table[var_matrix[i,j]]) > max_var else max_var
        max_list.append(max_var)

    head = '{0:>48s}-'.format('read_id \\ site idx')
    sys.stdout.write(head)
    idx = 1
    for index, j in enumerate(max_list):
        sys.stdout.write(addpadding(str(idx), j)+' ')
        idx = idx + 1
        idx = idx % 10
    print()

   #print matrix
    for i in range(var_matrix.shape[0]):
        head = '{0:>48s}-'.format(v_matrix_ids[i][0])
        sys.stdout.write(head)
        if(prange == ':'):
            for j in range(var_matrix.shape[1]):
                sys.stdout.write(addpadding(encoding_table[var_matrix[i,j]], max_list[j]) + ',')
        else:
            start, end = map(int, prange.split(":"))
            for j in range(start, end):
                sys.stdout.write(addpadding(encoding_table[var_matrix[i,j]], max_list[j]) + ',')
        print()

    return

def addpadding(seq, l):
    if(l > len(seq)):
        for i in range(l-len(seq)):
            seq = seq + ' '

    return seq


def split_vcf_by_chrom(variant_file, indel=False):

    current_chrom = None

    variant_chrom_dict = dict()

    locus_list = []
    locus_Mm_list = []
    locus_idx_list = []

    variants_vcf = open(variant_file,'r')
    for line in variants_vcf:
        if line[0] == '#':
            continue
        e = line.split("\t")
        c, i, _, m, n = e[0:5]

        if not indel and (len(m) > 1 or len(n) > 1):
            continue    

        if c != current_chrom:
            variant_chrom_dict[current_chrom] = [locus_list, locus_Mm_list, locus_idx_list]

            locus_list = []
            locus_Mm_list = []
            locus_idx_list = []

            current_chrom = c
            print(current_chrom)
            pass
        else:
            locus_list.append(int(i))
            locus_Mm_list.append([m, n])
            locus_idx_list.append(i)

    del variant_chrom_dict[None]

    return variant_chrom_dict
    ## miss last line

def add_arguments(parser):
    arg = parser.add_argument
    arg('variant_file', metavar='VCF', help='VCF file with variants needed to be phased')
    arg('input_file', metavar='BAM', help='BAM file')
    arg('--mms', dest='mms', default=0, type=int, help='mininal mapping score')

def main(args):
    variant_chrom_dict = split_vcf_by_chrom(args.variant_file)
    for key, (v1, v2, v3) in variant_chrom_dict.items():
        print(key, len(v1), len(v2), len(v3))
        connected_components, ref_alt_list, genomic_location_list = blocks_main(args, key, v1, v2, v3)
        print('connected_components,:', len(connected_components))
        main2(args, key, connected_components, ref_alt_list, genomic_location_list)

def main2(args, chrom, connected_components, ref_alt_list, genomic_location_list):

    ## search variants blocks
    ## return, (1). connected_conmponents, 
    ##         (2). list of ref/alt, 
    ##         (3). list of genomic location of varaint

    for cc in connected_components:
        ## remove sorted?
        cc_sorted = sorted(cc)
  
        ## output phase instance to file 
        cc_file = output_CC2VCF(cc_sorted, ref_alt_list, genomic_location_list)

        print("ref_alt_list")
        print(ref_alt_list)

        ## build phase instance from blocks search directly
        phase_input = []
        for c in cc_sorted:
            phase_input.append([chrom, int(genomic_location_list[c]), 0, 'snp', ref_alt_list[c][0], ref_alt_list[c][0], ref_alt_list[c][1]]) 

        ## build matrix creater
        reader = InputVcfReader(cc_file, args.input_file, 0, args.mms, cc_input=phase_input)

        ## read_matrix, phase_fragments_dict, phase_fragments_range_list, phase_variants_location_list, 
        ## InputVcfReader.encoding_table, phase_variants_called_list
        x, y, z, m, n, l = reader.get_readmtx()
        a = create_sv_freq_dict(z, m, x)
        b = create_pv_freq_dict(z, m, x)
        c = calc_score_matrix(b, l, n, 1)
        print_var_matrix(x, z, n, ':')

        vars_candidates_mapping = list()
        vars_candidates_pool = dict()

        print("start assembly")
        for i in range(len(m)):
            vars_candidates_mapping.append(str(i))
            ## variants idx and variants location
            node = HiMergeNode(str(i), m[i])
            ## node.name = str(i)
            vars_candidates_pool[node.name] = node

        print(vars_candidates_mapping)
        print(vars_candidates_pool)       

        r, s, t, u = hc_merge(vars_candidates_mapping, vars_candidates_pool, b, c, a, m, l, n, 'False', 0.9)
        rebuild_haplo(vars_candidates_pool, n, m, y, a)

    '''
    output_CCs2VCF(connected_component, ref_alt_list, genomic_idx_list)

    reader = InputVcfReader('1007.in.vcf', args.input_file, 0, 60)
    a, b, c ,d ,e ,f = reader.get_readmtx()
    #### read input_VCFs ####
    
    print_var_matrix(a, c, e, ':')
    '''
    print('in phase')
    pass
