#!/usr/bin/env python3
"""
Function phase
"""

import logging
import sys
import collections
import os
from GHaplo.blocks import main as blocks_main
from GHaplo.vcf import output_CCs2VCF, output_CC2VCF, output_phasing2VCF, split_vcf_by_chrom, output_cc2csv
from GHaplo.variants import InputVcfReader
from GHaplo.entropy import create_sv_freq_dict, create_pv_freq_dict, calc_score_matrix
from GHaplo.assembly import HiMergeNode, hc_merge, rebuild_haplo, save_phasing
from .timers import StageTimer
from .matrix import print_var_matrix

sys.setrecursionlimit(30000)

logger = logging.getLogger(__name__)

def add_arguments(parser):
    arg = parser.add_argument
    arg('variant_file', metavar='VCF', help='VCF file with variants needed to be phased')
    arg('input_file', metavar='BAM', help='BAM file')
    arg('output_file', metavar='OUT', help='Predicted file')
    arg('cc_file', metavar='CC', help='CC file')
    arg('--mms', dest='mms', default=0, type=int, help='mininal mapping score')

def main(args):
    if os.path.isfile(args.output_file):
        os.remove(args.output_file)

    if os.path.isfile(args.cc_file):
        os.remove(args.cc_file)

    timer = StageTimer()

    variant_chrom_dict = split_vcf_by_chrom(args.variant_file)
    variant_chrom_dict = collections.OrderedDict(sorted(variant_chrom_dict.items()))    

    for key, (v1, v2, v3) in variant_chrom_dict.items():
        logging.info(key)
        print(key, len(v1), len(v2), len(v3))
        #timer.start('blocks')
        connected_components, ref_alt_list, genomic_location_list = blocks_main(args, str(key), v1, v2, v3, timer)
        #timer.stop('blocks')
        #output_cc2csv(args.cc_file, key, connected_components, genomic_location_list)
        #print(connected_components)
        print('connected_components,:', len(connected_components))
        pipeline(args, str(key), connected_components, ref_alt_list, genomic_location_list, timer)
        output_cc2csv(args.cc_file, key, connected_components, genomic_location_list)

    print("timer.pre",timer.elapsed('pre'), round(timer.elapsed('pre')/timer.total(),3))
    print("timer.search",timer.elapsed('search'), round(timer.elapsed('search')/timer.total(),3))
    print("timer.putin",timer.elapsed('putin'), round(timer.elapsed('putin')/timer.total(),3))
    print("timer.walk",timer.elapsed('walk'), round(timer.elapsed('walk')/timer.total(),3))
    print("timer.prepare",timer.elapsed('prepare') ,round(timer.elapsed('prepare')/timer.total(),3))
    print("timer.prepare2",timer.elapsed('prepare2') ,round(timer.elapsed('prepare2')/timer.total(),3))
    print("timer.readmtx_init",timer.elapsed('readmtx_init'),round(timer.elapsed('readmtx_init')/timer.total(),3))
    #print("timer.read_variance",timer.elapsed('read_variance'),round(timer.elapsed('read_variance')/timer.total(),3))
    print("timer.fetch",timer.elapsed('fetch'),round(timer.elapsed('fetch')/timer.total(),3))
    print("timer.aligned_pairs",timer.elapsed('aligned_pairs'),round(timer.elapsed('aligned_pairs')/timer.total(),3))
    print("timer.phase_fragments_dict",timer.elapsed('phase_fragments_dict'),round(timer.elapsed('phase_fragments_dict')/timer.total(),3))
    print("timer.removed_set",timer.elapsed('removed_set'),round(timer.elapsed('removed_set')/timer.total(),3))

    print("timer.clear_matrix",timer.elapsed('clear_matrix'),round(timer.elapsed('clear_matrix')/timer.total(),3))
    print("timer.produce_var_matrix'",timer.elapsed('produce_var_matrix'),round(timer.elapsed('produce_var_matrix')/timer.total(),3))
        
    print("timer.sv_dict",timer.elapsed('sv_dict'),round(timer.elapsed('sv_dict')/timer.total(),3))
    print("timer.pv_dict",timer.elapsed('pv_dict'),round(timer.elapsed('pv_dict')/timer.total(),3))
    print("timer.calc_score_matrix",timer.elapsed('calc_score_matrix'),round(timer.elapsed('calc_score_matrix')/timer.total(),3))
    print("timer.create_pool",timer.elapsed('create_pool'),round(timer.elapsed('create_pool')/timer.total(),3))
    #print("timer.hc_merge",timer.elapsed('hc_merge'),round(timer.elapsed('hc_merge')/timer.total(),3)

    print("timer.copy_matrix_1",timer.elapsed('copy_matrix_1'),round(timer.elapsed('copy_matrix_1')/timer.total(),3))
    print("timer.get_maxscore",timer.elapsed('get_maxscore'),round(timer.elapsed('get_maxscore')/timer.total(),3))
    print("timer.get_maxscore_between_segments",timer.elapsed('get_maxscore_between_segments'),round(timer.elapsed('get_maxscore_between_segments')/timer.total(),3))
    print("timer.choice_haplos",timer.elapsed('choice_haplos'),round(timer.elapsed('choice_haplos')/timer.total(),3))
    print("timer.rearrange",timer.elapsed('rearrange'),round(timer.elapsed('rearrange')/timer.total(),3))
    print("timer.reduce_matrix",timer.elapsed('reduce_matrix'),round(timer.elapsed('reduce_matrix')/timer.total(),3))



    print("timer.check_null1",timer.elapsed('check_null'),round(timer.elapsed('check_null')/timer.total(),3))
    print("timer.check_null2",timer.elapsed('check_null2'),round(timer.elapsed('check_null2')/timer.total(),3))
    print("timer.output",timer.elapsed('output'),round(timer.elapsed('output')/timer.total(),3))

def pipeline(args, chrom, connected_components, ref_alt_list, genomic_location_list, timer):

    ## search variants blocks
    ## return, (1). connected_conmponents, 
    ##         (2). list of ref/alt, 
    ##         (3). list of genomic location of varaint

    if len(connected_components) == 0:
        return

    output_dict = dict()

    reader = None
    remove_dict = dict()

    for cc_idx, cc in enumerate(connected_components):
        ## remove sorted?
        cc_sorted = sorted(cc)

        cc_file = ''
        ## output phase instance to file 
        #cc_file = output_CC2VCF(cc_sorted, ref_alt_list, genomic_location_list)

        ## build phase instance from blocks search directly
        phase_input = []
        timer.start('prepare')
        for c in cc_sorted:
            phase_input.append([chrom, int(genomic_location_list[c]), 0, 'snp', ref_alt_list[c][0], ref_alt_list[c][0], ref_alt_list[c][1]]) 

        timer.stop('prepare')
        ## build matrix creater
        timer.start('prepare2')
        if reader == None:
            reader = InputVcfReader(cc_file, args.input_file, 0, args.mms, cc_input=phase_input)
        else:
            reader.reset(cc_file, args.input_file, 0, args.mms, cc_input=phase_input)
        timer.stop('prepare2')

        ## read_matrix, phase_fragments_dict, phase_fragments_range_list, phase_variants_location_list, 
        ## InputVcfReader.encoding_table, phase_variants_called_list
        #timer.start('get_readmtx')
        x, y, z, m, n, l = reader.get_readmtx(timer)
        #timer.stop('get_readmtx')
        if len(z) == 0:
            remove_dict[cc_idx] = None
            print("fragmentOUT")
            continue
        timer.start('sv_dict')
        a = create_sv_freq_dict(z, m, x)
        timer.stop('sv_dict')
        timer.start('pv_dict')
        b = create_pv_freq_dict(z, m, x)
        if len(b) == 0:
            remove_dict[cc_idx] = None
            print("pv OUT")
            continue
        timer.stop('pv_dict')
        timer.start('calc_score_matrix')
        c = calc_score_matrix(b, l, n, 1)
        timer.stop('calc_score_matrix')
        #print_var_matrix(x, z, n, ':')

        vars_candidates_mapping = list()
        vars_candidates_pool = dict()

        timer.start('create_pool')
        for i in range(len(m)):
            vars_candidates_mapping.append(str(i))
            ## variants idx and variants location
            node = HiMergeNode(str(i), m[i])
            ## node.name = str(i)
            vars_candidates_pool[node.name] = node
        timer.stop('create_pool')

        #timer.start('hc_merge')
        r, s, t, u = hc_merge(vars_candidates_mapping, vars_candidates_pool, b, c, a, m, l, n, 'False', 0.9, timer)
        #timer.stop('hc_merge')

        if len(vars_candidates_pool) > 1:
            timer.start('check_null')
            od = collections.OrderedDict(sorted(vars_candidates_pool.items(), reverse=True))
            r_list = []
            for k in od.keys():
                r_list.append(list(map(int, k.split('_'))))
            remove_dict[cc_idx] = r_list
            #print(vars_candidates_pool)
            #print("== ", str(genomic_location_list[cc_sorted[0]]), str(genomic_location_list[cc_sorted[-1]]))
            timer.stop('check_null')

        #rebuild_haplo(chrom, vars_candidates_pool, n, m, y, a)
        timer.start('output')
        phasing_dict = save_phasing(chrom, vars_candidates_pool, n, m, output_dict)
        timer.stop('output')

    timer.start('check_null2')
    rod = collections.OrderedDict(sorted(remove_dict.items(), reverse=True))
    for k,value in rod.items():
        item = connected_components[k]
        if value != None:
            for v in value:
                insert = [item[i] for i in v]        
                connected_components.insert(k+1,insert)
        del connected_components[k]
    timer.stop('check_null2')

    output_file = args.output_file
    timer.start('output')
    output_phasing2VCF(args.variant_file, output_file, output_dict, chrom, n)
    timer.stop('output')
