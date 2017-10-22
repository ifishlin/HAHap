#!/usr/bin/env python3
"""
Function phase
"""

import logging
import sys
import collections
import os
from HAHap.blocks import main as blocks_main
from HAHap.vcf import output_phasing2VCF, split_vcf_by_chrom
from HAHap.variants import InputVcfReader
from HAHap.entropy import create_pairs_sup, calc_cs_mx
from GHaplo.assembly import HiMergeNode, ha_phasing, save_phasing
from .timers import StageTimer
#from .matrix import print_var_matrix

sys.setrecursionlimit(30000)

logger = logging.getLogger(__name__)

def add_arguments(parser):
    arg = parser.add_argument
    arg('variant_file', metavar='VCF', help='VCF file with variants needed to be phased')
    arg('bam_file', metavar='BAM', help='BAM file')
    arg('output_file', metavar='OUT', help='Predicted file')
    # arg('cc_file', metavar='CC', help='CC file')
    arg('--mms', dest='mms', default=0, type=int, help='mininal mapping score')
    arg('--d', dest='distance', default=0, type=int, help='distance')
    arg('--lct', dest='lct', default=0, type=int, help='low coverage threshold')
    arg('--los_off', dest='los', action='store_false', help='off local optimal search')
    arg('--read4', dest='read4', action='store_true', help='read4')
    arg('--lastops', dest='lastops', action='store_true', help='lastops')


def main(args):
    if os.path.isfile(args.output_file):
        os.remove(args.output_file)

    # if os.path.isfile(args.cc_file):
    #     os.remove(args.cc_file)

    timer = StageTimer()

    # chrom: allele, str loc
    var_chrom_dict = split_vcf_by_chrom(args.variant_file)
    var_chrom_list = sorted(var_chrom_dict.items())

    print("mms", args.mms)
    print("distance", args.distance)
    print("lct", args.lct)
    print("los", args.los)
    print("read4", args.read4)
    print("lastops", args.lastops)

    for chrom, (var_allele, var_str_loc) in var_chrom_list:
        logging.info(chrom)
        logger.warning('a')
        logger.info('b')
        var_loc = list(map(int, var_str_loc))
        #print(chrom, len(var_loc), len(var_allele), len(var_str_loc))
        timer.start("00.blocks_main")
        connected_component, f, f_s, f_e = blocks_main(args, str(chrom), var_loc, timer, var_allele)
        timer.stop("00.blocks_main")
        print('connected_component,:', len(connected_component))
        #print(connected_component)
        f_s_name, f_s_idx = zip(*f_s)
        f_e_name, f_e_idx = zip(*f_e)
        pipeline(args, str(chrom), connected_component, var_allele, var_str_loc, timer, f, f_s_name, f_s_idx, f_e_name, f_e_idx)
        # output_cc2csv(args.cc_file, key, connected_component, var_loc)
    _timer_summary(timer)


def _timer_summary(timer):
    print("_timer_summary")
    total = timer.total()
    sorted_x = sorted(timer.tags())
    for k, v in sorted_x:
        print(k, v, round(v/total,3))


def pipeline(args, chrom, connected_component, var_allele, var_loc, timer, f, f_s_name, f_s_idx, f_e_name, f_e_idx):

    #  search vars blocks
    #  return, (1). connected_conmponents,
    #          (2). list of ref/alt,
    #          (3). list of genomic loc of varaint


    if len(connected_component) == 0:
        return

    output_dict = dict()
    remove_dict = dict()
    reader = None
    from bisect import bisect_left, bisect_right
    for cc_idx, cc_sorted in enumerate(connected_component):
        f_s_cc = cc_sorted[0]
        f_e_cc = cc_sorted[-1]

        timer.start('01.search')
        a = bisect_left(f_s_idx, f_s_cc)
        b = bisect_right(f_s_idx, f_e_cc)
        c = bisect_left(f_e_idx, f_s_cc)
        d = bisect_right(f_e_idx, f_e_cc) 
        timer.stop('01.search')
        timer.start('01.d')
        n = set(f_e_name[c:d]) | set(f_s_name[a:b])
        d = dict((k, f[k]) for k in n)
        timer.stop('01.d')
        #print(d)
        #print(len(d), len(n))

        #  output phase instance to file
        # cc_file = output_CC2VCF(cc_sorted, var_allele, var_loc)

        #  build phase instance from blocks search directly

        timer.start('01.prepare')
        phase_allele = [(var_allele[c][0], var_allele[c][1]) for c in cc_sorted]
        phase_loc = [int(var_loc[c]) for c in cc_sorted]
        phase_total = len(phase_loc)
        timer.stop('01.prepare')

        #print(phase_loc)

        #  build matrix creater
        timer.start('02.prepare2')
        if reader is None:
            reader = InputVcfReader(args.bam_file, 0, args.mms)
        #else:
        #    reader.reset(args.bam_file, 0, args.mms)
        timer.stop('02.prepare2')

        # SNP-fragment matrix
        timer.start('02.get_readmtx')
        sf_mx, fragments, fragment_se, codes = reader.get_sf_mx(chrom, phase_loc, phase_allele, timer, d)
        timer.stop('02.get_readmtx')
        if len(fragment_se) == 0:
            remove_dict[cc_idx] = None
            print("No qualified fragment")
            continue

        print(len(f))
        print(len(fragments))

        '''
        for x, y in f.items():
            if x not in fragments:
                print("NOT IN")
                pass
            elif y != fragments[x]:
                print("NO")
        '''

        timer.start('03.pv_dict')
        pairs_sup = create_pairs_sup(fragment_se, phase_total, sf_mx)

        if len(pairs_sup) == 0:
            remove_dict[cc_idx] = None
            print("No trusted allele")
            continue
        timer.stop('03.pv_dict')
        timer.start('04.calc_score_matrix')
        cs_mx = calc_cs_mx(pairs_sup, phase_loc, args.distance, args.lct)
        timer.stop('04.calc_score_matrix')
        # print_var_matrix(x, fragment_se, n, ':')

        vars_pool = dict()

        timer.start('05.create_pool')
        for i in range(len(phase_loc)):
            # vars idx and vars loc
            node = HiMergeNode(str(i), phase_loc[i])
            # node.name = str(i)
            vars_pool[node.name] = node
        timer.stop('05.create_pool')

        timer.start('06.ha_phasing')
        #print(args.los, args.lastops)
        ha_phasing(vars_pool, pairs_sup, cs_mx, phase_loc, phase_allele, codes, fragments, fragment_se, args.read4, args.los, args.lastops, timer)
        #ha_phasing(vars_pool, pairs_sup, cs_mx, phase_loc, phase_allele, codes, timer, args.los, fragments, fragment_se, args.read4, args.lastops)
        timer.stop('06.ha_phasing')

        if len(vars_pool) > 1:
            timer.start('07.check_null')
            od = collections.OrderedDict(sorted(vars_pool.items(), reverse=True))
            r_list = []
            for k in od.keys():
                r_list.append(list(map(int, k.split('_'))))
            remove_dict[cc_idx] = r_list
            timer.stop('07.check_null')

        timer.start('08.output')
        save_phasing(vars_pool, phase_loc, output_dict)
        timer.stop('08.output')

    timer.start('09.check_null2')
    rod = collections.OrderedDict(sorted(remove_dict.items(), reverse=True))
    for k, value in rod.items():
        item = connected_component[k]
        if value is not None:
            for v in value:
                insert = [item[i] for i in v]        
                connected_component.insert(k+1, insert)
        del connected_component[k]
    timer.stop('09.check_null2')

    timer.start('08.output')
    output_phasing2VCF(args.variant_file, args.output_file, output_dict, chrom, codes)
    timer.stop('08.output')
