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
from HAHap.confident import create_pairs_sup, calc_cs_mx
from HAHap.assembly import HiMergeNode, ha_phasing, save_phasing
from .timers import StageTimer
from .matrix import print_var_matrix

sys.setrecursionlimit(30000)

logger = logging.getLogger(__name__)

def add_arguments(parser):
    arg = parser.add_argument
    arg('variant_file', metavar='VCF', help='VCF file with variants needed to be phased')
    arg('bam_file', metavar='BAM', help='Read mapped file')
    arg('output_file', metavar='OUT', help='VCF file with predicted haplotype. (HP tags)')
    arg('--mms', dest='mms', default=0, type=int, help='Minimum read mapping quality (default:0)')
    arg('--lct', dest='lct', default=0, type=int, help='Threshold of low coverage pairs (default:median)')
    arg('--embed_disable', dest='embed_disable', action='store_false', help='Disable optimal search in embed case.')
    arg('--last_disable', dest='last_disable', action='store_false', help='Disable optimal search in ambiguous case.')


def main(args):
    if os.path.isfile(args.output_file):
        os.remove(args.output_file)

    logger.info("=== Start HAHap phasing ===")
    logger.info("Parameters: Minimum mapping quality = " + str(args.mms))
    logger.info("Parameters: Threshold of low coverage " + ("... Median" if args.lct == 0 else "= " + str(args.lct)))
    logger.info("Parameters: Embed optimal search ... " + ("Enable" if args.embed_disable == True else "Disable"))
    logger.info("Parameters: Last optimal search ... " + ("Enable" if args.last_disable == True else "Disable"))

    timer = StageTimer()

    logger.info("")
    logger.info("=== Read Heterozygous Data ===")
    var_chrom_dict = split_vcf_by_chrom(args.variant_file)
    var_chrom_list = sorted(var_chrom_dict.items())

    for chrom, (var_allele, var_str_loc) in var_chrom_list:
        logger.info("")
        logger.info("=== Build Connected Component ===")
        var_loc = list(map(int, var_str_loc))
        timer.start("00.blocks_main")
        connected_component = blocks_main(args, str(chrom), var_loc, timer)
        timer.stop("00.blocks_main")
        logger.info("Chromosome " + str(chrom))
        logger.info("Found variant block : " + str(len(connected_component)))
        logger.info("Phaing proceeding ... ")
        pipeline(args, str(chrom), connected_component, var_allele, var_str_loc, timer)
        logger.info("Phaing end and Output results")
    #_timer_summary(timer)

    logger.info("")
    logger.info("=== End HAHap phasing ===")
    logger.info("")
    _timer_summary(timer)

def _timer_summary(timer):
    print("_timer_summary")
    total = timer.total()
    sorted_x = sorted(timer.tags())
    for k, v in sorted_x:
        print(k, v, round(v/total,3))


def pipeline(args, chrom, connected_component, var_allele, var_loc, timer):

    #  search vars blocks
    #  return, (1). connected_conmponents,
    #          (2). list of ref/alt,
    #          (3). list of genomic loc of varaint


    if len(connected_component) == 0:
        return

    output_dict = dict()
    remove_dict = dict()
    reader = None

    for cc_idx, cc_sorted in enumerate(connected_component):
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
            reader = InputVcfReader(args.bam_file, args.mms)
        #else:
        #    reader.reset(args.bam_file, 0, args.mms)
        timer.stop('02.prepare2')

        # SNP-fragment matrix
        timer.start('02.get_readmtx')
        sf_mx, fragments, fragment_se, codes = reader.get_sf_mx(chrom, phase_loc, phase_allele, timer)
        timer.stop('02.get_readmtx')
        if len(fragment_se) == 0:
            remove_dict[cc_idx] = None
            #print("No qualified fragment")
            continue

        timer.start('03.pv_dict')
        pairs_sup = create_pairs_sup(fragment_se, phase_total, sf_mx)

        if len(pairs_sup) == 0:
            remove_dict[cc_idx] = None
            #print("No trusted allele")
            continue
        timer.stop('03.pv_dict')
        timer.start('04.calc_score_matrix')
        cs_mx = calc_cs_mx(pairs_sup, phase_loc, args.lct)
        timer.stop('04.calc_score_matrix')
        #print_var_matrix(sf_mx, fragment_se, codes, ':')

        vars_pool = dict()

        timer.start('05.create_pool')
        for i in range(len(phase_loc)):
            # vars idx and vars loc
            node = HiMergeNode(str(i), phase_loc[i])
            # node.name = str(i)
            vars_pool[node.name] = node
        timer.stop('05.create_pool')

        timer.start('06.ha_phasing')
        #print(args.embed_disable, args.last_disable)
        ha_phasing(vars_pool, pairs_sup, cs_mx, phase_loc, phase_allele, codes, fragments, fragment_se, args.embed_disable, args.last_disable, timer)
        #ha_phasing(vars_pool, pairs_sup, cs_mx, phase_loc, phase_allele, codes, timer, args.embed_disable, fragments, fragment_se, args.read4, args.last_disable)
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
