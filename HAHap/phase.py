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
    arg('--minj', dest='minj', default=4, type=int, help='Minimum junctions number (default:4)')
    arg('--pl', dest='pl', default=0.49, type=float, help='The likelihood of P1 and P2 (default:0.49)')
    #arg('--last_disable', dest='last_disable', action='store_false', help='Disable optimal search in ambiguous case.')


def main(args):
    if os.path.isfile(args.output_file):
        os.remove(args.output_file)

    logger.info("=== Start HAHap phasing ===")
    logger.info("Parameters: Minimum mapping quality = " + str(args.mms))
    logger.info("Parameters: Threshold of low coverage " + ("= Median" if args.lct == 0 else "= " + str(args.lct)))
    logger.info("Parameters: Minimum junction number = " + str(args.minj))
    logger.info("Parameters: Likelihood of P1 and P2 = " + str(args.pl))
    if args.pl*2 > 1:
        logger.info("Parameters incorrect, (P1 + P2) = " + str(args.pl*2) + " > 1, the sum of likelihood is less than or equal to one.")
        logger.info("Program exists")
        return

    #logger.info("Parameters: Embed optimal search ... " + ("Enable" if args.embed_disable == True else "Disable"))
    #logger.info("Parameters: Last optimal search ... " + ("Enable" if args.last_disable == True else "Disable"))


    timer = StageTimer()

    logger.info("")
    logger.info("=== Read Heterozygous Data ===")
    var_chrom_dict = split_vcf_by_chrom(args.variant_file)
    var_chrom_list = sorted(var_chrom_dict.items())

    for chrom, (var_allele, var_str_loc) in var_chrom_list:
        logger.info("")
        logger.info("=== Build Connected Component ===")
        var_loc = list(map(int, var_str_loc))
        connected_component = blocks_main(args, str(chrom), var_loc, timer)
        logger.info("Chromosome " + str(chrom))
        logger.info("Found variant block : " + str(len(connected_component)))
        logger.info("Phaing proceeding ... ")
        pipeline(args, str(chrom), connected_component, var_allele, var_str_loc, timer)
        logger.info("Phaing end and Output results")

    logger.info("")
    logger.info("=== End HAHap phasing ===")
    logger.info("")

def _timer_summary(timer):
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

        phase_allele = [(var_allele[c][0], var_allele[c][1]) for c in cc_sorted]
        phase_loc = [int(var_loc[c]) for c in cc_sorted]
        phase_total = len(phase_loc)

        #print(phase_loc)

        #  build matrix creater
        if reader is None:
            reader = InputVcfReader(args.bam_file, args.mms)
        #else:
        #    reader.reset(args.bam_file, 0, args.mms)

        # SNP-fragment matrix
        sf_mx, fragments, fragment_se, codes = reader.get_sf_mx(chrom, phase_loc, phase_allele, timer)
        if len(fragment_se) == 0:
            remove_dict[cc_idx] = None
            #print("No qualified fragment")
            continue

        pairs_sup = create_pairs_sup(fragment_se, phase_total, sf_mx)

        if len(pairs_sup) == 0:
            remove_dict[cc_idx] = None
            #print("No trusted allele")
            continue
        cs_mx = calc_cs_mx(pairs_sup, phase_loc, args.lct, args.pl)
        #print_var_matrix(sf_mx, fragment_se, codes, ':')

        vars_pool = dict()

        for i in range(len(phase_loc)):
            # vars idx and vars loc
            node = HiMergeNode(str(i), phase_loc[i])
            # node.name = str(i)
            vars_pool[node.name] = node

        ha_phasing(vars_pool, pairs_sup, cs_mx, phase_loc, phase_allele, codes, fragments, fragment_se, True, True, timer, args.minj)
        #ha_phasing(vars_pool, pairs_sup, cs_mx, phase_loc, phase_allele, codes, fragments, fragment_se, args.embed_disable, args.last_disable, timer)

        if len(vars_pool) > 1:
            od = collections.OrderedDict(sorted(vars_pool.items(), reverse=True))
            r_list = []
            for k in od.keys():
                r_list.append(list(map(int, k.split('_'))))
            remove_dict[cc_idx] = r_list

        save_phasing(vars_pool, phase_loc, output_dict)

    rod = collections.OrderedDict(sorted(remove_dict.items(), reverse=True))
    for k, value in rod.items():
        item = connected_component[k]
        if value is not None:
            for v in value:
                insert = [item[i] for i in v]        
                connected_component.insert(k+1, insert)
        del connected_component[k]

    output_phasing2VCF(args.variant_file, args.output_file, output_dict, chrom, codes)
