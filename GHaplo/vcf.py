import logging
import sys
import collections
import os

logger = logging.getLogger()

def output_CCs2VCF(connected_components, cc_allele_Mm, cc_idx_list, chrom=1):
    for cc in connected_components:
       vcf_file = open(str(cc[0]) + ".in.vcf", 'w')
       for locus in cc:
           vcf_file.write('\t'.join([str(chrom), cc_idx_list[locus], '.', '\t'.join(cc_allele_Mm[locus]), '.\t.\t.\t.\t.\t\n']))

def output_CC2VCF(connected_component, cc_allele_Mm, cc_idx_list, chrom=1):
    """
    """
    sorted_cc = sorted(connected_component)
    out_file = str(sorted_cc[0]) + ".in.vcf"
    vcf_file = open(out_file, 'w')
    for locus in sorted_cc:
        vcf_file.write('\t'.join([str(chrom), cc_idx_list[locus], '.', '\t'.join(cc_allele_Mm[locus]), '.\t.\t.\t.\t.\t\n']))
    
    vcf_file.close()
    return out_file

def output_phasing2VCF(input_vcf, output_file, output_dict, chrom, encoding_table):
    """
    """
    with open(output_file, "a") as myfile:
        input_file = open(input_vcf, 'r')
        for line in input_file:
            #line = line.strip()
            e = line.split('\t')
            if chrom == e[0] and e[1] in output_dict and len(e[3]) < 2 and len(e[4]) < 2:
                ps_id, h1, h2 = output_dict[e[1]]
                ps_id = str(ps_id)
                h1 = '1' if encoding_table[h1] == e[3] else '2'
                h2 = '1' if encoding_table[h2] == e[3] else '2'
                e[8] += ':HP'
                e[9] = e[9].strip()
                e[9] += ''.join([':', ps_id, '-', h1, ',', ps_id, '-', h2, '\n'])
                myfile.write('\t'.join(e))
            elif chrom == e[0]:
                myfile.write(line)

        input_file.close()

def output_cc2csv(output_file, chrom, connected_components, cc_idx, block_min=2):
    """
    """
    chrom = str(chrom)
    with open(output_file, "a") as myfile:
        serial = 0
        for c in connected_components:
            if len(c) < block_min:
                continue
            serial += 1
            myfile.write("---- chromesome " + chrom + " ---- serial " + str(serial) + " ----" + '\n')
            for i in c:
                myfile.write('\t'.join([chrom, cc_idx[i], '\n']))
            #myfile.write('\t'.join([chrom, c]))

def split_vcf_by_chrom(variant_file, indel=False):
    current_chrom = None
    pre_location = None

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
            if current_chrom != None:
                variant_chrom_dict[int(current_chrom)] = [locus_list, locus_Mm_list, locus_idx_list]

            locus_list = []
            locus_Mm_list = []
            locus_idx_list = []

            locus_list.append(int(i))
            locus_Mm_list.append([m, n])
            locus_idx_list.append(i)

            pre_location = i
            current_chrom = c
            logger.info("chromosome "+ current_chrom)
        else:
            if pre_location == i:
                logger.warning("duplicate location " + i)
                while locus_list[-1] == int(pre_location):
                    logger.warning("DELETE")
                    del locus_list[-1]
                    del locus_Mm_list[-1]
                    del locus_idx_list[-1]
                continue

            pre_location = i
            locus_list.append(int(i))
            locus_Mm_list.append([m, n])
            locus_idx_list.append(i)

    variant_chrom_dict[int(current_chrom)] = [locus_list, locus_Mm_list, locus_idx_list]
    #print(locus_list)

    return variant_chrom_dict

