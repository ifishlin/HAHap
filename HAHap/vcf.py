import logging

logger = logging.getLogger()


def output_CCs2VCF(connected_components, cc_allele_Mm, cc_idx_list, chrom=1):
    for cc in connected_components:
        vcf_file = open(str(cc[0]) + ".in.vcf", 'w')
        for locus in cc:
            vcf_file.write('\t'.join([str(chrom), cc_idx_list[locus], '.', '\t'.join(cc_allele_Mm[locus]),
                                      '.\t.\t.\t.\t.\t\n']))


def output_CC2VCF(connected_component, cc_allele_Mm, cc_idx_list, chrom=1):
    """
    """
    sorted_cc = sorted(connected_component)
    out_file = str(sorted_cc[0]) + ".in.vcf"
    vcf_file = open(out_file, 'w')
    for locus in sorted_cc:
        vcf_file.write('\t'.join([str(chrom), cc_idx_list[locus], '.', '\t'.join(cc_allele_Mm[locus]),
                                  '.\t.\t.\t.\t.\t\n']))
    
    vcf_file.close()
    return out_file


def output_phasing2VCF(input_vcf, output_file, output_dict, chrom, encoding_table):
    """
    """
    with open(output_file, "a") as myfile:
        input_file = open(input_vcf, 'r')
        for line in input_file:
            # line = line.strip()
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
            # myfile.write('\t'.join([chrom, c]))


def split_vcf_by_chrom(variant_file, indel=False):
    chrom = None
    pre = None
    var_chrom_dict = dict()

    vars_allele = []
    vars_str_loc = []

    variants_vcf = open(variant_file, 'r')
    for line in variants_vcf:
        if line[0] == '#':
            continue
        e = line.split("\t")
        c, i, _, m, n = e[0:5]

        if not indel and (len(m) > 1 or len(n) > 1):
            continue

        if c != chrom:
            if chrom is not None:
                var_chrom_dict[int(chrom)] = [vars_allele, vars_str_loc]

            vars_allele = []
            vars_str_loc = []

            vars_allele.append([m, n])
            vars_str_loc.append(i)

            pre = i
            chrom = c
            logger.info("chromosome " + chrom)
        else:
            if pre == i:
                logger.warning("duplicate location " + i)
                while vars_str_loc[-1] == pre:
                    logger.warning("DELETE")
                    del vars_allele[-1]
                    del vars_str_loc[-1]
                continue

            pre = i
            vars_allele.append([m, n])
            vars_str_loc.append(i)

    var_chrom_dict[int(chrom)] = [vars_allele, vars_str_loc]

    return var_chrom_dict
