import logging

logger = logging.getLogger()

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
                logger.info("Found variants : " + str(len(vars_str_loc)) + "\n")
                var_chrom_dict[int(chrom)] = [vars_allele, vars_str_loc]
                logger.info("Reading Chromosome " + c)
            else:
                logger.info("Reading Chromosome " + c)

            vars_allele = []
            vars_str_loc = []

            vars_allele.append([m, n])
            vars_str_loc.append(i)

            pre = i
            chrom = c
        else:
            if pre == i:
                logger.warning("Duplicate location " + i + ", retain last record")
                while vars_str_loc[-1] == pre:
                    del vars_allele[-1]
                    del vars_str_loc[-1]
                continue

            pre = i
            vars_allele.append([m, n])
            vars_str_loc.append(i)

    logger.info("Chromosome " + chrom + ". Found variants : " + str(len(vars_str_loc)))
    var_chrom_dict[int(chrom)] = [vars_allele, vars_str_loc]

    return var_chrom_dict
