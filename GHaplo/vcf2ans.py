#!/usr/bin/env python3
"""
"""

def add_arguments(parser):
    arg = parser.add_argument
    arg('phased_vcf', metavar='VCF', help='')
    arg('ans_csv', metavar='OUT', help='')
    arg('ans_tag', metavar='TYPE', default='HP', help='')

def main(args):
    output_file = args.ans_csv
    with open(output_file, "w") as myfile:
        vcf_file = open(args.phased_vcf)
        for line in vcf_file:
            if line[0] == '#':
                continue

            line = line.strip()
            e = line.split('\t')
            keys   = e[8].split(':')
            values = e[9].split(':')

            if "HP" in keys and args.ans_tag == "HP":
                hp_idx = keys.index("HP")
                hp_value = values[hp_idx]
                hp1, hp2 = hp_value.split(",")
                ps_id = hp1.split('-')[0]
                h1 = int(hp1.split('-')[1])-1
                h2 = int(hp2.split('-')[1])-1
                myfile.write('\t'.join([e[0], e[1], ps_id, str(h1), str(h2), '\n']))

            elif "PS" in keys and "GT" in keys and args.ans_tag == "PS":
                try:
                    gt_idx = keys.index("GT")
                    g1, g2 = values[gt_idx].split('|')
                    ps_idx = keys.index("PS")
                    ps_id = values[ps_idx]
                    myfile.write('\t'.join([e[0], e[1], ps_id, str(g1), str(g2), '\n']))
                except ValueError:
                    pass
