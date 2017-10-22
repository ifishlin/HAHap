#!/home/fish/anaconda3/bin/python
import random
import sys

nucletide = ['A','C','G','T','a','c','g','t']

def randomSNP(args):
    fasta_file = args.fasta
    fasta_file_f = fasta_file
    if isinstance(fasta_file, str):
        fasta_file_f = open(fasta_file, 'r')

    h1_fasta_f = open(args.prefix + '_1.fasta','w')
    h2_fasta_f = open(args.prefix + '_2.fasta','w')

    snps = genSNPs(args.min, args.max)
    print(args.min, args.max)
    #snps = [[10,], [20,], [15,], [16,], [50,], [32,], [300,]]
    #print(snps)

    cidx = 0
    snps_idx = 0
    snps_next = snps[snps_idx][0]
    next_flag = False
    for line in fasta_file_f:
        line = line.strip()
        if(line[0] == '>'):
            h1_fasta_f.write(line + '\n')
            h2_fasta_f.write(line + '\n')
            print(line)
            continue
 
        line_len = len(line)
        #print(snps_next, line_len)
        if(snps_next <= line_len):
            import copy
            line1 = copy.deepcopy(line)
            line2 = copy.deepcopy(line)

            bidx = 0
            while snps_next <= (line_len - bidx):
                cidx = bidx + snps_next
                l1_pre  = line1[:cidx - 1] if snps_next > 0 else ''
                l1_post = line1[cidx:] if snps_next > 0 else line1[1:]
                l2_pre  = line2[:cidx - 1] if snps_next > 0 else ''
                l2_post = line2[cidx:] if snps_next > 0 else line2[1:]
                snp  = makeSNP(line[cidx - 1])

                dice = random.randint(0,1)
                #[pos, ref, alt, h1, h2]
                if(line[cidx - 1] in nucletide):
                    #dict = 0 => h1 is mutated,
                    #dict = 1 => h2 is mutated
                    if(dice == 0): 
                        snps[snps_idx].append(line[cidx - 1].upper()) #ref
                        snps[snps_idx].append(snp)                    #alt
                        snps[snps_idx].append(snp)                    #mutated 
                        snps[snps_idx].append(line[cidx - 1].upper()) #not mutated
                        line1 = ''.join([l1_pre, snp, l1_post])
                        line2 = ''.join([l2_pre, line[cidx - 1], l2_post])
                    else:
                        snps[snps_idx].append(line[cidx - 1].upper())
                        snps[snps_idx].append(snp) 
                        snps[snps_idx].append(line[cidx - 1].upper())
                        snps[snps_idx].append(snp)
                        line2 = ''.join([l2_pre, snp, l2_post])
                        line1 = ''.join([l1_pre, line[cidx - 1], l1_post])
                else:
                    snps[snps_idx].append(line[cidx - 1])
                    snps[snps_idx].append(line[cidx - 1])
                    snps[snps_idx].append(line[cidx - 1])
                    snps[snps_idx].append(line[cidx - 1])
                    line1 = ''.join([l1_pre, line[cidx - 1], l1_post])
                    line2 = ''.join([l2_pre, line[cidx - 1], l2_post])
                    
                print(snps[snps_idx])

                bidx = cidx
                snps_idx = snps_idx + 1
                if(snps_idx >= len(snps)):
                    tmp = genSNPs(args.min, args.max)
                    snps = snps + tmp
                    #break
                snps_next = snps[snps_idx][0]
           
            snps_next = snps_next - (line_len - bidx)
        else:
            line1 = line
            line2 = line
            snps_next = snps_next - line_len 


        if(snps_next == 0):
            sys.exit()

        h1_fasta_f.write(line1 + '\n')
        h2_fasta_f.write(line2 + '\n')

    snps = snps[:snps_idx]
    #print(args.rgid, snps)

    make_whatshap_input(args.prefix, args.rgid, snps, args.chrid)
    #make_fishhap_input(args.prefix, snps)
    h1_fasta_f.close()
    h2_fasta_f.close()
 
def genSNPs(bmin, bmax):
    #distance of next snps, (>0)
    snps = []
    for x in range(10000):
        i = random.randint(bmin, bmax)
        #i = random.randint(1000,1500)
        snps.append([i,])
    return snps

def makeSNP(nucle):
    #idx = nucletide.index(nucle.upper()) + 1
    #idx = idx % len(nucletide)
    #return nucletide[idx]
    idx = random.randint(0,3)
    while nucletide[idx] == nucle.upper():
        idx = random.randint(0,3)
    return nucletide[idx] 

def make_whatshap_input(prefix, rgid, snps, chrid):
    whatshap_f = open(prefix + '_whatshap.vcf','w')
    whatshap_f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + rgid + '\n')
    whatshap2_f = open(prefix + '_whatshap_notATCG.vcf','w')
    whatshap_ans_f = open(prefix + '_ans.vcf','w')

    cs_pos = 0
    s_pos = 0
    e_pos = 0
    max_range = 0
    pos = 0

    bid = None
    for idx, x in enumerate(snps):
        pos = pos + x[0]
        if(x[1] in nucletide):
            if(cs_pos == 0):
                cs_pos = pos
            gt_tag = '0|1' if x[1] == x[3] else '1|0'
            whatshap_f.write('\t'.join([chrid, str(pos), 'snp_' + str(idx+1), x[1], x[2], '100', '.', '.', 'GT:PS', gt_tag + ':1\n']))
 
            if bid == None:
                bid = str(pos)

            if gt_tag == '0|1': 
                whatshap_ans_f.write('\t'.join([chrid, bid, str(pos), x[1], x[2], '.', x[1], x[2], '\n']))
            else:
                whatshap_ans_f.write('\t'.join([chrid, bid, str(pos), x[1], x[2], '.', x[2], x[1], '\n']))
        else:
            if(cs_pos != 0):
                if(pos - cs_pos) > max_range:
                    max_range = pos - cs_pos
                    s_pos = cs_pos
                    e_pos = pos - x[0]
                cs_pos = 0 

            gt_tag = '0|1' if x[1] == x[3] else '1|0'
            whatshap2_f.write('\t'.join([chrid, str(pos), 'snp_' + str(idx+1), x[1], x[2], '100', '.', '.', 'GT:PS', gt_tag + ':1\n']))         
    whatshap_f.close()
    whatshap2_f.close()
    whatshap_ans_f.close()
    print(s_pos, e_pos)

'''
def make_fishhap_input(prefix, snps):
    fishhap_f = open(prefix + '_fishhap.vcf','w')
    fishhap_f.write('#fishhap\n')
    pos = 0
    for idx, x in enumerate(snps):
        pos = pos + x[0]
        if(x[1] in nucletide):
            fishhap_f.write('\t'.join(['M', '1', '.', str(pos), x[1], x[2], '0', x[3], x[4]]))
            fishhap_f.write('\n') 
    fishhap_f.close()
'''
   

if __name__ == '__main__':
    import sys
    import argparse
    from textwrap import dedent

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="print test")
    subparsers = parser.add_subparsers(title='subcommands', description='valid subcommands', help='additional help')

    parser_a = subparsers.add_parser('randomSNP', help='')
    parser_a.add_argument('fasta', type=str, help='input fasta')
    parser_a.add_argument('prefix', type=str, help='output prefix')
    parser_a.add_argument('rgid', type=str, help='@RG id')
    parser_a.add_argument('chrid', type=str, help='chrom id')
    parser_a.add_argument('min', type=int)
    parser_a.add_argument('max',type=int)
    parser_a.set_defaults(func=randomSNP)

    args = parser.parse_args()
    args.func(args)
