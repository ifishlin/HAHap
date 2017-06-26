import logging
import sys
#import vcf #PyVCF

def output_CCs2VCF(connected_components, cc_allele_Mm, cc_idx_list, chrom=1):
    for cc in connected_components:
       vcf_file = open(str(cc[0]) + ".in.vcf", 'w')
       for locus in cc:
           vcf_file.write('\t'.join([str(chrom), cc_idx_list[locus], '.', '\t'.join(cc_allele_Mm[locus]), '.\t.\t.\t.\t.\t\n']))

def output_CC2VCF(connected_component, cc_allele_Mm, cc_idx_list, chrom=1):
    sorted_cc = sorted(connected_component)
    out_file = str(sorted_cc[0]) + ".in.vcf"
    vcf_file = open(out_file, 'w')
    for locus in sorted_cc:
        vcf_file.write('\t'.join([str(chrom), cc_idx_list[locus], '.', '\t'.join(cc_allele_Mm[locus]), '.\t.\t.\t.\t.\t\n']))
    
    vcf_file.close()
    return out_file
