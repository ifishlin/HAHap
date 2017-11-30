# HAHap: A read-based haplotyping method using hierarchical assembly

About
---
HAHap is a method to infer the haplotype using sequencing data. HAHap is not designed to solve the MEC problem. Instead, it attempts to eliminate the influence of noises through the process of assembly, though it remains the spirit of minimum error correc-tion in certain conditions. We develop an adjusted multinomial probabilistic metric for evaluating the reliability of a variant pair, and the derived scores guide the assembly process.

The method is validated on DNA-seq short reads (Illumina HiSeq), and it supposed to handle all kinds of BAM files.


Required
---
HAHap is a pure-python program. It required fellowing packages. 

* Python 3.x
* Numpy (version >= 1.10)
* Pysam (version >= 0.12)


Usage
---
Git clone and execute **bin/HAHap**.
```
git clone https://github.com/ifishlin/HAHap
cd HAHap/bin
python HAHap phase vcf bam out
```

```
usage: python HAHap phase [--mms MMS] [--lct LCT] [--embed] [--lastops] VCF BAM OUT

positional arguments
VCF          VCF file with heterozygous variants needed to be phased
BAM          Read mapped file
OUT          VCF file with predicted haplotype. (HP tags)

optional arguments:
--mms            Minimum read mapping quality (default:0)
--lct            Threshold of low coverage pairs (int, default:median)
--embed_disable  Enable optimal search in embed case (default:enabled)
--last_disable   Enable optimal search in ambiguous case (default:enabled) 
        
```

Data (Ashkenazim family)
---
The answer set used in real data experiment which created by taking intersection between 

* the haplotype prediction of 10xGenomics (ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/10XGenomics_ChromiumGenome_LongRanger2.1_09302016/README) and 
* the variants calling results on the read sets (ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/novoalign_bams/README) (caller, GATK HaplotypeCaller 3.6). 
* The input vcf files of real data experiment are located at /data 
* The BAM of real data experiment is located at Base URL: ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/AshkenazimTrio
    * Base_URL/HG002_NA24385_son/NIST_Illumina_2x250bps/novoalign_bams/HG002.hs37d5.2x250.bam
    * Base_URL/HG003_NA24149_father/NIST_Illumina_2x250bps/novoalign_bams/HG003.hs37d5.2x250.bam
    * Base_URL/HG004_NA24143_mother/NIST_Illumina_2x250bps/novoalign_bams/HG004.hs37d5.2x250.bam



Authors
---
Yu-Yu Lin, Pei-Lung Chen, Yen-Jen Oyang and Chien-Yu Chen. 
National Taiwan University, Taiwan.
