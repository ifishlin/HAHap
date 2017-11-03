# HAHap: A read-based haplotyping method using hierarchical assembly

About
---
HAHap is a method to infer the haplotype using sequencing data. This method does not focus on solving MEC problem. Instead, it attempts to eliminate the error on the process of assembly, though it remains the spirit of minimum error correction in certain conditions. We develop an adjusted multinomial probabilistic metric for evaluating the reliability of two variant pair, and the score guides the assembly process.  

The method is validated on NGS short reads (Illumina HiSeq), and it supposed to handle all kinds of BAM files.


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
usage: python HAHap phase [--mms MMS] [--lct LCT] [--embed] [--lastops] VCF BAM OUT

positional arguments
VCF          VCF file with heterozygous variants needed to be phased
BAM          Read mapped file
OUT          VCF file with predicted haplotype. (HP tags)

optional arguments:
--mms            Minimum read mapping quality (default:0)
--lct            Threshold of low coverage pairs (default:median)
--embed_disable  Enable optimal search in embed case (default:enabled)
--last_disable   Enable optimal search in ambiguous case (default:enabled) 
        
```

Data (Ashkenazim family)
---
The answer set used in real data experiment. It created by the intersection between 

* the [haplotype prediction of 10xGenomics](https://ftp-trace.ncbi.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/10XGenomics_ChromiumGenome_LongRanger2.1_09302016/README) and 
* the [variants calling result](https://ftp-trace.ncbi.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/novoalign_bams/README) by GATK HaplotypeCaller 3.6. 



Authors
---
Yu-Yu Lin, Yin-Hung Lin, Pei-Lung Chen, Yen-Jen Oyang and Chien-Yu Chen. 
National Taiwan University, Taiwan.
