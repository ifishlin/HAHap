#!/usr/bin/env python3
"""
Build SNP-fragment matrix


"""
import logging
import sys
from GHaplo.blocks import main as blocks_main

def add_arguments(parser):
    arg = parser.add_argument
    arg('variant_file', metavar='VCF', help='VCF file with variants needed to be phased')
    arg('input_file', metavar='BAM', help='BAM file')

def main(args):
    print(args)
    print(args.variant_file, args.input_file)
    blocks_main(args)
    print('in matirx')
    pass

