#!/usr/bin/env python3
"""
util function
"""


def add_arguments(parser):
    arg = parser.add_argument
    arg('blocks_file', metavar='BLK', help='')
    arg('answer_file', metavar='ANS', help='')
    arg('blocks_ans', metavar='BLK_ANS', help='')


def main(args):
    get_phased_ans(args.blocks_file, args.answer_file, args.blocks_ans)


def get_phased_ans(blocks_file, answer_file, blocks_ans):
    """
    """
    ans = dict()
    a_file = open(answer_file, 'r')
    for line in a_file:
        line = line.strip()
        e = line.split('\t')
        keys = e[8].split(':')
        values = e[9].split(':')

        if "PS" in keys and "GT" in keys:
            ps_idx = keys.index("PS")
            ps_value = values[ps_idx]
            gt_idx = keys.index("GT")
            gt_value = values[gt_idx]   
        else:
            continue

        ans[e[0]+"_"+e[1]] = (gt_value.split("|"), ps_value)

    output_file = blocks_ans
    with open(output_file, "w") as myfile:
        b_file = open(blocks_file, 'r')
        for line in b_file:
            e = line.strip().split('\t')
            if len(e) == 2:
                key = e[0]+"_"+e[1]
                myfile.write('\t'.join([e[0], e[1], ans[key][1], ans[key][0][0], ans[key][0][1], '\n']))
            else:
                myfile.write(line)
