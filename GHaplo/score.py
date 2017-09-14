#!/usr/bin/env python3

import os
import sys
from GHaplo.timers import StageTimer
from six.moves import cPickle as pickle

timer = StageTimer()


def check(block, predicted):
    # align first nucleotide

    b1, b2, b3, b4 = block[0]
    p1, p2, p3, p4 = predicted[0]

    first_alinged_switch = False

    if b3 == p3 and b4 == p4:
        pass
    elif b3 == p4 and b4 == p3:
        first_alinged_switch = True
        timer.start("align")
        # block = list(map(lambda b:(b[0],b[1],b[3],b[2]),block))
        timer.stop("align")
    else:
        print("==ERROR==") 
        sys.exit()

    last_is_identical = True
    last_time_error = False
    error = 0
    unphased = 0
    case = 0
    total = 0

    timer.start("loop")

    b_next = 0
    for pidx, (p1, p2, p3, p4) in enumerate(predicted):
        is_find = False
        if pidx != 0:
            total += 1
        for bidx, (b1, b2, b3, b4) in enumerate(block[b_next:]):
            if first_alinged_switch:
                b3, b4 = b4, b3

            # is_find = True
            if p2 != b2:
                continue
            else:
                b_next = bidx
                is_find = True
                if pidx != 0:
                    case += 1

            if b3 == p3 and b4 == p4:
                if not last_is_identical:
                    last_is_identical = True
                
                    if not last_time_error:
                        error += 1
                        last_time_error = True
                    else:
                        last_time_error = False
                else:
                    last_time_error = False

            elif b3 == p4 and b4 == p3:
                if last_is_identical:
                    last_is_identical = False

                    if not last_time_error:
                        error += 1
                        last_time_error = True
                    else:
                        last_time_error = False
                else:
                    last_time_error = False

            if first_alinged_switch:
                print(p2, b4, b3, p3, p4, error, last_is_identical, last_time_error)
            else:
                print(p2, b3, b4, p3, p4, error, last_is_identical, last_time_error)

        if not is_find:
            unphased += 1
            print("NOT IN SAME BLOCK ", predicted[0][0], predicted[0][1], p2, pidx)
            return (predicted[0][0], int(predicted[0][1]), error, case, unphased, total), pidx 

    timer.stop("loop")

    # print(predicted[0][0], predicted[0][1], error)

    return (predicted[0][0], int(predicted[0][1]), error, case, unphased, total), pidx + 1

predicted_fn = sys.argv[1]
blocks_ans_fn = sys.argv[2]

predicted_file = open(predicted_fn, 'r')
blocks_ans_file = open(blocks_ans_fn, 'r')

timer.start("block_ans")
blocks_dict = dict()
loc2blocks = dict()
# read blocks ans
for line in blocks_ans_file:
    line = line.strip()
    chrom, location, ps_id, h1, h2 = line.split("\t")
    key = chrom + "_" + ps_id
    if key not in blocks_dict:
        loc2blocks[chrom + "_" + location] = key
        blocks_dict[key] = [(chrom, location, h1, h2), ]
    else:
        loc2blocks[chrom + "_" + location] = key
        blocks_dict[key].append((chrom, location, h1, h2)) 

blocks_ans_file.close()

timer.stop("block_ans")

timer.start("predict")
predicted_dict = dict()
for line in predicted_file:
    line = line.strip()
    chrom, location, ps_id, h1, h2 = line.split("\t")
    key = chrom + "_" + ps_id
    if key not in predicted_dict:
        predicted_dict[key] = [(chrom, location, h1, h2), ]
    else:
        predicted_dict[key].append((chrom, location, h1, h2))

predicted = None
answer_block = None
   
timer.stop("predict")


def b_search_while(lst, target):
    left = 0
    right = len(lst)-1
    while left <= right:
        avg = (left + right)//2
        mid = int(lst[avg][1])

        if mid < target:
            left = avg + 1
        elif mid > target:
            right = avg - 1
        else:
            return avg

    return left


u_error = 0
u_case = 0
u_unphased = 0
u_total = 0

ps_data = []
s_total = 0
s_case = 0
s_unphased = 0
s_error = 0

repeat_count = 0 


break_in_predict = dict()

for k, vals in predicted_dict.items():
    v_idx = 0
    v_length = len(vals)
    first_ps = vals[0][1]
    while v_length != v_idx:
        f_row = vals[v_idx]
        key = f_row[0] + "_" + f_row[1]
        block_key = loc2blocks[key]
        block = blocks_dict[block_key]
        if block_key in break_in_predict:
            a, b = break_in_predict[block_key]
            break_in_predict[block_key] = [a+1, b]
        else:
            break_in_predict[block_key] = [0, len(block)]

        idxb = b_search_while(block, int(f_row[1]))
        d, v_step = check(block[idxb:], vals[v_idx:])
        v_idx += v_step        

        u_error += d[2]
        u_case += d[3]
        u_unphased += d[4]
        u_total += d[5]
        repeat_count += 1
        if repeat_count > 5:
            print("DEAD")
            u_error = -1
            break
            # sys.exit()
        print("==")

    if u_error != -1:

        s_error += u_error
        s_case += u_case
        s_unphased += u_unphased
        s_total += u_total
        ps_data.append((d[0], first_ps, u_error, u_case, u_unphased, u_total, repeat_count))

    repeat_count = 0
    u_error = 0
    u_case = 0
    u_unphased = 0
    u_total = 0

out_pickle_file = sys.argv[3] + '.pickle'


if os.path.isfile(out_pickle_file):
    os.remove(out_pickle_file)

print(len(ps_data))
fh = open(out_pickle_file, 'wb')
# pack = {'ps_data': ps_data, 'break_in_predict': break_in_predict}
pack = {'ps_data': ps_data}
pickle.dump(pack, fh)
fh.close()

print(s_error, s_case, s_unphased, s_total)
