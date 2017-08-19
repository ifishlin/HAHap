#!/usr/bin/env python3

import os
import sys
from six.moves import cPickle as pickle

def overlap(blocks, predicts):
    pass 

def check(blocks, predicts):

    chrom = blocks[0][3]
    b1, b2, b3, b4 = blocks[0]
    p1, p2, p3 = predicts[0]
    ps_id = b1

    if b2 == p2 and b3 == p3:
        pass
    elif b2 == p3 and b3 == p2:
        blocks = list(map(lambda b:(b[0],b[2],b[1],b[3]),blocks))
        pass
    else:
        print("== ERROR")
        sys.exit()

    bidx = 1
    pidx = 1
    last_is_identical = True
    last_time_error = False
    error = 0
    unphased = 0
    case  = 0

    for p_idx, (p1, p2, p3) in enumerate(predicts[pidx:]):
        is_find = False
        for b_idx, (b1, b2, b3, b4) in enumerate(blocks[bidx:]):
            is_find = True
            if p1 != b1:
                continue
            else:
                case += 1
                bidx = b_idx

            if b2 == p2 and b3 == p3:
                if not last_is_identical:
                    last_is_identical = True

                    if not last_time_error:
                        error += 1
                        last_time_error = True
                    else:
                        last_time_error = False

            elif b2 == p3 and b3 == p2:
                if last_is_identical:
                    last_is_identical = False

                    if not last_time_error:
                        error += 1
                        last_time_error = True
                    else:
                        last_time_error = False

            #print(b1, error, b2, b3, p2, p3, last_is_identical, last_time_error, 'chrom'+chrom)

        if is_find == False:
            ## non phased
            case += 1
            unphased += 1
            pass


    '''
    for b_idx, (b1, b2, b3, b4) in enumerate(blocks[bidx:]):
        is_find = False
        for p_idx, (p1, p2, p3) in enumerate(predicts[pidx:]):
            is_find = True
            if b1 != p1:
                continue
            else: 
                case += 1 
                pidx = p_idx

            if b2 == p2 and b3 == p3:
                if not last_is_identical:
                    last_is_identical = True

                    if not last_time_error:
                        error += 1
                        last_time_error = True
                    else:
                        last_time_error = False

            elif b2 == p3 and b3 == p2:
                if last_is_identical:
                    last_is_identical = False

                    if not last_time_error:
                        error += 1 
                        last_time_error = True
                    else:
                        last_time_error = False

            #print(b1, error, b2, b3, p2, p3, last_is_identical, last_time_error, 'chrom'+chrom)

        if is_find == False:
            ## non phased
            case += 1
            unphased += 1
            pass
    ''' 

    #print(error, unphased, case)
    return (error, unphased, case)


def check_weired(blocks, predicts):
    print(blocks)

#pickle

ps_data = dict()

#read who's answer 
predicted_fname = sys.argv[1]
predicted_file = open(predicted_fname, 'r')

predict_dict = dict()
predict_belongTo_dict = dict()

#place phased set in dict by key chrom_ps_id 
for line in predicted_file:
    line = line.strip()
    chrom, location, ps_id, h1, h2 = line.split('\t')
    key = chrom + "_" + ps_id
    ps_id = int(ps_id)
    if key not in predict_dict:
        predict_dict[key] = [(int(location), h1, h2),]
    else:
        predict_dict[key].append((int(location), h1, h2))

    #dict for reverse checking
    predict_belongTo_dict[chrom + '_' + location] = ps_id

#read blocks we concerned
blocks_fname = sys.argv[2]
blocks_file  = open(blocks_fname, 'r')

t_error = 0
t_unphased = 0
t_case = 0

case_id = None
pre_case_id = None

blocks = []
ps_id = None
skip_case = False
for line in blocks_file:
    line = line.strip()
    e = line.split('\t')
    if len(e) == 5:
        #block data
        if ps_id != None and ps_id != e[2]:
            #print("== NOT IN SAME BLOCK")    
            skip_case = True
        else:
            ps_id = e[2]
        
        blocks.append((int(e[1]), e[3], e[4], e[0]))
    else:
        case_id = line.split()[5]

        #block header
        if skip_case == True:
            predicts = predict_dict[chrom + '_' + str(pre_case_id)]
            #check_weired(blocks, predicts)
            #ps_data[chrom + "_" + pre_case_id] = [-1, -1, -1, -1, -1, -1, len(blocks)]
            print(chrom + "_" + pre_case_id, -1, -1, -1, -1, -1, -1, len(blocks))
            blocks = []
            skip_case = False
            ps_id = None
            pre_case_id =  case_id
            continue

        
        if len(blocks) == 0:
            ps_id = None
            pre_case_id =  case_id
            continue

        try:

            chrom = blocks[0][3]

            f_ps_id = None
            e_ps_id = None

            for fidx, (f1, f2, f3, f4) in enumerate(blocks):
                if f4 + '_' + str(f1) in predict_belongTo_dict:
                    f_ps_id = predict_belongTo_dict[f4 + '_' + str(f1)]
                    break

            for eidx , (e1, e2, e3, e4) in enumerate(blocks[::-1]):
                if e4 + '_' + str(e1) in predict_belongTo_dict: 
                    e_ps_id = predict_belongTo_dict[e4 + '_' + str(e1)]
                    break

            if f_ps_id == None or e_ps_id == None:
                ps_data[chrom + "_" + pre_case_id] = [f_ps_id, e_ps_id, fidx, eidx, -1, -1, -1]
                #print(chrom + "_" + pre_case_id, f_ps_id, e_ps_id, fidx, eidx, -1, -1, -1)
                print("== None")
            elif f_ps_id != e_ps_id:
                ps_data[chrom + "_" + pre_case_id] = [f_ps_id, e_ps_id, fidx, eidx, -1, -1, -1]
                #print(chrom + "_" + pre_case_id, f_ps_id, e_ps_id, fidx, eidx, -1, -1, -1)
                #print("== NOT IN SAME BLOCKS")
            #elif fidx == 0 and eidx == 0:
            else:
                predicts = predict_dict[chrom + '_' + str(f_ps_id)]
                predict_locations = list(map(lambda i:i[0], predicts))
                idx = predict_locations.index(f1)
                e, u, c  = check(blocks[fidx:], predicts[idx:])
                t_error += e
                t_unphased += u
                t_case += c
                ps_data[chrom + "_" + pre_case_id] = [f_ps_id, e_ps_id, fidx, eidx, e, u, c]
                #print(chrom + "_" + pre_case_id, f_ps_id, e_ps_id, fidx, eidx, e, u, c)
        except KeyError:
            print(f_ps_id, e_ps_id)
            print("== KeyError")


        pre_case_id =  case_id
        blocks = []
        ps_id = None
        pass

if os.path.isfile('test.pickle'):
   os.remove('test.pickle')

print(len(ps_data))
fh = open('test.pickle', 'wb')
pack = {'ps_data': ps_data}
pickle.dump(pack,fh)
fh.close()

print(t_error, t_unphased, t_case)
print('done')
