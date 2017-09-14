import sys
from collections import OrderedDict
import numpy as np
from math import pow, log, floor
import operator
from GHaplo.math import normalized_score_scale_to_max


def create_sv_freq_dict(read_lst, var_lst, read_mtx):
    """
    Collect single variant infomation

    """
    EMPTY_CODE = 0
    sv_freq_dict = OrderedDict()

    for ridx, read in enumerate(read_lst):
        for vidx in range(read[1], read[2] + 1):
            if read_mtx[ridx, vidx] == EMPTY_CODE:
                continue
            loc = var_lst[vidx]

            try:

                fdict = sv_freq_dict[loc]
                
                try:
                    fdict[read_mtx[ridx, vidx]] += 1
                except KeyError:
                    fdict[read_mtx[ridx, vidx]] = 1

            except KeyError:
                fdict = dict()
                fdict[read_mtx[ridx, vidx]] = 1
                sv_freq_dict[loc] = fdict

    return sv_freq_dict


def create_pv_freq_dict(read_lst, var_lst, read_mtx):
    """
    Collect all variants pair information

    """
    EMPTY_CODE = 0
    GAP_CODE = 6
    MISSING_CODE = [EMPTY_CODE, GAP_CODE]

    var_amount = len(var_lst)
    phase_pairs_called_list = OrderedDict()

    for ridx, read in enumerate(read_lst):
        for fidx in range(read[1], read[2] + 1):
            for sidx in range(fidx + 1, read[2] + 1):
                if fidx > sidx:
                    continue

                pairs_key = fidx * var_amount + sidx
                fvar = read_mtx[ridx, fidx]
                svar = read_mtx[ridx, sidx]

                if fvar not in MISSING_CODE and svar not in MISSING_CODE:
                    pv_id = '_'.join(map(str, [fvar, svar]))
                    if pairs_key not in phase_pairs_called_list:
                        pv_dict = {pv_id: 1}
                        phase_pairs_called_list[pairs_key] = pv_dict
                    else:
                        pv_dict = phase_pairs_called_list[pairs_key]
                        pv_dict[pv_id] = pv_dict[pv_id] + 1 if pv_id in pv_dict else 1

    return phase_pairs_called_list


def _calc_n1_n2(outcome_sum, sol_cmt_0, sol_cmt_1):
    """
    """
    base_size = 100
    scale_ratio = base_size / outcome_sum
    n = round(base_size)
    n1 = floor(sol_cmt_0 * scale_ratio)
    n2 = round(sol_cmt_1 * scale_ratio)
    n3 = n - n1 - n2
    return n, n1, n2, n3


def calc_score_matrix(phase_pairs_called_list, phase_variants_called_list, encoding_tb, alpha):
    """
    Calc_priority_of_solutions:
    Calculate based on solutions(two unconficted outcomes).
    """
    # create score_matrix
    variants_number = len(phase_variants_called_list)
    score_matrix = np.empty(shape=(variants_number, variants_number))
    score_matrix[:] = np.NINF

    # Quartile Q2
    pairs_values = map(lambda p: p.values(), phase_pairs_called_list.values())
    pairs_values_sorted = sorted([item for sublist in pairs_values for item in sublist])

    try:

        q2 = pairs_values_sorted[len(pairs_values_sorted)//4]
        pairs_count = sum(pairs_values_sorted)

    except:
        print(phase_pairs_called_list)
        print(phase_variants_called_list)
        print(pairs_values_sorted)
        print(len(pairs_values_sorted)//4)
        sys.exit()

    # calculate average ave_coverage
    pairs_count = sum(pairs_values_sorted)
    ave_coverage = floor(pairs_count / len(phase_pairs_called_list))

    for pairs_key, pairs_dict in phase_pairs_called_list.items():
        base = pairs_key // variants_number
        shift = pairs_key % variants_number

        pairs_sorted_dict = sorted(pairs_dict.items(), key=operator.itemgetter(1), reverse=True)
        pairs_number = sum([p[1] for p in pairs_sorted_dict])
 
        max_score = np.NINF
        for idx, p in enumerate(pairs_sorted_dict):
            singleton_observed = True
            for q in pairs_sorted_dict[idx+1:]:
                p_l, p_r = p[0].split("_")
                q_l, q_r = q[0].split("_")

                # heterozygous assumption
                if p_l != q_l and p_r != q_r:
                    singleton_observed = False
                    n, n1, n2, n3 = _calc_n1_n2(pairs_number, p[1], q[1])
                    score = normalized_score_scale_to_max(n, n1, n2, p[1]+q[1], ave_coverage, alpha)
                    if pairs_number < 5:
                        score -= 1e12
                    if score > max_score:
                        max_score = score

            if singleton_observed:
                n, n1, n2, n3 = _calc_n1_n2(pairs_number, p[1], 0)
                score = normalized_score_scale_to_max(n, n1, n2, p[1], ave_coverage, alpha) - 1e10
                if score > max_score:
                    max_score = score

        score_matrix[base, shift] = max_score
        score_matrix[shift, base] = max_score
 
    return score_matrix
