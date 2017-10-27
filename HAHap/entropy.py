from collections import OrderedDict
import numpy as np
from math import floor
import operator
from HAHap.math import normalized_score_scale_to_max

EMPTY_CODE = 0
GAP_CODE = 6
MISSING_CODE = [EMPTY_CODE, GAP_CODE]

def create_pairs_sup(fragment_se, phase_total, sf_mx):
    """
    Collect all variants pair information

    """
    pairs_sup = OrderedDict()

    for r_idx, (_, s, e) in enumerate(fragment_se):
        e += 1
        for f_idx in range(s, e):
            for s_idx in range(f_idx + 1, e):
                pairs_key = f_idx * phase_total + s_idx
                f_allele = sf_mx[r_idx, f_idx]
                s_allele = sf_mx[r_idx, s_idx]

                if f_allele not in MISSING_CODE and s_allele not in MISSING_CODE:
                    p_id = '_'.join(map(str, [f_allele, s_allele]))
                    if pairs_key not in pairs_sup:
                        pairs_sup[pairs_key] = {p_id: 1}
                    else:
                        sup = pairs_sup[pairs_key]
                        sup[p_id] = sup[p_id] + 1 if p_id in sup else 1

    return pairs_sup


def _calc_n1_n2(outcome_sum, sol_cnt_0, sol_cnt_1):
    """
    """
    base_size = 100
    scale_ratio = base_size / outcome_sum
    n = round(base_size)
    n1 = floor(sol_cnt_0 * scale_ratio)
    n2 = round(sol_cnt_1 * scale_ratio)
    n3 = n - n1 - n2
    return n, n1, n2, n3


def calc_cs_mx(pairs_sup, phase_loc, distance, lct, alpha=1):
    """
    Calc_priority_of_solutions:
    Calculate based on solutions(two unconficted outcomes).
    """

    phase_total = len(phase_loc)
    # create cs_mx
    cs_mx = np.empty(shape=(phase_total, phase_total))
    cs_mx[:] = np.NINF

    # Quartile Q2
    pairs_count_sorted = sorted([sum(p.values()) for p in pairs_sup.values()])
    q2 = pairs_count_sorted[len(pairs_count_sorted)//2]
    max_coverage = pairs_count_sorted[-1]

    print("q2:",q2)

    for pairs_key, sups in pairs_sup.items():
        base = pairs_key // phase_total
        shift = pairs_key % phase_total

        sup = sorted(sups.items(), key=operator.itemgetter(1), reverse=True)
        pairs_total = sum([p[1] for p in sup])

        #print(base, shift, sups)

        max_score = np.NINF
        for idx, p in enumerate(sup):
            singleton_observed = True
            for q in sup[idx+1:]:
                p_l, p_r = p[0].split("_")
                q_l, q_r = q[0].split("_")
                #print(p_l, p_r , q_l, q_r)
                # heterozygous assumption
                if p_l != q_l and p_r != q_r:
                    # print("C,",singleton_observed)
                    singleton_observed = False
                    #print(singleton_observed)
                    n, n1, n2, n3 = _calc_n1_n2(pairs_total, p[1], q[1])
                    score = normalized_score_scale_to_max(n, n1, n2, p[1]+q[1], max_coverage, alpha)
                    if distance == 5:
                        score -= abs(phase_loc[shift] - phase_loc[base]) * 1e-5
                    elif distance == 10:
                        score -= abs(phase_loc[shift] - phase_loc[base]) * 1e-10
                    else:
                        pass
                    # if pairs_total < q2
                    if lct == 0:
                        if pairs_total < q2:
                            score -= 1e12
                        if score > max_score:
                            max_score = score
                    else:
                        if pairs_total < lct:
                            score -= 1e12
                        if score > max_score:
                            max_score = score
                #print("A ",score)

            #print("D ",singleton_observed)
            if singleton_observed:
                n, n1, n2, n3 = _calc_n1_n2(p[1], p[1], 0)
                score = normalized_score_scale_to_max(n, n1, n2, p[1], max_coverage, alpha) - 1e10
                if lct == 0:
                    if pairs_total < q2:
                        score -= 1e12
                    if score > max_score:
                        max_score = score
                else:
                    if pairs_total < lct:
                        score -= 1e12
                    if score > max_score:
                        max_score = score

                if score > max_score:
                    max_score = score

                #print("B ",score)

        #print("cs_mx["+str(base)+","+str(shift)+"] = ",max_score)
        cs_mx[base, shift] = max_score
        cs_mx[shift, base] = max_score
 
    return cs_mx
