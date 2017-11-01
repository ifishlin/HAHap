import numpy as np
import operator
import copy
import logging
import sys
from math import ceil

logger = logging.getLogger(__name__)


class HiMergeNode:
    """
    """
    def __init__(self, name, value, left=None, right=None, h1=None, h2=None):
        self.name = name
        self.value = value
        self.left = left
        self.right = right
        self.h1 = h1
        self.h2 = h2
        self.name_split = list(map(int, name.split("_")))

    def set_name(self, name):
        self.name = name
        self.name_split = list(map(int, name.split("_")))


def get_max_cs(c1_name, c2_name, cs_mx, rank):
    """
    """

    cs_dict = dict()
    c1_ns = list(map(int, c1_name.split('_')))
    c2_ns = list(map(int, c2_name.split('_')))

    for i in c1_ns:
        for j in c2_ns:
            cs = cs_mx[i][j] if j > i else cs_mx[j][i]
            cs_dict[str(i) + "_" + str(j)] = cs

    sorted_cs = sorted(cs_dict.items(), key=operator.itemgetter(1, 0), reverse=True)

    a, b = sorted_cs[rank][0].split("_")

    return int(a), int(b), sorted_cs[rank][1]


def get_pv_key(c1_idx, c2_idx, phase_total):
    """
    """
    if c1_idx < c2_idx:
        return int(c1_idx) * phase_total + int(c2_idx)
    elif c1_idx > c2_idx:
        return int(c2_idx) * phase_total + int(c1_idx)
    else:
        sys.exit()


def _get_anchors(c1_node, c2_node, c1_idx, c2_idx):
    """
    """
    c1_node_name = c1_node.name.split('_')
    c2_node_name = c2_node.name.split('_')

    node1_idx = c1_node_name.index(str(c1_idx))
    node2_idx = c2_node_name.index(str(c2_idx))

    c1_anchor, c2_anchor = None, None

    if c1_node.h1 is not None:
        c1_anchor = [c1_node.h1[node1_idx], c1_node.h2[node1_idx]]

    if c2_node.h1 is not None:
        c2_anchor = [c2_node.h1[node2_idx], c2_node.h2[node2_idx]]

    return c1_anchor, c2_anchor


def _ha_merge_two_nodes(c1_node, c2_node):
    """
    """
    node = HiMergeNode('_'.join([c1_node.name, c2_node.name]), value='intermediate')
    node.left = c1_node
    node.right = c2_node
    return node


def _ha_phasing_main(sorted_sup, c1_anchor, c2_anchor, c1_idx, c2_idx, c1_node, c2_node,
                   phase_loc, phase_allele, codes, fragments,
                   fragment_se, embed_ops, last_ops, timer):
    """
    """

    timer.start("061.reasonable solution")
    #if c1_idx >= c2_idx:
    #    sys.exit()

    base, shift = (c1_idx, c2_idx)

    # prioritize sols

    exclude = []
    sols = []
    heter_assumption = True

    #  two pair sols first

    tmp_sols = []
    tmp_read_cnt = dict()

    for idx, x in enumerate(sorted_sup):
        for y in sorted_sup[idx + 1:]:
            x1, x2 = x[0].split("_")
            y1, y2 = y[0].split("_")

            if heter_assumption and x1 != y1 and x2 != y2:
                tmp_sols.append([[[int(x1), int(x2)], x[1]], [[int(y1), int(y2)], y[1]]])
                tmp_read_cnt[len(tmp_read_cnt)] = x[1] + y[1]
                exclude.append(x)
                exclude.append(y)

    read_cnt_sorted = sorted(tmp_read_cnt.items(), key=operator.itemgetter(1))
    for key, _ in read_cnt_sorted[::-1]:
        sols.append(tmp_sols[key])
    # one pair sols later

    tmp_sols = []
    tmp_read_cnt = dict()

    for x in sorted_sup:
        x1, x2 = x[0].split("_")
        if x not in exclude:
            tmp_sols.append([[[int(x1), int(x2)], x[1]]])
            tmp_read_cnt[len(tmp_read_cnt)] = x[1]

    read_cnt_sorted = sorted(tmp_read_cnt.items(), key=operator.itemgetter(1))
    for key, value in read_cnt_sorted[::-1]:
        sols.append(tmp_sols[key])

    result = []


    timer.stop("061.reasonable solution")

    # ambigeous pairs remain only
    if last_ops and not (c1_anchor is None and c2_anchor is None):
        if c1_anchor is None:
            f_maj, f_min = list(map(lambda h: codes[h], phase_allele[base]))
            #c1_node.set_name(str(base)) #
            c1_node.h1 = [f_maj, ]
            c1_node.h2 = [f_min, ]
        if c2_anchor is None:
            e_maj, e_min = list(map(lambda h: codes[h], phase_allele[shift]))
            #c2_node.set_name(str(shift)) #
            c2_node.h1 = [e_maj, ]
            c2_node.h2 = [e_min, ]

        result = _local_pre_MEC_search(c1_node, c2_node, fragments, fragment_se, phase_loc, codes, timer)

    elif c1_anchor is None and c2_anchor is None:
        f_maj, f_min = list(map(lambda h: codes[h], phase_allele[base]))
        e_maj, e_min = list(map(lambda h: codes[h], phase_allele[shift]))

        for s in sols:
            if len(s) > 1:
                result = [s[0][0], s[1][0]]
                break
            elif len(s) == 1:
                f_h1, e_h1 = s[0][0]

                f_h2 = f_maj if f_h1 == f_min else f_min
                e_h2 = e_maj if e_h1 == e_min else e_min

                '''
                # alleles of sinlge pair solution must be in alleles of variants calling
                if f_h1 in (f_maj, f_min) and e_h1 in (e_maj, e_min):
                    f_h2 = f_maj if f_h1 == f_min else f_min
                    e_h2 = e_maj if e_h1 == e_min else e_min
                else:
                    continue
                '''

                h1 = s[0][0]
                h2 = [f_h2, e_h2]
                result = [h1, h2]
                break

    elif c1_anchor is not None and c2_anchor is None:
        e_maj, e_min = list(map(lambda h: codes[h], phase_allele[shift]))

        for s in sols:
            if len(s) > 1:
                m1, m2 = s[0][0]
                m3, m4 = s[1][0]
                if m1 == c1_anchor[0] and m3 == c1_anchor[1]:
                    h1 = copy.copy(c1_node.h1)
                    h2 = copy.copy(c1_node.h2)
                    h1.append(m2)
                    h2.append(m4)
                    result = [h1, h2]
                    break
                elif m1 == c1_anchor[1] and m3 == c1_anchor[0]:
                    h1 = copy.copy(c1_node.h1)
                    h2 = copy.copy(c1_node.h2)
                    h1.append(m4)
                    h2.append(m2)
                    result = [h1, h2]
                    break

            elif len(s) == 1:
                m1, m2 = s[0][0]
                m4 = e_maj if m2 == e_min else e_min

                if m1 == c1_anchor[0]:
                    h1 = copy.copy(c1_node.h1)
                    h2 = copy.copy(c1_node.h2)
                    h1.append(m2)
                    h2.append(m4)
                    result = [h1, h2]
                    break
                elif m1 == c1_anchor[1]:
                    h1 = copy.copy(c1_node.h1)
                    h2 = copy.copy(c1_node.h2)
                    h1.append(m4)
                    h2.append(m2)
                    result = [h1, h2]
                    break

    elif c1_anchor is None and c2_anchor is not None:
        f_maj, f_min = list(map(lambda h: codes[h], phase_allele[base]))

        for s in sols:
            if len(s) > 1:
                m1, m2 = s[0][0]
                m3, m4 = s[1][0]
                if m2 == c2_anchor[0] and m4 == c2_anchor[1]:
                    h1 = [m1] + c2_node.h1
                    h2 = [m3] + c2_node.h2
                    result = [h1, h2]
                    break
                elif m2 == c2_anchor[1] and m4 == c2_anchor[0]:
                    h1 = [m1] + c2_node.h2
                    h2 = [m3] + c2_node.h1
                    result = [h1, h2]
                    break

            elif len(s) == 1:
                m1, m2 = s[0][0]
                m3 = f_maj if m1 == f_min else f_min

                if m2 == c2_anchor[0]:
                    h1 = [m1] + c2_node.h1
                    h2 = [m3] + c2_node.h2
                    result = [h1, h2]
                    break
                elif m2 == c2_anchor[1]:
                    h1 = [m1] + c2_node.h2
                    h2 = [m3] + c2_node.h1
                    result = [h1, h2]
                    break

    elif c1_anchor is not None and c2_anchor is not None:

        if embed_ops:
            junction_cnt = 0
            if len(c1_node.name_split) > len(c2_node.name_split):
                for c in c2_node.name_split:
                    junction_cnt = junction_cnt + 1 if c - 1 in c1_node.name_split else junction_cnt
                    junction_cnt = junction_cnt + 1 if c + 1 in c1_node.name_split else junction_cnt
            else:
                for c in c1_node.name_split:
                    junction_cnt = junction_cnt + 1 if c - 1 in c2_node.name_split else junction_cnt
                    junction_cnt = junction_cnt + 1 if c + 1 in c2_node.name_split else junction_cnt

        # embedded threshold
        embed_threshold = 3
        if embed_ops and junction_cnt > embed_threshold:

            result = _local_pre_MEC_search(c1_node, c2_node, fragments, fragment_se, phase_loc, codes, timer)  

        else:
            for s in sols:
                if len(s) > 1:

                    m1, m2 = s[0][0]
                    m3, m4 = s[1][0]

                    if [m1, m3] == c1_anchor and [m2, m4] == c2_anchor:
                        n1 = c1_node.h1 + c2_node.h1
                        n2 = c1_node.h2 + c2_node.h2
                        result = [n1, n2]
                        break
                    elif [m3, m1] == c1_anchor and [m2, m4] == c2_anchor:
                        n1 = c1_node.h1 + c2_node.h2
                        n2 = c1_node.h2 + c2_node.h1
                        result = [n1, n2]
                        break
                    elif [m1, m3] == c1_anchor and [m4, m2] == c2_anchor:
                        n1 = c1_node.h1 + c2_node.h2
                        n2 = c1_node.h2 + c2_node.h1
                        result = [n1, n2]
                        break
                    elif [m3, m1] == c1_anchor and [m4, m2] == c2_anchor:
                        n1 = c1_node.h1 + c2_node.h1
                        n2 = c1_node.h2 + c2_node.h2
                        result = [n1, n2]
                        break

                elif len(s) == 1:

                    m1, m2 = s[0][0]
                    if m1 == c1_anchor[0] and m2 == c2_anchor[0]:
                        n1 = c1_node.h1 + c2_node.h1
                        n2 = c1_node.h2 + c2_node.h2
                        result = [n1, n2]
                        break
                    elif m1 == c1_anchor[0] and m2 == c2_anchor[1]:
                        n1 = c1_node.h1 + c2_node.h2
                        n2 = c1_node.h2 + c2_node.h1
                        result = [n1, n2]
                        break
                    elif m1 == c1_anchor[1] and m2 == c2_anchor[0]:
                        n1 = c1_node.h1 + c2_node.h2
                        n2 = c1_node.h2 + c2_node.h1
                        result = [n1, n2]
                        break
                    elif m1 == c1_anchor[1] and m2 == c2_anchor[1]:
                        n1 = c1_node.h1 + c2_node.h1
                        n2 = c1_node.h2 + c2_node.h2
                        result = [n1, n2]
                        break

    # Merge Failure
    if len(result) < 2:
        logging.warning("Merge Fail, Can't find two haplo that connect both sides "
                        + str(c1_idx) + ' ' + str(c2_idx))
        return None, None

    return [r for r in result[:2]]

def _local_pre_MEC_search(c1_node, c2_node, fragments, fragment_se, phase_loc, codes, timer):
    """
    """
    sol1_h1, sol1_h2, sol2_h1, sol2_h2, sorted_name, t_sol1, t_sol2 \
        = _build_two_sols(c1_node, c2_node)

    penalty_1 = 0
    penalty_2 = 0

    for f_name, s, e in fragment_se:
        c1_f = False
        c2_f = False

        if e < sorted_name[0]:
            continue
        if s > sorted_name[-1]:
            break

        for i in range(s, e):
            if c1_f and c2_f:
                break
            elif i in c1_node.name_split:
                c1_f = True
            elif i in c2_node.name_split:
                c2_f = True

        if c1_f and c2_f:
            fragment_loc = []
            fragment_allele = []
            for p in fragments[f_name]:
                chrom, loc, allele, _ = p
                fragment_loc.append(phase_loc.index(loc))
                fragment_allele.append(codes[allele])
            timer.start("062._local_MEC_search")
            p1, p2 = _local_MEC_search(fragment_loc, fragment_allele, sol1_h1, sol1_h2, sol2_h1, sol2_h2, sorted_name)
            timer.stop("062._local_MEC_search")
            penalty_1 += p1
            penalty_2 += p2

    result = t_sol2 if penalty_1 > penalty_2 else t_sol1
    return result


def _local_MEC_search(fragment_loc, fragment_allele, sol1_h1, sol1_h2, sol2_h1, sol2_h2, sorted_name):
    """
    """
    p1, p2, p3, p4 = (0, 0, 0, 0)
    for loc, a in zip(fragment_loc, fragment_allele):
        # check candidate 1
        if loc in sorted_name:
            o1 = sol1_h1[sorted_name.index(loc)]
            o2 = sol1_h2[sorted_name.index(loc)]
            o3 = sol2_h1[sorted_name.index(loc)]
            o4 = sol2_h2[sorted_name.index(loc)]
            p1 = p1 if o1 == a else p1 + 1
            p2 = p2 if o2 == a else p2 + 1
            p3 = p3 if o3 == a else p3 + 1
            p4 = p4 if o4 == a else p4 + 1

    return min(p1, p2), min(p3, p4)


def _build_two_sols(c1_node, c2_node):
    """
    """
    pre_name = c1_node.name_split + c2_node.name_split
    sorted_name = sorted(pre_name)

    pre_sol1_h1 = c1_node.h1 + c2_node.h1
    pre_sol1_h2 = c1_node.h2 + c2_node.h2
    pre_sol2_h1 = c1_node.h1 + c2_node.h2
    pre_sol2_h2 = c1_node.h2 + c2_node.h1

    sol1_h1 = list()
    sol1_h2 = list()
    sol2_h1 = list()
    sol2_h2 = list()

    for s in sorted_name:
        idx = pre_name.index(s)
        sol1_h1.append(pre_sol1_h1[idx])
        sol1_h2.append(pre_sol1_h2[idx])
        sol2_h1.append(pre_sol2_h1[idx])
        sol2_h2.append(pre_sol2_h2[idx])

    return sol1_h1, sol1_h2, sol2_h1, sol2_h2, sorted_name, [pre_sol1_h1, pre_sol1_h2], \
                [pre_sol2_h1, pre_sol2_h2]


def argmax_single_link_vector(vector, timer):
    timer.start("argmax_single_link_vector")
    max_i = -1
    max_j = -1
    max_cs = -np.inf
    for idx, (_, _, closest, cs) in enumerate(vector):
        if cs > max_cs:
            max_i = idx
            max_j = closest
            max_cs = cs

    timer.stop("argmax_single_link_vector")
    return max_i, max_j


def merge_single_link_vector(reduced_mx, vector, i, j, timer):
    timer.start("merge_single_link_vector")
    max_idx = 0
    max_value = -np.inf
    i_row = reduced_mx[i]
    j_row = reduced_mx[j]

    # copy cs to data i
    for idx in range(len(vector)):
        a = i_row[idx]
        b = j_row[idx]

        t = a if a > b else b
        reduced_mx[idx][i] = t
        i_row[idx] = t

    # remove col and row of data j
    reduced_mx[:, j] = -np.inf
    reduced_mx[j] = -np.inf

    # keep m(i,i) none 
    reduced_mx[i][i] = -np.inf

    # update closest of data j to none
    vector[j][2:] = [-1, -np.inf]

    # search closest of data i
    i_row = reduced_mx[i]
    for idx in range(len(reduced_mx)):
        if i_row[idx] > max_value:
            max_idx = idx
            max_value = i_row[idx]

    # update closest of data i
    vector[i][0] = vector[i][0] + "_" + vector[j][0]
    vector[i][2:] = [max_idx, float(max_value)]
    vector[i][1] += 1

    # change cloest to j to i
    for idx in range(len(vector)):
        if vector[idx][2] == j:
            vector[idx][2] = i

    timer.stop("merge_single_link_vector")

def ha_phasing(vars_pool, pairs_sup, cs_mx, phase_loc, phase_allele, codes, fragments, fragment_se, embed_ops, last_ops, timer):
    """
    """

    phase_total = len(phase_loc)
    reduced_mx = np.copy(cs_mx)

    single_link_vector = []
    for i in range(phase_total):
        col = reduced_mx[i, :]
        idx = np.argmax(col)
        single_link_vector.append([str(i), 1, idx, float(col[idx])])

    merge_cnt = 0
    ops_threshold = 0.005
    ops_cutoff = ceil(phase_total * ops_threshold)
    
    last_ops_flag = False

    # vertex_1_idx, < vertex_2_idx
    vertex_1_idx, vertex_2_idx = argmax_single_link_vector(single_link_vector, timer)
    while vertex_1_idx != -1 and vertex_2_idx != -1:
        merge_cnt += 1
        v1_name, v1_cnt = single_link_vector[vertex_1_idx][:2]
        v2_name, v2_cnt = single_link_vector[vertex_2_idx][:2]

        for rank in range(0, v1_cnt * v2_cnt):
            x, y, cs = get_max_cs(v1_name, v2_name, cs_mx, rank)

            if (phase_total - merge_cnt) <= ops_cutoff and last_ops  and cs < -100:
                # if phase_total - merge_cnt == 1 and last_ops:
                last_ops_flag = True

            c1_node, c2_node, c1_idx, c2_idx = (vars_pool[v1_name], vars_pool[v2_name], x, y) \
                                    if y > x else (vars_pool[v2_name], vars_pool[v1_name], y, x)

            pv_key = get_pv_key(c1_idx, c2_idx, phase_total)

            if cs == np.NINF:
                break

            if pv_key in pairs_sup:
                node = _ha_merge_two_nodes(c1_node, c2_node)
                #print("name")
                #print(c1_node.name, c2_node.name)
                #print(node.name)
                c1_anchor, c2_anchor = _get_anchors(c1_node, c2_node, c1_idx, c2_idx)

                sorted_sup = sorted(pairs_sup[pv_key].items(), key=operator.itemgetter(1, 0), reverse=True)

                timer.start("060._ha_phasing_main")
                node.h1, node.h2 = _ha_phasing_main(sorted_sup, c1_anchor, c2_anchor, c1_idx, c2_idx, c1_node, c2_node,
                                        phase_loc, phase_allele, codes,
                                        fragments, fragment_se, embed_ops, last_ops_flag, timer)
                timer.stop("060._ha_phasing_main")

                # careful for end condition? only one None?
                if node.h1 is not None and node.h2 is not None:
                    break
            else:
                logging.error("can't find in pairs_sup")
                sys.exit()

        vars_pool[node.name] = node
        del vars_pool[c1_node.name]
        del vars_pool[c2_node.name]

        if x < y:
            merge_single_link_vector(reduced_mx, single_link_vector, vertex_1_idx, vertex_2_idx, timer)
        else:
            merge_single_link_vector(reduced_mx, single_link_vector, vertex_2_idx, vertex_1_idx, timer)

        vertex_1_idx, vertex_2_idx = argmax_single_link_vector(single_link_vector, timer)



def save_phasing(vars_pool, phase_loc, output_dict):
    """
    """

    for k, v in vars_pool.items():

        cluster_name = k
        cluster_node = v

        cluster_name_seq = list(map(int, cluster_name.split("_")))
        cluster_name_seq_sorted = sorted(cluster_name_seq)

        h1 = cluster_node.h1
        h2 = cluster_node.h2

        # single variant
        if h1 is None or h2 is None:
            continue

        ps_id = phase_loc[cluster_name_seq_sorted[0]]

        for s in cluster_name_seq_sorted:
            location = str(phase_loc[s])
            idx = cluster_name_seq.index(s)
            h1_allele = h1[idx]
            h2_allele = h2[idx]

            output_dict['_'.join([location])] = (ps_id, h1_allele, h2_allele)

    return output_dict
