import numpy as np
import operator
import copy
import logging
import sys

b_threshold=3
local_brute=True

class HiMergeNode:
    """
    Tree model of hierarchical assembly
    """
    def __init__(self, name, value, left = None, right= None, h1 = None, h2 = None):
        self.name  = name
        self.value = value
        self.left  = left
        self.right = right
        self.h1    = h1
        self.h2    = h2

    def set_left(self, left):
        self.left = left

    def set_right(self, right):
        self.right = right

    def set_value(self, value):
        self.value = value

    def set_h1(self, h1):
        self.value = h1

    def set_h2(self, h2):
        self.value = h2

def copy_score_matrix(mat):
    """
    Copy a new matrix
    Args:
    Return:
    """
    return np.copy(mat)

def get_maxscore(mat):
    """
    Get max score in the matrix

    Args:

    Return:
        base: row
        shift: col
        mut_info: max value
    """
    mlen = len(mat[0])
    idx = np.argmax(mat)
    base  = idx // mlen
    shift = idx % mlen
    return base, shift, mat[base][shift]

def reduce_matrix( mat, row, col):
    """
    Update and reduce matrix, update values in row 'row'
    delete column 'row' and 'col'

    Args:
    Return:
           np.delete create new matrix.
    """
    keep = row
    remove = col
    for i in range(len(mat)):
        if(i == remove):
            continue

        if(keep > i):
            item1 = mat[i][keep]
            item2 = mat[i][remove] if remove > i else mat[remove][i]
            mat[i][keep] = max(item1, item2)
        if(i > keep):
            item3 = mat[keep][i]
            item4 = mat[remove][i] if i > remove else mat[i][remove]
            mat[keep][i] = max(item3, item4)
    #mat[keep][keep] = np.NINF

    #print("DELETE row col ROW:", row, "COL:", col)
    mat = np.delete(mat, col, 0)
    mat = np.delete(mat, col, 1)
    return mat

def get_maxscore_between_segments(group1, group2, mut_info, rank):
    """
    Find the most related pairs in group1 and group2

    Args:
    Return:
          related idx in group1
          related idx in group2
          score
    """
    rank_dict = dict()
    g1 = list(map(int, group1.split('_')))
    g2 = list(map(int, group2.split('_')))
    for i in g1:
        for j in g2:
            tvalue = mut_info[i][j] if j > i else mut_info[j][i]
            rank_dict[str(i)+"_"+str(j)] = tvalue

    sorted_x = sorted(rank_dict.items(), key=operator.itemgetter(1), reverse=True)


    a ,b = sorted_x[rank][0].split("_")

    return int(a), int(b), sorted_x[rank][1]

def get_matrix_encode_key(idx1, idx2, var_num):
    """
    Get encode key of matrix cell.

    Args:
    Return:
    """
    if(idx1 < idx2):
        return int(idx1) * var_num + int(idx2)
    elif(idx1 > idx2):
        return int(idx2) * var_num + int(idx1)
    else:
        sys.exit()
        return 0

def get_matrix_idx(key, var_num):
    """
    """
    base  = key // var_num
    shift = key % var_num
    return base, shift


def get_phased_haplo(cluster1, cluster2, idx_in_c1, idx_in_c2):
    """
    Get determined haplotypes on idx_in_c1 of node1 and idx_in_c2 of node2 
    """

    cluster1_name = cluster1.name.split('_')
    cluster2_name = cluster2.name.split('_')

    # what condition not in name??
    # node1_idx = cluster1_name.index(str(idx_in_c1)) if str(idx_in_c1) in cluster1_name else None
    # node2_idx = cluster2_name.index(str(idx_in_c2)) if str(idx_in_c2) in cluster2_name else None

    node1_idx = cluster1_name.index(str(idx_in_c1))
    node2_idx = cluster2_name.index(str(idx_in_c2))

    if cluster1.h1 != None:
        c1_anchor = [cluster1.h1[node1_idx], cluster1.h2[node1_idx]]
    else:
        c1_anchor = None

    if(cluster2.h1 != None):
        c2_anchor = [cluster2.h1[node2_idx], cluster2.h2[node2_idx]]
    else:
        c2_anchor = None

    return c1_anchor, c2_anchor

def hc_create_parent_node(l_node, r_node):
    """
    Hirarchical Clustering
    Create parent node of l_node, r_node.
    Concatenate node1/node2 name as it's name.

    Args:

    Return:
        HiMergeNode
    """
    node = HiMergeNode('_'.join([l_node.name, r_node.name]), value='intermediate')
    node.left  = l_node
    node.right = r_node
    return node

def choice_haplos(observed_pairs_dict, cluster1, cluster2, idx_in_c1, idx_in_c2, sv_dict, phase_variants_location_list, phase_variants_called_list, encoding_tb, heter, hratio, pv_dict, variant_num):
    """
    g1_idx, g2_idx, the closet two nodes in each group.
    """
    c1_anchor, c2_anchor = get_phased_haplo(cluster1, cluster2, idx_in_c1, idx_in_c2)

    sorted_x = sorted(observed_pairs_dict.items(), key=operator.itemgetter(1,0), reverse=True)

    h1, h2 = hc_merge_haplo(sorted_x, 
                                  c1_anchor, 
                                  c2_anchor, 
                                  idx_in_c1, idx_in_c2, 
                                  cluster1, cluster2, 
                                  sv_dict,
                                  phase_variants_location_list, 
                                  phase_variants_called_list, encoding_tb,
                                  heter, hratio, pv_dict, variant_num)

    return h1, h2


def hc_merge_haplo(sorted_x, c1_anchor, c2_anchor, idx_in_c1, idx_in_c2, cluster1, cluster2, sv_dict, phase_variants_location_list, reference, encoding_tb, heter, hratio, pv_dict, variant_num):
    """
    """

    if idx_in_c1 >= idx_in_c2:
        sys.exit()

    swap = lambda x, i, j: (x[0], x[1]) if i < j else (x[1], x[0])
    base, shift = (idx_in_c1, idx_in_c2)

    ## prioritize solutions

    exclude = dict()
    solutions = []
    heter_assumption = True

    ## two pair solutions first

    tmp_solutions = []
    tmp_read_count = dict()

    for idx, x in enumerate(sorted_x):
        for y in sorted_x[idx+1:]:
            x1, x2 = x[0].split("_")
            y1, y2 = y[0].split("_")
            if heter_assumption and x1 != y1 and x2 != y2:
                tmp_solutions.append([[[int(x1),int(x2)],x[1]],[[int(y1),int(y2)],y[1]]])
                tmp_read_count[len(tmp_read_count)] = x[1] + y[1]
                exclude[x] = 0
                exclude[y] = 0

    read_count_sorted = sorted(tmp_read_count.items(), key=operator.itemgetter(1))
    for key, value in read_count_sorted[::-1]:
        solutions.append(tmp_solutions[key])
  

    ## one pair solutions later

    tmp_solutions = []
    tmp_read_count = dict()

    for x in sorted_x:
        x1, x2 = x[0].split("_")
        if x not in exclude:
            tmp_solutions.append([[[int(x1),int(x2)], x[1]]])
            tmp_read_count[len(tmp_read_count)] = x[1]

    read_count_sorted = sorted(tmp_read_count.items(), key=operator.itemgetter(1))
    for key, value in read_count_sorted[::-1]:
        solutions.append(tmp_solutions[key])

    sol = solutions

    ## i think it's necessary
    result = []


    ## test

    if c1_anchor == None and c2_anchor == None:
        f_maj, f_min = list(map(lambda h: encoding_tb[h], reference[base][1:]))
        e_maj, e_min = list(map(lambda h: encoding_tb[h], reference[shift][1:]))
 
        for s in solutions:
            if len(s) > 1:
                result = [s[0][0], s[1][0]]
                break
            elif len(s) == 1:
                f_h1, e_h1 = s[0][0]
           
                ## alleles of sinlge pair solution must be in alleles of variants calling
                if f_h1 in (f_maj, f_min) and e_h1 in (e_maj, e_min):               
                    f_h2 = f_maj if f_h1 == f_min else f_min
                    e_h2 = e_maj if e_h1 == e_min else e_min
                else:
                    continue

                n1 = s[0][0]
                n2 = [f_h2, e_h2]
                result = [n1, n2] 
                break
            else:
                logging.error("merge error")
                sys.exit()

    elif c1_anchor != None and c2_anchor == None:

        e_maj, e_min = list(map(lambda h: encoding_tb[h], reference[shift][1:]))

        for s in solutions:

            if len(s) > 1:
                m1, m2 = s[0][0]
                m3, m4 = s[1][0]
                if m1 == c1_anchor[0] and m3 == c1_anchor[1]:
                    n1 = copy.copy(cluster1.h1)
                    n2 = copy.copy(cluster1.h2)
                    n1.append(m2)
                    n2.append(m4)
                    result = [n1, n2]
                    break
                elif m1 == c1_anchor[1] and m3 == c1_anchor[0]:
                    n1 = copy.copy(cluster1.h1)
                    n2 = copy.copy(cluster1.h2)
                    n1.append(m4)
                    n2.append(m2)
                    result = [n1, n2]
                    break

            elif len(s) == 1:
                m1, m2 = s[0][0]
                m4 = e_maj if m2 == e_min else e_min

                if m1 == c1_anchor[0]:
                    n1 = copy.copy(cluster1.h1)
                    n2 = copy.copy(cluster1.h2)
                    n1.append(m2)
                    n2.append(m4)
                    result = [n1, n2]
                    break
                elif m1 == c1_anchor[1]:
                    n1 = copy.copy(cluster1.h1)
                    n2 = copy.copy(cluster1.h2)
                    n1.append(m4)
                    n2.append(m2)
                    result = [n1, n2] 
                    break
            else:
                logging.error("merge error")
                sys.exit() 

    elif c1_anchor == None and c2_anchor != None:

        f_maj, f_min = list(map(lambda h: encoding_tb[h], reference[base][1:]))

        for s in solutions:
            if len(s) > 1:

                m1, m2 = s[0][0]
                m3, m4 = s[1][0]
                if m2 == c2_anchor[0] and m4 == c2_anchor[1]:
                    n1 = [m1] + cluster2.h1
                    n2 = [m3] + cluster2.h2
                    result = [n1, n2]
                    break
                elif m2 == c2_anchor[1] and m4 == c2_anchor[0]:
                    n1 = [m1] + cluster2.h2
                    n2 = [m3] + cluster2.h1
                    result = [n1, n2] 
                    break

            elif len(s) == 1:

                m1, m2 = s[0][0]
                m3 = f_maj if m1 == f_min else f_min

                if m2 == c2_anchor[0]:
                    n1 = [m1] + cluster2.h1
                    n2 = [m3] + cluster2.h2
                    result = [n1, n2]
                    break
                elif m2 == c2_anchor[1]:
                    n1 = [m1] + cluster2.h2
                    n2 = [m3] + cluster2.h1
                    result = [n1, n2] 
                    break

            else:
                logging.error("merge error")
                sys.exit() 

    elif c1_anchor != None and c2_anchor != None:

        c1_name_list = list(map(int,cluster1.name.split("_")))
        c2_name_list = list(map(int,cluster2.name.split("_")))

        seam_point_count = 0
        if len(c1_name_list) > len(c2_name_list):
            for c in c2_name_list:
                seam_point_count = seam_point_count + 1 if c-1 in c1_name_list else seam_point_count
                seam_point_count = seam_point_count + 1 if c+1 in c1_name_list else seam_point_count
        else:
            for c in c1_name_list:
                seam_point_count = seam_point_count + 1 if c-1 in c2_name_list else seam_point_count
                seam_point_count = seam_point_count + 1 if c+1 in c2_name_list else seam_point_count
            
 
        if seam_point_count > b_threshold and local_brute: 
            # first four for voting(sorted), last four for return (non-sorted)
            sol1_h1, sol1_h2, sol2_h1, sol2_h2, sorted_name , t_sol1_h1, t_sol1_h2, t_sol2_h1, t_sol2_h2 = build_two_possible_sol(cluster1, cluster2)

            penalty_1 = 0
            penalty_2 = 0

            if len(c1_name_list) > len(c2_name_list):
                for c2 in c2_name_list:
                    if (c2 - 1) in c1_name_list:
                        key = (c2 - 1) * variant_num + c2
                        if key in pv_dict:
                            pv = pv_dict[key]
                            penalty_1 += vote_haplotype(sol1_h1, sol1_h2, pv, c2 - 1, c2, sorted_name)
                            penalty_2 += vote_haplotype(sol2_h1, sol2_h2, pv, c2 - 1, c2, sorted_name)
                    if (c2 + 1) in c1_name_list:
                        key = c2 * variant_num + c2 + 1
                        if key in pv_dict:
                            pv = pv_dict[key]
                            penalty_1 += vote_haplotype(sol1_h1, sol1_h2, pv, c2, c2 + 1, sorted_name)
                            penalty_2 += vote_haplotype(sol2_h1, sol2_h2, pv, c2, c2 + 1, sorted_name) 
            else:
                for c1 in c1_name_list:
                    if c1-1 in c2_name_list:
                        key = (c1-1) * variant_num +c1
                        if key in pv_dict:
                            pv = pv_dict[key]
                            penalty_1 += vote_haplotype(sol1_h1, sol1_h2, pv, c1-1, c1, sorted_name)
                            penalty_2 += vote_haplotype(sol2_h1, sol2_h2, pv, c1-1, c1, sorted_name) 
                    if c1+1 in c2_name_list:
                        key = c1 * variant_num + c1 + 1
                        if key in pv_dict:
                            pv = pv_dict[key]
                            penalty_1 += vote_haplotype(sol1_h1, sol1_h2, pv, c1, c1+1, sorted_name)
                            penalty_2 += vote_haplotype(sol2_h1, sol2_h2, pv, c1, c1+1, sorted_name)

            '''
            penalty_1 = 0
            penalty_2 = 0

            print("d1, d2")
            for d1 in c1_name_list:
                for d2 in c2_name_list:
                    l, h = (d1, d2) if d1 < d2 else (d2, d1)
                    k = l * variant_num + h
                    if k in pv_dict:
                        pv = pv_dict[k]
                        penalty_1 += vote_haplotype(sol1_h1, sol1_h2, pv, l, h, sorted_name)
                        #print("=", l, h, penalty_1)
                        penalty_2 += vote_haplotype(sol2_h1, sol2_h2, pv, l, h, sorted_name)
                        #print("-", l, h, penalty_2)

            print("penalty_1:", penalty_1)
            print("penalty_2:", penalty_2)  
            '''

            if penalty_1 > penalty_2:
                result = [t_sol2_h1, t_sol2_h2]
            else:
                result = [t_sol1_h1, t_sol1_h2]
        else:
            for s in sol:
                if len(s) > 1:

                    m1, m2 = s[0][0]
                    m3, m4 = s[1][0]

                    if [m1, m3] == c1_anchor and [m2, m4] == c2_anchor:
                        n1 = cluster1.h1 + cluster2.h1
                        n2 = cluster1.h2 + cluster2.h2
                        result = [n1, n2]
                        break 
                    elif [m3, m1] == c1_anchor and [m2, m4] == c2_anchor:
                        n1 = cluster1.h1 + cluster2.h2
                        n2 = cluster1.h2 + cluster2.h1
                        result = [n1, n2]
                        break 
                    elif [m1, m3] == c1_anchor and [m4, m2] == c2_anchor:
                        n1 = cluster1.h1 + cluster2.h2
                        n2 = cluster1.h2 + cluster2.h1
                        result = [n1, n2]
                        break 
                    elif [m3, m1] == c1_anchor and [m4, m2] == c2_anchor:
                        n1 = cluster1.h1 + cluster2.h1
                        n2 = cluster1.h2 + cluster2.h2
                        result = [n1, n2]
                        break 

                elif len(s) == 1:

                    m1, m2 = s[0][0]
                    if m1 == c1_anchor[0] and m2 == c2_anchor[0]:
                        n1 = cluster1.h1 + cluster2.h1
                        n2 = cluster1.h2 + cluster2.h2
                        result = [n1, n2]
                        break
                    elif m1 == c1_anchor[0] and m2 == c2_anchor[1]:
                        n1 = cluster1.h1 + cluster2.h2
                        n2 = cluster1.h2 + cluster2.h1
                        result = [n1, n2]
                        break
                    elif m1 == c1_anchor[1] and m2 == c2_anchor[0]:
                        n1 = cluster1.h1 + cluster2.h2
                        n2 = cluster1.h2 + cluster2.h1
                        result = [n1, n2]
                        break
                    elif m1 == c1_anchor[1] and m2 == c2_anchor[1]:
                        n1 = cluster1.h1 + cluster2.h1
                        n2 = cluster1.h2 + cluster2.h2
                        result = [n1, n2]
                        break

                else:
                    logging.error("merge error")
                    sys.exit()                    

    else:
        logging.error("merge error")
        sys.exit()

    # Merge Failure
    if(len(result) < 2):
        logging.warning("Merge Fail, Can't find two haplo that connect both sides " + str(idx_in_c1) + ' ' + str(idx_in_c2))
        return None, None

    return [r for r in result[:2]]


def vote_haplotype(h1, h2, pv, sidx, bidx, sorted_name):
    """
    """
    sol1 = (h1[sorted_name.index(sidx)], h1[sorted_name.index(bidx)])
    sol2 = (h2[sorted_name.index(sidx)], h2[sorted_name.index(bidx)])

    penalty = 0
    for k, v in pv.items():
        m1, m2 = list(map(int,k.split('_')))
        if (m1, m2) == sol1 or (m1, m2) == sol2:
            pass
        elif m1 in sol1 or m1 in sol2 or m2 in sol1 or m2 in sol2:
            penalty = penalty + v
        else:
            penalty = penalty + 2*v
    return penalty


def build_two_possible_sol(cluster1, cluster2):
    template_name = "_".join([cluster1.name, cluster2.name])
    template_name = list(map(int, template_name.split("_")))
    template_sol1_h1 = cluster1.h1 + cluster2.h1
    template_sol1_h2 = cluster1.h2 + cluster2.h2
    template_sol2_h1 = cluster1.h1 + cluster2.h2
    template_sol2_h2 = cluster1.h2 + cluster2.h1

    sol1_h1 = list()
    sol1_h2 = list()
    sol2_h1 = list()
    sol2_h2 = list()

    sorted_name = sorted(template_name)

    for s in sorted_name:
        idx = template_name.index(s)
        sol1_h1.append(template_sol1_h1[idx])
        sol1_h2.append(template_sol1_h2[idx])
        sol2_h1.append(template_sol2_h1[idx])
        sol2_h2.append(template_sol2_h2[idx])

    return sol1_h1, sol1_h2, sol2_h1, sol2_h2, sorted_name, template_sol1_h1, template_sol1_h2, template_sol2_h1, template_sol2_h2





def hc_merge(vars_cddt_mp, vars_cddt_pl, pv_dict, score_matrix, sv_dict, phase_variants_location_list, phase_variants_called_list, encoding_tb, heter, hratio, timer):
    """
    Args:
    Return:
    """
    ##
    # score_matrix_reduced : score matrix of segment relationship, 
    #                     it reduced along mergence process.
    # score_matrix_lookup     : copy of original matrix. change value to np.NINF,
    #                     when the pair merge fail in process 
    ##
    timer.start('copy_matrix_1')
    score_matrix_reduced = copy_score_matrix(score_matrix)
    score_matrix_lookup  = copy_score_matrix(score_matrix)
    timer.stop('copy_matrix_1')
    variant_num = len(score_matrix)

    stop_flag = False

    ## evaluation
    loss_list = []
    threshold_list = []
    connected_list = []
    observed_list = []

    while not stop_flag:
        timer.start('get_maxscore')
        base, shift, score = get_maxscore(score_matrix_reduced)
        if(score == np.NINF):
            break

        ## base < shift

        group1_id = vars_cddt_mp[base]
        group2_id = vars_cddt_mp[shift]
        group1_num = len(group1_id.split("_"))
        group2_num = len(group2_id.split("_"))
        timer.stop('get_maxscore')

        ##
        # Try each combination from high score to low score.
        # If find non-conflict major and minor haplotypes, then break the loop
        # else find next until no possible pair.
        ##
        for rank in range(0, group1_num * group2_num):
            timer.start('get_maxscore_between_segments') 
            idx_in_c1, idx_in_c2, score = get_maxscore_between_segments(group1_id, group2_id, score_matrix_lookup, rank)
            timer.stop('get_maxscore_between_segments')

            if idx_in_c1 < idx_in_c2:
                cluster1 = vars_cddt_pl[group1_id]
                cluster2 = vars_cddt_pl[group2_id]
            else:
                cluster1 = vars_cddt_pl[group2_id]
                cluster2 = vars_cddt_pl[group1_id]
                idx_in_c2, idx_in_c1 = (idx_in_c1, idx_in_c2)
                

            pv_key = get_matrix_encode_key(idx_in_c1, idx_in_c2, variant_num)

            # no connected-pair exist
            if(score == np.NINF):
                break

            timer.start('choice_haplos')
            if pv_key in pv_dict:
                observed_pairs_dict = pv_dict[pv_key]
           
                ## switch clusters

                #cluster1 = vars_cddt_pl[group1_id]
                #cluster2 = vars_cddt_pl[group2_id]
                node = hc_create_parent_node(cluster1, cluster2)
                node.h1, node.h2 = choice_haplos(observed_pairs_dict, cluster1, cluster2, idx_in_c1, idx_in_c2, sv_dict, phase_variants_location_list, phase_variants_called_list, encoding_tb, heter, hratio, pv_dict, variant_num)

                #careful for end condition? only one None?
                if(node.h1 != None and node.h2 != None):
                    break
      
            else:
                logging.error("can't find in pv_dict")
                sys.exit()

            timer.stop('choice_haplos')
        ##
        # If two segment can't merge together in previous step,
        # remove it from lookup matrix, then this pair wouldn't be considered again.
        # remove it from reduced matrix. (segment-wised)
        ##
        timer.start('rearrange')
        if(node.h1 == None and node.h2 == None):
            score_matrix_lookup[idx_in_c1][idx_in_c2] = np.NINF
            score_matrix_lookup[idx_in_c2][idx_in_c1] = np.NINF
            score_matrix_reduced[base][shift] = np.NINF
            continue

        '''
        ## seems useless?
        if cluster1.h1 == None:
            cluster1.h1, cluster1.h2 = (node.h1[0], node.h2[0])
        if cluster2.h1 == None:
            cluster2.h1, cluster2.h2 = (node.h1[1], node.h2[1])
        '''

        ##
        #  re-arrange?
        ##
        t = node.name.split("_")
        tmp = sorted(map(int, node.name.split("_")))
        tmp_name = "_".join(map(str, tmp))
        h1 = list()
        h2 = list()
        for i in tmp:
            h1.append(node.h1[t.index(str(i))])
            h2.append(node.h2[t.index(str(i))])

        # ? node.name should always be in pl?
        #delete old name of the node
        if node.name in vars_cddt_pl: 
            del vars_cddt_pl[node.name]
        node.h1 = h1
        node.h2 = h2
        node.name = tmp_name

        ##
        # Put new node into pool, and remove old two from it.
        ##
        vars_cddt_pl[node.name] = node
        del vars_cddt_pl[cluster1.name]
        del vars_cddt_pl[cluster2.name]

        '''
        print("start checking")
        loss = checking_multi_answer(vars_cddt_pl, self.ans_dict, self.code_tb)
        print("loss_total ", loss, score, observed_pairs_dict)
        loss_list.append(loss)
        threshold_list.append(score)
        connected_list.append(str(idx_in_c1)+"_"+str(idx_in_c2))
        observed_list.append(pv_dict[pv_key])
        '''

        ##
        # update mapping (var_lst)
        # keep and upate group1 entry as new segment, remove group2.
        ##
        #vars_cddt_mp[base] = '_'.join([group1_id, group2_id])
        #del vars_cddt_mp[shift]

        vars_cddt_mp[base] = tmp_name
        del vars_cddt_mp[shift]

        timer.stop('rearrange')
        ##
        # Reduce process, keep base, remove shift
        ##
        timer.start('reduce_matrix')
        score_matrix_reduced = reduce_matrix(score_matrix_reduced, base, shift)
        timer.stop('reduce_matrix')

    return loss_list, threshold_list, connected_list, observed_list


def save_phasing(chrom, vars_cddt_pl, encoding_tb, phase_variants_location_list, output_dict):
    """
    """

    for k, v in vars_cddt_pl.items():

        cluster_name = k
        cluster_node = v
       
        cluster_name_seq = list(map(int, cluster_name.split("_")))
        cluster_name_seq_sorted = sorted(cluster_name_seq)

        h1 = cluster_node.h1
        h2 = cluster_node.h2

        # single variant
        if h1 == None or h2 == None:
            continue

        ps_id = phase_variants_location_list[cluster_name_seq_sorted[0]]

        for s in cluster_name_seq_sorted:
            location = str(phase_variants_location_list[s])
            idx = cluster_name_seq.index(s)
            h1_allele = h1[idx]
            h2_allele = h2[idx]

            output_dict['_'.join([location])] = (ps_id, h1_allele, h2_allele)

    return output_dict

def rebuild_haplo(chrom, vars_cddt_pl, encoding_tb, phase_variants_location_list, vars_dict, sv_dict):
    """
    """
    print('======', len(vars_cddt_pl), '=====')

    print(vars_cddt_pl)

    sorted_dict = dict()
    for key, haplo in vars_cddt_pl.items():
        nkey = int(key.split('_')[0])
        sorted_dict[nkey] = (key, haplo)

    import collections
    import operator

    ps_name = None
    a = collections.OrderedDict(sorted(sorted_dict.items()))
    for key, content in a.items():
        print(key, content)
        key_o = content[0].split('_')
        key_s = sorted(map(int, key_o))
        haplo = content[1]

        if(haplo.h1 == None):
            '''
            key_out = list(map(lambda y: phase_variants_location_list[y], key_s))
            print('\t'.join(map(str,[key_s[0], key_out[0], 'N', 'N'])))
            print('==============')
            '''
            continue

        haplo1 = list(map(lambda h: encoding_tb[int(h)], haplo.h1))
        haplo2 = list(map(lambda h: encoding_tb[int(h)], haplo.h2))
        haplo1_out = []
        haplo2_out = []
        for i in key_s:
            haplo1_out.append(haplo1[key_o.index(str(i))])
            haplo2_out.append(haplo2[key_o.index(str(i))])

        key_out = list(map(lambda y: phase_variants_location_list[y], key_s))

        if(ps_name == None):
            ps_name = key_out[0]

        for i in range(len(key_s)):
            
            b = sv_dict[key_out[i]]
            h = "homo " if haplo1_out[i] == haplo2_out[i] else "heter"
            s = [(encoding_tb[k], v) for k, v  in b.items()]
            print('\t'.join(map(str,[key_s[i], ps_name, key_out[i], haplo1_out[i], haplo2_out[i], h, s])))

        ps_name = None


