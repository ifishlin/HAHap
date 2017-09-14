#!/usr/bin/env python3
import sys

def print_var_matrix(var_matrix, v_matrix_ids, encoding_table, prange):
    """
    """

    max_list = []

    # max len in each variance
    if var_matrix.shape[0] == 0:
        sys.exit()
    for j in range(var_matrix.shape[1]):
        max_var = 0
        for i in range(var_matrix.shape[0]):
            max_var = len(encoding_table[var_matrix[i, ]]) \
                if len(encoding_table[var_matrix[i, j]]) > max_var else max_var
        max_list.append(max_var)

    head = '{0:>48s}-'.format('read_id \\ site idx')
    sys.stdout.write(head)
    idx = 1
    for index, j in enumerate(max_list):
        sys.stdout.write(addpadding(str(idx), j)+' ')
        idx += 1
        idx %= 10
    print()

    for i in range(var_matrix.shape[0]):
        head = '{0:>48s}-'.format(v_matrix_ids[i][0])
        sys.stdout.write(head)
        if prange == ':':
            for j in range(var_matrix.shape[1]):
                sys.stdout.write(addpadding(encoding_table[var_matrix[i, j]], max_list[j]) + ',')
        else:
            start, end = map(int, prange.split(":"))
            for j in range(start, end):
                sys.stdout.write(addpadding(encoding_table[var_matrix[i, j]], max_list[j]) + ',')
        print()

    return


def addpadding(seq, l):
    if l > len(seq):
        for i in range(l-len(seq)):
            seq += ' '
    return seq
