__author__ = 'ifish'

import logging
from collections import OrderedDict
import pysam
import numpy as np
from itertools import chain
from scipy.sparse import csr_matrix, lil_matrix
from HAHap.blocks import trans_cigartuples_to_region, trans_loc_to_cigaridx, trans_loc_to_alignidx

class InputVcfReader(object):
    """
    """
    coding = {'A': 1, 'T': 2, 'C': 3, 'G': 4, 1: 'A', 2: 'T', 3: 'C', 4: 'G', 0: ' ', ' ': 0, 6: '-', '-': 6}

    def __init__(self, bam_file, min_read_quality,
                 indels=False, trust=True):
        self.bam_file = bam_file
        self.indels = indels
        self.min_read_quality = min_read_quality
        self.trust = trust
        self.reads_proxy = None
        self.logging = logging


    def get_sf_mx(self, chrom, phase_loc, phase_allele, timer):
        """"""

        fragments = OrderedDict()

        if self.reads_proxy is None:
            self.reads_proxy = pysam.AlignmentFile(self.bam_file, "rb")

        timer.start('021._read_variance_faster') 
        self._read_variance_faster(chrom, fragments, phase_loc, phase_allele, self.trust, timer)
        timer.stop('021._read_variance_faster')

        timer.start('022.clear_matrix')
        self.clear_matrix(fragments) #?
        timer.stop('022.clear_matrix')

        sf_mx = lil_matrix((len(fragments), len(phase_loc)), dtype=np.int8)

        timer.start('023.produce_var_matrix')
        fragment_se = self._fill_sf_mx(sf_mx, fragments, phase_loc, self.coding)
        timer.stop('023.produce_var_matrix')

        return sf_mx, fragments, fragment_se, self.coding

    @staticmethod
    def clear_matrix(fragments):
        """"""
        noninfo = []
        for x, y in fragments.items():
            a = [z[2] for z in y]
            if '-' in a:
                a.remove('-')
            if len(a) < 2:
                noninfo.append(x)

        for x in noninfo:
            del fragments[x]

    def _read_variance_faster(self, chrom, fragments, phase_loc, phase_allele, trust, timer):
        """"""

        # 0-base, end value not included
        f_loc = phase_loc[0] - 1
        l_loc = phase_loc[-1]

        visited_set = set()
        removed_set = set()

        #print("==", chrom, f_loc, l_loc)
        timer.start('024.fetch')
        for read in self.reads_proxy.fetch(chrom, f_loc, l_loc):
            timer.stop('024.fetch')

            timer.start('024.if')
            if read.mapping_quality < self.min_read_quality:
                continue

            if not read.is_proper_pair:
                continue
            timer.stop('024.if')

            timer.start('024.id')
            read_id = read.query_name
            r_st = read.reference_start
            # read.reference_end does not involved
            r_ed = read.reference_end - 1
            r_cigar = read.cigartuples

            '''
            if read.is_read1:
                read_id = read.query_name + "_" + str(read.next_reference_start)
            else:
                read_id = read.query_name + "_" + str(read.reference_start)
            '''
            timer.stop('024.id') 

            timer.start('024.align0')
            allele = ''
            #aligned_pairs = read.get_aligned_pairs()
            timer.stop('024.align0')
            timer.start('024.align1')
            #aligned_pairs = read.get_aligned_pairs(matches_only=True)
            timer.stop('024.align1')
            timer.start('024.align2')
            #aligned = [i[1] for i in aligned_pairs]
            timer.stop('024.align2')
            timer.start('024.align3')
            #region = trans_cigartuples_to_region(r_st, read.cigartuples)   
            timer.stop('024.align3')
            timer.start('024.align4')
            #chn = [i for j in [range(r[0], r[1]+1) for r in region] for i in j]
            timer.stop('024.align4')
            timer.start('024.align5')
            chn2 = range(r_st, r_ed + 1)
            timer.stop('024.align5')


            timer.start('024.loop')
            for var_idx, v in enumerate(phase_loc):
                v_0_base = v - 1

                if v_0_base < r_st:
                    continue
                if v_0_base > r_ed:
                    break

                #if v_0_base in chn:
                #if v_0_base in range(read.reference_start, read.reference_end):
                if v_0_base in chn2:
                    timer.start('024.loop_pre')
                    #t = trans_loc_to_cigaridx(v_0_base, region)
                    r = trans_loc_to_alignidx(v_0_base, r_st, read.cigartuples)
                    if r is None:
                        observed = '-'
                    else:
                        observed = read.query_sequence[r]
                    '''
                    if t is None and r is not None:
                        print("---")
                        print(v_0_base, t, r)
                        print(aligned_pairs)
                        print(read.get_aligned_pairs())
                        print(read.cigartuples)
                        print("---")
                    if t is None:
                        continue
                    aligned_idx = aligned.index(v_0_base)
                    #if t != aligned_idx:
                    if r != aligned_pairs[t][0]:
                        print(v_0_base)
                        print(t, aligned_idx, aligned_pairs[t][0], r)
                        print(aligned_pairs)
                        print(read.get_aligned_pairs())
                    observed = read.query_sequence[aligned_pairs[t][0]]

                    if v_0_base not in aligned:
                       continue
                    aligned_idx = aligned.index(v_0_base)

                    if aligned_pairs[aligned_idx][0] is None:
                        # gap in read
                        observed = '-'
                    else:
                        observed = read.seq[aligned_pairs[aligned_idx][0]]
                    '''

                    if observed not in phase_allele[var_idx] and trust:
                        allele = '-'
                    else:
                        allele = observed

                    if read_id in visited_set:
                        is_pair_read = True
                    else:
                        is_pair_read = False
                        visited_set.add(read_id)

                    timer.stop('024.loop_pre')
                    timer.start('024.overlap')
                    if read_id in fragments:
                        if not is_pair_read:
                            fragments[read_id].append((chrom, v, allele, ''))
                        else:
                            fragment = fragments[read_id]
                            exist = [i[1] for i in fragment]
                            observeds = [i[2] for i in fragment]
                            if v in exist and allele != observeds[exist.index(v)]:
                                removed_set.add(read_id)
                                break
                            elif v in exist:
                                continue
                            else:
                                fragment.append((chrom, v, allele, '')) 
                    else:
                        fragments[read_id] = [(chrom, v, allele, '')]
                    timer.stop('024.overlap')

            timer.stop('024.loop')
            timer.start('024.fetch')
        timer.stop('024.fetch')
        timer.start('025.removed_set')

        for i in removed_set:
            del fragments[i]
        timer.stop('025.removed_set')


        return


    @staticmethod
    def _fill_sf_mx(sf_mx, fragments, phase_loc, coding):
        """"""
        # vertical
        fragment_se = []
        v_idx = 0
        code_num = len(coding)/2

        loc_idx_dict = {y: x for x, y in dict(enumerate(phase_loc, 0)).items()}

        for read_id, variances in fragments.items():
            min_loc = loc_idx_dict[variances[0][1]]
            max_loc = loc_idx_dict[variances[-1][1]]

            fragment_se.append((read_id, min_loc, max_loc))
            for chrom, loci, allele, atype in variances:
                h_idx = loc_idx_dict[loci]
                if allele in coding:
                    code = coding[allele]
                else:
                    code = code_num
                    code_num += 1
                    coding[allele] = code
                    coding[code] = allele
                sf_mx[v_idx, h_idx] = code
            v_idx += 1

        return fragment_se
