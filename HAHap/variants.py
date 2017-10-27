__author__ = 'ifish'

import logging
from collections import OrderedDict
import pysam
import numpy as np
from scipy.sparse import csr_matrix, lil_matrix


class InputVcfReader(object):
    """
    """
    coding = {'A': 1, 'T': 2, 'C': 3, 'G': 4, 1: 'A', 2: 'T', 3: 'C', 4: 'G', 0: ' ', ' ': 0, 6: '-', '-': 6}

    def __init__(self, bam_file, min_base_quality, min_read_quality,
                 indels=False, trust=True):
        self.bam_file = bam_file
        self.indels = indels
        self.min_base_quality = min_base_quality
        self.min_read_quality = min_read_quality
        self.trust = trust
        self.reads_proxy = None
        self.logging = logging

    # @remove 
    def reset(self, bam_file, min_base_quality, min_read_quality, indels=False,
              trust=True):
        self.bam_file = bam_file
        self.indels = indels
        self.min_base_quality = min_base_quality
        self.min_read_quality = min_read_quality
        self.trust = trust
        self.logging = logging

    def get_sf_mx(self, chrom, phase_loc, phase_allele, timer, f):
        """
        Args:

        Return:
        """

        fragments = OrderedDict()

        '''
        if self.reads_proxy is None:
            self.reads_proxy = pysam.AlignmentFile(self.bam_file, "rb")

        timer.start('021._read_variance_faster') 
        self._read_variance_faster(chrom, fragments, phase_loc, phase_allele, self.trust, timer)
        timer.stop('021._read_variance_faster')

        timer.start('022.clear_matrix')
        self.clear_matrix(fragments) #?
        timer.stop('022.clear_matrix')
          
        print(fragments)
        print(f)

        for k, v in fragments.items():
            if k in f and v != f[k]:
                print("old:",v)
                print("new:",f[k])
                f[k] = v

        print("====")
        for k, v in fragments.items():
            if k in f and v != f[k]:
                print("old:",v)
                print("new:",f[k])
        '''
        fragments = f
        # test
        sf_mx = lil_matrix((len(fragments), len(phase_loc)), dtype=np.int8)
        fragment_se = self._fill_sf_mx(sf_mx, fragments, phase_loc, self.coding)
        # old
        '''
        sf_mx = lil_matrix((len(fragments), len(phase_loc)), dtype=np.int8)

        timer.start('023.produce_var_matrix')
        fragment_se = self._fill_sf_mx(sf_mx, fragments, phase_loc, self.coding)
        timer.stop('023.produce_var_matrix')
        '''
        
        #fragment_se3 = sorted(fragment_se2, key=lambda tup: tup[1])
        #print(fragment_se3)
        #return sf_mx2, f, fragment_se2, self.coding
        return sf_mx, fragments, fragment_se, self.coding

    # can we don't read this in
    @staticmethod
    def clear_matrix(fragments):
        """
        Remove unqualified read from v_desc_dict. (1). remove reads that cover less than two variants.
        """
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
        """
        Read raw NGS data of each variants from sam file.

        """

        # 0-base, end value not included
        f_loc = phase_loc[0] - 1
        l_loc = phase_loc[-1]

        visited_set = set()
        removed_set = set()

        timer.start('024.fetch')
        for read in self.reads_proxy.fetch(chrom, f_loc, l_loc):
            timer.stop('024.fetch')
            if read.mapping_quality < self.min_read_quality:
                continue

            if not read.is_proper_pair:
                continue

            if read.is_read1:
                read_id = read.query_name + "_" + str(read.next_reference_start)
            else:
                read_id = read.query_name + "_" + str(read.reference_start)

            if read_id == 'HISEQ1:93:H2YHMBCXX:1:1108:4879:55953_81952627':
                print("***", read_id)
 
            allele = ''
            aligned_pairs = read.get_aligned_pairs()
            aligned = [i[1] for i in aligned_pairs]

            for var_idx, v in enumerate(phase_loc):
                v_0_base = v-1
                # switch to binary search
                if v_0_base > read.reference_end - 1:  # one past last alignment
                    break
                if v_0_base < read.reference_start:
                    continue

                if v_0_base in range(read.reference_start, read.reference_end):
                    aligned_idx = aligned.index(v_0_base)
                    if aligned_pairs[aligned_idx][0] is None:
                        # gap in read
                        observed = '-'
                    else:
                        observed = read.seq[aligned_pairs[aligned_idx][0]]

                    if observed not in phase_allele[var_idx] and trust:
                        print("none trust", read_id)
                        allele = '-'
                    else:
                        allele = observed

                    if read_id in visited_set:
                        is_pair_read = True
                    else:
                        is_pair_read = False
                        visited_set.add(read_id)

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

                    '''
                    if(read_id in fragments):
                        exist     = list(map(lambda i:i[1],fragments[read_id]))
                        observeds = list(map(lambda i:i[2],fragments[read_id]))
                        if v in exist and allele != observeds[exist.index(v)]:
                            removed_set.add(read_id) 
                        elif v in exist:
                            continue
                        else:
                            fragments[read_id].append((chrom, v, allele, ''))
                    else:
                        fragments[read_id] = [(chrom, v, allele, '')]
                    timer.stop('fragments')
                    '''

            timer.start('024.fetch')
        timer.stop('024.fetch')
        timer.start('025.removed_set')

        print("** remvoe", removed_set)
        for i in removed_set:
            del fragments[i]
        timer.stop('025.removed_set')

        return


    @staticmethod
    def _fill_sf_mx(sf_mx, fragments, phase_loc, coding):
        """

        """
        # vertical
        fragment_se = []
        v_idx = 0
        code_num = len(coding)/2

        loc_idx_dict = {y: x for x, y in dict(enumerate(phase_loc, 0)).items()}

        for read_id, variances in fragments.items():

            min_loc = -1
            max_loc = -1

            for i in range(len(variances)):
                if variances[i][1] in loc_idx_dict:
                    min_loc = loc_idx_dict[variances[i][1]]

            for i in range(len(variances) -1, -1, -1):
                if variances[i][1] in loc_idx_dict:
                    max_loc = loc_idx_dict[variances[i][1]]


            max_loc = loc_idx_dict[variances[-1][1]]

            for chrom, loci, allele, atype in variances:
                if loci not in phase_loc:
                    continue

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

            fragment_se.append((read_id, min_loc, max_loc))

        return fragment_se
    '''
    @staticmethod
    def _fill_sf_mx(sf_mx, fragments, phase_loc, coding):
        """

        """
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
     '''
