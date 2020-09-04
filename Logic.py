from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import blosum62

import LogicPrep
class Logics:
    def __init__(self):
        pass

    """
    checkSeqByChar : match sequences by char
    :param
        seq_char :
        target_char : 
    :return
        boolean
    """
    def checkSeqByChar(self, seq_char, target_char):
        flag = False
        if target_char == 'N':
            return True
        elif target_char in 'ACGTU':
            if seq_char == target_char:
                return True
        elif target_char == 'R':
            if seq_char in 'AG':
                return True
        """
        add more rules of "ACTGU"
        """

        return flag

    """
    match : match sequence with "same length" strings
    :param
        i : index of seq
        seq_str : targeted DNA/RNA sequence 
        rule_str : rules with "ACGTU", "N", "R",...
    :return
        boolean
    """
    def match(self,i, seq_str, rule_str):
        if len(seq_str) == i:
            return True
        if self.checkSeqByChar(seq_str[i], rule_str[i]):
            return self.match(i + 1, seq_str, rule_str)
        else:
            return False

    def get_idx_of_matching_seq(self, full_seq, trgt_seq):
        plus_strand_list = []
        minus_strand_list = []

        minus_seq = full_seq.complement()

        for idx in range(len(full_seq) - len(trgt_seq)):

            plus_seq_str = str(full_seq[idx: idx + len(trgt_seq)])
            minus_seq_str = str(minus_seq[idx: idx + len(trgt_seq)])

            if self.match(0, plus_seq_str, trgt_seq):
                # add the first index of trgt_seq
                plus_strand_list.append(idx)

            if self.match(0, minus_seq_str, trgt_seq[::-1]):
                # add the last index of trgt_seq
                minus_strand_list.append(idx + len(trgt_seq))

        return plus_strand_list, minus_strand_list

    def get_idx_in_list(self, idx_list, trgt_list, plus=True, cleave=3):
        result_list = []
        for tmp_idx in idx_list:
            idx = tmp_idx
            if plus:
                idx -= cleave
            else:
                idx += cleave
            for trgt_arr in trgt_list:
                if idx in trgt_arr:
                    result_list.append(tmp_idx)

        return result_list
    def is_too_short(self, trgt_seq, len_trgt):
        if len(trgt_seq) < len_trgt:
            return True
        return False

    def get_trgt_seq_in_idx_list(self, full_seq, idx_list, init, plus=True):
        result_list = []
        len_spacer = init[0]
        cleavage = init[1]
        pam_seq = init[2]
        len_pam = len(pam_seq)
        len_aftr_pam = init[3]

        for idx in idx_list:
            if plus:
                trgt_seq = str(full_seq[idx - len_spacer: idx + len_pam + len_aftr_pam])

                if self.is_too_short(trgt_seq, len_spacer + len_pam + len_aftr_pam):
                    continue

                result_list.append([trgt_seq, "+"])
            else:
                trgt_seq = str(full_seq.complement()[idx - len_pam - len_aftr_pam: idx + len_spacer])

                if self.is_too_short(trgt_seq, len_spacer + len_pam + len_aftr_pam):
                    continue

                result_list.append([trgt_seq, "-"])

        return result_list





