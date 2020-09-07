import time
import os

import Util
import Logic
import LogicPrep
############### start to set env ################
WORK_DIR = os.getcwd() + "/"
PROJECT_NAME = WORK_DIR.split("/")[-2]
NCBI = "NCBI/"
FILE_NAME = ["Apob_NCBI", "Pcsk9_NCBI", "Hpd_NCBI"]
LEN_SPACER = 43
CLEAVAGE = 3
PAM = 'NNGRRT'
LEN_AFTR_PAM = 3
INIT_0 = [LEN_SPACER, CLEAVAGE, PAM, LEN_AFTR_PAM]
INIT_1 = [41, CLEAVAGE, "NNRGAA", LEN_AFTR_PAM]
INIT_RULE = [INIT_0, INIT_0, INIT_1]
############### end setting env #################

def main():
    util = Util.Utils()
    logic_prep = LogicPrep.LogicPreps()
    logic = Logic.Logics()

    for file_idx in range(len(INIT_RULE)):
        seq_record = util.get_seq_record_from_genbank(WORK_DIR + NCBI + FILE_NAME[file_idx] + ".gb")
        cds_idx_list = logic_prep.get_cds_idx_arr_to_list(seq_record)
        # add type 'mRNA' 20200907
        cds_idx_list.extend(logic_prep.get_cds_idx_arr_to_list(seq_record, 'mRNA'))

        init_rule = INIT_RULE[file_idx]
        pam_seq = init_rule[2]
        plus_strand_list, minus_strand_list = logic.get_idx_of_matching_seq(seq_record.seq, pam_seq)

        plus_idx_list = logic.get_idx_in_list(plus_strand_list, cds_idx_list)
        minus_idx_list = logic.get_idx_in_list(minus_strand_list, cds_idx_list, False)

        filtered_plus_idx_list = logic_prep.filter_out_dupl(plus_idx_list)
        filtered_minus_idx_list = logic_prep.filter_out_dupl(minus_idx_list)

        plus_seq_list = logic.get_trgt_seq_in_idx_list(seq_record.seq, filtered_plus_idx_list, init_rule)
        minus_seq_list = logic.get_trgt_seq_in_idx_list(seq_record.seq, filtered_minus_idx_list, init_rule, False)

        merge_list = logic_prep.merge_list([plus_seq_list, minus_seq_list])
        tot_list = logic_prep.sort_list_by_ele(merge_list, 0)

        header = ["sequence", "strand"]
        util.make_excel(WORK_DIR + "output/result_mRNA_" + FILE_NAME[file_idx], header, tot_list)
        # util.make_excel(WORK_DIR + "output/result_" + FILE_NAME[file_idx], header, tot_list)


if __name__ == '__main__':
    start_time = time.perf_counter()
    print("start [ " + PROJECT_NAME + " ]>>>>>>>>>>>>>>>>>>")
    main()
    print("::::::::::: %.2f seconds ::::::::::::::" % (time.perf_counter() - start_time))
