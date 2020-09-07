
class LogicPreps:

    def __init__(self):
        self.ext_fa = ".fa"
        self.ext_dat = ".dat"
        self.ext_gtf = ".gtf"

    def get_cds_idx_arr_to_list(self, seq_rec):
        result_list = []

        if seq_rec.features:
            for feature in seq_rec.features:
                if feature.type == "CDS":
                    result_list.append([idx for idx in feature.location])

                # add type 'mRNA' 20200907
                if feature.type == "mRNA":
                    result_list.append([idx for idx in feature.location])

                    # TODO make dict by transcript_id
                    # print("feature.qualifiers : ", feature.qualifiers)
                    # result_dict.update({feature.qualifiers["protein_id"][0]: [idx for idx in feature.location]})
        return result_list

    def filter_out_dupl(self, idx_list):
        result_list = []
        result_set = set()

        for idx in idx_list:
            result_set.add(idx)

        for flt_idx in result_set:
            result_list.append(flt_idx)

        return result_list

    def merge_list(self, arr_list):
        result_list = []
        for tmp_arr in arr_list:
            result_list.extend(tmp_arr)
        return result_list

    def sort_list_by_ele(self, data_list, ele_idx, up_down_flag=True):
        result_list = []
        for tmp_arr in sorted(data_list, key=lambda tmp_arr: tmp_arr[ele_idx], reverse=up_down_flag):
            result_list.append(tmp_arr)
        return result_list