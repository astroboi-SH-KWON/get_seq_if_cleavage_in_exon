
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

                    # TODO make dict by transcript_id
                    # print("feature.qualifiers : ", feature.qualifiers)
                    # result_dict.update({feature.qualifiers["protein_id"][0]: [idx for idx in feature.location]})
        return result_list