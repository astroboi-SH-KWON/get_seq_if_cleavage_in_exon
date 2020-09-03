import glob
from Bio import SeqIO
import openpyxl
import os

class Utils:
    def __init__(self):
        self.ext_txt = ".txt"
        self.ext_dat = ".dat"
        self.ext_xlsx = ".xlsx"

    def get_seq_record_from_genbank(self, path):
        return SeqIO.read(path, "genbank")

