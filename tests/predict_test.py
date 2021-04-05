# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 09:52:17 2021

@author: ananthansadagopan
"""

import os
import pandas as pd
from neop.main import predict

path = os.getcwd()
# file Upload
OUTPUT_FOLDER = os.path.join(path, 'output_files')

if not os.path.isdir(OUTPUT_FOLDER):
    os.mkdir(OUTPUT_FOLDER)

maf_filename_input = os.path.join(path, "test_maf.maf")

allele_annot = os.path.join(path, "allele_ref.csv")

output_json, output_csv = predict(allele_annot, maf_filename_input, OUTPUT_FOLDER, "predict_test")