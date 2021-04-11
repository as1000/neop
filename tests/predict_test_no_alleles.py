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

with open(os.path.join(OUTPUT_FOLDER, "predict_test_no_alleles_output.csv"), "r") as f:
    original_csv_content = f.read()

with open(os.path.join(OUTPUT_FOLDER, "predict_test_no_alleles_output.json"), "r") as f:
    original_json_content = f.read()

output_json, output_csv = predict("", maf_filename_input, OUTPUT_FOLDER, "predict_test_no_alleles")

with open(output_csv, "r") as f:
    generated_csv_content = f.read()

with open(output_json, "r") as f:
    generated_json_content = f.read()

assert original_csv_content == generated_csv_content and original_json_content == generated_json_content
