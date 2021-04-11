# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 09:52:17 2021

@author: ananthansadagopan
"""

import os
from neop.main import annotate_maf

path = os.getcwd()
# file Upload
OUTPUT_FOLDER = os.path.join(path, 'output_files')

if not os.path.isdir(OUTPUT_FOLDER):
    os.mkdir(OUTPUT_FOLDER)

maf_filename_input = os.path.join(path, "test_maf.maf")
log_file_input = os.path.join(path, "log_file.log")

with open(os.path.join(OUTPUT_FOLDER, "annotate_maf_test.maf"), "r") as f:
    original_file_content = f.read()

output_maf = annotate_maf(maf_filename_input, OUTPUT_FOLDER, "annotate_maf_test")

with open(output_maf, "r") as f:
    generated_file_content = f.read()

assert generated_file_content == original_file_content
