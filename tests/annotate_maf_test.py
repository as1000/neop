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

output_maf = annotate_maf(maf_filename_input, OUTPUT_FOLDER, "annotate_maf_test")