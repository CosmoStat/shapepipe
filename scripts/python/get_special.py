#!/usr/bin/env python

import re

path = "summary/special_job_128_ngmix_runner_6.txt"
path = "summary/special_job_16_sextractor_runner_7.txt"

with open(path) as f:
    lines = f.readlines()

for line in lines:
    # Regular expression to match the pattern
    match = re.search(r'process-(\d+)-(\d+)\.log', line)

    if match:
        result = f"{match.group(1)}.{match.group(2)}"
        print(result)  # Output: 439.250
    else:
        print("Pattern not found.")

