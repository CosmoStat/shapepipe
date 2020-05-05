Scripts directory
=================

This directory contain scripts that are used on pipeline outputs or to prepare
files. They are not run throught the pipeline framework. For more details on how
to run each them see below.

Python scripts
==============

1. `create_log_exp_headers`_
2. `create_star_cat`_

create_log_exp_headers
======================

This as to run after the module `split_exp_runner` it will create a master log
file containing all the WCS information for each CCDs of each single exposures.
To run the script :
`python create_log_exp_headers.py path/to/split_exp_runner/output path/to/srcipt/output_dir`

create_star_cat
===============

This script create all the star catalogs required to run the mask module for a
computational node without internet access.
To run the script :
`python create_star_cat.py path/to/image_dir path/to/script/output_dir`
