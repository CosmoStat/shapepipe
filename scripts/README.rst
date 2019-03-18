Scripts directory
=================

This directory contain scripts that are used on pipeline outputs or to prepare
files. They are not run throught the pipeline framework. For more details on how
to run each them see below.

Python scripts
==============

1. `create_log_exp_headers`_

create_log_exp_headers
======================

This as to run after the module `split_exp_runner` it will create a master log
file containing all the WCS information for each CCDs of each single exposures.
To run the script :
`create_log_exp_headers path/to/split_exp_runner/output path/to/srcipt/output`
