#!/usr/bin/env python3

"""Script canfar_run_analyse.py

Analyse a run (directory) on canfar using the job log files.

:Author: Martin Kilbinger

:Date: 05/2020
"""

import re, os, sys
from collections import Counter

# Results
res_ok = 0

res_unk = -1
res_wait = -2
res_subm = -3
res_abort = -4

res_noout = 1


def get_status(tile_num):

    base_name = 'log_canfar_sp_'
    #base_name = 'log_sp_tile_'

    log_name = '{}{}.log'.format(base_name, tile_num)
    out_name = '{}{}.out'.format(base_name, tile_num)
    err_name = '{}{}.err'.format(base_name, tile_num)

    status = res_unk, 'unknown status'

    if not os.path.exists(log_name):
        status = res_wait,  'waiting for submission'
    else:
        if os.path.exists(out_name):
            final_cat_found = False
            with open(out_name) as out_file:
                for line in out_file:
                    m = re.match('Upload final_cat results, (\d+) files', line)
                    if m:
                        final_cat_found = True
                        if m[1] == '0':
                            status = res_noout, '0 final_cat output files'
                        else:
                            status = res_ok, 'success, {} final_cat output files'.format(m[1])
                        break
            if final_cat_found == False:
                status = res_unk, 'Failed before final_cat'
        else:
            log_file = open(log_name)
            lines = log_file.readlines()
            log_file.close()
            for line in lines:
                if re.match('.*aborted', line):
                    status = res_abort, 'job aborted by user'
                    return status
            for line in lines:
                if re.match('.*executing', line):
                    status = res_wait, 'job running'
                    return status
            status = res_subm, 'job submitted, not started yet'

    return status


def output(status):

    print('## Issues')
    n_issue = 0
    for tile_num in sorted(status.keys()):
        if status[tile_num][0] == res_noout or status[tile_num][0] == res_noout:
            print('   ', tile_num, status[tile_num], res_wait)
            n_issue = n_issue + 1
    if n_issue == 0:
        print('   none')

    hist = Counter(status.values())
    print('## Summary')
    print('# Nb: status (code)') 
    for s in hist:
        #print('{:2d}, {}: {}'.format(int(s[0]), s[1], hist[s]))
        print('{:4d}: {} ({})'.format(hist[s], s[1], int(s[0])))
        #print(hist[s], s)


def main(argv=None):

    job_tile_path = "job_tile.sh"

    status = {}
    with open(job_tile_path) as job_file:
        for line in job_file:
            m = re.match('(\d{3}\.\d{3})', line)
            if m:
                tile_num = m[0]

                if len(argv) == 2 and argv[1] != tile_num:
                    continue
                status[tile_num] = get_status(tile_num)

    output(status)


if __name__ == "__main__":
    sys.exit(main(sys.argv))

