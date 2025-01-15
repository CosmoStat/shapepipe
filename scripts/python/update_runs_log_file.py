#!/usr/bin/env python

import re
import sys
import os

# TODO: move to cs_util
def matching_subdirs(base_dir, pattern, tail=False):                                         
                                                                                 
    # Find all matching subdirectories                                           
    subdirs = []                                                                 
    for entry in os.listdir(base_dir):                                           
        full_path = os.path.join(base_dir, entry)                                
        if os.path.isdir(full_path):
            found = False

            # Look for pattern at start or end
            if pattern in entry:               

                # Get full path or last part ("tail")
                if not tail:
                    path = full_path
                else:
                    head, tail = os.path.split(full_path) 
                    path = tail

                # Remove postfix in case of multiple runs of same module
                path = re.sub("_run_.\d?", "", path)

                # Append to result
                subdirs.append(path)
                                                                                 
    # Sort according to creation date                                            
    if not tail:
        subdirs.sort(key=os.path.getctime)                                           
                                                                                 
    return subdirs


def get_module_runs(subdirs):

    all_runs = {}
    for subdir in subdirs:
        runs = matching_subdirs(subdir, "_runner", tail=True)
        if len(runs) > 0:
            all_runs[subdir] = runs

    return all_runs


def update_log_file(module_runs, log_name):

    with open(log_name, "w") as f_out: 
        for key in module_runs:
            print(key, file=f_out, end=" ")
            print(",".join(module_runs[key]), file=f_out)


def main(argv=None):                                                             
                                                                                 
    # Set default parameters                                                     
    #p_def = params_default()                                                     
                                                                                 
    # Command line options                                                       
    #options, args = parse_options(p_def)                                         
                                                                                 
    #if check_options(options) is False:                                          
        #return 1                                                                 
                                                                                 
    #param = update_param(p_def, options)                                         
                                                                                 
    base_dir = "./output"
    pattern = "run_sp_"
    log_name = f"{base_dir}/log_run_sp.txt"

    subdirs = matching_subdirs(base_dir, pattern)
    module_runs = get_module_runs(subdirs)
    #save_prev(log_name)
    update_log_file(module_runs, log_name)
                                                                                 
    return 0                                                                     
                                                                                 
                                                                                 
if __name__ == "__main__":                                                       
    sys.exit(main(sys.argv))
