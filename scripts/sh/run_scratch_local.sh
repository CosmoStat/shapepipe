#!/usr/bin/bash

# Command line arguments                                                         
## Default values                                                                
job=-1                                                                           
ID=-1                                                                            
N_SMP=1                                                                          
dry_run=0                                                                        
dir=`pwd`                                                                        
debug_out=-1                                                                     
scratch=/n17data/`whoami`/scratch
exec_path=$HOME/shapepipe/scripts/sh
slurm=1
                                                                                 
# mh_local is 0 (1) if merge_header_runner is run on all exposures,              
# which is standard so far (run on exposures of given tile only; new)            
mh_local=0                                                                       
                                                                                 
# sp_local is 0 (1) is split_headers_runner and mask_runner is run               
# on all exposures (locally). Not 100% automatic yet.                            
sp_local=1                                                                       
VERBOSE=1                                                                        
                                                                                 
pat="-- "

# Help string                                                                   
usage="Usage: $(basename "$0") -j JOB -e ID -k KIND [OPTIONS]                    
\n\nOptions:\n                                                                   
   -h\tthis message\n                                                            
   -j, --job JOB\tRUnning JOB, bit-coded\n                                       
   -e, --exclusive ID                                                            
    \timage ID\n                                                                 
   -p, --psf MODEL\n                                                             
    \tPSF model, one in ['psfex'|'mccd'], default='$psf'\n                       
   -m, --mh_local MH\n                                                           
   \tmerged header file local (MH=0) or global (MH=1); default is $mh_local\n    
   -N, --N_SMP N_SMOp\n                                                          
    \tnumber of jobs (SMP mode only), default from original config files\n       
   -d, --directory\n                                                             
    \trun directory, default is pwd ($dir)\n                                     
   -n, --dry_run\n                                                               
    \tdry run, no actuall processing\n                                           
   --debug_out PATH\n                                                            
   \tdebug output file PATH, default not used\n                                  
"                                                                                
                                                                                 
## Help if no arguments                                                          
if [ -z $1 ]; then                                                               
  echo -ne $usage                                                                
  exit 1                                                                         
fi                                                                               
                                                                                 
## Parse command line                                                            
while [ $# -gt 0 ]; do                                                           
  case "$1" in                                                                   
    -h)                                                                          
      echo -ne $usage                                                            
      exit 0                                                                     
      ;;                                                                         
    -j|--job)                                                                    
      job="$2"                                                                   
      shift                                                                      
      ;;                                                                         
    -e|--exclusive)                                                              
      ID="$2"                                                                    
      shift                                                                      
      ;;                                                                         
    -p|--psf)                                                                    
      psf="$2"                                                                   
      shift                                                                      
      ;;                                                                         
    -m|--mh_local)                                                               
      mh_local="$2"                                                              
      shift                                                                      
      ;;                                                                         
    -N|--N_SMP)                                                                  
      N_SMP="$2"
      shift                                                                      
      ;;                                                                         
    -d|--directory)                                                              
      dir="$2"                                                                   
      shift                                                                      
      ;;                                                                         
    -n|--dry_run)                                                                
      dry_run=1                                                                  
      ;;                                                                         
    --debug_out)                                                                 
      debug_out="$2"                                                             
      shift                                                                      
      ;;                                                                         
  esac                                                                           
  shift                                                                          
done                                                                             
                                                                                 
## Check options                                                                 
if [ "$job" == "-1" ]; then                                                      
  echo "No job indicated, use option -j"                                         
  exit 2                                                                         
fi                                                                               
                                                                                 
if [ "$exclusive" == "-1" ]; then                                                
  echo "No image ID indicated, use option -e"                                    
  exit 3                                                                         
fi                                                                               
                                                                                 
if [ "$psf" != "psfex" ] && [ "$psf" != "mccd" ]; then                           
  echo "PSF (option -p) needs to be 'psfex' or 'mccd'"                           
  exit 4                                                                         
fi


source $HOME/shapepipe/scripts/sh/functions.sh                                   


kind=$(get_kind_from_job $job)


# Load common functions
source $HOME/shapepipe/scripts/sh/functions.sh


# Start script

if [ "$scratch" != "-1" ]; then

  command "mkdir -p $scratch/${kind}_runs" $dry_run
  command "cp -R ${kind}_runs/$ID $scratch/${kind}_runs" $dry_run                     
  command "cd $scratch" $dry_run

  if [ "$slurm" == "0" ]; then
    command "init_run_exclusive_canfar.sh -j $job -p $psf -m $mh_local -N $N_SMP -e $ID" $dry_run
  else
    STATUS=$(sbatch --output=./sbatch-$ID.out --partition=comp --job-name="j${job}_${ID}" --ntasks-per-node=$N_SMP --time=32:00:00 --mem=64G $exec_path/init_run_exclusive_canfar.sh -j $job -p $psf -m $mh_local -N $N_SMP -e $ID)

    JOB_ID=$(echo $STATUS | cut -d ' ' -f 4)
    echo "JOB_ID=$JOB_ID"

    # Wait for the job to finish
    while true; do
      STATUS=$(squeue -j "$JOB_ID" -h -o "%T")
      if [[ -z "$STATUS" ]]; then
        echo "job $JOB_ID no longer in the queue"
        break
      fi

      echo "Waiting for job $JOB_ID in state '$STATUS' to complete..."
      sleep 10
    done

    echo "Job $JOB_ID has completed. Proceeding with the script..."
  fi

  if [ "$job" == "32" ]; then
    command "mv ${kind}_runs/$ID/output/run_sp_exp_SxSe* $dir/${kind}_runs/$ID/output" $dry_run
  elif [ "$job" == "64" ]; then
    command "mv ${kind}_runs/$ID/output/run_sp_tile_PsViSm** $dir/${kind}_runs/$ID/output" $dry_run
  fi

  command "rm -rf ${kind}_runs/$ID" $dry_run
  command "cd $dir/${kind}_runs/$ID" $dry_run
  # Gave Input/Output python error
  #command "update_runs_log_file.py" $dry_run
  command "cd $dir" $dry_run

fi

exit 0
