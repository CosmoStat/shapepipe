#!/bin/bash

# init_run_exclusive_canfar.sh

# Command line arguments
## Default values
job=-1
ID=-1
N_SMP=1
dry_run=0
dir=`pwd`
debug_out=-1
scratch=-1
fix=0
test_only=0
sm=1

# mh_local is 0 (1) if merge_header_runner is run on all exposures,
# which is standard so far (run on exposures of given tile only; new)
mh_local=0

# sp_local is 0 (1) is split_headers_runner and mask_runner is run
# on all exposures (locally). Not 100% automatic yet. 
sp_local=0
VERBOSE=1

pat="-- "


## Help string
usage="Usage: $(basename "$0") -j JOB -e ID -k KIND [OPTIONS]
\n\nOptions:\n
   -h\tthis message\n
   -j, --job JOB\tRUnning JOB, bit-coded\n
   -e, --exclusive ID
    \timage ID\n
   -p, --psf MODEL\n
    \tPSF model, one in ['psfex'|'mccd'], default='$psf'\n
   -m, --mh_local MH\n
   \tmerge header file local (MH=1) or global (MH=0); default is $mh_local\n
   -s, --sp_local SP\n
   \tsplit local run local (SP=1) or global (SP=r0wwdefault is $sp_local\n
   --sm SM\n
   \tWith (SM=1; default) or without (SM=0) spread model input\n
   -N, --N_SMP N_SMOp\n
    \tnumber of jobs (SMP mode only), default from original config files\n
   -d, --directory\n
    \trun directory, default is pwd ($dir)\n
   -S, --scratch\n
    \tprocessing scratch directory, default is None ($scratch)\n
   -F, --fix FIX\n
    \tfix missing data (re-download tile, unzip) for FIX=1; default is $FIX\n
   -n, --dry_run LEVEL\n
    \tdry run (LEVEL=1), no actual processing; default is $dry_run\n
   --debug_out PATH\n
   \tdebug output file PATH, default not used\n
   --test\n
   \ttest mode, no processind\n
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
    -s|--sp_local)
      sp_local="$2"
      shift
      ;;
    --sm)
      sm="$2"
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
    -S|--scratch)
      scratch="$2"
      shift
      ;;
    -n|--dry_run)
      dry_run="$2"
      shift
      ;;
    -F|--fix)
      fix="$2"
      shift
      ;;
    --debug_out)
      debug_out="$2"
      shift
      ;;
    --test)
      test_only=1
      ;;
  esac                                                                          
  shift                                                                         
done


# functions
function message() {
  msg=$1
  my_debug_out=$2
  my_exit=$3

  echo $msg
  if [ "$my_debug_out" != "-1" ]; then
    echo ${pat}$msg >> $my_debug_out
  fi

  if [ "$my_exit" != "-1" ]; then
    echo "${pat}exiting with code $my_exit" >> $my_debug_out
    exit $my_exit
  fi
}


# Init message
message "test=$test_only" $debug_out -1
if [ "$test_only" == "1" ]; then
  msg="init_run_exclusive.py script test mode, exiting."
  ex=0
else
  msg="init_run_exclusive.py script processing mode, starting."
  ex=-1
fi
message "$msg" $debug_out $ex


## Check options
message "checking options" $debug_out -1

if [ "$job" == "-1" ]; then
  message "No job indicated, use option -j" $debug_out 2
fi

if [ "$exclusive" == "-1" ]; then
  message "No image ID indicated, use option -e" $debug_out 3
fi

if [ "$psf" != "psfex" ] && [ "$psf" != "mccd" ]; then
  message "PSF (option -p) needs to be 'psfex' or 'mccd'" $debug_out 4
fi

if [ "$mh_local" != "0" ] && [ "$mh_local" != "1" ]; then
  message "mh_local (option -m) needs to be 0 or 1" $debug_out 5
fi

if [ "$sp_local" != "0" ] && [ "$sp_local" != "1" ]; then
  message "sp_local (option -m) needs to be 0 or 1" $debug_out 6
fi

if [ "$dry_run" != "0" ] && [ "$dry_run" != 1 ]; then
  message "dry_run must be 0 or 1, not $dry_run" $debug_out 8
fi


# Start script

source $HOME/shapepipe/scripts/sh/functions.sh

msg="Starting $(basename "$0")"
message "$msg" $debug_out -1
message "`date`" $debug_out -1
message "ID=$ID" $debug_out -1

# Set kind
kind=$(get_kind_from_job $job)

if [ "$kind" == "none" ]; then
  message "Error: invalid job $job" $debug_out 5
fi

message "kind=$kind" $debug_out -1


if [ "$dry_run" == "1" ]; then
  message "running in dry run mode" $debug_out -1
else
  message "not running in dry run mode" $debug_out -1
fi

CONDA_PREFIX=$HOME/.conda/envs/shapepipe
PATH=$PATH:$CONDA_PREFIX/bin

cd $dir

if [ ! -d ${kind}_runs ]; then
  command "mkdir ${kind}_runs" $dry_run
fi

if [ "$fix" == "1" ]; then
  message "Fixing missing data" $debug_out -1

  message "Download tile images" $debug_out -1
  cd data_tiles
  file_name=CFIS.$ID.r.fits
  if [ ! -e $file_name ]; then
    message "Downloading $file_name" $debug_out -1
    vcp vos:cfis/tiles_DR5/$file_name .
  else
    message "File $file_name exists, skipping" $debug_out -1
  fi
  file_name=CFIS.$ID.r.weight.fits.fz
  if [ ! -e $file_name ]; then
    message "Downloading $file_name" $debug_out -1
    vcp vos:cfis/tiles_DR5/$file_name .
  else
    message "File $file_name exists, skipping" $debug_out -1
  fi
  cd ..

  message "link to tiles..." $debug_out -1
  cd output/run_sp_GitFeGie_202*/get_images_runner_run_1/output
  IDt=`echo $ID | tr "." "-"`
  path_tiles=/arc/home/kilbinger/cosmostat/v2/pre_v2/$psf/$patch/data_tiles
  link_name=CFIS_image-$IDt.fits
  if [ ! -e $link_name ]; then
    message "Creating link $link_name" $debug_out -1
    # force to overwrite broken link
    ln -sf $path_tiles/CFIS.$ID.r.fits $link_name
  else
    message "Link $link_name exists, skipping" $debug_out -1
  fi

  link_name=CFIS_weight-$IDt.fitsfz
  if [ ! -e $link_name ]; then
    message "Creating link $link_name" $debug_out -1
    ln -sf $path_tiles/CFIS.$ID.r.weight.fits.fz $link_name
  else
    message "Link $link_name exists, skipping" $debug_out -1
  fi

  command "cd ../../../.." $dry_run

  message "Unzip weight ($dry_run)" $debug_out -1
  command "cd tile_runs/$ID" $dry_run
  export SP_RUN=`pwd`
  command "shapepipe_run -c cfis/config_tile_Uz.ini -e $ID" $dry_run

  cd $dir
else
  message "Not fixing missing data" $debug_out -1
fi


cd ${kind}_runs

if [ ! -d "$ID" ]; then
  command "mkdir $ID" $dry_run
fi

cd $ID
pwd

# Point cfis to local link, to be independent of platform
ln -sf ~/shapepipe/example/cfis


if [ ! -d "output" ]; then
  command "mkdir output" $dry_run
fi

cd output


# Update links to global run directories (GiFeGie)
# New: 27/11/2024: Remove link to Uz, conflict with fix
for my_dir in $dir/output/run_sp_[G]*; do
  command "ln -sf $my_dir" $dry_run
done

# The following could be done also with fix=1
if [ "$fix" == "0" ] && [ "$kind" == "tile" ]; then
  if  [ ! -e "run_sp_Uz*" ]; then

    message "previous uncompress run not found..." $debug_out -1

    # Remove broken links
    for my_dir in run_sp_Uz*; do
      if  [ -L $my_dir ]; then
        command "rm $my_dir" $dry_run
        message "remove broken link $my_dir" $debug_out -1
      fi
    done

    # Create working link
    for my_dir in $dir/output/run_sp_Uz*; do
      command "ln -sf $my_dir" $dry_run                                             
      message "create valid link $my_dir" $debug_out -1
    done
    cd ..
    command "update_runs_log_file.py" $dry_run
    cd output
  fi
fi


# Combined flags

## Tiles
command "ln -sf $dir/output/run_sp_Ma_tile" $dry_run

if [ "$sp_local" == "0" ]; then
  # Exposures
  command "ln -sf $dir/output/run_sp_Ma_exp" $dry_run
#else
  #command "rm -f $dir/output/run_sp_Ma_exp" $dry_run
fi

# Check for existing exp SpMh dir
dir_SpMh="run_sp_exp_SpMh"
if [[ -d "$dir/output/$dir_SpMh" ]]; then
  # create link
  command "ln -sf $dir/output/$dir_SpMh" $dry_run
else
  message "No global SpMh directory found" $debug_out -1
  if [ -L "$dir_SpMh" ] && [ ! -e "$dir_SpMh" ]; then
      message "Remove broken link $dir_SpMh" $debug_out -1
      rm $dir_SpMh
  fi
fi

(( do_job = $job & 2 ))
if [ $do_job != 0 ] && [ "$sp_local" == "1" ]; then
#if [ ! -e "run_exp_Sp_shdu" ] && [ "$sp_local" == "1" ]; then
  # run local Sp if not done already; works only with mh_local=1; this step needs to be done
  # before following mh_local=1 steps 
  message "run local sp" $debug_out -1
  #command "rm -rf run_sp_GitFeGie*/get_images_runner_run_2" $dry_run
  command "rm -rf run_sp_Gie*" $dry_run
  command "rm -rf run_sp_exp_Sp*" $dry_run

  # Create new get_image dir
  new_dir="run_sp_Gie/get_images_runner/output"
  command "rm -rf $new_dir" $dry_run
  command "mkdir -p $new_dir" $dry_run

  # Link to image, weight, and flag file of current exposure

  # Remove HDU extension
  command "cd $new_dir" $dry_run
  exp_ID=$(echo "$ID" | sed 's/-[0-9]\{1,2\}//')
  for file in $dir/output//run_sp_GitFeGie_20*/get_images_runner_run_2/output/*${exp_ID}.* ; do
    echo $file
    command "ln -s $file" $dry_run
  done
  command "cd ../../.." $dry_run

  # Run Sp
  command "cd .." $dry_run
  if [ "$scratch" != "-1" ]; then
    command "cd ../.." $dry_run
    command "mkdir -p $scratch/exp_runs" $dry_run
    command "cp -R exp_runs/$ID $scratch/exp_runs" $dry_run
    command "cd $scratch/exp_runs/$ID" $dry_run
  fi
  command "update_runs_log_file.py" $dry_run
  export SP_RUN=`pwd`
  command "shapepipe_run -c cfis/config_exp_Sp.ini -e $exp_ID" $dry_run

  # Only keep CCD of this ID 
  command "mkdir -p output/run_sp_exp_Sp_shdu/split_exp_runner/output" $dry_run
  command "mv output/run_sp_exp_Sp/split_exp_runner/output/*$ID.* output/run_sp_exp_Sp_shdu/split_exp_runner/output" $dry_run
  command "mv output/run_sp_exp_Sp/split_exp_runner/output/headers* output/run_sp_exp_Sp_shdu/split_exp_runner/output" $dry_run
  command "rm -rf output/run_sp_exp_Sp" $dry_run
  if [ "$scratch" != "-1" ]; then
    command "mv output/run_sp_exp_Sp_shdu $dir/exp_runs/$ID/output" $dry_run
    command "cd .." $dry_run
    command "rm -rf $ID" $dry_run
    command "cd $dir" $dry_run
  fi
  command "update_runs_log_file.py" $dry_run
  cd output

  if [ "$job" == "2" ]; then
    exit 0
  fi

fi

if [ "$kind" == "tile" ] && [ "$sp_local" == "1" ]; then
  cd ../../..
  command "link_to_exp_for_tile.py -t $ID -i tile_runs -I exp_runs -s $sp_local" $dry_run
  cd tile_runs/$ID
  command "combine_runs.bash -p psfex -c shdu" $dry_run
  cd output
fi


if [ "$mh_local" == "0" ]; then
  if [ ! -f log_exp_headers.sqlite ]; then
    # Global Mh and file does not exist -> symlink to
    # gllobal mh file
    message "creating global link to exp headers (mh_local=0)" $debug_out -1
    command "ln -s $dir/output/log_exp_headers.sqlite" $dry_run
  else
    message "global link to exp headers (mh_local=0) exists" $debug_out -1
  fi
else
  # Local Mh
  message "not creating global link to exp headers (mh_local=1)" $debug_out -1
  if [ "$ID" == "-1" ]; then
    message "ID needs to be given (option -e) for mh_local" $debug_out 6
  fi

  # Check and remove symbolic (global) mh file link
  if [ -L log_exp_headers.sqlite ]; then
    # Local Mh and symlink -> remove previous link to
    # (potentially incomplete) global mh file
    message "Removing previous mh sym link" $debug_out -1
    command "rm log_exp_headers.sqlite" $dry_run
  else
    message "no mh link found" $debug_out -1
  fi

  # Check size of existing header file
  if [ -e log_exp_headers.sqlite ]; then
    size=$(stat -c %s log_exp_headers.sqlite)
    if (( size > 15000 )); then
      message "Found valid local mh file, continuing" $debug_out -1
    else
      message "Existing local mh file looks invalid, deleting" $debug_out -1
      rm -f log_exp_headers.sqlite
    fi
  fi

  if [ ! -e log_exp_headers.sqlite ]; then
    message "Creating local mh file" $debug_out -1

    cd ..
    command "update_runs_log_file.py" $dry_run
    export SP_RUN=`pwd`
    command "shapepipe_run -c cfis/config_exp_Mh.ini" $dry_run
    cd output
  else
    message "Found local mh file, continuing" $debug_out -1
  fi
fi


(( do_job = $job & 8 ))
if [ $do_job != 0 ] && [ "$sp_local" == "1" ]; then
  # Remove previous local Ma runs
  message "cdsclient = $CDSCLIENT" $debug_out -1
  command "rm -rf run_sp_exp_Ma*" $dry_run
fi

(( do_job = $job & 16 ))
if [[ $do_job != 0 ]]; then
  # Remove previous Sx runs
  command "rm -rf run_sp_tile_Sx_*" $dry_run
fi

# Update links to exposure run directories, which were created in job 32
(( do_job = $job & 64 ))
if [[ $do_job != 0 ]]; then
  if [ "$kind" == "tile" ]; then
    cd ../../..
    command "link_to_exp_for_tile.py -t $ID -i tile_runs -I exp_runs -s $sp_local" $dry_run
    cd ${kind}_runs/$ID/output

    # Remove duplicate job-16 runs (tile detection)
    # New (P8) commented
    #n_16=`ls -rt1d run_sp_tile_Sx_* | wc -l`
    #if [ "$n_16" != "1" ]; then
      #n_remove="$(($n_16-1))"
      #echo "removing $n_remove duplicate old job-16 runs"
      #command "rm -rf `ls -rt1d run_sp_tile_Sx_* | head -$n_remove`" $dry_run
    #fi

    # Remove previous runs of this job
    rm -rf run_sp_tile_PsViSmVi*
  fi
fi

(( do_job = $job & 128 ))
if [[ $do_job != 0 ]]; then

    echo

    cat_ngmix="run_sp_tile_ngmix_Ng1u/ngmix_runner/output/ngmix-*.fits"
    dir_ngmix_prev="run_sp_tile_ngmix_Ng1u_prev/ngmix_runner/output"
    cat_ngmix_prev="$dir_ngmix_prev/ngmix-*.fits"

    # Check whether ngmix output exists
    if [ -e $cat_ngmix ]; then
        message "ngmix output catalogue exists" $debug_out -1

        # Check whether previous ngmix directory and output cat exist
        exists="1"
        if [ ! -d $dir_ngmix_prev ]; then
            exists="0"
        elif [ ! -e "$cat_ngmix_prev" ]; then
            exists="1"
        fi
        if [ "$exists" == "0" ]; then
            message "Moving to previous batch-save dir (does not exist yet)" $debug_out -1
            command "mkdir -p $dir_ngmix_prev" $dry_run
            command "mv $cat_ngmix $dir_ngmix_prev" $dry_run
        else
            # Compare file sizes
            size_cat_ngmix=$(stat -c%s $cat_ngmix)
            size_cat_ngmix_prev=$(stat -c%s $cat_ngmix_prev)
            if [ "$size_cat_ngmix" -gt "$size_cat_ngmix_prev" ]; then
                message "Moving to batch-save dir, overwriting smaller batch-save cat" $debug_out -1
                command "mv $cat_ngmix $dir_ngmix_prev" $dry_run
            else
                message "Previous batch-save dir not smaller, removing ngmix output" $debug_out -1
                command "rm $cat_ngmix" $dry_run
            fi
        fi
    else
        # Whether or not previous ngmix exists, job_sp_canfar will handle it
        message "No ngmix output exists, continuing..." $debug_out -1
    fi

    echo
fi

(( do_job = $job & 256 ))
if [[ $do_job != 0 ]]; then

  # Remove previous runs of this job
  rm -rf run_sp_Ms_20??_*

fi

(( do_job = $job & 512 ))
if [[ $do_job != 0 ]]; then

  # Remove previous runs of this job
  rm -rf run_sp_Mc_20??_*

fi


cd ..

# Update log file
command update_runs_log_file.py $dry_run

echo -n "pwd: "
pwd

echo -n "environment: "
echo $CONDA_PREFIX

# To avoid (new?) qt error with setools (-j 32)
export DISPLAY=:1.0


if [ "$scratch" != "-1" ]; then
  # Copy inputs to scratch
  command "mkdir -p $scratch/${kind}_runs" $dry_run
  cd ../..
  command "pwd" $dry_run
  command "cp -R ${kind}_runs/$ID $scratch/${kind}_runs" $dry_run
  command "cd $scratch/${kind}_runs/$ID" $dry_run
fi

command "job_sp_canfar.bash -p psfex -j $job -e $ID --n_smp $N_SMP --nsh_jobs $N_SMP --debug_out $debug_out --sm $sm " $dry_run

if [ "$scratch" != "-1" ]; then
  cd ../..
  if [ "$job" == "16" ]; then                                                   
    command "mv ${kind}_runs/$ID/output/run_sp_Sx_* $dir/${kind}_runs/$ID/output" $dry_run
  elif [ "$job" == "32" ]; then                                                   
    command "mv ${kind}_runs/$ID/output/run_sp_exp_SxSe* $dir/${kind}_runs/$ID/output" $dry_run
  elif [ "$job" == "64" ]; then                                                 
    command "mv ${kind}_runs/$ID/output/run_sp_tile_PsViSm** $dir/${kind}_runs/$ID/output" $dry_run
  elif [ "$job" == "128" ]; then
    command "mv ${kind}_runs/$ID/output/run_sp_tile_ngmix* $dir/${kind}_runs/$ID/output" $dry_run
  else
    echo "scratch mode with job=$job not implemented yet, exiting..."
    exit 8
  fi 

  command "rm -rf ${kind}_runs/$ID" $dry_run
  command "cd $dir/${kind}_runs/$ID" $dry_run
  command "update_runs_log_file.py" $dry_run
  command "cd $dir" $dry_run
fi

cd $dir

msg="End $(basename "$0")"
echo $msg
if [ "$debug_out" != "-1" ]; then
  echo $pat$msg >> $debug_out-
fi
