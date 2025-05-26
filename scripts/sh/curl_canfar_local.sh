# Global variables
SSL=~/.ssl/cadcproxy.pem
SESSION=https://ws-uv.canfar.net/skaha/v0/session
IMAGE=images.canfar.net/unions/shapepipe
NAME=shapepipe

# :TODO: Not working
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $HOME/shapepipe/scripts/sh/functions.sh

# Command line arguments

## Default values
job=-1
psf="psfex"
ID=-1
file_IDs=-1
N_SMP=1
fix=0
version="1.1"
cmd_remote="$HOME/shapepipe/scripts/sh/init_run_exclusive_canfar.sh"
batch=30
batch_max=200
dry_run=0
mh_local=0
sp_local=0
test_only=0
debug_out="-1"
scratch="-1"
sm=1

pat="- "

## Help string
usage="Usage: $(basename "$0") -j JOB -[e ID |-f file_IDs] [OPTIONS]
\n\nOptions:\n
   -h\tthis message\n
   -j, --job JOB\tRunning JOB, bit-coded\n
   -e, --exclusive ID
    \timage ID\n
   -f, --file_IDs path
    \tfile containing IDs\n
   -p, --psf MODEL\n
    \tPSF model, one in ['psfex'|'mccd'], default='$psf'\n
   -m, --mh_local MH\n
    \tmerged header file local (MH=0) or global (MH=1); default is $mh_local\n
   -s, --sp_local SP\n
    \tsplit local run local (SP=1) or global (SP=0); default is SP=$sp_local\n
   --sm SM\n
    \tWith (SM=1; default) or without (SM=0) spread model input\n
   -N, --N_SMP N_SMOp\n
    \tnumber of jobs (SMP mode only), default=$N_SMP\n
   -F, --fix FIX\n
    \tfix missing data (re-download tile, unzip) for FIX=1; default is $fix\n
   -V, --version\n
    \tversion of docker image, default='$version'\n
   -C, --command_remote\n
    \tremote command to run on canfar, default='$cmd_remote'\n
   -S, --scratch\n
    \tprocessing scratch directory, default is None ($scratch)\n
   -B, --batch\n
    \tbatch size = number of jobs per iteration, default=$batch\n
   -b, --batch_max\n
    \tmaximum batch size = number of jobs run simultaneously, default=$batch_max\n
   --debug_out PATH\n
    \tdebug output file PATH, default not used\n
   -n, --dry_run LEVEL\n
    \tdry run, from LEVEL=2 (no processing) to 0 (full run; default)\n
   --test\n
    \ttest mode, no processing\n
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
    -e|--exclusive)
      ID="$2"
      shift
      ;;
    -f|--file_IDs)
      file_IDs="$2"
      shift
      ;;
    -N|--N_SMP)
      N_SMP="$2"
      shift
      ;;
    -F|--fix)
      fix="$2"
      shift
      ;;
    -S|--scratch)
      scratch="$2"
      shift
      ;;
    -V|--version)
      version="$2"
      shift
      ;;
    -B|--batch)
      batch="$2"
      shift
      ;;
    -b|--batch_max)
      batch_max="$2"
      shift
      ;;
    --debug_out)
      debug_out="$2"
      shift
      ;;
    -n|--dry_run)
      dry_run="$2"
      shift
      ;;
    --test)
      test_only=1
      ;;
  esac
  shift
done


## Check options                                                                 

if [ "$test_only" == "1" ]; then
  test_arg="--test"
else
  test_arg=""
fi

if [ "$job" == "-1" ]; then                                                     
  echo "No job indicated, use option -j"                                        
  exit 2                                                                        
fi                                                                              
                                                                                
if [ "$ID" == "-1" ] && [ "$file_IDs" == "-1" ]; then                                               
  echo "No image ID(s) indicated, use option -e ID or -f file_IDs"                                   
  exit 3                                                                        
fi                                                                              

if [ "$psf" != "psfex" ] && [ "$psf" != "mccd" ]; then
  echo "PSF (option -p) needs to be 'psfex' or 'mccd'"
  exit 4
fi

if [ "$dry_run" != 0 ] && [ "$dry_run" != 1 ] && [ "$dry_run" != 2 ]; then
  echo "Invalid dry_run option, allowed are 0, 1, and 2"
  exit 5
fi

if [ "$debug_out" != "-1" ]; then
  echo "${pat}Starting $(basename "$0") $test_arg" >> $debug_out
  echo "${pat}curl ID=$ID" >> $debug_out
  echo ${pat}`date` >> $debug_out
fi

source activate shapepipe
if [ "$debug_out"  != "-1" ]; then
    echo "${pat}conda prefix = ${CONDA_PREFIX}" >> $debug_out
    echo "${pat}script version = ${script_version}" >> $debug_out
fi

# command line arguments for remote script:
# collect into string


RESOURCES="ram=4&cores=$N_SMP"
dir=`pwd`


function submit_batch() {
  path=$1

  for ID in `cat $path`; do
    IDt=`echo $ID | tr "." "-"`
    my_name="SP-${patch}-J${job}-${IDt}"
    call_curl $my_name $job $psf $ID $N_SMP $dry_run $dir $mh_local $sp_local $sm $debug_out $fix $scratch $test_arg
  done
}

batch=50
if [ "$batch" -ge "$batch_max" ]; then
  ((batch=batch_max/2))
  echo "Reducing batch size to $batch"
fi
sleep=75

((n_thresh=batch_max-batch))


if [ "$dry_run" == 2 ]; then

  # Do not call curl (dry run = 2)
  echo "Running command dry run:"

  if [ "$ID" == "-1" ]; then


    # Submit file (dry run = 2)
    for ID in `cat $file_IDs`; do
      IDt=`echo $ID | tr "." "-"`
      my_name="SP-${patch}-J${job}-${IDt}"
      call_curl $my_name $job  $psf $ID $N_SMP $dry_run $dir $mh_local $sp_local $sm $debug_out $fix $scratch $test_arg
    done

  else

    # Submit image (dry run = 2)
    IDt=`echo $ID | tr "." "-"`
    my_name="SP-${patch}-J${job}-${IDt}"
    call_curl $my_name $job  $psf $ID $N_SMP $dry_run $dir $mh_local $sp_local $sm $debug_out $fix $scratch $test_arg

  fi

else

  # Call curl
  rm -rf session_IDs.txt session_image_IDs.txt

  if [ "$ID" == "-1" ]; then

    # Submit file
    n_jobs=`cat $file_IDs | wc -l`
    if [ "$n_jobs" -gt "$batch_max" ]; then

      # Split into batches 
      prefix="${file_IDs}_split_"
      split -d -l $batch $file_IDs $prefix
      n_split=`ls -l $prefix* | wc -l`
      echo "Split '$file_IDs' into $n_split batches of size $batch"

      count=1
      n_queued=`stats_jobs_canfar.sh -w all`
      for batch in $prefix*; do
        echo "Number of queued jobs = $n_queued"
        echo "Submitting batch $batch ($count/$n_split)"
        echo -ne "\033]0;curl patch=$patch job=$job $count/$n_split\007"
        submit_batch $batch
        ((count=count+1))

        n_queued=`stats_jobs_canfar.sh -w all`

        while [ "$n_queued" -gt "$n_thresh" ]; do
          echo "Wait for #jobs = $n_queued jobs to go < $n_thresh ..."
          sleep $sleep
          n_queued=`stats_jobs_canfar.sh -w all`
        done

      done

    else

      # Submit entire file (single batch)
      echo "Submit '$file_IDs' in single batch"
      submit_batch $file_IDs

    fi

  else

    # Submit image
    IDt=`echo $ID | tr "." "-"`
    my_name="SP-${patch}-J${job}-${IDt}"
    call_curl $my_name $job  $psf $ID $N_SMP $dry_run $dir $mh_local $sp_local $sm $debug_out $fix $scratch $test_arg

  fi

fi

echo "Done $(basename "$0")" 

if [ "$debug_out" != "-1" ]; then
  echo "${pat}End $(basename "$0") $test_arg" >> $debug_out
fi
