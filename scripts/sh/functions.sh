# Global variables
SSL=~/.ssl/cadcproxy.pem
SESSION=https://ws-uv.canfar.net/skaha/v0/session
IMAGE=images.canfar.net/unions/shapepipe
NAME=shapepipe

version="1.1"
cmd_remote="$HOME/shapepipe/scripts/sh/init_run_exclusive_canfar.sh"
pat="- "
STOP=0


# Return argument for local script to be called via curl
function set_arg() {
  my_arg="-j $job -p $psf -e $ID -N $N_SMP $arg_dry_run -d $dir -m $mh_local --debug_out $debug_out"
  echo $my_arg
}


# Add session and image IDs to log files
function update_session_logs() {
  echo $my_session >> session_IDs.txt
  echo "$my_session $ID" >> session_image_IDs.txt
}

function call_curl() {
  my_name=$1
  dry_run=$2
  debug_out=$3

  my_arg=$(set_arg)

  if [ "$dry_run" == "0" ]; then

    #my_session=`curl -E $SSL "$SESSION?$RESOURCES" -d "image=$IMAGE:$version" -d "name=${my_name}" -d "cmd=$cmd_remote" --data-urlencode "args=${my_arg[@]}" &> /dev/null`
    my_session=`curl -E $SSL "$SESSION?$RESOURCES" -d "image=$IMAGE:$version" -d "name=${my_name}" -d "cmd=$cmd_remote" --data-urlencode "args=${my_arg[@]}"`
  fi


  cmd=("curl" "-E" "$SSL" "$SESSION?$RESOURCES" "-d" "image=$IMAGE:$version" "-d" "name=${my_name}" "-d" "cmd=$cmd_remote" "--data-urlencode" "args=\"${my_arg}\"")

  if [ "$debug_out" != "-1" ]; then
    echo "${pat}call_curl $my_name $my_arg" >> $debug_out
    echo "${pat}Running ${cmd[@]} (dry_run=$dry_run)" >> $debug_out
  fi
  echo "${cmd[@]} (dry_run=$dry_run)"


  # Running $cmd does not work due to unknown problems with passing of args

  update_session_logs
}


## Print string, executes command, and prints return value.
function command () {
   cmd=$1
   dry_run=$2

   RED='\033[0;31m'
   GREEN='\033[0;32m'
   NC='\033[0m' # No Color
   # Color escape characters show up in log files
   #RED=''
   #GREEN=''
   #NC=''

   msg="running '$cmd' (dry run=$dry_run)"
   if [ $VERBOSE == 1 ]; then
        echo $msg
   fi
   if [ "$debug_out" != "-1" ]; then
        echo ${pat}$msg >> $debug_out
   fi

   if [ "$dry_run" == "0" ]; then
        $cmd
        res=$?

        if [ "$debug_out" != "-1" ]; then
          echo "${pat}exit code = $res" >> $debug_out
        fi

        if [ $VERBOSE == 1 ]; then
            if [ $res == 0 ]; then
              echo -e "${GREEN}success, return value = $res${NC}"
            else
              echo -e "${RED}error, return value = $res${NC}"
              if [ $STOP == 1 ]; then
                  echo "${RED}exiting  $(basename "$0")', error in command '$cmd'${NC}"
                  exit $res
              else
                  echo "${RED}continuing '$(basename "$0")', error in command '$cmd'${NC}"
              fi
            fi
        fi
   fi
}


function get_kind_from_job() {
    my_job=$1

    job_to_test=2
    kind="none"

    # loop over possible job numbers
    while  [ $job_to_test -le 1024 ]; do

      (( do_job = $job & $job_to_test ))
      if [[ $do_job != 0 ]]; then

        if [ $job_to_test == 32 ]; then
          if [ "$kind" == "tile" ]; then
            echo "Error: Invalid job $job. mixing tile and exp kinds"
            exit 6
          fi

          # job=32 -> set kind to exp
          kind="exp"
        elif [ $job_to_test == 2 ]; then
          if [ "$kind" == "tile" ]; then
            echo "Error: Invalid job $job. mixing tile and exp kinds"
            exit 6
          fi

          kind="exp"
        elif [ $job_to_test == 8 ]; then
          if [ "$kind" == "tile" ]; then
            echo "Error: Invalid job $job. mixing tile and exp kinds"
            exit 6
          fi

          kind="exp"
        else
          if [ "$kind" == "exp" ]; then
            echo "Error: Invalid job $job. mixing tile and exp kinds"
            exit 6
          fi

          # job != 32 -> set kind to tile
          kind="tile"
        fi

    fi

    # Multiply job number by two to get next bitwise number
    job_to_test=$((job_to_test * 2))
  done

  echo $kind
}