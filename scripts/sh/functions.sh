# Shell helper sourced by run_job_sp_canfar_v2.0.bash.
#
# command() reads $VERBOSE, $debug_out, $pat and $STOP as free variables set
# by its caller; the two defaults below cover callers that set neither.

pat="---- "
STOP=0


## Print string, executes command, and prints return value.
function command () {
   cmd=$1
   my_dry_run=$2

   RED='\033[0;31m'
   GREEN='\033[0;32m'
   NC='\033[0m' # No Color
   # Color escape characters show up in log files
   #RED=''
   #GREEN=''
   #NC=''

   msg="running '$cmd' (dry run=$my_dry_run)"
   if [ $VERBOSE == 1 ]; then
        echo $msg
   fi
   if [ -n "$debug_out" ]; then
        echo ${pat}$msg >> $debug_out
   fi

   if [ "$my_dry_run" == "0" ]; then
        $cmd
        res=$?

        if [ -n "$debug_out" ]; then
          echo "${pat}exit code=$res" >> $debug_out
        fi

        if [ $VERBOSE == 1 ]; then
            if [ $res == 0 ]; then
              echo -e "${GREEN}success, return value = $res${NC}"
            else
              echo -e "${RED}error, return value = $res${NC}"
              if [ $STOP == 1 ]; then
                  echo "${RED}exiting  $(basename "$0")', error in command '$cmd'${NC}"
                  if [ -n "$debug_out" ]; then
                      echo "${pat}${RED}exiting  $(basename "$0")', error in command '$cmd'${NC}" >> $debug_out
                  fi
                  if [ -n "$debug_out" ]; then
                      echo "${pat}${RED}exiting  $(basename "$0")', error in command '$cmd'${NC}" >> $debug_out
                  fi
                  exit $res
              else
                  echo "${RED}continuing '$(basename "$0")', error in command '$cmd'${NC}"
              fi
            fi
        fi
   fi
}
