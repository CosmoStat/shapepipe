# Global variables                                                              
SSL=~/.ssl/cadcproxy.pem                                                        
SESSION=https://ws-uv.canfar.net/skaha/v0/session                               
IMAGE=images.canfar.net/unions/shapepipe                                        
NAME=shapepipe

version="1.1"
cmd_remote="$HOME/shapepipe/scripts/sh/init_run_exclusive_canfar.sh"
pat="- "


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
                                                                                
  #if [ "$debug_out" != "-1" ]; then                                             
    #echo "${pat}call_curl $my_name $my_arg" >> $debug_out                       
    #echo "${pat}Running ${cmd[@]} (dry_run=$dry_run)" >> $debug_out             
  #fi                                                                            
  echo "${cmd[@]} (dry_run=$dry_run)"                                           
                                                                                
                                                                                
  # Running $cmd does not work due to unknown problems with passing of args     
                                                                                
  update_session_logs                                                           
}

