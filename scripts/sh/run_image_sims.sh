base=$HOME/v2.0/image_sims

type="grid"

tile_ID="233.293"

if [ "$type" == "grid" ];
  str_type="_grid"
else
  str_type=""

job=4091

num=1

dgs=("m" "z" "p")

mkdir -p $base/$type
cd $base/$type

for dg in dgs; do

  name="1${dg}2${dx}${str_type}_$num"

  init_run_v2.0.sh -t image_sims -s $name

  cd $name

  apptainer exec --bind /n09data,/home \
    run_job_sp_canfar_v2.0.bash -e ${tile_ID} -t image_sims -j $job

done

#--env PATH=/home/mkilbing/.local/bin:/home/mkilbing/astro/repositories/github/shapepipe/bin:$PATH --env PYTHONPATH=/home/mkilbing/astro/repositories/github/shapepipe/src:$PYTHONPATH /n17data/mkilbing/shapepipe_ngmix_v2.0.sif
