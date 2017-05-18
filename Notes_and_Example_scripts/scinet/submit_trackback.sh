#!/bin/bash
#PBS -l nodes=1:ppn=8,walltime=10:00:00
#PBS -q batch
#PBS -N Ramses_all_cell_reduction
#PBS -M dwmurray@uwm.edu
#PBS -m abe

cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=1
module load gcc
source ~/yt-x86_64/bin/activate

cd
cd test-ramses/jeans_32
python ~/ramses-scripts/all_profiles_yt3.py 100 111 1 1 32 --shellsphere --backwards &
python ~/ramses-scripts/all_profiles_yt3.py 120 141 1 1 32 --shellsphere &
python ~/ramses-scripts/all_profiles_yt3.py 112 120 1 1 32 --shellsphere &
#cd ../jeans_16; python ~/ramses-scripts/all_profiles_yt3.py 100 194 1 1 16 --shellsphere --backwards &
#cd ../jeans_8/; python ~/ramses-scripts/all_profiles_yt3.py 100 266 1 1 8 --shellsphere --backwards &
#cd ../jeans_4/; python ~/ramses-scripts/all_profiles_yt3.py 100 255 1 1 4 --shellsphere --backwards &
wait

# 32 got 6 done in ~5 hours
# 16 got 18 done in ~5 hours
# 08 got 6 done in ~5 hours
# 04 got 6 done in ~5 hours
