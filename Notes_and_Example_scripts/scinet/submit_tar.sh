#!/bin/bash
#PBS -l nodes=1:ppn=8,walltime=03:00:00
#PBS -q batch
#PBS -N hd_tarball
#PBS -M dwmurray@uwm.edu
#PBS -m abe

cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=1
#module load gcc
#source ~/yt-x86_64/bin/activate

cd
cd scratch/test-ramses/hd_32jeans_1024_v2

for f in output_0006*; do tar -zcvf tar_$f.tar.gz $f; echo "Processing $f file.."; done > ~/scratch/test-ramses/hd_32jeans_1024_v2/tar_outputs0006.txt 2>&1 &

for f in output_0007*; do tar -zcvf tar_$f.tar.gz $f; echo "Processing $f file.."; done > ~/scratch/test-ramses/hd_32jeans_1024_v2/tar_outputs0007.txt 2>&1 &
for f in output_0008*; do tar -zcvf tar_$f.tar.gz $f; echo "Processing $f file.."; done > ~/scratch/test-ramses/hd_32jeans_1024_v2/tar_outputs0008.txt 2>&1 &

#for f in output_0009*; do tar -zcvf tar_$f.tar.gz $f; echo "Processing $f file.."; done > ~/scratch/test-ramses/hd_8jeans_1024_v1/tar_outputs0009.txt 2>&1 &
#for f in output_0005*; do tar -zcvf tar_$f.tar.gz $f; echo "Processing $f file.."; done > ~/scratch/test-ramses/mhd_32jeans_1024/tar_outputs0005.txt 2>&1 &

wait


# for loop example
#for f in output_0015*; do tar -zcvf tar_$f.tar.gz $f; echo "Processing $f file.."; done > ~/scratch/test-ramses/mhd_32jeans_1024/tar_outputs0015.txt 2>&1 &
# Individual file example
#tar -zcvf tar_output_00079.tar.gz output_00079 &

# Tar multiple files into one tarball.
#tar -zcvf tar_BB_hdf5_plt_cnt_01.tar.gz BB_hdf5_plt_cnt_01*\
#tar -zcvf tar_BB_hdf5_plt_cnt_014.tar.gz BB_hdf5_plt_cnt_014* &