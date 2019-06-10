#!script that run ".in" files in LAMMPS one after one
#!enter run directory name
cd nvt_BB_real_1
echo "\n\nRUN with 1 THREAD\n\n"
OMP_NUM_THREADS=1 /home/ofekba/Desktop/original_lammps/lammps/src/lmp_omp -sf omp < in.nvt
cd -
cd nvt_BB_real_2
echo "\n\nRUN with 2 THREAD\n\n"
OMP_NUM_THREADS=2 /home/ofekba/Desktop/original_lammps/lammps/src/lmp_omp -sf omp < in.nvt
cd -
cd nvt_BB_real_4
echo "\n\nRUN with 4 THREAD\n\n"
OMP_NUM_THREADS=4 /home/ofekba/Desktop/original_lammps/lammps/src/lmp_omp -sf omp < in.nvt
cd -
cd nvt_BB_real_8
echo "\n\nRUN with 8 THREAD\n\n"
OMP_NUM_THREADS=8 /home/ofekba/Desktop/original_lammps/lammps/src/lmp_omp -sf omp < in.nvt
cd -
cd nvt_BB_real_12
echo "\n\nRUN with 12 THREAD\n\n"
OMP_NUM_THREADS=12 /home/ofekba/Desktop/original_lammps/lammps/src/lmp_omp -sf omp < in.nvt
cd -
cd nvt_BB_real_serial
echo "\n\nRUN with serial\n\n"
/home/ofekba/Desktop/original_lammps/lammps/src/lmp_serial  < in.nvt
cd -
echo "\n\n\n*~*~*~ finito kilililiiiii ~*~*~* \n\n\n"
