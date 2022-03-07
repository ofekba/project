### The code of lammps I changed

##### fix_reaxc_checkFourset:
The code that implement the algorithem in the article

##### pair_reaxc: (or pair_reaxc_omp for multithreading)
The code that calculate the Boost Bond extra potential with the ReaxFF force field potential.

##### nvt_run:
The script code and input files to test the code in lammps


##### code.py:
The fix_reaxc_checkFourset code creates "dists.reax" text file that follows the distance between each 2 atoms.
code.py file reads this file and display Graphs of the distance between each 2 atoms as a function of time, reads the "log.lammps" output file to create temperature, pressure, total energy and potential energy Graphs as a function of time, and reads "extra_energy.txt" output file from pair_reaxc code that follows the energy that the extra potential adds to the system to create Graph of extra energy as a function of time
to run this code:
```
python3 code.py
```

##### conv_xyz_to_dat.py:
To convert XYZ file into DAT file that lammps support as an input file, run this python code.
while running the code, you will be asked to type the file name, the box borders, and the atom's mases.

#### how to run it?
1. Download lammps source code from [here](https://github.com/lammps/lammps.git)
2. Copy fix_reaxc_checkFourset.h, fix_reaxc_checkFourset.cpp to lammps/src directory.
3. Copy the files pair_reaxc.h, pair_reaxc.cpp to lammps/src/USER-REAXC directory
3. Open terminal in the lammps/src directory
5. Use those commands to compile the lammps code:
##### Serial run
```
make yes-USER-REAXC
make serial
  ```
  
6. Download nvt_run folder and open terminal in it
7. Run the command:
  >> path/to/lammps/src/lmp_serial < in.nvt \n ofek
 
 ##### multiThreading run (with OpenMP)
 first, make sure that OpenMP is a feature of your compiler.
  ```
  make yes-USER-REAXC 
  make yes-USER-OMP
  make omp
   ```
  
6. Download nvt_run folder and open terminal in it
7. Run the command:
```
OMP_NUM_THREADS=x path/to/lammps/src/lmp_omp -sf omp < in.nvt
```




