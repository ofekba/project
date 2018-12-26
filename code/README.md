### The code of lammps I changed

##### fix_reaxc_ofek:
The code that implement the algorithem in the article

##### nvt_run:
The script code and input files to test the code in lammps

#### how to run it?
1. download lammps source code from [here](https://github.com/lammps/lammps.git)
2. copy fix_reaxc_ofek.h, fix_reaxc_ofek.cpp to lammps/src directory.
3. open terminal in the lammps/src directory
5. use those commands:
  >> make yes-USER-REAXC
  
  >> make serial //compile the lammps code for serial use.
6. download nvt_run folder and open terminal in it
7. run the command:
  >>path/to/lammps/src/lmp_serial < in.nvt




