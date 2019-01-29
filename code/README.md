### The code of lammps I changed

##### fix_reaxc_ofek:
The code that implement the algorithem in the article

##### pair_reaxc:
The code that calculate the Boost Bond extra potential with the ReaxFF force field potential.

##### nvt_run:
The script code and input files to test the code in lammps


##### code.py:
to use it please add the line:
>> followDistFunc();

as a second line in the function end_of_step() on file fix_reaxc_ofek.cpp.
this function create text file that follows the distance between each 2 atoms.
code.py file reads the created dist file and display Graph of the distance between each 2 atoms.
to run this file:
>> python3 code.py

#### how to run it?
1. Download lammps source code from [here](https://github.com/lammps/lammps.git)
2. Copy fix_reaxc_ofek.h, fix_reaxc_ofek.cpp to lammps/src directory.
3. Copy the files pair_reaxc.h, pair_reaxc.cpp to lammps/src/USER-REAXC directory
3. Open terminal in the lammps/src directory
5. Use those commands to compile the lammps code:
  >> make yes-USER-REAXC
  
  >> make serial
  
6. Download nvt_run folder and open terminal in it
7. Run the command:
  >>path/to/lammps/src/lmp_serial < in.nvt




