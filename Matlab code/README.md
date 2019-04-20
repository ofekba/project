## Visualize Lammps Dump files using Matlab
To visulize dump files using Matlab's "molviewer" command you need to convert
dump files into xyz files with this attached code. then, using "molviewer" command.
#### How to do that?
open folder with this code and dump file
run command:
>> xyz_from_lammps_dump(DUMP_FILE_NAME, NUM_OF_ATOMS)

then run:
>> molviewer(DUMP_FILE_NAME.xyz)
