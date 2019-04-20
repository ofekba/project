function [x,y,z,atom_name]=xyz_from_lammps_dump(fname, natoms)
% Convert LAMMPS dump output file with xyz coordinates, but not in XYZ format (for MATERIALS STUDIO ETC)
%  into FORMAL XYZ FILE
% Example - input file is    fname='min_end';
%                      ouput file is   fname='min_end.xyz';
dir_out='';
fname_out=[dir_out,fname,'.xyz'];
% data=dlmread([fname,'.xyz'],'',7,0);
 data1=dlmread(fname,'',9,0);
 data = sortrows(data1);
x=data(:,3);y=data(:,4);z=data(:,5);
atom_num=data(:,2);
for k=1:length(atom_num)
    j=atom_num(k);
    if(j==1)
        atom_name(k)='C';
    elseif (j==2)
        atom_name(k)='H';
    elseif (j==3)
        atom_name(k)='O';
    elseif (j==4)
        atom_name(k)='N';
    end
end
disp(fname_out)
fname_out = fname+".xyz"
disp(fname_out)
fid = fopen(fname_out, 'wt');
fprintf(fid, '%6.0f \n \n',natoms );
for i=1:natoms
    fprintf(fid, '%c %8.5f %8.5f %8.5f \n',atom_name(i),x(i),y(i),z(i) );
end
fclose (fid);
return
end