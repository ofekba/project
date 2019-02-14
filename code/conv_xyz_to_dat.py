import sys

typedict =	{
  "C": 1,
  "H": 2,
  "O": 3,
  "N": 4
}

fp = open("epon2_detda_a_14_94.xyz","r") 
text=fp.read() 

text_list = text.split("\n")
index=0
atom_type=0
atom_x=0
atom_y=0
atom_z=0

num_atoms=0

dat_file="#epon2 and detda molecules\n"
text_list = [x for x in text_list if x != ""]



num_atoms=int(text_list[0])
text_list=text_list[1:]
dat_file+=str(num_atoms)+" atoms\n"
dat_file+="4 atom types\n\n"

box_x = float (input("Enter a box x value: ") )
box_y = float (input("Enter a box y value: ") )
box_z = float (input("Enter a box z value: ") )
dat_file+="0 "+str(box_x)+" xlo xhi\n"
dat_file+="0 "+str(box_y)+" ylo yhi\n"
dat_file+="0 "+str(box_z)+" zlo zhi\n\n"

dat_file+="Masses\n\n"
C_mass = float (input("Enter a C mass: ") )
H_mass = float (input("Enter a H mass: ") )
O_mass = float (input("Enter a O mass: ") )
N_mass = float (input("Enter a N mass: ") )
dat_file+="1 "+str(C_mass)+"\n"
dat_file+="2 "+str(H_mass)+"\n"
dat_file+="3 "+str(O_mass)+"\n"
dat_file+="4 "+str(N_mass)+"\n\n"


dat_file+="Atoms\n\n"
atom_tag=1
for line in text_list:
	ln=line.split(" ")
	ln = [x for x in ln if x != ""]
	atom_type=typedict[ln[0]]
	atom_x=float(ln[1])
	atom_y=float(ln[2])
	atom_z=float(ln[3])
	dat_file+=str(atom_tag)+" "+str(atom_type)+" 0 "+str(atom_x)+" "+str(atom_y)+" "+str(atom_z)
	atom_tag+=1
	if atom_tag<=num_atoms:
		dat_file+="\n"

new_file= open("epon2_detda_a_14_94.dat","w+")
new_file.write(dat_file)




