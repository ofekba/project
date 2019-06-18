"""this script convert XYZ file to DAT file using inputs from the user:
XYZ file name and box boundaries (high and low x,y,z values) """
import sys

def conv_to_dat(text, file_name):
	text_list = text.split("\n")
	index=0
	atom_type=0
	atom_x=0
	atom_y=0
	atom_z=0
	num_atoms=0
	#get the file name from the user
	file_name=file_name[:-4]
	dat_file="#"+file_name+" molecules\n"
	file_name+=".dat"
	text_list = [x for x in text_list if x != ""]

	num_atoms=int(text_list[0])
	text_list=text_list[1:]
	dat_file+=str(num_atoms)+" atoms\n"
	dat_file+="4 atom types\n\n"

	#get the box boundaries value from the user
	box_xlo = float (input("Enter a box xlo value: ") )
	box_xhi = float (input("Enter a box xhi value: ") )
	box_ylo = float (input("Enter a box ylo value: ") )
	box_yhi = float (input("Enter a box yhi value: ") )
	box_zlo = float (input("Enter a box zlo value: ") )
	box_zhi = float (input("Enter a box zhi value: ") )
	dat_file+=str(box_xlo)+" "+str(box_xhi)+" xlo xhi\n"
	dat_file+=str(box_ylo)+" "+str(box_yhi)+" ylo yhi\n"
	dat_file+=str(box_zlo)+" "+str(box_zhi)+" zlo zhi\n\n"

	dat_file+="Masses\n\n"
	dat_file+="1 12.0\n"
	dat_file+="2 1.008\n"
	dat_file+="3 15.999\n"
	dat_file+="4 14.0\n\n"

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
	#create new DAT file
	new_file= open(file_name,"w+")
	new_file.write(dat_file)


typedict =	{
  "C": 1,
  "H": 2,
  "O": 3,
  "N": 4
}

file_name = str (input("Enter an XYZ file name: ") )
try:
    with open(file_name,"r") as fp:
        text = fp.read()
        conv_to_dat(text, file_name)
except IOError as e:
	print("ERROR-No such a file.\nUSAGE: FILENAME.xyz or FILENAME.XYZ")


