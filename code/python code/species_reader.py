"""this script reads the species.out LAMMPS output file and creates a
graph of num of species for each species as a function of time """
import sys
import math
import numpy as np
from matplotlib import pyplot as plt


"""create graph using mathplot.
	input- axis x, axis y, label for each axis and title for the graph"""
def make_graph(axis_x, axis_y, xlabel="x", ylabel="y", title="Graph"):
	plt.plot(axis_x,axis_y)
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	x = species_dict.get(ylabel)
	if x==0 or x==1: plt.gca().invert_yaxis()
	plt.show()
	

"""remove tabs from strings """	
def remove(string): 
    return string.replace("\t", "")  


"""count the number of atom in fragment"""
def count_shows(freg, f):
	toNum=0
	while(f<len(freg)):
		if 48<=ord(freg[f])<=57:
			toNum=toNum*10+int(freg[f])
			f+=1
		else:
			if toNum==0: toNum=1
			return f,toNum
	if toNum==0: toNum=1
	return f,toNum


"""main func"""
def species_reader():
	i=0
	num_species=4
	time_Step=[]
	No_Moles=[]
	No_Specs=[]
	trucker=[[],[],[],[]]
	fp = open("species.out","r") 
	text=fp.read() 
	text_list = text.split("\n")
	species_at_end=[]
	species_at_start=[]
	while i<len(text_list)-1:
		ln = text_list[i].split(" ")
		ln=[x for x in ln if x]
		species_list=ln[4].split("\t")
		species_list=[x for x in species_list if x]
		i+=1
		ln = text_list[i].split(" ")
		ln=[remove(x) for x in ln if x]
		time_Step.append(int(ln[0]))
		No_Moles.append(int(ln[1]))
		No_Specs.append(int(ln[2]))
		for j in range(len(species_list)):
			x = species_dict.get(species_list[j])
			if x is None:
				species_dict.update( {species_list[j] : num_species} )
				trucker.append([])
				trucker[num_species]=[0 for ns in range(len(time_Step)-1)]
				trucker[num_species].append(int(ln[3+j]))
				num_species+=1
			else:
				x=int(x)
				trucker[x].append(int(ln[3+j]))
		for j in range(len(trucker)):
			if len(trucker[j])<len(time_Step): trucker[j].append(0)
			else:
				if i>=len(text_list)-2: species_at_end.append(j)
				if i<=2: species_at_start.append(j)
		i+=1

	#create list of all species at the end and start of the run without duplicates
	result_list=list(set(species_at_end)|set(species_at_start))
	spec_list=[list(species_dict.keys())[i] for i in result_list]
	
	#count the num of each atom in each species
	c_num=[]
	h_num=[]
	o_num=[]
	n_num=[]
	
	for freg in spec_list:
		f=0
		freg_c=0
		freg_h=0
		freg_o=0
		freg_n=0
		if(freg=="H"):
			print("1")
		if freg[f]=="c" or freg[f]=="C":
			f,freg_c=count_shows(freg, f+1)
		c_num.append(freg_c)
		if(f==len(freg)):
			h_num.append(freg_h)
			o_num.append(freg_o)
			n_num.append(freg_n)
			continue
		if(freg=="H"):
			print("2")		
		if freg[f]=="h" or freg[f]=="H":
			f,freg_h=count_shows(freg, f+1)
		if(f==len(freg)):
			h_num.append(freg_h)
			o_num.append(freg_o)
			n_num.append(freg_n)
			continue
		h_num.append(freg_h)
		if freg[f]=="o" or freg[f]=="O":
			f,freg_o=count_shows(freg, f+1)
		o_num.append(freg_o)
		if(f==len(freg)):
			n_num.append(freg_n)
			continue
		if freg[f]=="n" or freg[f]=="N":
			f,freg_n=count_shows(freg, f+1)
		n_num.append(freg_n)
		
	
	#translate each species into DEDTA and EPON 
	new_label=["D","E"]
	for _ff in range(len(spec_list)):
		if(_ff<2): continue
		num_d=math.floor((n_num[_ff]/n_num[0]))
		num_e=math.floor((o_num[_ff]/o_num[1]))
		bf=""
		sub=c_num[_ff]-num_d*c_num[0]-num_e*c_num[1]
		if (sub>0):
			bf+="C"
			if(sub!=1): bf+=str(sub)
		sub=h_num[_ff]-num_d*h_num[0]-num_e*h_num[1]
		if (sub>0):
			bf+="H"
			if(sub!=1): bf+=str(sub)
		sub=o_num[_ff]-num_d*o_num[0]-num_e*o_num[1]
		if (sub>0):
			bf+="O"
			if(sub!=1): bf+=str(sub)
		sub=n_num[_ff]-num_d*n_num[0]-num_e*n_num[1]
		if (sub>0):
			bf+="N"
			if(sub!=1): bf+=str(sub)
		
		freg_new_label=""
		if (num_d>0):
			freg_new_label+="D"
			if(num_d>1):
				freg_new_label+=str(num_d)
		if (num_e>0):
			freg_new_label+="E"
			if(num_e>1):
				freg_new_label+=str(num_e)
		if(bf!=""):
			if freg_new_label!="":
				freg_new_label+="+"
			freg_new_label+=bf
		new_label.append(freg_new_label)
	print(new_label)


	#for i in range(len(trucker)):
		#make_graph(time_Step, trucker[i], xlabel="time_Step", ylabel=list(species_dict.keys())[i], title="Graph")
	#print(species_dict)
	_freg_name_index=0
	for i in result_list:
		if i==len(text_list)-2: species_at_end.append(j)
		i+=1

	print(species_at_end)
	print(species_dict)
	for i in species_at_end:
		spec_label=list(species_dict.keys())[i]
		plt.plot(time_Step,trucker[i], label=new_label[_freg_name_index])
		x = species_dict.get(spec_label)
		_freg_name_index+=1
	plt.xlabel("timeStep")
	plt.ylabel("Number of molecules of each species")
	plt.legend(loc='upper right')
	plt.title("Number of molecules of each species as a function of time")
	plt.show()
			
species_dict =	{
	  "C11H18N2": 0,
	  "C19H20O4": 1,
	  "C30H38O4N2": 2,
	  "C49H58O8N2": 3
}		
species_reader()
