import sys
import math
import numpy as np
from matplotlib import pyplot as plt


def make_graph(axis_x, axis_y, xlabel="x", ylabel="y", title="Graph"):
	plt.plot(axis_x,axis_y)
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.title(title)
	plt.show()
	

def dist_graph():
	fp = open("dists.reax","r") 
	text=fp.read() 

	text_list = text.split("\n")
	timeStep=0
	timeStep_index=0
	atom_tag=0
	atom_type=0
	neigh_tag=0
	neigh_dist=0

	ln=text_list[0].split(" ")
	if ln[1] == "totalTimesteps":
		total_timesteps=int(ln[2])

	neverty_fix_check=10

	size=int(total_timesteps/neverty_fix_check) + 1

	arr=np.zeros((size, 117, 117))
	axis_x=[]

	fourset=[]
	fourset_timestep=[]
	num_of_foursets=0
	for line in text_list:
		ln=line.split(" ")
		i=0 #iterator on the line
		if ln[i]=="#": i+=1
		if(i>=len(ln)): continue
		elif ln[i]=="Timestep":
			timeStep=int(ln[i+1])
			timeStep_index=int(timeStep/neverty_fix_check)
			axis_x.append(timeStep)
			continue
		elif ln[i]=="atom":
			atom_tag=int(ln[i+1])
			atom_type=int(ln[i+3])
			i+=4
			j=i #iterator for atoms line
			while j < len(ln)-1:
				neigh_tag=int(ln[j])
				j+=1
				neigh_dist=float(ln[j])
				j+=1
				arr[timeStep_index][atom_tag-1][neigh_tag-1]=neigh_dist
		elif ln[i]=="fourset":
			i+=7
			fs_ts=int(ln[i]) #current fourset timestep
			i+=2
			fourset.append([int(x) for x in ln[i+1:i+5]])
			num_of_foursets+=1
			fourset_timestep.append(fs_ts)
			if ln[i]=="1/2-":
				i+=5
				fourset.append([int(x) for x in ln[i+1:i+5]])
				fourset_timestep.append(fs_ts)
				num_of_foursets+=1


	o_list, h_list, n_list, c_list= [],[],[],[]
	for list in fourset:
		o_list.append(list[0])
		h_list.append(list[1])
		n_list.append(list[2])
		c_list.append(list[3])


	typedict =	{
	  "C": 1,
	  "H": 2,
	  "O": 3,
	  "N": 4
	}

	axis_y=[]
	for j in range(4*num_of_foursets):
		axis_y.append([])
	print(axis_y)
	for i in range(len(fourset_timestep)):
		#fourset 1
		for j in range(timeStep_index+1):
			index=i*4
			temp_dist=arr[j][c_list[i]-1][o_list[i]-1]
			if temp_dist==0:
				temp_dist=arr[j][o_list[i]-1][c_list[i]-1]
			if temp_dist==0: temp_dist=10
			axis_y[index].append(temp_dist) #C-O dist
			
			temp_dist=arr[j][o_list[i]-1][h_list[i]-1]
			if temp_dist==0:
				temp_dist=arr[j][h_list[i]-1][o_list[i]-1]
			if temp_dist==0: temp_dist=10
			axis_y[index+1].append(temp_dist) #O-H dist
			
			temp_dist=arr[j][n_list[i]-1][h_list[i]-1]
			if temp_dist==0:
				temp_dist=arr[j][h_list[i]-1][n_list[i]-1]
			if temp_dist==0: temp_dist=10
			axis_y[index+2].append(temp_dist) #N-H dist
			
			temp_dist=arr[j][n_list[i]-1][c_list[i]-1]
			if temp_dist==0:
				temp_dist=arr[j][c_list[i]-1][n_list[i]-1]
			if temp_dist==0: temp_dist=10
			axis_y[index+3].append(temp_dist) #N-C dist
			


	for i in range(len(fourset_timestep)):
		index=i*4
		plt.plot(axis_x,axis_y[index], label="C("+str(c_list[i])+")-O("+str(o_list[i])+")")
		plt.plot(axis_x,axis_y[index+1], label="H("+str(h_list[i])+")-O("+str(o_list[i])+")")
		plt.plot(axis_x,axis_y[index+2], label="H("+str(h_list[i])+")-N("+str(n_list[i])+")")
		plt.plot(axis_x,axis_y[index+3], label="C("+str(c_list[i])+")-N("+str(n_list[i])+")")
		plt.legend(loc='upper right')
		plt.xlabel("timeStep")
		plt.ylabel("distance")
		plt.title("Distance Between Atoms As A Function Of Time, Extra potential started at "+str(fourset_timestep[i]))
		plt.show()
	#, linestyle='--'


def extraE_graph():
	e_fp = open("energy.reax","r") 
	text=e_fp.read() 
	text_list = text.split("\n")

	text_list=text_list[1:]
	#print(text_list)
	y=[]
	for e in text_list:
		if e == "finish": y.append(0.0)
		elif e == "start": continue
		else: y.append(float(e))

	#x=list(range(len(y)))
	#axis_x=list(range(1, len(y)))
	#print(axis_x)

	x= [i for i in range(len(y))]
	make_graph(x,y, "TimeStep", "added energy", "Additional Energy As A Function Of Time")

	
def log_graphs():
	fp = open("log.lammps","r") 
	text=fp.read() 

	text_list = text.split("\n")

	i=0
	for line in text_list:
		ln=line.split(" ")
		i+=1
		if(ln[0]=="Step"): break

	timeStep_arr=[]
	potE_arr=[]
	totalE_arr=[]
	temp_arr=[]
	press_arr=[]
	for line in text_list[i:]:
		ln=line.split(" ")
		ln = [x for x in ln if x != ""]
		if ln[0]=="Loop": break
		timeStep_arr.append(int(ln[0]))
		temp_arr.append(float(ln[1]))
		potE_arr.append(float(ln[2]))
		totalE_arr.append(float(ln[3]))
		press_arr.append(float(ln[4]))

	make_graph(timeStep_arr, temp_arr, "TimeStep", "Temprature", "Temprature As A Function Of Time")
	make_graph(timeStep_arr, potE_arr, "TimeStep", "Potential Energy", "Potential Energy As A Function Of Time")
	make_graph(timeStep_arr, press_arr, "TimeStep", "Pressure", "Pressure As A Function Of Time")
	make_graph(timeStep_arr, totalE_arr, "TimeStep", "Total Energy", "Total Energy As A Function Of Time")
	
def remove(string): 
    return string.replace("\t", "")  
	
def species_reader():
	i=0
	species_dict =	{
	  "C11H18N2": 0,
	  "C19H20O4": 1,
	  "C30H38O4N2": 2,
	  "C49H58O8N2": 3
	}
	time_Step=[]
	No_Moles=[]
	No_Specs=[]
	species=[]
	trucker=[[],[],[],[]]
	fp = open("species.out","r") 
	text=fp.read() 
	text_list = text.split("\n")
	while i<len(text_list)-1:
		ln = text_list[i].split(" ")
		ln=[x for x in ln if x]
		species_list=ln[4].split("\t")
		species_list=[x for x in species_list if x]
		i+=1
		ln = text_list[i].split(" ")
		ln=[remove(x) for x in ln if x]
		time_Step.append(int(ln[0]))
		No_Moles.append(ln[1])
		No_Specs.append(ln[2])
		for j in range(len(species_list)):
			x = species_dict.get(species_list[j])
			if x is not None:
				x=int(x)
				trucker[x].append(ln[3+j])
		for j in range(4):
			if len(trucker[j])<len(time_Step): trucker[j].append(0)
		i+=1
	make_graph(time_Step, trucker[0], xlabel="time_Step", ylabel="C11H18N2", title="Graph")
	make_graph(time_Step, trucker[1], xlabel="time_Step", ylabel="C19H20O4", title="Graph")
	make_graph(time_Step, trucker[2], xlabel="time_Step", ylabel="C30H38O4N2", title="Graph")
	make_graph(time_Step, trucker[3], xlabel="time_Step", ylabel="C49H58O8N2", title="Graph")
		
			
			
			
	
	
		
species_reader()	
#dist_graph()
#extraE_graph()
#log_graphs()
species_reader()

