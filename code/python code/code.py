"""python code for short runs (with only one dists file and small amount of atoms
that create distance Graphs between N-C, O-H, O-C, H-N pairs for each
	fourset we apply the extra potential on during simulation.
	and by reading log.lammps file and energy.reax file, create 
	tempeture, preassure, potential energy, total energy, added energy by
	the extra potential graphs as a function of time"""

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
	plt.title(title)
	plt.show()
	
"""create distance graph between each N-C, O-H, O-C, H-N  pairs for each
	fourset we apply the extra potential on during simulation"""
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

	ln=text_list[1].split(" ")
	if ln[1] == "totalAtomNum":
		total_atomNum=int(ln[2])
		print("total_atomNum ", total_atomNum)
	
	
	neverty_fix_check=10 #default val
	ln=text_list[2].split(" ")
	if ln[1] == "fix_nevery":
		neverty_fix_check=int(ln[2])
		print("neverty_fix_check ", neverty_fix_check)


	size=int(total_timesteps/neverty_fix_check) + 1
	#array that document the distance between any 2 atoms at any timestep
	#arr[timestpe][i][j]=the distance between atom i, atom j at time step.
	arr=np.zeros((size, total_atomNum, total_atomNum))
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
			#create timesteps axis- axis_x
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
			#create foursets array that document the foursets we apply the extra potential on
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

	#seperate the foursets array by type
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
	#create the axis_y array, each cell is an axis y.
	#each 4 axis_y are distances between 4 pairs of each fourset.
	axis_y=[]
	for j in range(4*num_of_foursets):
		axis_y.append([])
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
			

	#print graphs of distances between 4 pairs for each fourset.
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


"""create Graph of the added energy by the extra potential as a function of time"""
def extraE_graph():
	e_fp = open("energy.reax","r") 
	text=e_fp.read() 
	text_list = text.split("\n")
	text_list=text_list[1:]
	y=[]
	for e in text_list:
		if e == "finish": y.append(0.0)
		elif e == "start": continue
		else: y.append(float(e))
	x= [i for i in range(len(y))]
	make_graph(x,y, "TimeStep", "added energy", "Additional Energy As A Function Of Time")


"""create Graphs of the potential energy, total energy, tempeture and pressure by the extra potential as a function of time"""	
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
	
		
#operate functions as you want to create the graphs you want.	
dist_graph()
#extraE_graph()
#log_graphs()

