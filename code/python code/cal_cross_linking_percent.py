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


"""this code calculate the cross linking percent of a simulation at a timestep
	that the user types. the calculation:
	cross linking % = number of created N-C bonds / 4* num of DETDA molecoles.""" 
def cal_cross_linking_percent(time_step_to_cal):
	#gets from the user the timestep he wants to check the cross linking %
	#time_step_to_cal = int (input("Enter timeStep to calculate cross linking % on: ") )
	count_n_c_bond=0 #the number of N-C bonds created
	#open the bonds LAMMPS output file that contain all the bonds of each atom any 500 timesteps
	fp = open("bonds.reax","r") 
	text=fp.read() 
	text_list = text.split("\n")
	text_list=[x for x in text_list if x]
	ln=text_list[2].split(" ")
	ln=[x for x in ln if x]
	num_of_atoms=int(ln[4]) #num of atoms at the simulation
	tag_to_type = [0] * (num_of_atoms+1) #get a type of atom by his tag
	for i in range(num_of_atoms):
		ln=text_list[7+i].split(" ")
		ln=[x for x in ln if x]
		tag_to_type[int(ln[0])]=int(ln[1])

	start_n_c_bond=0 #number of N-C bonds at the start of the run

	#check timestep 0 at the bonds file the number of N-C bonds at the start of the run
	for i in range(num_of_atoms):
		ln=text_list[i+7].split(" ")
		ln=[x for x in ln if x]
		if(len(ln)<7): break
		_id=int(ln[0])
		_type=int(ln[1])
		if _type!=4: continue
		_num_neigh=int(ln[2])
		for n in range(_num_neigh):
			neigh_id=int(ln[3+n])
			if(tag_to_type[neigh_id]==1): start_n_c_bond+=1
		
	#find the bonds table at the chosen timestep	
	for i in range(len(text_list)):
		ln=text_list[i].split(" ")
		ln=[x for x in ln if x]
		if len(ln)>2:
			if ln[1]=="Timestep" and int(ln[2])==time_step_to_cal: break
	if i>=len(text_list)-1:
		print("ERROR, CANNOT FIND TIMESTEP ",time_step_to_cal)
		return
	i+=7
	#check the number of N-C bonds at the chosen timestep
	for j in range(num_of_atoms):
		ln=text_list[i+j].split(" ")
		ln=[x for x in ln if x]
		if(len(ln)<7): break
		_id=int(ln[0])
		_type=int(ln[1])
		if _type!=4: continue
		_num_neigh=int(ln[2])
		for n in range(_num_neigh):
			neigh_id=int(ln[3+n])
			if(tag_to_type[neigh_id]==1): count_n_c_bond+=1
	#print the results
	print("\n=============\nRESULTS:")
	print("num of N-C bond: ",count_n_c_bond)
	print("num of N-C bond at the start: ",start_n_c_bond)
	num_of_detda=num_of_atoms/117
	cross_linking_percent=(count_n_c_bond-start_n_c_bond)/(4*num_of_detda)
	print("num of DETDA molecules: ",num_of_detda)
	print("\nRESULT-cross linking percent: ",cross_linking_percent)
	return cross_linking_percent
	

axis_x=[]
axis_y=[]
for i in range(11):
    x=i*50000
    axis_x.append(x)
    y=cal_cross_linking_percent(x)
    axis_y.append(y)
make_graph(axis_x, axis_y, "timeStep", "cross linking %", "Graph of cross linking % as a function of time")


	
	
