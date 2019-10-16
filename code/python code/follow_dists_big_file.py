"""python code for long runs (with many dists files and large amount of atoms)
	that create distance Graphs between N-C, O-H, O-C, H-N pairs for each
	fourset we apply the extra potential on during simulation.
	and by reading log.lammps file and energy.reax file, create 
	tempeture, preassure, potential energy, total energy, added energy by
	the extra potential graphs as a function of time"""

import sys
import math
import numpy as np
from matplotlib import pyplot as plt



"""this code find the timesteps when a reaction happend succefully""" 
def find_reaction_ts(time_step_to_cal, ts_list):
	#gets from the user the timestep he wants to check the cross linking %
	#time_step_to_cal = int (input("Enter timeStep to calculate cross linking % on: ") )
	count_n_c_bonds=0 #the number of N-C bonds created
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

	start_n_c_bonds=0 #number of N-C bonds at the start of the run

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
			if(tag_to_type[neigh_id]==1): start_n_c_bonds+=1
				
	old_count_n_c_bonds=start_n_c_bonds
		
	#find the bonds table at the chosen timestep	
	for i in range(len(text_list)):
		ln=text_list[i].split(" ")
		ln=[x for x in ln if x]
		if len(ln)>2:
			if ln[1]=="Timestep" and (int(ln[2])%time_step_to_cal)==0:
				time_step=int(ln[2])
				i+=7
				count_n_c_bonds=0
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
						if(tag_to_type[neigh_id]==1): count_n_c_bonds+=1
				if(count_n_c_bonds>old_count_n_c_bonds):
					ts_list.append(time_step)
				if(count_n_c_bonds<old_count_n_c_bonds):
					if len(ts_list)>0: del ts_list[-1]
				old_count_n_c_bonds=count_n_c_bonds
				
	
ts_list=[]
find_reaction_ts(1000, ts_list)
print(ts_list)


"""create Energy graph using mathplot. mark the timesteps when reaction happend.
	input- axis x, axis y, label for each axis and title for the graph"""
def make_E_graph(axis_x, axis_y, xlabel="x", ylabel="y", title="Graph"):
	plt.plot(axis_x,axis_y)
	for ts in ts_list: plt.axvline(x=ts, color='#d62728')
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.title(title)
	plt.show()
	
	
"""create Energy graph using mathplot
	input- axis x, axis y, label for each axis and title for the graph"""
def make_graph(axis_x, axis_y, xlabel="x", ylabel="y", title="Graph"):
	plt.plot(axis_x,axis_y)
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.title(title)
	plt.show()
	
"""this function recieve fourset tags to create array that document the distance
	between any 2 atoms at any timestep to create the axis_y array, each cell is an axis y.
	each 4 axis_y are distances between 4 pairs of this fourset.
	at the end create a graph of the distances between 4 pairs as a function of time"""	
def create_dist_nparray(o_tag,h_tag,n_tag,c_tag, fourset_ts, text_list, total_timesteps, nevery_file_dists, fourset_index, save_flag):
	atom_tag=0
	neigh_tag=0
	neigh_dist=0
	axis_x=[]
	axis_y=[[],[],[],[]]
	timestep_counter=0
	i_on_textlist=0
	end_of_files_flag=1
	curr_file=0
	while end_of_files_flag==1:
		i_on_textlist+=1
		if i_on_textlist==len(text_list):
			curr_file+=nevery_file_dists
			file_to_open="dists.reax."+str(curr_file)
			fp = open(file_to_open,"r") 
			text_list = fp.read().split("\n")
			i_on_textlist=0
		line=text_list[i_on_textlist]
		ln=line.split(" ")
		i=0 #iterator on the line
		if ln[i]=="#": i+=1
		if(i>=len(ln)): continue
		elif ln[i]=="Timestep":
			if(int(ln[i+1])==total_timesteps): end_of_files_flag=0
			else:
				axis_x.append(int(ln[i+1]))
				timestep_counter+=1
				if len(axis_y[0])>0:
				    axis_y[0].append(axis_y[0][-1])
				    axis_y[1].append(axis_y[1][-1])
				    axis_y[2].append(axis_y[2][-1])
				    axis_y[3].append(axis_y[3][-1])
				else:
				    axis_y[0].append(10)
				    axis_y[1].append(10)
				    axis_y[2].append(10)
				    axis_y[3].append(10)
			continue
		elif ln[i]=="atom":
			atom_tag=int(ln[i+1])
			i+=4
			j=i #iterator for atoms line
			if(atom_tag==c_tag):
				while j < len(ln)-1:
					neigh_tag=int(ln[j])
					j+=1
					neigh_dist=float(ln[j])
					j+=1
					if(neigh_tag==o_tag): axis_y[0][timestep_counter-1]=(neigh_dist)
					if(neigh_tag==n_tag): axis_y[3][timestep_counter-1]=(neigh_dist)
						
			elif(atom_tag==n_tag):
				while j < len(ln)-1:
					neigh_tag=int(ln[j])
					j+=1
					neigh_dist=float(ln[j])
					j+=1
					if(neigh_tag==h_tag):
						axis_y[2][timestep_counter-1]=(neigh_dist)
					if(neigh_tag==c_tag):
						axis_y[3][timestep_counter-1]=(neigh_dist)
			
			elif(atom_tag==o_tag):
				while j < len(ln)-1:
					neigh_tag=int(ln[j])
					j+=1
					neigh_dist=float(ln[j])
					j+=1
					if(neigh_tag==h_tag): axis_y[1][timestep_counter-1]=(neigh_dist)

	plt.plot(axis_x,axis_y[0], label="C("+str(c_tag)+")-O("+str(o_tag)+")")
	plt.plot(axis_x,axis_y[1], label="H("+str(h_tag)+")-O("+str(o_tag)+")")
	plt.plot(axis_x,axis_y[2], label="H("+str(h_tag)+")-N("+str(n_tag)+")")
	plt.plot(axis_x,axis_y[3], label="C("+str(c_tag)+")-N("+str(n_tag)+")")
	plt.legend(loc='upper right')
	plt.xlabel("timeStep")
	plt.ylabel("distance")
	plt.title("Distance Between Atoms As A Function Of Time, Extra potential started at "+str(fourset_ts))
	if save_flag==1:
		pic_name="fourset_"+str(fourset_index)+".png"
		fig = plt.gcf()
		fig.set_size_inches((15, 8), forward=False)
		fig.savefig(pic_name)
		plt.clf()
	else:
		plt.show()
	return axis_y
	
""" creates distance graph for each fourset that we apply the extra potential on duting the simulation run """
def dist_graph():
	fp = open("dists.reax","r") 
	text=fp.read() 
	text_list = text.split("\n")

	ln=text_list[0].split(" ")
	if ln[1] == "totalTimesteps":
		total_timesteps=int(ln[2])

	ln=text_list[1].split(" ")
	if ln[1] == "totalAtomNum":
		total_atomNum=int(ln[2])
		print("total_atomNum ", total_atomNum)
	
	nevery_dists_follow=10
	ln=text_list[2].split(" ")
	if ln[1] == "nevery_dists_follow":
		nevery_dists_follow=int(ln[2])
		print("nevery_dists_follow ", nevery_dists_follow)

	nevery_file_dists=100000
	ln=text_list[2].split(" ")
	if ln[1] == "nevery_file_dists":
		nevery_file_dists=int(ln[2])
		print("nevery_file_dists ", nevery_file_dists)

	fourset=[]
	fourset_timestep=[]
	num_of_foursets=0
	i_on_textlist=0
	end_of_files_flag=1
	curr_file=0
	while end_of_files_flag==1:
		i_on_textlist+=1
		if i_on_textlist==len(text_list):
			curr_file+=nevery_file_dists
			file_to_open="dists.reax."+str(curr_file)
			fp = open(file_to_open,"r") 
			text_list = fp.read().split("\n")
			i_on_textlist=0
		line=text_list[i_on_textlist]
		ln=line.split(" ")
		i=0 #iterator on the line
		if ln[i]=="#": i+=1
		if(i>=len(ln)): continue
		if ln[i]=="Timestep":
			#if int(ln[i+1])<fff_start_point: fff_start+=1
			#if int(ln[i+1])<fff_end_point: fff_end+=1
			if int(ln[i+1])==total_timesteps: end_of_files_flag=0
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

	print(num_of_foursets)
	#show_only_success(fourset,fourset_timestep,total_timesteps,nevery_file_dists)
	#create graph for each fourset
	for i in range(len(fourset)):
		fp = open("dists.reax","r") 
		text_list = fp.read().split("\n")
		create_dist_nparray(fourset[i][0],fourset[i][1],fourset[i][2],fourset[i][3], fourset_timestep[i], text_list, total_timesteps, nevery_file_dists,i+1,1)


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
	make_E_graph(x,y, "TimeStep", "added energy", "Additional Energy As A Function Of Time")


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
	make_E_graph(timeStep_arr, potE_arr, "TimeStep", "Potential Energy", "Potential Energy As A Function Of Time")
	make_graph(timeStep_arr, press_arr, "TimeStep", "Pressure", "Pressure As A Function Of Time")
	make_E_graph(timeStep_arr, totalE_arr, "TimeStep", "Total Energy", "Total Energy As A Function Of Time")


def show_only_success(fourset,fourset_timestep,total_timesteps,nevery_file_dists):
    arr=[3,5,14,16,17,18,22,38,40,49,68,77,84,94,99,106,108,110,137,140,144,145,153,161,166,174,176]
    for item in arr:
        fp = open("dists.reax","r") 
        text_list = fp.read().split("\n")
        i=item-1
        create_dist_nparray(fourset[i][0],fourset[i][1],fourset[i][2],fourset[i][3], fourset_timestep[i], text_list, total_timesteps,nevery_file_dists,i+1,0)


		
#operate functions as you want to create the graphs you want.	
#dist_graph()
extraE_graph()
log_graphs()
