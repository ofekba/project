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
	
	
def create_dist_nparray(filename):
	fp = open(filename,"r") 
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
		total_timesteps=time_steps_to_check

	ln=text_list[1].split(" ")
	if ln[1] == "totalAtomNum":
		total_atomNum=int(ln[2])
	
	
	ln=text_list[2].split(" ")
	if ln[1] == "fix_nevery":
		neverty_fix_check=int(ln[2])
	neverty_fix_check=10

	size=int(total_timesteps/neverty_fix_check) + 1


	axis_y=[[],[],[],[]]
	fourset=[]
	fourset_timestep=[]
	timestep_counter=0;
	num_of_foursets=0
	for line in text_list:
		ln=line.split(" ")
		i=0 #iterator on the line
		if ln[i]=="#": i+=1
		if(i>=len(ln)): continue
		elif ln[i]=="Timestep":
			timeStep=int(ln[i+1])
			timestep_counter+=1
			if timeStep>time_steps_to_check: break
			axis_y[0].append(10)
			axis_y[1].append(10)
			axis_y[2].append(10)
			axis_y[3].append(10)
			continue
		elif ln[i]=="atom":
			atom_tag=int(ln[i+1])
			atom_type=int(ln[i+3])
			i+=4
			j=i #iterator for atoms line
			if(atom_tag==1263):
				while j < len(ln)-1:
					neigh_tag=int(ln[j])
					j+=1
					neigh_dist=float(ln[j])
					j+=1
					if(neigh_tag==1255): axis_y[0][timestep_counter-1]=(neigh_dist)
					if(neigh_tag==1726): axis_y[3][timestep_counter-1]=(neigh_dist)
						
			
			elif(atom_tag==1726):
				while j < len(ln)-1:
					neigh_tag=int(ln[j])
					j+=1
					neigh_dist=float(ln[j])
					j+=1
					if(neigh_tag==1748):
						axis_y[2][timestep_counter-1]=(neigh_dist)
					if(neigh_tag==1263):
						axis_y[3][timestep_counter-1]=(neigh_dist)
			
			elif(atom_tag==1255):
				while j < len(ln)-1:
					neigh_tag=int(ln[j])
					j+=1
					neigh_dist=float(ln[j])
					j+=1
					if(neigh_tag==1748): axis_y[1][timestep_counter-1]=(neigh_dist)
			
		elif ln[i]=="fourset":
			i+=7
			fs_ts=int(ln[i]) #current fourset timestep
			print("found fourset of ", filename, " at ", fs_ts)
			print(ln[10:])
			i+=2
			fourset.append([int(x) for x in ln[i+1:i+5]])
			num_of_foursets+=1
			fourset_timestep.append(fs_ts)
			if ln[i]=="1/2-":
				i+=5
				fourset.append([int(x) for x in ln[i+1:i+5]])
				fourset_timestep.append(fs_ts)
				num_of_foursets+=1
	print("\n")
	return axis_y
	
#=================================================================

time_steps_to_check=22000
axis_x=[x*10 for x in range(2201)]
axis_y_1=create_dist_nparray("dists1.reax")
axis_y_2=create_dist_nparray("dists2.reax")
axis_y_4=create_dist_nparray("dists4.reax")
axis_y_8=create_dist_nparray("dists8.reax")
axis_y_12=create_dist_nparray("dists12.reax")
axis_y_s=create_dist_nparray("distss.reax")
title=["C-O dist","O-H dist","H-N dist","C-N dist"]
for i in range(4):
	plt.plot(axis_x,axis_y_1[i], label="1 threads")
	#plt.plot(axis_x,axis_y_2[i], label="2 threads")
	#plt.plot(axis_x,axis_y_4[i], label="4 threads")
	#plt.plot(axis_x,axis_y_8[i], label="8 threads")
	#plt.plot(axis_x,axis_y_12[i], label="12 threads")
	plt.plot(axis_x,axis_y_s[i], label="serial")
	plt.legend(loc='upper right')
	plt.xlabel("timeStep")
	plt.ylabel("distance")
	plt.title(title[i])
	plt.show()

