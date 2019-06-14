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
	

def create_dist_nparray(o_tag,h_tag,n_tag,c_tag, fourset_ts, text_list):
	atom_tag=0
	atom_type=0
	neigh_tag=0
	neigh_dist=0
	axis_x=[]
	axis_y=[[],[],[],[]]
	timestep_counter=0;
	num_of_foursets=0
	for line in text_list:
		ln=line.split(" ")
		i=0 #iterator on the line
		if ln[i]=="#": i+=1
		if(i>=len(ln)): continue
		elif ln[i]=="Timestep":
			axis_x.append(int(ln[i+1]))
			timestep_counter+=1
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
	plt.show()
	return axis_y
	

def dist_graph():
	fp = open("dists.reax","r") 
	text=fp.read() 
	
	fff = open("demofile2.txt", "a")
	fff_start=0
	fff_end=0
	fff.write("Now the file has more content!\n")

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
	
	
	ln=text_list[2].split(" ")
	if ln[1] == "fix_nevery":
		neverty_fix_check=int(ln[2])
		print("neverty_fix_check ", neverty_fix_check)
	neverty_fix_check=10

	size=int(total_timesteps/neverty_fix_check) + 1
	fourset=[]
	fourset_timestep=[]
	num_of_foursets=0
	for line in text_list:
		ln=line.split(" ")
		i=0 #iterator on the line
		if ln[i]=="#": i+=1
		if(i>=len(ln)): continue
		if ln[i]=="Timestep":
			if int(ln[i+1])<90000: fff_start+=1
			if int(ln[i+1])<110000: fff_end+=1
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
	
	axis_y=[]
	for ln in text_list[fff_start:fff_end]:
		fff.write(ln)
	
	for i in range(len(fourset)):
		create_dist_nparray(fourset[i][0],fourset[i][1],fourset[i][2],fourset[i][3], fourset_timestep[i], text_list)


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
	
		
	
dist_graph()
#extraE_graph()
#log_graphs()

