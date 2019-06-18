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


def species_reader():
	i=0
	num_species=4
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
			if x is None:
				species_dict.update( {species_list[j] : num_species} )
				trucker.append([])
				trucker[num_species]=[0 for ns in range(len(time_Step)-1)]
				trucker[num_species].append(ln[3+j])
				num_species+=1
			else:
				x=int(x)
				trucker[x].append(ln[3+j])
		for j in range(len(trucker)):
			if len(trucker[j])<len(time_Step): trucker[j].append(0)
		i+=1

		
	for i in range(len(trucker)):
		make_graph(time_Step, trucker[i], xlabel="time_Step", ylabel=list(species_dict.keys())[i], title="Graph")
		
			
species_dict =	{
	  "C11H18N2": 0,
	  "C19H20O4": 1,
	  "C30H38O4N2": 2,
	  "C49H58O8N2": 3
}		
species_reader()