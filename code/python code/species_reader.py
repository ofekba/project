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