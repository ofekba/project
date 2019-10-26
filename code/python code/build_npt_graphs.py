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

	
	
"""create Energy graph using mathplot
	input- axis x, axis y, label for each axis and title for the graph"""
def make_graph(axis_x, axis_y, xlabel="x", ylabel="y", title="Graph", name="graph.png"):
	plt.plot(axis_x,axis_y)
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.title(title)
	plt.savefig(name)
	plt.show()
	

	

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
	volume_arr=[]
	density_arr=[]
	for line in text_list[i:]:
		ln=line.split(" ")
		ln = [x for x in ln if x != ""]
		if ln[0]=="Loop": break
		timeStep_arr.append(int(ln[0]))
		temp_arr.append(float(ln[1]))
		potE_arr.append(float(ln[2]))
		totalE_arr.append(float(ln[3]))
		press_arr.append(float(ln[4]))
		volume_arr.append(float(ln[5]))
		density_arr.append(float(ln[6]))

	make_graph(timeStep_arr, temp_arr, "TimeStep", "Temprature", "Temprature As A Function Of Time","temp.png")
	make_graph(timeStep_arr, press_arr, "TimeStep", "Pressure", "Pressure As A Function Of Time","press.png")
	make_graph(timeStep_arr, potE_arr, "TimeStep", "Potential Energy", "Potential Energy As A Function Of Time","potE.png")
	make_graph(timeStep_arr, totalE_arr, "TimeStep", "Total Energy", "Total Energy As A Function Of Time","totE.png")
	make_graph(timeStep_arr, volume_arr, "TimeStep", "Volume", "Volume As A Function Of Time","volume.png")
	make_graph(timeStep_arr, density_arr, "TimeStep", "Density", "Density As A Function Of Time","density.png")




		
#operate functions as you want to create the graphs you want.	
log_graphs()