import sys
import math
import numpy as np
from matplotlib import pyplot as plt


def cal_cross_linking_percent():
	time_step_to_cal = int (input("Enter timeStep to calculate cross linking % on: ") )
	count_n_c_bond=0
	fp = open("bonds.reax","r") 
	text=fp.read() 
	text_list = text.split("\n")
	text_list=[x for x in text_list if x]
	ln=text_list[2].split(" ")
	ln=[x for x in ln if x]
	num_of_atoms=int(ln[4])
	tag_to_type = [0] * (num_of_atoms+1)
	for i in range(num_of_atoms):
		ln=text_list[7+i].split(" ")
		ln=[x for x in ln if x]
		tag_to_type[int(ln[0])]=int(ln[1])

	start_n_c_bond=0
	for i in range(num_of_atoms):
		ln=text_list[i+7].split(" ")
		ln=[x for x in ln if x]
		if(len(ln)<7): break
		#id type nb id_1...id_nb mol bo_1...bo_nb abo nlp q
		_id=int(ln[0])
		_type=int(ln[1])
		if _type!=4: continue
		_num_neigh=int(ln[2])
		for n in range(_num_neigh):
			neigh_id=int(ln[3+n])
			if(tag_to_type[neigh_id]==1): start_n_c_bond+=1
		
		
	for i in range(len(text_list)):
		ln=text_list[i].split(" ")
		ln=[x for x in ln if x]
		if len(ln)>2:
			if ln[1]=="Timestep" and int(ln[2])==time_step_to_cal: break
	if i>=len(text_list)-1:
		print("ERROR, CANNOT FIND TIMESTEP ",time_step_to_cal)
		return
	i+=7
	for j in range(num_of_atoms):
		ln=text_list[i+j].split(" ")
		ln=[x for x in ln if x]
		if(len(ln)<7): break
		#id type nb id_1...id_nb mol bo_1...bo_nb abo nlp q
		_id=int(ln[0])
		_type=int(ln[1])
		if _type!=4: continue
		_num_neigh=int(ln[2])
		for n in range(_num_neigh):
			neigh_id=int(ln[3+n])
			if(tag_to_type[neigh_id]==1): count_n_c_bond+=1
	print("=============\nRESULTS:")
	print("num of N-C bond: ",count_n_c_bond)
	print("num of N-C bond at the start: ",start_n_c_bond)
	num_of_detda=num_of_atoms/117
	cross_linking_percent=(count_n_c_bond-start_n_c_bond)/(4*num_of_detda)
	print("num of DETDA molecules: ",num_of_detda)
	print("cross linking percent: ",cross_linking_percent)
	

cal_cross_linking_percent()
	
	
