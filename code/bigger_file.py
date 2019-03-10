import sys
import numpy as np
from matplotlib import pyplot as plt

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
for line in text_list:
	ln=line.split(" ")
	i=0 #iterator on the line
	if ln[i]=="#":
		i+=1
	if(i>=len(ln)):
		continue;
	elif ln[i]=="Timestep":
		timeStep=int(ln[i+1])
		timeStep_index=int(timeStep/neverty_fix_check)
		axis_x.append(timeStep)
		continue;
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
		fourset_timestep.append(fs_ts)
		if ln[i]=="1/2-":
			i+=5
			fourset.append([int(x) for x in ln[i+1:i+5]])
			fourset_timestep.append(fs_ts)
		
		
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

axis_y=[ [], [], [], [], [], [], [], [],[], [], [], []]
for i in range(len(fourset_timestep)):
	#fourset 1
	for j in range(timeStep_index+1):
		index=i*4
		axis_y[index].append(arr[j][c_list[i]-1][o_list[i]-1]) #C-O dist
		axis_y[index+1].append(arr[j][o_list[i]-1][h_list[i]-1]) #O-H dist
		axis_y[index+2].append(arr[j][n_list[i]-1][h_list[i]-1]) #N-H dist
		axis_y[index+3].append(arr[j][n_list[i]-1][c_list[i]-1]) #N-C dist
	

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


e_fp = open("energy.reax","r") 
text=e_fp.read() 
text_list = text.split("\n")

text_list=text_list[1:]
#print(text_list)
y=[]
for e in text_list:
	if e == "finish":
		y.append(0.0)
	elif e == "start":
		continue;
	else:
		y.append(float(e))
	
#x=list(range(len(y)))
#axis_x=list(range(1, len(y)))
#print(axis_x)

x= [i for i in range(len(y))]
plt.plot(x,y)
plt.xlabel("timeStep")
plt.ylabel("energy")
plt.title("Additional Energy As A Function Of Time")
plt.show()

	
