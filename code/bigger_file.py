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

neverty_fix_check=1
total_timesteps=7000

size=int(total_timesteps/neverty_fix_check) + 1

arr=np.zeros((size, 117, 117))
axis_x=[]


#print(text_list)

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

#print(list(range(11)),"\n")
#print(arr)

typedict =	{
  "C": 1,
  "H": 2,
  "O": 3,
  "N": 4
}

axis_y=[ [], [], [], [], []]
for j in range(timeStep_index+1):
	axis_y[0].append(arr[j][91][95]) #C-O dist
	axis_y[1].append(arr[j][95][27]) #O-H dist
	axis_y[2].append(arr[j][27][7]) #N-H dist
	axis_y[3].append(arr[j][7][91]) #N-C dist
	axis_y[4].append(arr[j][7][90]) #N-C91 dist

plt.plot(axis_x[:int(5500/neverty_fix_check)],axis_y[2][:int(5500/neverty_fix_check)], label="N-H")
plt.plot(axis_x[:int(5500/neverty_fix_check)],axis_y[1][:int(5500/neverty_fix_check)], label="H-O")
plt.plot(axis_x[:int(5500/neverty_fix_check)],axis_y[0][:int(5500/neverty_fix_check)], label="C-O")
plt.plot(axis_x[:int(5500/neverty_fix_check)],axis_y[3][:int(5500/neverty_fix_check)], label="N-C")
plt.plot(axis_x[:int(5500/neverty_fix_check)],axis_y[4][:int(5500/neverty_fix_check)], label="N-C91")



plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.xlabel("timeStep")
plt.ylabel("distance")
plt.title("Distance Between Atoms As A Function Of Time")
plt.show()		
		
