import sys
import numpy as np
from matplotlib import pyplot as plt

fp = open("dist.reax","r") 
text=fp.read() 

text_list = text.split("\n")
timeStep=0
timeStep_index=0
atom_tag=0
atom_type=0
neigh_tag=0
neigh_dist=0

arr=np.zeros((81, 11, 11))
listAtoms=np.zeros(11)
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
		timeStep_index=int(timeStep/50)
		axis_x.append(timeStep)
		continue;
	elif ln[i]=="atom":
		atom_tag=int(ln[i+1])
		atom_type=int(ln[i+3])
		if(timeStep_index==0):
			listAtoms[atom_tag-1]=atom_type
		i+=4
		j=i #iterator for atoms line
		while j < len(ln)-1:
			neigh_tag=int(ln[j])
			j+=1
			neigh_dist=float(ln[j])
			j+=1
			arr[timeStep_index][atom_tag-1][neigh_tag-1]=neigh_dist

print(list(range(11)),"\n", listAtoms)
print(arr)

typedict =	{
  "C": 1,
  "H": 2,
  "O": 3,
  "N": 4
}

axis_y=[ [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [] ]
for j in range(timeStep_index+1):
	axis_y[0].append(arr[j][0][1]) #1-2 dist
	axis_y[1].append(arr[j][0][2]) #1-3 dist
	axis_y[2].append(arr[j][0][3]) #1-4 dist
	axis_y[3].append(arr[j][0][4]) #1-5 dist
	axis_y[4].append(arr[j][0][6]) #1-7 dist
	axis_y[5].append(arr[j][1][3]) #2-4 dist
	axis_y[6].append(arr[j][2][4]) #3-5 dist
	axis_y[7].append(arr[j][2][6]) #3-7 dist
	axis_y[8].append(arr[j][4][5]) #5-6 dist
	axis_y[9].append(arr[j][4][6]) #5-7 dist
	axis_y[10].append(arr[j][4][8]) #5-9 dist
	axis_y[11].append(arr[j][4][10]) #5-11 dist
	axis_y[12].append(arr[j][5][6]) #6-7 dist
	axis_y[13].append(arr[j][5][7]) #6-8 dist
	axis_y[14].append(arr[j][5][9]) #6-10 dist
	axis_y[15].append(arr[j][8][10]) #9-11 dist
	
	
plt.plot(axis_x,axis_y[1], label="1-3")
plt.plot(axis_x,axis_y[3], label="1-5")
plt.plot(axis_x,axis_y[7], label="3-7")
plt.plot(axis_x,axis_y[9], label="5-7", linestyle='--')

plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.xlabel("timeStep")
plt.ylabel("distance")
plt.title("Distance Between Atoms As A Function Of Time")
plt.show()





		
		
		
			
			
		
