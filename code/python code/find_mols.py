import sys
import math
import numpy as np

fp = open("bond.txt","r") 
text=fp.read() 
text_list = text.split("\n")
i=0
mol_list=[1137,1133,1777,1756]
neigh_list=list(range(1873))
type_list=list(range(1873))

while i<len(text_list):
		ln = text_list[i].split(" ")
		ln=[x for x in ln if x]
		if len(ln)==0: break
		if i==0: print(ln)
		i+=1
		atom_tag=int(ln[0])
		atom_type=int(ln[1])
		type_list[atom_tag]=atom_type
		atom_num_neigh=int(ln[2])
		neigh_list[atom_tag]=[int(x) for x in ln[3:3+atom_num_neigh]]


def create_mol_list(mol_list):
	mol_type_list=[0,0,0,0,0]
	i=0
	found_atom=list(range(1873))
	while i<len(mol_list):
		atom=mol_list[i]
		i+=1
		if(found_atom[atom]==-1): continue
		atom_neigh_list=neigh_list[atom]
		mol_type_list[type_list[atom]]+=1
		found_atom[atom]=-1
		for x in atom_neigh_list: mol_list.append(x)
	mol_list=set(mol_list)
	mol_list=list(mol_list)
	mol_list.sort()
	return mol_list,mol_type_list

mol_list_1=[1137,1133]
mol_list_1, mol_type_1=create_mol_list(mol_list_1)
mol_list_2=[1777,1756]
mol_list_2, mol_type_2=create_mol_list(mol_list_2)
print(mol_type_1)
print(mol_list_1)
_sum=0
for i in mol_type_1: _sum+=i
print(_sum, len(mol_list_1))
print("===============")
print(mol_type_2)
print(mol_list_2)
_sum=0
for i in mol_type_2: _sum+=i
print(_sum, len(mol_list_2))

mol_list=mol_list_1+mol_list_2
mol_list.sort()
print("len", len(mol_list))
for index in range(len(mol_list)):
	if mol_list[index]==1137: print ("c new tag is",(index+1))
	if mol_list[index]==1133: print ("o new tag is",(index+1))
	if mol_list[index]==1777: print ("h new tag is",(index+1))
	if mol_list[index]==1756: print ("n new tag is",(index+1))
iter_mol=0
iter_file=0
fp = open("dump.nvt.1170000.xyz","r")
text=fp.read() 
text_list = text.split("\n")
item_to_del=[]

while iter_file<len(text_list):
	if iter_mol>=len(mol_list):
	   item_to_del+=range(iter_file,1873)
	   break
	if iter_file+1==mol_list[iter_mol]:
		iter_file+=1
		iter_mol+=1
		continue
	item_to_del.append(iter_file)
	iter_file+=1
print(len(item_to_del))
item_to_del.sort(reverse=True)
for i in item_to_del: del text_list[i]
print("text_list len", len(text_list))
fp.close()
fp = open("dump.nvt.1170000.xyz","w")
for tx in text_list: fp.write(tx+"\n") 
fp.close()
		




