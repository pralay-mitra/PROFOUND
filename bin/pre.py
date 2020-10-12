import numpy
import math
import re
import os
import sys 

argumentList = sys.argv 
exc_dir=sys.argv[0]
exc_dir_list = exc_dir.strip().split("/")
exc_dir_list = exc_dir_list[:-2]
exc_dir1 = "/".join(str(a) for a in exc_dir_list)


with open (exc_dir1+"/runinfo.txt") as f:
    lines = f.readlines()
    for line in lines:
        name,path = line.strip().split()
        if (name=="Workdir"):
            WorkdirPath = path
        elif(name=="Inpdir"):
            InpdirPath=path
        elif(name=="Bindir"):
            BindirPath=path
        elif(name=="Stride"):
            StridePath=path
        elif(name=="HBplus"):
            HBplusPath=path
        elif(name=="Foldx"):
            FoldxPath=path

f.close()

listPDB = []
pattern = "*.pdb"
contents = os.listdir(InpdirPath)
for content in contents:
    if(content.endswith(".pdb") or content.endswith(".PDB") ):
        listPDB.append(content)
listPDB.sort()
renumPdbPath=WorkdirPath+"/Updated_INPUT"
renum_folder_name="mkdir "+renumPdbPath
os.system(renum_folder_name)
f1=open(WorkdirPath+"/"+"pre_data.txt","w")

for pdb_file_name in range(len(listPDB)):
    pdbid,ext = listPDB[pdb_file_name].strip().split('.')
    pdbid_list= pdbid.strip().split('_')
    chain_name=pdbid_list[0][-1:]

    f=open(InpdirPath+"/"+listPDB[pdb_file_name],"r")
    lines=f.readlines()
    f.close()

    f=open(renumPdbPath+"/"+listPDB[pdb_file_name],"w")
    for line in lines:
        if(line[:4]=="ATOM" and chain_name == line[21:22].strip()):
            f.write(line)
    f.close()

    f=open(renumPdbPath+"/"+listPDB[pdb_file_name],"r")
    lines=f.readlines()
    f.close()


    temp_count=int(lines[0][22:26].strip())
    count=int(lines[0][22:26].strip())-1
    flag=1
    for line_index in lines:
        temp_count1=int(line_index[22:26].strip())
        if(not(temp_count==temp_count1 or temp_count+1==temp_count1)):
            flag=0
            break
        temp_count=temp_count1
    count1 = int(lines[len(lines)-1][22:26].strip())-count
    f1.write(pdbid+".pdb"+"\t"+str(flag)+"\t"+str(count)+"\t"+str(count1)+"\n")

    renum_command1="perl "+BindirPath+"/occupancyrect2a.pl " +renumPdbPath+"/"+pdbid
    rename_command1="mv "+renumPdbPath+"/"+ pdbid+".rect"+" "+renumPdbPath+"/"+pdbid+".pdb"
    renum_command2="perl "+BindirPath+"/resnumchange.pl "+renumPdbPath+"/"+pdbid+ " 1 "+chain_name+" " +renumPdbPath+"/"+ pdbid+".rect"
    rename_command2="mv "+renumPdbPath+"/"+pdbid+".rect"+" "+renumPdbPath+"/"+pdbid+".pdb"
    os.system(renum_command1)
    os.system(rename_command1)
    os.system(renum_command2)
    os.system(rename_command2)


