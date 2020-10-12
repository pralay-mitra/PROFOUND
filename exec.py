import numpy
import math
import re
import os
from sklearn.externals import joblib
import sys 

argumentList = sys.argv 
exc_dir=sys.argv[0]
exc_dir_list = exc_dir.strip().split("/")
exc_dir_list = exc_dir_list[:-1]
exc_dir1 = "/".join(str(a) for a in exc_dir_list)

if(exc_dir=="exec.py"):
    exc_dir1="."

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

os.chdir(WorkdirPath)

loop_model = BindirPath+"/Model/loop.mod"
clf = joblib.load(loop_model)   

nonloop_model = BindirPath+"/Model/nonloop.mod"
clf1 = joblib.load(nonloop_model)   

os.system("python3 "+BindirPath+"/pre.py")

PDBpath=WorkdirPath+"/Updated_INPUT"

listPDB = []
pattern = "*.pdb" 
contents = os.listdir(PDBpath)
for content in contents:
    if(content.endswith(".pdb") or content.endswith(".PDB") ):
        listPDB.append(content)

listPDB.sort()
f_pre=open(WorkdirPath+"/pre_data.txt","r")
lines_pre=f_pre.readlines()
pre_list=[]

for i in lines_pre:
    i=i[:-1]
    l = i.strip().split('\t')
    l[1]=int(l[1])
    l[2]=int(l[2])
    l[3]=int(l[3])
    pre_list.append(l)


pre_count=-1

f_res=open(exc_dir1+"/output.txt","w")

res = "Index | Protein Name | Start deletion position | End deletion position | Result"
f_res.write(res+"\n")

res = "==========================================================================================================================================="
f_res.write(res+"\n")

xx=0
for pdb_file in listPDB:
    xx=xx+1
    pre_count=pre_count+1
    pdb_info_list=re.split('_',pdb_file)
    pdb_info_list1=pdb_info_list[2].split('.')
    start_deletion_position=int(pdb_info_list[1])-pre_list[pre_count][2]
    end_deletion_position=int(pdb_info_list1[0])-pre_list[pre_count][2]

    if(end_deletion_position<start_deletion_position):
        res1 =  "Starting deletion position must smaller than ending deletion position"
        res = str(xx)+"\t"+pdb_info_list[0]+"\t\t\t"+pdb_info_list[1]+"\t\t\t"+pdb_info_list1[0]+"\t\t"+res1
        f_res.write(res+"\n")
        continue


    if(end_deletion_position-start_deletion_position+1>15):
        res1 =  "Deletion length more than 15 is not allowed"
        res = str(xx)+"\t"+pdb_info_list[0]+"\t\t\t"+pdb_info_list[1]+"\t\t\t"+pdb_info_list1[0]+"\t\t"+res1
        f_res.write(res+"\n")
        continue


    if(start_deletion_position==1):
        res1 =  "You can not delete first amino acid of PDB"
        res = str(xx)+"\t"+pdb_info_list[0]+"\t\t\t"+pdb_info_list[1]+"\t\t\t"+pdb_info_list1[0]+"\t\t"+res1
        f_res.write(res+"\n")
        continue

    if(end_deletion_position==pre_list[pre_count][3]):
        res1 =  "You can not delete last amino acid of PDB"
        res = str(xx)+"\t"+pdb_info_list[0]+"\t\t\t"+pdb_info_list[1]+"\t\t\t"+pdb_info_list1[0]+"\t\t"+res1
        f_res.write(res+"\n")
        continue

    if(pre_list[pre_count][1]==0):
        res1 =  "Chain break in your PDB"
        res = str(xx)+"\t"+pdb_info_list[0]+"\t\t\t"+pdb_info_list[1]+"\t\t\t"+pdb_info_list1[0]+"\t\t"+res1
        f_res.write(res+"\n")
        continue

    if(start_deletion_position>=end_deletion_position):
        res1 =  "Starting deletion position should be less than ending deletion position."
        res = str(xx)+"\t"+pdb_info_list[0]+"\t\t\t"+pdb_info_list[1]+"\t\t\t"+pdb_info_list1[0]+"\t\t"+res1
        f_res.write(res+"\n")
        continue

    pdb_file_path=PDBpath+"/"+pdb_file
    f_pdb = open(pdb_file_path, "r")
    lines=f_pdb.readlines()
    final_list=[]
    stride_output_path= StridePath+" "+ pdb_file_path + " > "+ WorkdirPath+"/"+pdb_file +".ent"
    print (stride_output_path)
    os.system(stride_output_path)
    stride_output_name=WorkdirPath+"/"+pdb_file +".ent"
    f = open(stride_output_name, "r")
    lines=f.readlines()
    print (lines)
    f.close()
    os.system("rm "+stride_output_name)
    line_no=0
    final_list_ent=[]
    for line in lines:
        line_no+=1
        if(line[:3]=='REM'):
            rem_line=line_no
    line_no=0
    for line in lines:
        line_no+=1
        if(line_no>=rem_line+1):
            line.strip
            l1=line.split(" ")
            l1 = [x for x in l1 if x]
            final_list_ent.append(l1)

# #     #########################################################################

    x=(int)(final_list_ent[0][3])
    len_ent=len(final_list_ent)
    loop_position=[]
    count=0
    for i  in range(x,len_ent+x):
        if(final_list_ent[i-x][5]=='b' or final_list_ent[i-x][5]=='B' or final_list_ent[i-x][5]=='T' or final_list_ent[i-x][5]=='C'):
            count+=1
        else:
            if(count>1):
                loop_position_entry=[]
                loop_position_entry.append(i-count)
                loop_position_entry.append(i-1)
                loop_position.append(loop_position_entry)
                count=0
            if(count==1):
                count=0


    if(count>1):
        loop_position_entry=[]
        loop_position_entry.append(i-count+1)
        loop_position_entry.append(i-1+1)
        loop_position.append(loop_position_entry)
        count=0

    flag=0
    for loop_itration in range(len(loop_position)):
        if(loop_position[loop_itration][0]==start_deletion_position and end_deletion_position==loop_position[loop_itration][1]):
            flag=2
            break
        elif(loop_position[loop_itration][0]<=start_deletion_position and end_deletion_position<=loop_position[loop_itration][1]):
            flag=1
            break


    if(flag==2):
        res1 =  "You can not delete entire loop."
        res = str(xx)+"\t"+pdb_info_list[0]+"\t\t\t"+pdb_info_list[1]+"\t\t\t"+pdb_info_list1[0]+"\t\t"+res1
        f_res.write(res+"\n")
        continue
    elif(flag==1):
        f=open(WorkdirPath+"/"+"deletion_detail.txt","w")
        f.write(pdb_file+" "+str(start_deletion_position)+" "+str(end_deletion_position))
        f.close()
        feature_calculation_command="python3 "+BindirPath+"/feature_loop.py"
        os.system(feature_calculation_command)

        f=open(WorkdirPath+"/"+"new_feature.txt","r")
        lines=f.readlines()
        line=lines[0][:-1]
        new_feature_array=line.split("\t")
        new_feature_list=[]
        new_feature_array_list=[]
        for i in new_feature_array:
            new_feature_list.append(float(i))
        new_feature_array=numpy.array(new_feature_list)
        new_feature_array_list.append(new_feature_array)
        input_array=numpy.array(new_feature_array_list)
        f.close()
        result = clf.predict_proba(input_array)
        res1 =  "The probability of the protien folding to its native fold after deletion is "+str(round(result[0][1],4)*100)+"%"
        res = str(xx)+"\t"+pdb_info_list[0]+"\t\t\t"+pdb_info_list[1]+"\t\t\t"+pdb_info_list1[0]+"\t\t"+res1
        f_res.write(res+"\n")
    else:
        f=open(WorkdirPath+"/"+"deletion_detail.txt","w")
        f.write(pdb_file+" "+str(start_deletion_position)+" "+str(end_deletion_position))
        f.close()
        feature_calculation_command="python3 "+BindirPath+"/feature_nonloop.py"
        os.system(feature_calculation_command)

        f=open(WorkdirPath+"/"+"new_feature.txt","r")
        lines=f.readlines()
        line=lines[0][:-1]
        new_feature_array=line.split("\t")
        new_feature_list=[]
        new_feature_array_list=[]
        for i in new_feature_array:
            new_feature_list.append(float(i))
        new_feature_array=numpy.array(new_feature_list)
        new_feature_array_list.append(new_feature_array)
        input_array=numpy.array(new_feature_array_list)
        f.close()
        result = clf1.predict_proba(input_array)
        res1 =  "The probability of the protien folding to its native fold after deletion is "+str(round(result[0][1],4)*100)+"%"
        res = str(xx)+"\t"+pdb_info_list[0]+"\t\t\t"+pdb_info_list[1]+"\t\t\t"+pdb_info_list1[0]+"\t\t"+res1
        f_res.write(res+"\n")

os.system("rm "+WorkdirPath+"/new_feature.txt")
os.system("rm "+WorkdirPath+"/deletion_detail.txt")
os.system("rm "+WorkdirPath+"/pre_data.txt")
os.system("rm Unrecognized_molecules.txt")
os.system("rm -rf "+WorkdirPath+"/Updated_INPUT")
os.system("rm "+WorkdirPath+"/hbdebug.dat")
