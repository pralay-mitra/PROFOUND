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

SaltBridgePath = BindirPath
projectPath = WorkdirPath
foldxPath = FoldxPath


f.close()

f=open(WorkdirPath+"/deletion_detail.txt","r")
lines=f.readlines()
line=lines[0]
pair_loop_position_list=line.strip().split(' ')
listPDB = [pair_loop_position_list[0]]

feature_file=open(WorkdirPath+"/new_feature.txt","w")


for pdb_file_name in range(0,len(listPDB)):
        hydro_scale={'ILE':4.5,'VAL':4.2,'LEU':3.8,'PHE':2.8,'CYS':2.5,'MET':1.9,'ALA':1.8,'GLY':-0.9,'THR':-0.7,'SER':-0.8,'TRP':-0.9,'TYR':-1.3,'PRO':-1.6,'HIS':-3.2,'GLU':-3.5,'ASP':-3.5,'GLN':-3.5,'ASN':-3.5,'LYS':-3.9,'ARG':-4.5}
        area_scale={'ILE':181.0,'VAL':164.5,'LEU':193.1,'PHE':222.8,'CYS':146.1,'MET':203.4,'ALA':118.1,'GLY':88.1,'THR':152.5,'SER':129.8,'TRP':266.3,'TYR':236.8,'PRO':146.8,'HIS':202.5,'GLU':186.2,'ASP':158.7,'GLN':193.2,'ASN':165.5,'LYS':225.8,'ARG':256.0}    
        prop_scale={'ILE':0.47,'VAL':0.50,'LEU':0.59,'PHE':0.60,'CYS':1.19,'MET':0.60,'ALA':0.66,'GLY':1.56,'THR':0.96,'SER':1.43,'TRP':0.96,'TYR':1.14,'PRO':1.52,'HIS':0.95,'GLU':0.74,'ASP':1.46,'GLN':0.98,'ASN':1.56,'LYS':1.01,'ARG':0.95}
        pdb_path_file=WorkdirPath+"/Updated_INPUT/"+listPDB[pdb_file_name]
        f = open(pdb_path_file, "r")
        lines=f.readlines()
        final_list=[]

        for line_index in lines:
            l=[]
            l.append(line_index[0:4].strip())
            l.append(line_index[7:11].strip())
            l.append(line_index[13:17].strip())
            l.append(line_index[17:20].strip())
            l.append(line_index[21:22].strip())
            l.append(line_index[22:26].strip())
            l.append(line_index[26:38].strip())
            l.append(line_index[38:46].strip())
            l.append(line_index[46:54].strip())
            l.append(line_index[54:60].strip())
            l.append(line_index[60:66].strip())
            l.append(line_index[66:81].strip())
            final_list.append(l)

        ##########################################

        count =0
        xValue=0
        yValue=0
        zValue=0

        newList=[]
        previousResidueNum = final_list[0][5]
        flagg=0
        for entry in final_list:
            if(entry[5] == previousResidueNum):
                count +=1
                xValue+=float(entry[6])
                yValue+=float(entry[7])
                zValue+=float(entry[8])
                n=entry[3]
            else:
                myTuple = (n,previousResidueNum ,xValue/count,yValue/count,zValue/count)
                newList.append(myTuple)
                previousResidueNum=entry[5]
                xValue = float(entry[6])
                yValue = float(entry[7])
                zValue = float(entry[8])
                count =1
        myTuple = (n,previousResidueNum ,xValue/count,yValue/count,zValue/count)
        newList.append(myTuple)
        previousResidueNum=entry[5]
        xValue = float(entry[6])
        yValue = float(entry[7])
        zValue = float(entry[8])
        count =1

        ##############################################

        pdbid,ext = listPDB[pdb_file_name].strip().split('.')
        stride_output_path= StridePath+" "+ pdb_path_file + " > "+ WorkdirPath+"/"+pdbid +".ent"
        os.system(stride_output_path)
        stride_output_name=WorkdirPath+"/"+pdbid +".ent"
        f = open(stride_output_name, "r")
        lines=f.readlines()
        line_no=0
        final_list_ent=[]
        pair_loop_entry=[[int(pair_loop_position_list[1]),int(pair_loop_position_list[2])]]
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

        #########################################################################

        x=(int)(final_list_ent[0][3])
        # x=0
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
          

        ###########################################################################
        def phi_psi(start_loop,end_loop):
            x=(int)(final_list_ent[0][3])
            len_ent=len(final_list_ent)
            phi_list=[]
            psi_list=[]
            for i  in range(start_loop,end_loop+1):
                phi_list.append(float(final_list_ent[i-x][7]))
                psi_list.append(float(final_list_ent[i-x][8]))
            max_phi=max(phi_list)
            min_phi=min(phi_list)
            max_psi=max(psi_list)
            min_psi=min(psi_list)
            l=[max_phi,min_phi,max_psi,min_psi]
            return l

    ####################################################################################

        def calculateFoldx(start_loop,end_loop):
            destination = projectPath + "/"

            foldxPath_list = foldxPath.strip().split("/")
            foldxPath_list = foldxPath_list[:-1]
            foldxPath1 = "/".join(str(a) for a in foldxPath_list)

            foldxPath1_rota=foldxPath1+"/rotabase.txt"
            copyCommand = 'cp ' + foldxPath1_rota + ' ' + projectPath
            os.system(copyCommand)



            getpdb = open(pdb_path_file,'r')
            pdblines = getpdb.readlines()
            pdbout1name = destination + "temp1.pdb"
            pdbout2name = destination + "temp2.pdb"
            pdbout1 = open(pdbout1name,'w')
            pdbout2 = open(pdbout2name,'w')

            for line in pdblines:
                l=[]
                l.append(line[0:4].strip())
                l.append(line[7:11].strip())
                l.append(line[13:17].strip())
                l.append(line[17:20].strip())
                l.append(line[21:22].strip())
                l.append(line[22:26].strip())
                l.append(line[26:38].strip())
                l.append(line[38:46].strip())
                l.append(line[46:54].strip())
                l.append(line[54:60].strip())
                l.append(line[60:66].strip())
                l.append(line[66:81].strip())
                if(int(l[5]) == start_loop and int(l[5]) <= end_loop):
                    continue
                if(int(l[5])< start_loop):
                    pdbout1.write(line )
                if(int(l[5])> end_loop):
                    pdbout2.write(line )
            pdbout1.close()
            pdbout2.close() 

            
            copyCommand = 'cp ' + pdb_path_file + ' ' + destination



            os.system(copyCommand)

            
            statinfo = os.stat(pdbout1name)
            temp1_size = statinfo.st_size
            statinfo = os.stat(pdbout2name)
            temp2_size = statinfo.st_size


            foldxCommand1 = foldxPath +" --command=Stability --pdb="+pdbid+".pdb"+" >"+pdbid+".foldx"
            foldxCommand2 = foldxPath +" --command=Stability --pdb="+"temp1" +".pdb"+" > temp1.foldx"
            foldxCommand3 = foldxPath +" --command=Stability --pdb="+"temp2.pdb"+" > temp2.foldx"

            zero_ind1 = 0
            zero_ind2 = 0
            os.system(foldxCommand1)
            if(not temp1_size == 0):
                os.system(foldxCommand2)
            else:
                list12 = [0] * 16
                zero_ind1 = 1 

            if(not temp2_size == 0):
                os.system(foldxCommand3)
            else:
                list13 = [0] * 16
                zero_ind2 = 1

            foldxCommand1_output=pdbid+".foldx"
            foldxCommand2_output="temp1.foldx"
            foldxCommand3_output="temp2.foldx"                
            list1 =[]
            f = open(foldxCommand1_output,'r')
            while(True):
                line = f.readline()
                if(not line):
                    break

                line = line.strip()

                try:
                    if(line[:6]=="number"):
                        continue
                    name,value= line.split('=',1)
                    name = name.strip()
                    value = value.strip()
                except ValueError:
                    continue

                foldx_tuple = name,value
                list1.append(foldx_tuple)

            list11=[]
            for i in range(len(list1)):
                if(i==11 or i==13 or i==16 or i==17 or i==19):
                    continue
                else:
                    list11.append(float(list1[i][1]))



            if(zero_ind1 == 0):
                list2 =[]
                f1 = open(foldxCommand2_output,'r')
                while(True):
                    line = f1.readline()
                    if(not line):
                        break

                    line = line.strip()

                    try:
                        if(line[:6]=="number"):
                            continue
                        name,value= line.split('=',1)
                        name = name.strip()
                        value = value.strip()
                    except ValueError:
                        continue

                    foldx_tuple = name,value
                    list2.append(foldx_tuple)


                list12=[]
                for i in range(len(list2)):
                    if(i==11 or i==13 or i==16 or i==17 or i==19):
                        continue
                    else:
                        list12.append(float(list2[i][1]))

            if(zero_ind2 == 0):
                list3 =[]
                f2 = open(foldxCommand3_output,'r')
                while(True):
                    line = f2.readline()
                    if(not line):
                        break

                    line = line.strip()

                    try:
                        if(line[:6]=="number"):
                            continue
                        name,value= line.split('=',1)
                        name = name.strip()
                        value = value.strip()
                    except ValueError:
                        continue

                    foldx_tuple = name,value
                    list3.append(foldx_tuple)


                list13=[]
                for i in range(len(list3)):
                    if(i==11 or i==13 or i==16 or i==17 or i==19):
                        continue
                    else:
                        list13.append(float(list3[i][1]))


            resultArray1 = numpy.array(list11,dtype=float)
            resultArray2 = numpy.array(list12,dtype=float)
            resultArray3 = numpy.array(list13,dtype = float)

            resultArray4 = resultArray1 - resultArray2 -resultArray3
            list_temp=[]
            # for x in range(len(resultArray4)):
            list_temp.append(resultArray4[len(resultArray4)-1])
            os.system("rm "+ foldxCommand1_output)

            os.system("rm "+ foldxCommand2_output)
            os.system("rm "+ foldxCommand3_output)

            os.system("rm "+ pdbout1name)
            os.system("rm "+ pdbout2name)
            pdbfilepath_foldx = destination + pdbid + ".pdb"
            os.system("rm "+ pdbfilepath_foldx)
            os.system("rm rotabase.txt")
            return list_temp

    ####################################################################################

        def long_range(start_loop,end_loop,last):
            long_range_list=[]
            salt_bridge_list=list()
            x=(int)(final_list_ent[0][3])
        
            salt_bridge_output_path= "perl "+ SaltBridgePath+"/salt_bridges_new.pl"+ " "+ pdb_path_file
            os.system(salt_bridge_output_path)
            disulphide_output_path= "perl "+ SaltBridgePath+"/disulphide_new.pl"+ " "+ pdb_path_file
            os.system(disulphide_output_path)
            hb_output_path=HBplusPath +" "+ pdb_path_file
            os.system(hb_output_path)
            # mvcommand = "mv "+ SaltBridgePath +"/" + pdbid+".hb2 " + SaltBridgePath+"/HB_plus/"+pdbid+".hb2"
            # os.system(mvcommand)
            path1="salt_bridge.txt"
            path2="disulphide.txt"
            path3=pdbid+".hb2"
            f1 = open(path1, "r")
            lines=f1.readlines()
            sb=0
            sb1=0
            for line in range(len(lines)-1):
                lines[line].strip()
                l=re.split(' |,\t+',lines[line])
                l = [x1 for x1 in l if x1]
                v1=int(l[3])
                v2=int(l[5][:-9])
                temp_list=[v1,v2]
                salt_bridge_list.append(temp_list)

            flag=[0] * len(salt_bridge_list)
            for i in range(start_loop,end_loop+1):
                for j in range(len(salt_bridge_list)):
                    if(i==salt_bridge_list[j][0] and flag[j]==0):
                            sb+=1
                            flag[j]=1
                    elif(i==salt_bridge_list[j][1] and flag[j]==0):
                            sb+=1
                            flag[j]=1


            for i in salt_bridge_list:
                if (i[0]>x and i[0]<start_loop-1 and i[1]>end_loop+1 and i[1]<last):
                    sb1+=1
            
            ######################################
            Disulphide_list=list()
            f2 = open(path2, "r")
            lines1=f2.readlines()
            db=0
            db1=0


            for line in range(len(lines1)-1):
                lines1[line].strip()
                l1=re.split(' |,\t+',lines1[line])
                l1 = [x1 for x1 in l1 if x1]
                v11=int(l1[3])
                v22=int(l1[5][:-9])
                temp_list1=[v11,v22]
                Disulphide_list.append(temp_list1)

            flag=[0] * len(Disulphide_list)
            for i in range(start_loop,end_loop+1):
                for j in range(len(Disulphide_list)):
                    if(i==Disulphide_list[j][0] and flag[j]==0):
                            db+=1
                            flag[j]=1
                    elif(i==Disulphide_list[j][1]  and flag[j]==0):
                            db+=1
                            flag[j]=1



            for i in Disulphide_list:
                if (i[0]>x and i[0]<start_loop-1 and i[1]>end_loop+1 and i[1]<last):
                    db1+=1

            # ######################################
            hb_list=list()
            f3 = open(path3, "r")
            lines212=f3.readlines()
            hb=0
            hb1=0
            for line in range(len(lines212)):
                if(lines212[line][0:6]=='n    s'):
                    break

            for line11 in range(line+1,len(lines212)):
                lines212[line11] = lines212[line11].strip()
                l2=re.split(' |,\t+',lines212[line11])
                l2 = [x1 for x1 in l2 if x1]
                v12=int(l2[0][1:5])
                v23=int(l2[2][1:5])
                temp_list=[v12,v23]
                temp_list.sort()
                hb_list.append(temp_list)

            flag=[0] * len(hb_list)

            for i in range(start_loop,end_loop+1):
                for j in range(len(hb_list)):
                    if(i==hb_list[j][0] and flag[j]==0):
                            hb+=1
                            flag[j]=1
                    if(i==hb_list[j][1] and flag[j]==0):
                            hb+=1
                            flag[j]=1

            for i in hb_list:
                if(i[0]>x and i[0]<start_loop-1 and i[1]>end_loop+1 and i[1]<last):
                    hb1+=1

            long_range_list.append(sb)
            long_range_list.append(db)
            long_range_list.append(hb)
            long_range_list.append(sb1)
            long_range_list.append(db1)
            # long_range_list.append(hb1)
            return long_range_list

    #####################################################################################


        a=c=d=e=f=g=h=i1=k=l=m=n=p=q=r=s=t=v=w=y=0
        amino_dict_pdb=dict()
        for i  in range(len(final_list_ent)):
            if(final_list_ent[i][1]=='ALA'):
                a+=1
            elif(final_list_ent[i][1]=='CYS'):
                c+=1
            elif(final_list_ent[i][1]=='ASP'):
                d+=1
            elif(final_list_ent[i][1]=='GLU'):
                e+=1
            elif(final_list_ent[i][1]=='PHE'):
                f+=1
            elif(final_list_ent[i][1]=='GLY'):
                g+=1
            elif(final_list_ent[i][1]=='HIS'):
                h+=1
            elif(final_list_ent[i][1]=='ILE'):
                i1+=1
            elif(final_list_ent[i][1]=='LYS'):
                k+=1
            elif(final_list_ent[i][1]=='LEU'):
                l+=1
            elif(final_list_ent[i][1]=='MET'):
                m+=1
            elif(final_list_ent[i][1]=='ASN'):
                n+=1
            elif(final_list_ent[i][1]=='PRO'):
                p+=1
            elif(final_list_ent[i][1]=='GLN'):
                q+=1
            elif(final_list_ent[i][1]=='ARG'):
                r+=1
            elif(final_list_ent[i][1]=='SER'):
                s+=1
            elif(final_list_ent[i][1]=='THR'):
                t+=1
            elif(final_list_ent[i][1]=='VAL'):
                v+=1
            elif(final_list_ent[i][1]=='TRP'):
                w+=1
            elif(final_list_ent[i][1]=='TYR'):
                y+=1
        amino_dict_pdb[0]=a
        amino_dict_pdb[1]=c
        amino_dict_pdb[2]=d
        amino_dict_pdb[3]=e
        amino_dict_pdb[4]=f
        amino_dict_pdb[5]=g
        amino_dict_pdb[6]=h
        amino_dict_pdb[7]=i1
        amino_dict_pdb[8]=k
        amino_dict_pdb[9]=l
        amino_dict_pdb[10]=m
        amino_dict_pdb[11]=n
        amino_dict_pdb[12]=p
        amino_dict_pdb[13]=q
        amino_dict_pdb[14]=r
        amino_dict_pdb[15]=s
        amino_dict_pdb[16]=t
        amino_dict_pdb[17]=v
        amino_dict_pdb[18]=w
        amino_dict_pdb[19]=y


        def amino_acid_propensity(start_loop,end_loop):
            x=(int)(final_list_ent[0][3])
            len_ent=len(final_list_ent)
            a=c=d=e=f=g=h=i1=k=l=m=n=p=q=r=s=t=v=w=y=0
            amino_dict=dict()
            for i  in range(start_loop,end_loop+1):
                if(final_list_ent[i-x][1]=='ALA'):
                    a+=1
                elif(final_list_ent[i-x][1]=='CYS'):
                    c+=1
                elif(final_list_ent[i-x][1]=='ASP'):
                    d+=1
                elif(final_list_ent[i-x][1]=='GLU'):
                    e+=1
                elif(final_list_ent[i-x][1]=='PHE'):
                    f+=1
                elif(final_list_ent[i-x][1]=='GLY'):
                    g+=1
                elif(final_list_ent[i-x][1]=='HIS'):
                    h+=1
                elif(final_list_ent[i-x][1]=='ILE'):
                    i1+=1
                elif(final_list_ent[i-x][1]=='LYS'):
                    k+=1
                elif(final_list_ent[i-x][1]=='LEU'):
                    l+=1
                elif(final_list_ent[i-x][1]=='MET'):
                    m+=1
                elif(final_list_ent[i-x][1]=='ASN'):
                    n+=1
                elif(final_list_ent[i-x][1]=='PRO'):
                    p+=1
                elif(final_list_ent[i-x][1]=='GLN'):
                    q+=1
                elif(final_list_ent[i-x][1]=='ARG'):
                    r+=1
                elif(final_list_ent[i-x][1]=='SER'):
                    s+=1
                elif(final_list_ent[i-x][1]=='THR'):
                    t+=1
                elif(final_list_ent[i-x][1]=='VAL'):
                    v+=1
                elif(final_list_ent[i-x][1]=='TRP'):
                    w+=1
                elif(final_list_ent[i-x][1]=='TYR'):
                    y+=1
            amino_dict[0]=a
            amino_dict[1]=c
            amino_dict[2]=d
            amino_dict[3]=e
            amino_dict[4]=f
            amino_dict[5]=g
            amino_dict[6]=h
            amino_dict[7]=i1
            amino_dict[8]=k
            amino_dict[9]=l
            amino_dict[10]=m
            amino_dict[11]=n
            amino_dict[12]=p
            amino_dict[13]=q
            amino_dict[14]=r
            amino_dict[15]=s
            amino_dict[16]=t
            amino_dict[17]=v
            amino_dict[18]=w
            amino_dict[19]=y
            amino_propensity_list=list()
            for i in range(20):
                if (amino_dict_pdb[i]==0):
                    x=0.0
                else:
                    x=float((float(amino_dict[i])/float(end_loop-start_loop+1)))/(float((amino_dict_pdb[i]))/float(len_ent))
                amino_propensity_list.append(x)
            return amino_propensity_list

 

    #####################################################################################

        def surfacearea(start_loop,end_loop):
            temp_list12=[]
            len_loop=end_loop-start_loop+1
            for i in range(start_loop,end_loop+1):
                x123=float (final_list_ent[i-x][9])/float (float(area_scale[final_list_ent[i-x][1]]))
                temp_list12.append(x123)
            net_sum=0
            for i in temp_list12:
                net_sum+=i
            avg=net_sum/len_loop

            net_deviation_sum=0
            for i in temp_list12:
                net_deviation_sum+=(i-avg)*(i-avg)
            net_deviation=math.sqrt(net_deviation_sum/len_loop)

            temp_list22=[avg,net_deviation]
            return temp_list22




    #####################################################################################
        result=[]
        def entry(start_loop,end_loop,loop_number):
            row=[]
            long_range_output_list=long_range(start_loop,end_loop,len_ent)
            phi_psi_output_list=phi_psi(start_loop,end_loop)
            amino_acid_propensity_list=amino_acid_propensity(start_loop,end_loop)
            surfacearea_list=surfacearea(start_loop,end_loop)
            foldx_list=calculateFoldx(start_loop,end_loop)
            score_list=[]
            len_loop=end_loop-start_loop+1
            # row.append(len_loop)
            last=len_ent


            
            a1 = numpy.array(((float)(newList[start_loop-x][2]),(float)(newList[start_loop-x][3]),(float)(newList[start_loop-x][4])))
            b1 = numpy.array(((float)(newList[end_loop-x][2]),(float)(newList[end_loop-x][3]),(float)(newList[end_loop-x][4])))
            end_end_distance = numpy.linalg.norm(a1-b1)

            walk_distance= 0
            for i in range(start_loop,end_loop):
                a = numpy.array(((float)(newList[i-x][2]),(float)(newList[i-x][3]),(float)(newList[i-x][4])))
                b = numpy.array(((float)(newList[i-x+1][2]),(float)(newList[i-x+1][3]),(float)(newList[i-x+1][4])))
                dist1 = numpy.linalg.norm(a-b)
                walk_distance+=dist1

            for i in range(start_loop,end_loop+1):
                s=0
                s1=0
                for j in range(start_loop,end_loop+1):
                    if(i!=j):
                        a = numpy.array(((float)(newList[i-x][2]),(float)(newList[i-x][3]),(float)(newList[i-x][4])))
                        b = numpy.array(((float)(newList[j-x][2]),(float)(newList[j-x][3]),(float)(newList[j-x][4])))
                        dist = numpy.linalg.norm(a-b)
                        s+=(1/((dist)*(dist)))
                        s1+=(hydro_scale[newList[j-x][0]]/((dist)*(dist)))
                score=[]
                score.append(i)
                score.append(s)
                score.append(s1)
                score_list.append(score)
            net_sum=0
            net_hydro_sum=0
            for i in score_list:
                net_sum+=i[1]
                net_hydro_sum+=i[2]
            avg=net_sum/len_loop
            avg1=net_hydro_sum/len_loop

            net_deviation_sum=0
            net_deviation_hydro_sum=0
            for i in score_list:
                net_deviation_sum+=(i[1]-avg)*(i[1]-avg)
                net_deviation_hydro_sum+=(i[2]-avg1)*(i[2]-avg1)
            net_deviation=math.sqrt(net_deviation_sum/len_loop)
            net_deviation_hydro=math.sqrt(net_deviation_hydro_sum/len_loop)

            row.append(round(end_end_distance,1))
            # row.append(round(walk_distance,1))
            row.append(round(avg,1))
            row.append(round(net_deviation,1))
            row.append(round(avg1,1))
            row.append(round(net_deviation_hydro,1))
            for var in long_range_output_list:
                row.append(var)

            for var in phi_psi_output_list:
                row.append(round(float(var),1))
            
            row.append(start_loop-x)
            row.append(last-end_loop)

            row.append(round(surfacearea_list[0],1))
            row.append(round(surfacearea_list[1],1))
            for var in amino_acid_propensity_list:
                row.append(round(var,1))
            for var in foldx_list:
                row.append(round(var,1))
            result.append(row)




        for i in range(len(pair_loop_entry)):
            entry(pair_loop_entry[i][0],pair_loop_entry[i][1],i+1)
    
        string1 = "\t".join(str(a) for a in result[0])
        feature_file.write(string1+"\n")
        os.system("rm "+WorkdirPath+"/"+pdbid +".ent")
        os.system("rm "+pdbid+".hb2")
        os.system("rm "+pdbid+"_0_ST.fxout")
        os.system("rm temp1_0_ST.fxout")
        os.system("rm temp2_0_ST.fxout")
        os.system("rm disulphide.txt")
        os.system("rm salt_bridge.txt")
        os.system("rm temp_salt.txt")
        os.system("rm temp_dsp.txt")


feature_file.close()
