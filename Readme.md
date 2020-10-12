PROFOUND: predicting change in PROtein FOldability associated with mUlti poiNt Deletions in a newly defined protein structure database

For any doubts and suggestions contact-
Computations in Structural, Molecular, and Systems Biology Lab
Department of Computer Science and Engineering
Indian Institute of Technology, Kharagpur
West Bengal, India
Email: anupambanerjee@iitkgp.ac.in / amit20031996@iitkgp.ac.in /  pralay@cse.iitkgp.ac.in



This folder (unzipped PROFOUND.tar.gz) contains the codes and the input files necessary to change in protein foldability associated with Multi Point Deletions in a 
newly defined protein structure database.

Software requirement:
	a.	Stride- You can get the Stride executable from- http://webclu.bio.wzw.tum.de/stride/		
	b.	HBPlus- You can get the HBPlus executable from- http://www.ebi.ac.uk/thornton-srv/software/HBPLUS/
	c.	FoldX- You can get the FoldX executable from- http://foldxsuite.crg.eu/academic-license-info

Download and install them in your system if they are not already available. Make sure that the programs are running standalone and can be integrated with the given code. Make sure that the FoldX folder has the rotabase.txt file present in it.

Provide the absolute path of the softwares and relevant folders in runinfo.txt file. An example file is provided as runinfo_example.txt. 
Please note that runinfo_example.txt is prepared in our system. Your system specific location should be in runinfo.txt.
DO NO ALTER THE SEQEUNCE AND THE FORMAT OF the runinfo.txt file.

bin/Model/loop.tar.gz and bin/Model/nonloop.tar.gz files need to be extracted.

The exec.py code and runinfo.txt file should be in the same location.

Place the PDB files whose foldability you want to determine in the INPUT folder using the following naming convention- name_StartingDeletionPosition_EndingDeletionPosition.pdb.
The name is the desired name of the pdb file without any space or special characters.The StartingDeletionPosition is the starting deletion position and EndingDeletionPosition is ending deletion of the desired deletion will take place.

Ensure the following in the input PDB files-
a. There should be no chain break (backbone discontinuity) in the PDB files. You can fix it using ITASSER or MODELLER softwares.
b. The residues shouldn't have multiple occupancies.
c. The atom and residue index should be contiguous.
d. There should be a chain ID in the relevant location in the PDB files.
e. A minimum of 2 and a maximum of 15 contiguous residues can be deleted.
f. Deletions of N-terminal and C-termal residue (akin to truncation of the protein) is not allowed.

Ensure that your python (version Python 3.5.2 used) has the Scikit-learn (version 0.20.2. used) and numpy (version 1.15.4 used) packages installed in it.

Finally run exec.py as:
Run- python absolute_path_to_exec.py (Run exec.py using its absolute path)

The foldability predition will be present in the file output.txt in the workdir folder mentioned in the runinfo.txt file. 
