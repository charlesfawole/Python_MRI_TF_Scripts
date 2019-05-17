#file 5 in the lead voltage heat script series
#the final scrip

import cmath
import csv
import numpy as np
import os
import scipy.io as sio
from scipy.interpolate import interp1d
import collections as col


#the overview below applies to the matlab script and not to this python script
'''
%% overview
% This script is developed for calculating the RF induced heating/voltage
% by transfer function method.
% Following input files should be ready before running this script;
%   1) Tangential E field along the lead
%   2) Transfer function
%   3) simulation_list
% The simulation_list should following this format:
%   1) Amplitude of cos excitation
%   2) Amplitude of sin excitation
%   3) Folder of Etan files
%   4) directory of TF file
%   5) Folder of normalization files
%   6) Result type
%   7) Normalization type
%   8) directory of result files
% After running this script, one or more excel files will be generated
% which contains the final result.
'''
# build up a map contains coils diameter information

coil_key = ['1.5T bodycoil Male','1.5T bodycoil Female',
    '1.5T bodycoil Boy','1.5T bodycoil Girl','1.5T bodycoil Fat',
    '1.5T headcoil Male','1.5T headcoil Female',
    '1.5T headcoil Boy','1.5T headcoil Girl','1.5T headcoil Fat',
    '1.5T 82Dbodycoil Male','1.5T 82Dbodycoil Fat',
    '3T bodycoil Male','3T bodycoil Female',
    '3T bodycoil Boy','3T bodycoil Girl','3T bodycoil Fat',
    '3T headcoil Male','3T headcoil Female',
    '3T headcoil Boy','3T headcoil Girl','3T headcoil Fat',
    '3T 82Dbodycoil Male','3T 82Dbodycoil Fat'];
	
coil_diameter_list = [630,630,630,630,758,240,240,240,240,240,820,820,
    630,630,630,630,758,240,240,240,240,240,820,820];
	
coil_map =col.OrderedDict(zip(coil_key,coil_diameter_list))

#read the simulation list line and get info about each simulation
#simulation_list = 'E:\human model\script\simulation_list.txt';  # deprecating simulation_list in this python implementation
#f_sim_list = fopen(simulation_list);
#sim_set = fgetl(f_sim_list);

sim_index = 1; # index of simulation
result_number = 1; # index of line in result excel sheet
line_index = 1; # index of line in result summary excel sheet
pre_result_file = 'none'; # flag of result file


#cos_mag = 1                           # Hardcode this for now, user should be promted in updated code or should be obtained automatically
#sin_mag = complex(0,-1.0)		   # Hardcode this for now, user should be promted in updated code 
Etan_folder = 'E:\\human model\\Etan\\1.5T bodycoil Fat'
tf_file = 'E:\\human model\\TF\\1.5T_heating_ring_106_304.csv'   # the matalab version uses xls but this version uses CSV. User needs to convert to transfer function file to CSV before use
norm_folder = 'E:\\human model\\Matlab human model'
result_type = 'heating'
norm_type = 'B1_wbSAR'
result_file = 'E:\\human model\\result\\example3.xlsx'


cos_mags = [1,complex(0,-1)]
sin_mags = [1,complex(0,1)]
sources = np.vstack((cos_mags, sin_mags))
Etan_folders = [Etan_folder]*2
tf_files = [tf_file]*2
norm_folders = [norm_folder]*2
Etans_folders = [Etan_folder]*2 
tf_files = [tf_file] *2
norm_folders = [norm_folder]*2
result_types = [result_type]*2
norm_types = [norm_type]*2
result_files = [result_file]*2

for ii in range (len(cos_mags)): # use length of cos_mags list to determine how many entries need to be manipulated

	source = [sources[ii][0] , sources[ii][1]]#complex(cos_mags[i],sin_mags[i])
	polarization = None

	if cmath.phase(sources[ii][0] + sources[ii][1]) < 0:
		polarization = 'CCW'
	else:
		polarization = 'CW'
		
	Etan_folder = Etan_folders[ii]
	tf_file = tf_files[ii]
	norm_folder = norm_folders[ii]
	result_type = result_types[ii]
	norm_type = norm_types[ii]
	result_file = result_files[ii]
	
	if pre_result_file != result_file:
		result_number = 1;
		line_index = 1;
	
	with open(tf_file,'rb') as f:
		reader = csv.reader(f)
		tf_data = list(reader)
		
	tf_data = np.transpose(tf_data)[:],[0]
	tf_mag = tf_data[0][0]
	
	tf_mag = [float(ti) for ti in tf_mag]
	
	tf_phase =(tf_data[0][1] ) 
	tf_phase = [float(tp) for tp in tf_phase]
	
	tf_phase = np.multiply(tf_phase,  (1.0/180*np.pi))
	tf = np.multiply(tf_mag,np.exp(np.multiply(complex(0,1),tf_phase)))
	
	sim_type = Etan_folder.split('\\')[-1]
	
	B1_field_folder = norm_folder+'\\B1 file\\'
	
	B1_cos_file = B1_field_folder+sim_type+' cos.txt'
	
	B1_sin_file = B1_field_folder+sim_type+' sin.txt'
	
	
	#importfile(B1_cos_file,polarization)
	
	for B1_file in [B1_cos_file,B1_sin_file]: # this IDE is not allowing me to do functions well, so i do this

		with open(B1_file,'rb') as f:  # read the B1 cos txt file as a csv
			reader = csv.reader(f)
			B1_table = list(reader)	
			
		for  i in range(len(B1_table)): 
			B1_table[i] = B1_table[:][i][0].split()
		
		n_columns = len(B1_table[0])  # number of columns in B1 table
		B1_columns = [[]]*n_columns
		
		for i in range(n_columns): # extract each column
			 B1_columns[i] = [column[i] for column in B1_table]
		
		B1_rowNames = B1_columns[:][0]
		
		if polarization == 'CCW':
			B1_real = B1_columns[:][1]
			B1_imag = B1_columns[:][2]
			
		if polarization == 'CW':
			B1_real = B1_columns[:][3]
			B1_imag = B1_columns[:][4]
		
		B1_real = [float(i) for i in B1_real]
		
		B1_imag = [(float(i)) for i in B1_imag]
		B1_imag = np.multiply(B1_imag, complex(0,1))
		B1_cosOrSin = np.add(B1_real, B1_imag)

		if B1_file == B1_cos_file:
			B1_cos = B1_cosOrSin 
			
		if B1_file == B1_sin_file:
			B1_sin = B1_cosOrSin 
	B_len = len(B1_cos)
	B1 = [0]*B_len
	
	if polarization == 'CCW':
		source1_mod = source[1]
	if polarization == 'CW':
		source1_mod = source[1]*-1   # transpose 
	for i in range(B_len):
		B1_stack = (np.vstack((B1_cos, B1_sin)).T[i])   # doing the follow because it seems python does not allow for complex arithmetic
		a = B1_stack[0]*source[0]
		b = B1_stack[1]*source1_mod  #conjugate
		
		B1_real = a.real + b.real
		B1_imag = a.imag + b.imag
		
		B1_complex =np.array( complex(B1_real,B1_imag))
		B1[i] = np.abs(B1_complex)
	B1_table_unsorted = dict(zip(B1_rowNames,B1))
	B1_table = col.OrderedDict(sorted(B1_table_unsorted.items(), key = lambda t: t[0]))
	B1 = list(B1_table.values())  # THE ACTUAL B1, SORTED!
	# load M matrix and calculate norm factor for body SAR
	
	M_file_folder = norm_folder+'\\'+sim_type+'\\M'
	M_file_list = [f for f in os.listdir(M_file_folder) if os.path.isfile(os.path.join(M_file_folder, f))]
	

	M_N = len(M_file_list)
	M_name = M_file_list   # COMMENTED NOW UNCOMMETED, FORGOT
	
	norm = [0]*M_N
	
	if 'body' not in sim_type:
		standard = 3.2   # 3.2W/kg headSAR
		coil_type = 'Head'
	else:
		standard = 2
		coil_type = 'Body'   # 2W/kg wbSAR
		
	if '3T' not in sim_type:
			peak_standard = 30e-6
			RMS_standard = 6e-6
	else:
		peak_standard = 15e-6
		RMS_standard = 2.8e-6
		
	for M_i in range(M_N):
		M_file = M_file_list[M_i]
		M = sio.loadmat(os.path.join(M_file_folder,M_file)) ['M']
		
		# again numpy didnt seem to support complex matrix operation
		M_temp = [[]]*2
		for i in range(len(M)):
			a = np.transpose(M)[i][0]*source[0]
			b = np.transpose(M)[i][1]*source[1] 			
			
		
			M_temp_real = a.real + b.real
			M_temp_imag = a.imag + b.imag
		
			M_temp[i] =np.array( [complex(M_temp_real,M_temp_imag)])

		a = M_temp[0]*source[0]
		b = M_temp[1]*source[1]*-1
		
		norm_real = a.real + b.real
		norm_imag = a.imag + b.imag
		
		norm[M_i] = standard/complex(norm_real,norm_imag)
		M_name[M_i] = M_file_list[M_i][0:-4]
		
	norm_table_unsorted = dict(zip(M_name,norm))
	norm_table = col.OrderedDict(sorted(norm_table_unsorted.items(), key = lambda t: t[0]))
	norm_SAR = list(norm_table.values())   # Already sorted, no need to sort later
	sort_key = list(norm_table.keys())

	
	if norm_type == 'B1_wbSAR':
		B1_wbSAR = np.multiply(B1,np.sqrt(norm_SAR))
		#sort_key(norm_SAR==min(norm_SAR));  # DONT KNOW WHAT THIS DOES FROM MATLAB. IT DOES NOT SEEM NEEDED
		B1_ref = min(B1_wbSAR);
		norm_B1 = np.multiply(norm_SAR,np.square(np.divide(B1_ref,B1_wbSAR)))
		M_map = col.OrderedDict(zip(sort_key,norm_B1))
		B1_map = np.multiply(B1_ref,np.ones(M_N))
		B1_map = col.OrderedDict(zip(sort_key,B1_map))

	if norm_type == 'SAR':
		M_map = col.OrderedDict(zip(sort_key,norm_SAR))
		B1_wbSAR = np.multiply(B1,np.sqrt(norm_SAR))
		B1_map = col.OrderedDict(zip(sort_key,B1_wbSAR))
	if norm_type == 'B1_peak':
		M_map = np.square(np.divide(peak_standard,B1))
		M_map = col.OrderedDict(zip(sort_key,M_map))
		B1_map = np.multiply(peak_standard,np.ones(M_N))
		B1_map = col.OrderedDict(zip(sort_key,B1_map))
		
	if norm_type == 'B1_RMS':
		norm_temp = np.minimum(np.square(np.divide(np.multiply(RMS_standard,np.sqrt([2])),B1)),norm_SAR)
		M_map = col.OrderedDict(zip(sort_key,norm_temp))
		B1_map = np.divide(np.multiply(B1,np.sqrt(norm_temp)),np.sqrt([2]))
		B1_map = col.OrderedDict(zip(sort_key,B1_map))
	

		
		
		
		
	#go through the E field folder and calculate heating/voltage
	
	file_folder = Etan_folder+' cos'
	filelist = [f for f in os.listdir(file_folder) if os.path.isfile(os.path.join(file_folder, f))]
	N = len(filelist)	
	
	temp_list = [0]*N   # result ValueError
	full_result = ['']*N
	for q in range(N):
		e_file = filelist[q]
		index2 = e_file.find('Ele');  # using 'Ele' as signature to find sim type. A deviation from matlab script
		E_sim_type = e_file[0:index2-1]
		Edata_cos = np.transpose(np.loadtxt(Etan_folder+' cos\\'+e_file, skiprows=1))
		Efield_cos = np.add(Edata_cos[:][3] , np.multiply(Edata_cos[:][4], complex(0,1)))
		
		
		Edata_sin = np.transpose(np.loadtxt(Etan_folder+' sin\\'+e_file, skiprows=1))
		Efield_sin = np.add(Edata_sin[:][3] , np.multiply(Edata_sin[:][4], complex(0,1)))
		
		# Etan after normalization
		Efield = [0]*len(Efield_cos)
		for j in range(len(Efield_cos)):
			Efield_stack = (np.vstack((Efield_cos, Efield_sin)).T[j])   # doing the follow because again it seems python does not allow for complex arithmetic
			
			a = Efield_stack[0]*source[0]
			b = Efield_stack[1]*source[1]  #conjugate
			
			Efield_real = a.real + b.real
			Efield_imag = a.imag + b.imag
			
			
			cc = (M_map[E_sim_type].real)**(1.0/2)
			Efield_complex =  complex(Efield_real*cc , Efield_imag*cc)
		
			Efield[j] =Efield_complex
			#print str(Efield_complex) + "  "+str(cc)+"    "+E_sim_type
		
		
		#interplate Etan to same resolution of tf
		
		a = np.transpose(np.diff(((Edata_cos[0:3][:]))))
		b = np.square(a)
		c = np.sum(b,axis = 1)
		d = [0] + (np.sqrt(c)).tolist()

		distance = [0]*len(d)
		distance[0] = 0
		
		for iii in range (len(d)-1):
			distance[iii+1] = distance[iii] + d [iii +1]
			
		leadlength = 0.01*len(tf)
		
		dtf = np.arange(0.01,leadlength+0.01,0.01)
		
		e_tf = interp1d(distance,Efield,kind = 'nearest',fill_value = 'extrapolate')(dtf)
		
		
		#calcululate final ValueError
		
		
		if result_type == 'heating':
			temp = (abs(np.matmul(e_tf,tf)))**2
			#temp = (np.abs(np.multiply(e_tf,tf)))**2
		if result_type == 'voltage_RMS':
			temp = (abs(np.matmul(e_tf,tf)))
			
		if result_type == 'voltage_peak':
			temp = (abs(np.matmul(e_tf,tf)))*sqrt(2)
			
		if result_type == 'worst heating':
			temp = (np.matmul(abs(e_tf),abs(tf))/100)**2
			
		result_tab = 	e_file.split('_')
		tf_index = tf_file.rfind('\\')
		tf_name = tf_file[tf_index+1:-4]
		tf_info = tf_name.split('_')
		B1_res = B1_map[E_sim_type]
		full_result[q] =[sim_type,result_tab,result_type,norm_type,
                temp,tf_info,coil_map[sim_type],polarization,
                B1_res]
		
		file_o = open('E:\\c_test\\py_MRI_result_file2.csv', "ab")
		csv_o = csv.writer(file_o, delimiter=',')
		csv_o.writerow(full_result[q])
	file_o.close()	