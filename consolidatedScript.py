

from s4l_v1 import *
import XCoreModeling as XCore
import re
import numpy as np
import os
import scipy.io as io



pathway_folder = None
field_export_folder = None
Etan_folder = None  # folder for field along lead path

EXTRACT_LEAD_PATH = True  # set to true if you want to run the lead path extraction stub of the script
EXPORT_FIELD = True
CALC_M = True
FIELD_ALONG_PATH = True
CALC_HEAT_VOLTAGE = True


os.chdir(sys.path[0])  # set the working directory for scripts as the directory in which this script is stored. Necssary for relative file path access

'''
SECTION 0 OF SCRIPT

to calculate lead voltage/temperature
this script expects the active model to have the lead paths as an model entity. This script then indiviualizes the the paths, and store  each 
path into a file

the script generates the points on the path by converting the path in the model into a curve, and then iterates
over the curve in 'point_step" increment to evaluate the Vec at each point

if the user can generate individual lead paths or if the paths pre-exist, then this script section may be skipped


'''
print ("starting lead path extraction from model")

if (EXTRACT_LEAD_PATH):
	# Build up the segment dictionary
	entity_dict = model.AllEntities();  #all the entities in the active model

	# Input parameters for pathway
	pathway_folder_list = ['pathway_Fat'];    # the name of folder where the indiviidual lead paths will be stored for this virtual family member MAY PROMPT USER FOR THIS
	point_step = 2;	# Sample resolution for points on pathway.

	for pathway_folder in pathway_folder_list:
		 
		#result_folder = pathway_folder;
		
		

		# Build up the output directory
		if (not os.path.exists(pathway_folder)):
			os.mkdir(pathway_folder);
		group_list = entity_dict[pathway_folder]  # the pathways as a group
		pathway_list = group_list.Entities;  #all the different pathway in the model as individuals
		for pathway in pathway_list:  #for each individual pathway
			output_list = [];  #output list
			pathway_wire = model.GetWires(pathway);
			pathway_curve = pathway_wire[0].GetGeometry(True);  # not sure what the getgeometry(true) does. It is like it converts the path to an equation-based curve
			start = pathway_curve.ParameterRange.Start;  # in test, start is 0 but end the length of the wire which may not be integer
			end = pathway_curve.ParameterRange.End;  
			point = start
			while point<end:
				output_list.append(pathway_curve.Eval(point));
				point += point_step;
			output_list.append(pathway_curve.Eval(end));
			output_file_name = pathway_folder + '\\' + pathway.Name + '.txt';  #file name is extracted from the model name
			output_file = open(output_file_name,'w'); #file to which a file is writter
			for data in output_list:
				str_out = str(data[0]/1e3) + '\t' + str(data[1]/1e3) + '\t' + str(data[2]/1e3) + '\n';
				output_file.write(str_out);
			output_file.close();
			
			
		
	pathway_folder = pathway_folder  # this is done so that pathway folder can be accessible in later part of this script, if there multiple pathways, the last one will be used

	print pathway_folder
print ("end of lead path extraction from model")


if (EXPORT_FIELD):
	
	
	print ("starting field export..." )
	'''

	SECTION 1 OF SCRIPT

	export the E-field and B-1 field of the active model
	script 2 of the temperature voltage script series
	modified by cfawole


	saves the real and imaginary parts of the B1 puluse and B1 minuse fields at the center(in the head region) of the model to  txt file
	Still do not understand why the B1 field is calculated at the same point (center) irrespective of the relative location of the coil to body
	saves the e fields, SAR, etc to mat file
	'''


	import s4l_v1
	import numpy
	import re
	import os
	import scipy.io
	import XCoreModeling as XCore

	field_export_folder = 'full EM fields\\'   # top directory where the  EM fields will be saved  PROMPT USER FOR THIS

	if (not os.path.exists(field_export_folder)):
		os.mkdir(field_export_folder); 
		
		
	B1_file_folder = field_export_folder+'B1 file\\'      #folder for the B1 field
	if (not os.path.exists(B1_file_folder)):
		os.mkdir(B1_file_folder); 
	sim_list = list(s4l_v1.document.AllSimulations);  #simulation list
	project_name = sim_list[0].GetOutputFileName()
	project_name = project_name[0:project_name.find('.smash')].split('/')[-1]
	result_file = field_export_folder + 'B1 file\\'+project_name +'.txt'
	result_file_id = open(result_file,'w')
	result_folder = field_export_folder + project_name[0:-4]+'\\'


	if (not os.path.exists(result_folder)):
		os.mkdir(result_folder);
		
	result_folder = result_folder+project_name[-3:]+'\\'

	if (not os.path.exists(result_folder)):
		os.mkdir(result_folder);

	segment_dict = s4l_v1.model.AllEntities();
	head_box = segment_dict['Skull_cortical'];
	center = segment_dict['center'].Position   #center point of the head


	if 'head' in project_name:
		head_box_p0 = XCore.GetBoundingBox([head_box])[0]
		head_box_p1 = XCore.GetBoundingBox([head_box])[1]
		
	for sim in sim_list:
		
		sim_name = sim.Name;
		
		sim_info = re.split('_',sim_name);
		x_shift = float(sim_info[2]);  #get the relative shift of the model to the coil from the simulation name
		y_shift = float(sim_info[3]);
		z_shift = float(sim_info[4]);

		results = sim.Results();
		field_sensor = results['Overall Field'];

		efield_sensor = field_sensor['EM E(x,y,z,f0)']
		energy_sensor = field_sensor['El. Loss Density(x,y,z,f0)']
		SAR_sensor = field_sensor['SAR(x,y,z,f0)']
		B1_sensor = field_sensor['B1(x,y,z,f0)']

		efield_sensor.Update()
		energy_sensor.Update()
		SAR_sensor.Update()
		B1_sensor.Update()

		efield_array = efield_sensor.Data.Field(0)
		energy_array = energy_sensor.Data.Field(0)
		SAR_array = SAR_sensor.Data.Field(0)
		B1_array = B1_sensor.Data.Field(0)
		grid = SAR_sensor.Data.Grid;
		xList = grid.XAxis;
		yList = grid.YAxis;
		zList = grid.ZAxis;
		[center_x,center_y,center_z] = grid.FindCell(center[0]/1000.0,center[1]/1000.0,center[2]/1000.0)
		center_index = grid.ComputeCellIndex(center_x,center_y,center_z)
		B1_plus_real = B1_array[center_index][0].real
		B1_plus_imag = B1_array[center_index][0].imag
		B1_minus_real = B1_array[center_index][1].real
		B1_minus_imag = B1_array[center_index][1].imag

		if 'head' in project_name:  # this applies if dealing with head coil
			bx_p0 = numpy.divide(head_box_p0+s4l_v1.Vec3(x_shift,y_shift,z_shift)-s4l_v1.Vec3(15,30,60),1000.0)
			bx_p1 = numpy.divide(head_box_p1+s4l_v1.Vec3(x_shift,y_shift,z_shift)+s4l_v1.Vec3(35,30,35),1000.0)
			[low_x,low_y,low_z] = grid.FindCell(bx_p0[0],bx_p0[1],bx_p0[2])
			[up_x, up_y, up_z]  = grid.FindCell(bx_p1[0],bx_p1[1],bx_p1[2])
			Xn = xList.size
			Yn = yList.size
			Zn = zList.size
			efield_array = numpy.reshape(efield_array,(Zn-1,Yn-1,Xn-1,3))
			energy_array = numpy.reshape(energy_array,(Zn-1,Yn-1,Xn-1))
			SAR_array = numpy.reshape(SAR_array,(Zn-1,Yn-1,Xn-1))
			xList = xList[low_x:up_x+1]
			yList = yList[low_y:up_y+1]
			zList = zList[low_z:up_z+1]
			efield_array = efield_array[low_z:up_z,low_y:up_y,low_x:up_x,:]
			energy_array = energy_array[low_z:up_z,low_y:up_y,low_x:up_x]
			SAR_array = SAR_array[low_z:up_z,low_y:up_y,low_x:up_x]
			new_size = numpy.size(energy_array)
			efield_array = numpy.reshape(efield_array,(new_size,3))
			energy_array = numpy.reshape(energy_array,new_size)
			SAR_array = numpy.reshape(SAR_array,new_size)
			
		result_file_id.write(sim.Name+'\t' + str(B1_plus_real) + '\t' + str(B1_plus_imag)+'\t'+\
		str(B1_minus_real) + '\t' + str(B1_minus_imag)+'\n')    # write B1 field to file
		
		out_value = {'x':xList,'y':yList,'z':zList,'E':efield_array,'Loss':energy_array,'SAR':SAR_array}   #note the conbination of x,y,z do not form a tupule, each is to be manipulated independently becase the grid may have different sizes in the x, y and z directions
		scipy.io.savemat(result_folder+sim.Name,out_value,True)  # write model fields to file
	result_file_id.close()

	field_export_folder = result_folder  # so that other part of the stub can use this
	
	
	print ("End of field export ")
if (CALC_M):
	#scipt stub number 3 in the voltage temperature script series
	#script calculates the M-matrix
	
	print ("Starting M calculation...")

	import numpy as np
	import os
	import scipy.io as sio
	import matplotlib.pyplot as plt

	c = 299792458;   # speed of light
	pi = np.pi;   # pi
	mu = 4*pi*1e-7;  #  vacuum permeability
	eps0 = 1.0/(c*c)/mu; #Vacuum permittivity
	


	 # field_export_folder is name of directory that contains sin and cosin efields for different simulations
	
	if field_export_folder is None:  # this is necessary in case we run this stub before running other stubs, we dont want none folder
		field_export_folder = 'E:\\c_test\\Sasis Project\\Python MRI TF Scripts\\full EM fields\\3T bodycoil\\Fat'
		print (" HARD CODING ALERT:::the directory in which to look for the E field has been set to 'E:\\c_test\\Sasis Project\\Python MRI TF Scripts\\full EM fields\\3T bodycoil\\Fat'")

	#sub_folder1 = folder+'\\cos\\';  # name of directory containing  the human body fields that was obtained by cosine excitation of the MRI coil
	#sub_folder2 = folder+'\\sin\\';  # name of directory containing  the human body fields that was obtained by sine excitation of the MRI coil


	M_folder = field_export_folder+'\\M\\';  # the folder in which the M-Matrix will be stored

	if not (os.path.isdir(M_folder)): # if M folder directory does not exist, then create it 	
		 os.makedirs(M_folder)
		 print "made M folder"



	all_files_cos = [ f for f in os.listdir(field_export_folder) if 'Cos' in str(f) ] #create a list of all files in with cosine exciation in the directory
	all_files_sin = [ f for f in os.listdir(field_export_folder) if 'Sin' in str(f) ] #create a list of all files in with cosine exciation in the directory
	
	all_files = all_files_cos[:] # this will be used to decouple cos and sin labels in the final M file
	
	all_file_index = 0
	for afile in all_files_cos:
		splitFile = afile.split('_')
		all_files[all_file_index ] = splitFile[0]+'_'+splitFile[2]+'_'+splitFile[3]+'_'+splitFile[4]
		all_file_index = all_file_index+1

	N = len(all_files_cos);  # there should be an equal number of sine and cos  files

	for i in range(N):
		cosE_mat = sio.loadmat(os.path.join(field_export_folder,all_files_cos[i]))   #load the .mat file containing the cosine EM fields of the body model
		#cosE_mat['a']
		
		'''
		Load all the parameters in the mat file
		'''
		x =  cosE_mat['x']
		y = cosE_mat['y']
		z = cosE_mat['z']
		E = cosE_mat['E']
		Loss = cosE_mat['Loss']
		SAR = cosE_mat['SAR']
		
		
		Xn = len(x[0])
		Yn = len(y[0])
		Zn = len(z[0])

		
		Ex_cos = np.reshape(E[:,0],(Xn-1,Yn-1,Zn-1),order = "F");  # reshape the Efield, stored in 1D form by Sim4Life (probably for covenience) into 3D form
		Ey_cos = np.reshape(E[:,1],(Xn-1,Yn-1,Zn-1),order = "F");
		Ez_cos = np.reshape(E[:,2],(Xn-1,Yn-1,Zn-1),order = "F");
		
		
		
		E_rms = (np.power(np.absolute(Ex_cos),2) + np.power(np.absolute(Ey_cos),2) + np.power(np.absolute(Ez_cos),2))/2  
		Sigma = np.divide(np.reshape(Loss,(Xn-1,Yn-1,Zn-1),order = "F"),E_rms);  # divide loss by Erms to find conductivity
		Rho = np.reshape(np.divide(Loss,SAR),(Xn-1,Yn-1,Zn-1),order = "F");  # divide loss by SAR to find density
		
		Rho = np.nan_to_num(Rho)  # Replace NaN with zero and infinity with large finite numbers. This corrects for NaN in the  'undefined regions' in the model
		
		mg = np.meshgrid(x,y,z);  
		X = mg[0]
		Y = mg[1]
		Z = mg[2]
		
		print Xn
		print Yn
		print Zn
		
		print "X:"
		print(X.shape)
		
		print "Y:"
		print(Y.shape)
		
		print "Z:"
		print(Z.shape)
		
		
		'''
		calculating volume
		'''
		vx = np.diff(X[0:Yn-1,0:Xn,0:Zn-1],1,1)
		vy = np.diff(Y[0:Yn,0:Xn-1,0:Zn-1],1,0)
		vz = np.diff(Z[0:Yn-1,0:Xn-1,0:Zn],1,2)
		
		v = np.multiply(np.multiply(vx,vy),vz)  # np.multipy :)

		V = np.transpose(v,(1,0,2))
		
		mass = (np.sum(np.multiply(V,Rho)))   # calculating mass, becuase we have replaced NaN values in rho with zero, the undefined regions in the model will not contribute to the mass
		
		print mass   # sanity check. The value of mass should be a reasonable human body mass
		
		Sigma_V = np.multiply(Sigma,V)
		
		
		
		sinE_mat = sio.loadmat(os.path.join(field_export_folder,all_files_sin[i]))   #load the .mat file containing the sin E fields of the body model
		E = sinE_mat['E']
		Ex_sin = np.reshape(E[:,0],(Xn-1,Yn-1,Zn-1),order = "F");
		Ey_sin = np.reshape(E[:,1],(Xn-1,Yn-1,Zn-1),order = "F");
		Ez_sin = np.reshape(E[:,2],(Xn-1,Yn-1,Zn-1),order = "F");
		
		M  = np.array([[0+0j,0+0j],[0+0j,0+0j]],dtype = complex)
		
		#calculate the M matrix
		M[0,0] = np.sum(np.multiply(np.multiply(Ex_cos,np.conj(Ex_cos))+np.multiply(Ey_cos,np.conj(Ey_cos))+np.multiply(Ez_cos,np.conj(Ez_cos)),Sigma_V))
		M[0,1] = np.sum(np.multiply(np.multiply(Ex_cos,np.conj(Ex_sin))+np.multiply(Ey_cos,np.conj(Ey_sin))+np.multiply(Ez_cos,np.conj(Ez_sin)),Sigma_V))
		M[1,0] = np.conj(M[0,1])
		M[1,1] = np.sum(np.multiply(np.multiply(Ex_sin,np.conj(Ex_sin))+np.multiply(Ey_sin,np.conj(Ey_sin))+np.multiply(Ez_sin,np.conj(Ez_sin)),Sigma_V))
		M = M/2/mass
		print M
		
		sio.savemat(os.path.join(M_folder,all_files[i]), {'M':M})
			

			

	print ("End of M Calculation")		
if (FIELD_ALONG_PATH):

	print ("starting extraction of field along lead paths")

	'''
	#PYTHON SCRIPT STUB 4, SCRIPT NEEDED FOR VOLTAGE TEMPERATURE CALCULATION OF THE IPG LEAD
	#THE MODEL SHOULD BE ACTIVE AND SIMULATION RESULTS OF THE MODEL SHOULD HAVE BEEN COMPLETED BEFORE THE SCRIPT IS RUN
	# original code modified by cfawole
	#owning U of H code

	#extract human body electric fields along the lead paths specified in the text files in the pathway_xxx directory
	#save the fields in a folder with the same name as the short name of simulation file
	'''

	import s4l_v1
	import numpy
	import os
	import re

	sim_list = list(s4l_v1.document.AllSimulations);
	project_name = sim_list[0].GetOutputFileName()   #long name of simulation file
	result_subfolder = project_name[0:project_name.find('.smash')].split('/')[-1]   #short name of simulation file

	
	#pathway_folder = 'E:\\c_test\\Transalating Juliannas Code\\pathway_Fat\\'  # directory of the different  lead pathway files
	
	pathway_folder = pathway_folder  
	
	if (pathway_folder is None):
		raise Exception ("Pathway folder is None, allow code for generating pathway folder to run")
	

	Etan_folder = 'Etan6\\';  # top level directory of where the fields will be saved  PROMT USER FOR THIS

	if (not os.path.exists(Etan_folder)):
		os.mkdir(Etan_folder);                       # make top directory folder if it does not exist
		
	Etan_folder = Etan_folder+result_subfolder;		# actual folder where field files will be stored
	file_list = os.listdir(pathway_folder);           
		
	if (not os.path.exists(Etan_folder)):
		os.mkdir(Etan_folder);

	
	sim_index = -1
	for sim in sim_list:      #for each simulation in simulation list
		sim_index = sim_index+1
		sim_name = sim.Name;
		sim_info = re.split('_',sim_name);
		body_type = sim_info[0];
		x_shift = float(sim_info[2]);
		y_shift = float(sim_info[3]);
		z_shift = float(sim_info[4]);
		
		results = sim.Results()
		field_sensor = results[ 'Overall Field' ]
		efield_sensor = field_sensor[ 'EM E(x,y,z,f0)' ]
		efield_sensor.Update()	
		efield_data = efield_sensor.Data;
		efield_grid = efield_data.Grid;
		efield_array = efield_data.Field(0)   #the e-field
		
		for coo_file in file_list:   #coo_file is each lead path file
			dot_index = coo_file.find('.');
			if coo_file[dot_index:] == '.txt':
				InputFileName = pathway_folder + '/'+coo_file;
				path_coo_file = open(InputFileName,'r');
				coo_list_str = path_coo_file.readlines();
				path_coo_file.close();
				result_file_name = Etan_folder + '/' + sim_name + '_' + coo_file;
				field_result_file = open(result_file_name,'w');
				str_out = str('%coord_x\tcoord_y\tcoord_z\tEtan_real\tEtan_imag\n');
				field_result_file.write(str_out);
				j = 0;
				coo_list = [];
				for a_coordinat in coo_list_str:
					cord = re.split('\s+', a_coordinat);
					str_x = cord[0];
					coo_x = float(str_x);  # x-cooridinate in lead path file
					str_y = cord[1];
					coo_y = float(str_y);
					str_z = cord[2];
					coo_z = float(str_z);
					
					coo_value = coo_x+x_shift/1000.0,coo_y+y_shift/1000.0,coo_z+z_shift/1000.0;  #if the body model is shifted, shift the each lead path correspondingly
					temp_point = s4l_v1.Vec3(coo_value[0],coo_value[1],coo_value[2]);
					coo_list.append(temp_point);
				
				IndexNum = len(coo_list);
				Mid_Point = [];
				Point_Vec = [];
				for index1 in range(IndexNum-1):
					temp_point2 = coo_list[index1+1] + coo_list[index1]; # add the  next point to the current point
					Mid_Point.append(s4l_v1.Vec3([0.5*temp_point2[i] for i in range(0,3)]));  # find the mid point of the point obtained above 
					Point_Vec.append(coo_list[index1+1] - coo_list[index1]); #subtract the next point from the current point. This will be used to calculate the length of a segment
			
				ETan = [];
				Mag_Etan = [];
				TotalLength = 0;
				for index2 in range(len(Mid_Point)):
					cell_index = efield_grid.FindCell(Mid_Point[index2][0],Mid_Point[index2][1],Mid_Point[index2][2]);
					point_index = efield_grid.ComputeCellIndex(cell_index[0],cell_index[1],cell_index[2]);
					efield = efield_array[point_index]
					e_x = efield[0];  # x-component of the e-field
					e_y = efield[1];
					e_z = efield[2];
					Segment_Length = numpy.sqrt((Point_Vec[index2][0])**2+(Point_Vec[index2][1])**2\
					+(Point_Vec[index2][2])**2);	#length of a segment of the lead path
					TotalLength = TotalLength + Segment_Length; # total length of lead path
					UnitVec = Point_Vec[index2]*(1/Segment_Length);   #unit vector of a segment of the lead path
					ETanTmp = (UnitVec[0]*e_x+UnitVec[1]*e_y+UnitVec[2]*e_z);  # dot product of E-field and unit vector
					ETan.append(ETanTmp);
					str_out = str(Mid_Point[index2][0])+'\t'+str(Mid_Point[index2][1]) +'\t'\
					+str(Mid_Point[index2][2])+'\t'+ str(ETanTmp.real) + '\t' + str(ETanTmp.imag)+'\n';  #each output line is mid point and the real/imaginary path of E-field
					field_result_file.write(str_out);
				field_result_file.close();
				
				
	print("End of field extraction")
	
	
if (CALC_HEAT_VOLTAGE):	

	print("starting heat voltage calculation...")
	
	#stub 5 script in the lead voltage heat script series
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
	
	if (Etan_folder ==  None):
		Etan_folder = 'E:\\c_test\\Sasis Project\\Python MRI TF Scripts\\Etan6\\3T bodycoil Fat'
		print ("etan folder not prexisting in this script scope. It will be set to "+Etan_folder)
	
	tf_file = '1.5T_heating_ring_106_304.csv'  # prompt user for this
	print ("the transfer function file that will be used is this "+tf_file)
	#tf_file = 'E:\\human model\\TF\\1.5T_heating_ring_106_304.csv'   # the matalab version uses xls but this version uses CSV. User needs to convert to transfer function file to CSV before use
	norm_folder = 'full EM fields'
	print ("ALERT:: norm folder hardcoded as "+norm_folder)
	result_type = 'heating'
	norm_type = 'B1_wbSAR'
	heat_voltage_result_file = 'example6.csv'
	print ("result file is hardcoded as: "+result_file)
	
	


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
		
		
		all_B1_file = B1_field_folder+sim_type+'.txt'
		
		f_B1 = open(all_B1_file)
		
		
		
		B1_cos_file = B1_field_folder+sim_type+' cos.txt'
		
		B1_sin_file = B1_field_folder+sim_type+' sin.txt'
		
		cosB  = open(B1_cos_file,'w')
		sinB  = open(B1_sin_file,'w')
		
		B1_lines = f_B1.readlines()
		
		
		
		
		
		
		for B1_line in B1_lines:
			
			if 'Cos' in B1_line:
				cosB.write(B1_line)
			if 'Sin' in B1_line:
				sinB.write(B1_line)
		cosB.close()
		sinB.close()
				
		
		
		
		#B1_cos_file = B1_field_folder+sim_type+' cos.txt'
		
		#B1_sin_file = B1_field_folder+sim_type+' sin.txt'
		
		
		
		
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
		
		if field_export_folder is None:  # this is necessary in case we run this stub before running other stubs, we dont want none folder
			field_export_folder = 'E:\\c_test\\Sasis Project\\Python MRI TF Scripts\\full EM fields\\3T bodycoil\\Fat'
			print (" HARD CODING ALERT:::the directory in which to look for the M folder has been set to 'E:\\c_test\\Sasis Project\\Python MRI TF Scripts\\full EM fields\\3T bodycoil\\Fat'")

		
		M_file_folder = field_export_folder+'\\M'
		
		
		
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
		
		#partition Etan_folder into cos and sin folders
		
		import shutil
		
		Etan_folder_cos = Etan_folder+' cos'
		Etan_folder_sin = Etan_folder+' sin'
		
		if os.path.exists(Etan_folder_cos):
			shutil.rmtree(Etan_folder_cos)
			os.makedirs(Etan_folder_cos)
			
			
		if os.path.exists(Etan_folder_sin):
			shutil.rmtree(Etan_folder_sin)
			os.makedirs(Etan_folder_sin)

		
		
		src_files = os.listdir(Etan_folder)
		for file_name in src_files:
			full_file_name = os.path.join(Etan_folder, file_name)
			if (os.path.isfile(full_file_name)):
				if "Cos" in str(full_file_name): 
					shutil.copy(full_file_name, Etan_folder_cos)
					#c_file = os.path.join(Etan_folder_cos,full_file_name)
					#os.rename( c_file, c_file[0:c_file.find('Cos')]+c_file[c_file.find('Cos')+4:])
					
				if "Sin" in str(full_file_name): 
					shutil.copy(full_file_name, Etan_folder_sin)
		
		
		try:
			#strip off excitation type from file name for backward compactiblity
			cos_files = os.listdir(Etan_folder_cos)
			for cos_file in cos_files:
				full_cos_file_name = os.path.join(Etan_folder_cos, cos_file)
				os.rename( full_cos_file_name, full_cos_file_name[0:full_cos_file_name.find('Cos')]+full_cos_file_name[full_cos_file_name.find('Cos')+4:])
			
			
			sin_files = os.listdir(Etan_folder_sin)
			for sin_file in sin_files:
				full_sin_file_name = os.path.join(Etan_folder_sin, sin_file)
				os.rename( full_sin_file_name, full_sin_file_name[0:full_sin_file_name.find('Sin')]+full_sin_file_name[full_sin_file_name.find('Sin')+4:])

		except Exception as e:
			raise e
				
		
		
		file_folder = Etan_folder_cos
		
		
		
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
			file_o = open(heat_voltage_result_file, "ab")
			csv_o = csv.writer(file_o, delimiter=',')
			csv_o.writerow(full_result[q])
		file_o.close()	
	
	print("end of heat voltage calculation")
	
	