

from s4l_v1 import *
import XCoreModeling as XCore
import re
import numpy as np
import os
import scipy.io as io



result_folder = None  # give gloabal access to the  generic result folder used in stubs of script
pathway_folder = None

EXTRACT_LEAD_PATH = True  # set to true if you want to run the lead path extraction stub of the script
EXPORT_FIELD = False
CALC_M = True
FIELD_ALONG_PATH =True

'''
SECTION 0 OF SCRIPT

to calculate lead voltage/temperature
this script expects the active model to have the lead paths as an model entity. This script then indiviualizes the the paths, and store  each 
path into a file

the script generates the points on the path by converting the path in the model into a curve, and then iterates
over the curve in 'point_step" increment to evaluate the Vec at each point

if the user can generate individual lead paths or if the paths pre-exist, then this script section may be skipped


'''

if (EXTRACT_LEAD_PATH):
	# Build up the segment dictionary
	entity_dict = model.AllEntities();  #all the entities in the active model

	# Input parameters for pathway
	pathway_folder_list = ['pathway_Fat'];    # the name of folder where the indiviidual lead paths will be stored for this virtual family member MAY PROMPT USER FOR THIS
	point_step = 2;	# Sample resolution for points on pathway.

	for pathway_folder in pathway_folder_list:
		 
		result_folder = pathway_folder;
		
		

		# Build up the output directory
		if (not os.path.exists(result_folder)):
			os.mkdir(result_folder);
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
			output_file_name = result_folder + '\\' + pathway.Name + '.txt';  #file name is extracted from the model name
			output_file = open(output_file_name,'w'); #file to which a file is writter
			for data in output_list:
				str_out = str(data[0]/1e3) + '\t' + str(data[1]/1e3) + '\t' + str(data[2]/1e3) + '\n';
				output_file.write(str_out);
			output_file.close();
			
			
		
	pathway_folder = result_folder  # this is done so that pathway folder can be accessible in later part of this script, if there multiple pathways, the last one will be used

	
if (EXPORT_FIELD):
	
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

	folder = 'full EM fields\\'   # top directory where the  EM fields will be saved  PROMPT USER FOR THIS

	if (not os.path.exists(folder)):
		os.mkdir(folder); 
		
		
	B1_file_folder = folder+'B1 file\\'      #folder for the B1 field
	if (not os.path.exists(B1_file_folder)):
		os.mkdir(B1_file_folder); 
	sim_list = list(s4l_v1.document.AllSimulations);  #simulation list
	project_name = sim_list[0].GetOutputFileName()
	project_name = project_name[0:project_name.find('.smash')].split('/')[-1]
	result_file = folder + 'B1 file\\'+project_name +'.txt'
	result_file_id = open(result_file,'w')
	result_folder = folder + project_name[0:-4]+'\\'


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
	


	folder = result_folder;    # name of directory that contains sin and cosin efields for different simulations
	
	if folder is None:  # this is necessary in case we run this stub before running other stubs, we dont want none folder
		folder = 'E:\\c_test\\Sasis Project\\full EM fields\\3T bodycoil\\Fat\\'
		print (" HARD CODING ALERT:::the directory in which to look for the E field has been set to 'E:\\c_test\\Sasis Project\\full EM fields\\3T bodycoil\\Fat\\'")

	#sub_folder1 = folder+'\\cos\\';  # name of directory containing  the human body fields that was obtained by cosine excitation of the MRI coil
	#sub_folder2 = folder+'\\sin\\';  # name of directory containing  the human body fields that was obtained by sine excitation of the MRI coil


	M_folder = folder+'\\M\\';  # the folder in which the M-Matrix will be stored

	if not (os.path.isdir(M_folder)): # if M folder directory does not exist, then create it 	
		 os.makedirs(M_folder)



	all_files_cos = [ f for f in os.listdir(folder) if 'Cos' in str(f) ] #create a list of all files in with cosine exciation in the directory
	all_files_sin = [ f for f in os.listdir(folder) if 'Sin' in str(f) ] #create a list of all files in with cosine exciation in the directory
	
	all_files = all_files_cos[:] # this will be used to decouple cos and sin labels in the final M file
	
	all_file_index = 0
	for afile in all_files_cos:
		splitFile = afile.split('_')
		all_files[all_file_index ] = splitFile[0]+'_'+splitFile[2]+'_'+splitFile[3]+'_'+splitFile[4]
		all_file_index = all_file_index+1

	N = len(all_files_cos);  # there should be an equal number of sine and cos  files

	for i in range(N):
		cosE_mat = sio.loadmat(os.path.join(folder,all_files_cos[i]))   #load the .mat file containing the cosine EM fields of the body model
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
		
		
		
		sinE_mat = sio.loadmat(os.path.join(folder,all_files_sin[i]))   #load the .mat file containing the sin E fields of the body model
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
			

			

			
if (FIELD_ALONG_PATH):

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
	

	result_folder = '\\Etan5\\';  # top level directory of where the fields will be saved  PROMT USER FOR THIS

	if (not os.path.exists(result_folder)):
		os.mkdir(result_folder);                       # make top directory folder if it does not exist
		
	result_folder = result_folder+result_subfolder;		# actual folder where field files will be stored
	file_list = os.listdir(pathway_folder);           
		
	if (not os.path.exists(result_folder)):
		os.mkdir(result_folder);


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
				result_file_name = result_folder + '/' + sim_name + '_' + coo_file;
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