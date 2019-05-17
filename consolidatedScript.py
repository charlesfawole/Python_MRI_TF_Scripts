

from s4l_v1 import *
import XCoreModeling as XCore
import re
import numpy as np
import os
import scipy.io as io





EXTRACT_LEAD_PATH = False  # set to true if you want to run the lead path extraction stub of the script
EXPORT_FIELD = True

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

	
