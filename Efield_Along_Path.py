# -*- coding: utf-8 -*-
#PYTHON SCRIPT #X   X SCRIPT NEEDED FOR VOLTAGE TEMPERATURE CALCULATION OF THE IPG LEAD
#THE MODEL SHOULD BE ACTIVE AND SIMULATION RESULTS OF THE MODEL SHOULD HAVE BEEN COMPLETED BEFORE THE SCRIPT IS RUN
# original code modified by cfawole
#owning U of H code

#extract human body electric fields along the lead paths specified in the text files in the pathway_xxx directory
#save the fields in a folder with the same name as the short name of simulation file


import s4l_v1
import numpy
import os
import re

sim_list = list(s4l_v1.document.AllSimulations);
project_name = sim_list[0].GetOutputFileName()   #long name of simulation file
result_subfolder = project_name[0:project_name.find('.smash')].split('/')[-1]   #short name of simulation file


pathway_folder = 'E:\\c_test\\Transalating Juliannas Code\\pathway_Fat\\'  # directory of the different  lead pathway files
result_folder = 'E:\\c_test\Transalating Juliannas Code\\Etan2\\';  # top level directory of where the fields will be saved

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
	x_shift = float(sim_info[1]);
	y_shift = float(sim_info[2]);
	z_shift = float(sim_info[3]);
	
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