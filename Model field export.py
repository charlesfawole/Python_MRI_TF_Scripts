#export the E-field and B-1 field of the active model
#script 2 of the temperature voltage script series
#modified by cfawole
#owing U of H script

#saves the real and imaginary parts of the B1 puluse and B1 minuse fields at the center(in the head region) of the model to  txt file
#Still do not understand why the B1 field is calculated at the same point (center) irrespective of the relative location of the coil to body
#saves the e fields, SAR, etc to mat file

import s4l_v1
import numpy
import re
import os
import scipy.io
import XCoreModeling as XCore

folder = 'E:\c_test\Transalating Juliannas Code\\full EM fields\\'   # top directory where the  EM fields will be saved

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
	x_shift = float(sim_info[1]);  #get the relative shift of the model to the coil from the filename
	y_shift = float(sim_info[2]);
	z_shift = float(sim_info[3]);

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