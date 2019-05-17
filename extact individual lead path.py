#Script number 4 or 0 in the script series to calculate lead voltage/temperature
# this script expects the active model to have the lead paths as an model entity. This script then indiviualizes the the paths, and store  each 
#path into a file

#the script generates the points on the path by converting the path in the model into a curve, and then iterates
#over the curve in 'point_step" increment to evaluate the Vec at each point

#if the user can generate individual lead paths or if the paths pre-exist, then this script may be skipped

#Example of model that contain the lead path for the FATS virtual family member is pathways_FAT_0_0_0.smash

from s4l_v1 import *
import XCoreModeling as XCore
import re
import numpy as np
import os
import scipy.io as io


# Build up the segment dictionary
entity_dict = model.AllEntities();  #all the entities in the active model

# Input parameters for pathway
pathway_folder_list = ['pathway_Fat'];    # the name of folder where the indiviidual lead paths will be stored for this virtual family member 
point_step = 2;	# Sample resolution for points on pathway.

for pathway_folder in pathway_folder_list:
	result_folder = 'E:\\c_test\\Transalating Juliannas Code\\';  # customize result folder 
	result_folder = result_folder + pathway_folder;
	
	

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