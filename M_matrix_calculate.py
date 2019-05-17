# -*- coding: utf-8 -*-
#scipt number 3 in the voltage temperature script series
#script calculates the M-matrix

import numpy as np
import os
import scipy.io as sio
import matplotlib.pyplot as plt

c = 299792458;   # speed of light
pi = np.pi;   # pi
mu = 4*pi*1e-7;  #  vacuum permeability
eps0 = 1.0/(c*c)/mu; #Vacuum permittivity

folder = 'E:\\c_test\\Transalating Juliannas Code\\full EM fields';    # name of directory that contains the directrory of sin and cosin efields for different simulations

sub_folder1 = folder+'\\cos\\';  # name of directory containing  the human body fields that was obtained by cosine excitation of the MRI coil
sub_folder2 = folder+'\\sin\\';  # name of directory containing  the human body fields that was obtained by sine excitation of the MRI coil


M_folder = folder+'\\M\\';  # the folder in which the M-Matrix will be stored

if not (os.path.isdir(M_folder)): # if M folder directory does not exist, then create it 	
	 os.makedirs(M_folder)



all_files = [f for f in os.listdir(sub_folder1) if os.path.isfile(os.path.join(sub_folder1, f))]  # create a list of all files in the cosine directory

N = len(all_files);

for i in range(N):
	cosE_mat = sio.loadmat(os.path.join(sub_folder1,all_files[i]))   #load the .mat file containing the cosine EM fields of the body model
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
	
	
	
	sinE_mat = sio.loadmat(os.path.join(sub_folder2,all_files[i]))   #load the .mat file containing the sin E fields of the body model
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
	
	