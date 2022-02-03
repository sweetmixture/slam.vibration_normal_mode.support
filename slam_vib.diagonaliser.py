#!/bin/python

import numpy as np
import sys
import datetime

kd = lambda i,j : 1. if (i==j) else 0.

class slam_vib_dig:

	def __init__(self,*argv):
		#DEBUG?
		self.number_of_atoms = int(argv[0])			# GET NUMBER OF CORES IN THE CLUSTER
		self.delta	     = float(argv[1])*2.		# Mult "2" > */2dx
		#self.number_of_atoms = int(argv[1])			# GET NUMBER OF CORES IN THE CLUSTER
		#self.delta	     = float(argv[2])
		self.main_config_file= "geo.txt"
		self.main_config     = None
		self.main_config_type= None		

		#self.wavenumber_conv=5.8921E-5				# COV FACTOR INTO cm-1	... 2.71918.E+05? .. or 1./3.67758.E-06
		self.wavenumber_conv=3.67758E-06			# COV FACTOR INTO cm-1	(factor itself is in 'cm^2' unit : used for division

		# conversion factor  ...
		'''
			Unit of internal hessian elements are : eV / amu / angs^2 (also the eigenvalues)
			it has to be converted into (cm-1) unit.

			using relation 4 pi^2 c^2 v'^2 = (eval)

			v'^2 = (eval) / (4 pi^2 c^2 )	: unit of eval is eV / angs / angs /amum, unit of denominator (cm-2)(s^2)

			eV / angs / angs / amu <=> const * N (eV/angs) m-1 (/angs) kg-1 (/amu)

						<=> kg m / s^2 * m-1 * kg-1 * const <=> const (1/s^2)

		'''


		self.cla_cnt = 0;	self.sp_cnt = 0			# ATOM SPECIES COUNT

		self.x_ind = 0						#
		self.y_ind = 1						#
		self.z_ind = 2						# INTERNAL INDICATORS

		self.symm_checker_tol = 0.005				# HESSIAN SYMMETRIC CHECK TOLERANCE
	
		self.atomic_mass_dic = {'B': 10.811    , 'C' : 12.0107   ,'F'  : 18.9984032,'H'  : 1.00794   , 'D' : 2.014102 ,'T' : 3.01605
		,'I'  : 126.90447 , 'K' : 39.0983   , 'N' : 14.0067   ,'O'  : 15.9994   ,'P'  : 30.973762 , 'S' : 32.065   ,'U' : 238.02891 
		,'V'  : 50.9415   , 'W' : 183.84    , 'Y' : 88.90585  ,'He' : 4.002602  ,'Li' : 6.941     ,'Be' : 9.012182 ,'Ne': 20.1797   
		,'Na' : 22.989769 ,'Mg' : 24.3050   ,'Al' : 26.9815386,'Si' : 28.0855   ,'Cl' : 35.453    ,'Ar' : 39.948    
		,'Ca' : 40.078    ,'Sc' : 44.955912 ,'Ti' : 47.867    ,'Cr' : 51.9961   ,'Mn' : 54.938045 ,'Fe' : 55.845    
		,'Co' : 58.933195 ,'Ni' : 58.6934   ,'Cu' : 63.546    ,'Zn' : 65.409    ,'Ga' : 69.723    ,'Ge' : 72.64     
		,'As' : 74.92160  ,'Se' : 78.96     ,'Br' : 79.904    ,'Kr' : 83.798    ,'Rb' : 85.4678   ,'Sr' : 87.62     
		,'Zr' : 91.224    ,'Nb' : 92.90638  ,'Mo' : 95.94     ,'Tc' : 98        ,'Ru' : 101.07    ,'Rh' : 102.90550 
		,'Pd' : 106.42    ,'Ag' : 107.8682  ,'Cd' : 112.411   ,'In' : 114.818   ,'Sn' : 118.710   ,'Sb' : 121.760   
		,'Te' : 127.60    ,'Xe' : 131.293   ,'Cs' : 132.905451,'Ba' : 137.327   ,'La' : 138.90547 ,'Ce' : 140.116   
		,'Pr' : 140.90765 ,'Nd' : 144.242   ,'Pm' : 145       ,'Sm' : 150.36    ,'Eu' : 151.964   ,'Gd' : 157.25    
		,'Tb' : 158.92535 ,'Dy' : 162.500   ,'Ho' : 164.93032 ,'Er' : 167.259   ,'Tm' : 168.93421 ,'Yb' : 173.04    
		,'Lu' : 174.967   ,'Hf' : 178.49    ,'Ta' : 180.94788 ,'Re' : 186.207   ,'Os' : 190.23    ,'Ir' : 192.217   
		,'Pt' : 195.084   ,'Au' : 196.966569,'Hg' : 200.59    ,'Tl' : 204.3833  ,'Pb' : 207.2     ,'Bi' : 208.98040 
		,'Po' : 209       ,'At' : 210       ,'Rn' : 222       ,'Fr' : 223       ,'Ra' : 226       ,'Ac' : 227       
		,'Th' : 232.03806 ,'Pa' : 231.03588 ,'Np' : 237       ,'Pu' : 244       ,'Am' : 243       ,'Cm' : 247       
		,'Bk' : 247       ,'Cf' : 251       ,'Es' : 252       ,'Fm' : 257       ,'Md' : 258       ,'No' : 259       
		,'Lr' : 262       ,'Rf' : 267       ,'Db' : 268       ,'Sg' : 271       ,'Bh' : 272       ,'Hs' : 270       
		,'Mt' : 276       ,'Ds' : 281       ,'Rg' : 280       ,'Uub': 285       ,'Uut': 284       ,'Uuq': 289      
		,'Uup': 288       ,'Uuh' : 293      ,'Uuo' : 294}      
	
	def read_main_config(self):					# BUILD INTERTIA TENSOR ACCORDING TO THE GIVEN CONFIG
		try:							
			with open(self.main_config_file) as f:
				cm = f.readline()
				rl = f.readline()
				spl = rl.split()
				
				self.cla_cnt = int(spl[0])	
				self.sp_cnt  = int(spl[1])

				self.main_config      = np.matrix( [[ 0. for i in range(3) ] for j in range(self.number_of_atoms)] )	# CREAT ARR TO SAVE MAIN CONFIG WITH ATOM NAMES
				self.main_config_type = [ None for i in range(self.number_of_atoms) ]					# SAVE ELEMENT TYPE


				cnt=0

				for i in range(self.cla_cnt + self.sp_cnt):
					if i < self.cla_cnt:
						rl = f.readline()
						spl= rl.split()

						if spl[1] == 'c' :
							self.main_config_type[cnt] = spl[0]
							self.main_config.itemset((cnt,0),spl[2])		# X
							self.main_config.itemset((cnt,1),spl[3])		# Y
							self.main_config.itemset((cnt,2),spl[4])		# Z
							cnt += 1
					else:

						rl = f.readline()
						spl= rl.split()
						self.main_config_type[cnt] = spl[0]
						self.main_config.itemset((cnt,0),spl[1])		# X
						self.main_config.itemset((cnt,1),spl[2])		# Y
						self.main_config.itemset((cnt,2),spl[3])		# Z
						cnt += 1
				#print(self.main_config)
				'''
				for i in range(self.number_of_atoms):
					if i < self.cla_cnt :
						rl = f.readline()	
						spl= rl.split()
						self.main_config_type[i] = spl[0]
						self.main_config.itemset((i,0),spl[2])		# X
						self.main_config.itemset((i,1),spl[3])		# Y
						self.main_config.itemset((i,2),spl[4])		# Z
					else:
						rl = f.readline()
						spl= rl.split()
						self.main_config_type[i] = spl[0]
						self.main_config.itemset((i,0),spl[1])		# X
						self.main_config.itemset((i,1),spl[2])		# Y
						self.main_config.itemset((i,2),spl[3])		# Z
				'''
			
		except FileNotFoundError:
			print(FileNotFoundError);	sys.exit(1)


	def build_inertia_tensor(self):
		self.inertia_tensor    = np.matrix( [[ 0. for i in range(3) ] for j in range(3) ] )			# EMPTY 3 X 3 MATRIX FOR INERTIA TENSOR
		self.com	       = np.array( [ 0. for i in range(3) ] )						# EMPTY 3D VECTOR FOR COM
		self.main_config_shift = np.matrix( [[ 0. for i in range(3) ] for j in range(self.number_of_atoms) ] )	# FOR SHIFTED (RECENTRED CONFIG)
		self.mass_sum          = np.array([0.])

		for i in range(self.number_of_atoms):												# GET COM
			self.mass_sum = self.mass_sum + self.atomic_mass_dic[ self.main_config_type[i] ]					#
			self.com.itemset(0, self.com.item(0) + self.atomic_mass_dic[ self.main_config_type[i] ]*self.main_config.item((i,0)))	#
			self.com.itemset(1, self.com.item(1) + self.atomic_mass_dic[ self.main_config_type[i] ]*self.main_config.item((i,1)))	#
			self.com.itemset(2, self.com.item(2) + self.atomic_mass_dic[ self.main_config_type[i] ]*self.main_config.item((i,2)))	#
		self.com = self.com/self.mass_sum												# COM IN 'self.com'

		for i in range(self.number_of_atoms):												# RECENTRE
			self.main_config_shift.itemset((i,0), self.main_config.item((i,0)) - self.com.item(0))					#
			self.main_config_shift.itemset((i,1), self.main_config.item((i,1)) - self.com.item(1))					# 'r_com = r - R_com'
			self.main_config_shift.itemset((i,2), self.main_config.item((i,2)) - self.com.item(2))					# RECENTRED CONFIG IN 'self.main_conifig_shift'

		# COMPUTE INERTIA TENSOR
		for i in range(self.number_of_atoms):
			
			rsq = np.array( self.main_config_shift.item((i,0))**2. + self.main_config_shift.item((i,1))**2.	# GET R^2
						+ self.main_config_shift.item((i,2))**2. ) 				#
			m   = self.atomic_mass_dic[ self.main_config_type[i] ]						# GET ATOMIC MASS
			#print(self.main_config_type[i],m,rsq)

			#self.intertial_tensor.itemset((0,0), self.inertia_tensor.item((0,0)) + m*( rsq - self.main_config_shift.item((i,0))*self.main_config_shift.item((i,0)) ))

			# I_jk = sum_i ( rsq_i * kd_jk - r_i;j * r_i;k )*m_i			

			for j in range(3):
				for k in range(3):
					self.inertia_tensor.itemset((j,k), self.inertia_tensor.item((j,k)) 
						+ m*( float(rsq)*kd(j,k) - self.main_config_shift.item((i,j))*self.main_config_shift.item((i,k)) ))
		# COMPUTE INERTIA TENSOR DONE
		self.inertia_tensor_eval, self.inertia_tensor_evec = np.linalg.eig( self.inertia_tensor )	

		''' KEY ITEMS
		
			self.inertia_tensor_eval
			self.inertia_tensor_evec
		
		'''
	
	# END BUILD INERTIA TENSOR

	def build_tf_matrix(self):

		''' 
			Build a transformation matrix, transforming weighted hessian matrix for 
			projecting out trans/rotas motions
		'''
		
		stride = self.number_of_atoms * 3
		self.tf_matrix_ws = np.matrix( [[ 0. for i in range(stride) ] for j in range(stride) ] )
	
		for i in range(self.number_of_atoms):
			self.tf_matrix_ws.itemset((i*3 + 0, 0), np.sqrt( self.atomic_mass_dic[ self.main_config_type[i] ] )) 	# TRANSLATION - X
			self.tf_matrix_ws.itemset((i*3 + 1, 1), np.sqrt( self.atomic_mass_dic[ self.main_config_type[i] ] ))	# TRANSLATION - Y
			self.tf_matrix_ws.itemset((i*3 + 2, 2), np.sqrt( self.atomic_mass_dic[ self.main_config_type[i] ] ))	# TRANSLATION - Z

		# ROTATIONAL MOTIONS
		
		X = self.inertia_tensor_evec
		
		X_x = np.array( [ self.inertia_tensor_evec.item((0,0)), self.inertia_tensor_evec.item((0,1)), self.inertia_tensor_evec.item((0,2)) ] )
		X_y = np.array( [ self.inertia_tensor_evec.item((1,0)), self.inertia_tensor_evec.item((1,1)), self.inertia_tensor_evec.item((1,2)) ] )
		X_z = np.array( [ self.inertia_tensor_evec.item((2,0)), self.inertia_tensor_evec.item((2,1)), self.inertia_tensor_evec.item((2,2)) ] )

		P_x = np.array( [ 0. for i in range(self.number_of_atoms) ] )
		P_y = np.array( [ 0. for i in range(self.number_of_atoms) ] )
		P_z = np.array( [ 0. for i in range(self.number_of_atoms) ] )

		for i in range(self.number_of_atoms):
			R_i = np.array( [ self.main_config_shift.item((i,0)), self.main_config_shift.item((i,1)), self.main_config_shift.item((i,2)) ] )
			P_x.itemset(i, np.dot( R_i, X_x )) 
			P_y.itemset(i, np.dot( R_i, X_y )) 
			P_z.itemset(i, np.dot( R_i, X_z )) 
		
		for i in range(self.number_of_atoms):

			# D4
			self.tf_matrix_ws.itemset((i*3 + 0, 3) , ((P_y.item(i)*X_x.item(2) - P_z.item(i)*X_x.item(1)))*np.sqrt( self.atomic_mass_dic[ self.main_config_type[i] ] ))
			self.tf_matrix_ws.itemset((i*3 + 1, 3) , ((P_y.item(i)*X_y.item(2) - P_z.item(i)*X_y.item(1)))*np.sqrt( self.atomic_mass_dic[ self.main_config_type[i] ] ))
			self.tf_matrix_ws.itemset((i*3 + 2, 3) , ((P_y.item(i)*X_z.item(2) - P_z.item(i)*X_z.item(1)))*np.sqrt( self.atomic_mass_dic[ self.main_config_type[i] ] ))
			# D5
			self.tf_matrix_ws.itemset((i*3 + 0, 4) , ((P_z.item(i)*X_x.item(0) - P_x.item(i)*X_x.item(2)))*np.sqrt( self.atomic_mass_dic[ self.main_config_type[i] ] ))
			self.tf_matrix_ws.itemset((i*3 + 1, 4) , ((P_z.item(i)*X_y.item(0) - P_x.item(i)*X_y.item(2)))*np.sqrt( self.atomic_mass_dic[ self.main_config_type[i] ] ))
			self.tf_matrix_ws.itemset((i*3 + 2, 4) , ((P_z.item(i)*X_z.item(0) - P_x.item(i)*X_z.item(2)))*np.sqrt( self.atomic_mass_dic[ self.main_config_type[i] ] ))
			# D6
			self.tf_matrix_ws.itemset((i*3 + 0, 5) , ((P_x.item(i)*X_x.item(1) - P_y.item(i)*X_x.item(0)))*np.sqrt( self.atomic_mass_dic[ self.main_config_type[i] ] ))
			self.tf_matrix_ws.itemset((i*3 + 1, 5) , ((P_x.item(i)*X_y.item(1) - P_y.item(i)*X_y.item(0)))*np.sqrt( self.atomic_mass_dic[ self.main_config_type[i] ] ))
			self.tf_matrix_ws.itemset((i*3 + 2, 5) , ((P_x.item(i)*X_z.item(1) - P_y.item(i)*X_z.item(0)))*np.sqrt( self.atomic_mass_dic[ self.main_config_type[i] ] ))


		self.D, self.R = np.linalg.qr(self.tf_matrix_ws)	# ORTHONORMALISATION OF THE MATRIX 'D' (FOR TRANSFORMATION)

	def build_hessian(self):

		self.stride = self.number_of_atoms * 3
		
		self.ws_ford = np.matrix( [[ 0. for i in range(3) ] for j in range(self.number_of_atoms) ] )
		self.ws_back = np.matrix( [[ 0. for i in range(3) ] for j in range(self.number_of_atoms) ] )
		self.hessian = np.matrix( [[ 0. for i in range(self.stride) ] for j in range(self.stride) ] )

		for i in range( self.number_of_atoms ):

			file_x = "atom_" + str(i+1) + "_x"		# SET FILE NAMES TO READ
			file_y = "atom_" + str(i+1) + "_y"		#
			file_z = "atom_" + str(i+1) + "_z"		#

			try:										# READ FDM X 
				with open(file_x,'r') as f:						#
					cm = f.readline()						# READ COMMENT LINE
					for j in range(self.number_of_atoms):				# END READ FORWARD  + dx
						rl  = f.readline()					#
						spl = rl.split()					#
						self.ws_ford.itemset((j,0),spl[1])			#
						self.ws_ford.itemset((j,1),spl[2])			#
						self.ws_ford.itemset((j,2),spl[3])			#
					cm = f.readline()						# READ COMMENT LINE
					for j in range(self.number_of_atoms):				# END READ BACKWARD - dx
						rl  = f.readline()					#
						spl = rl.split()					#
						self.ws_back.itemset((j,0),spl[1])			#
						self.ws_back.itemset((j,1),spl[2])			#
						self.ws_back.itemset((j,2),spl[3])			#
				
				res = np.subtract(self.ws_ford,self.ws_back)				# CALCULATE NUMERICAL 2ND DERIVATIVES
				res = res/self.delta							#
				
				for j in range(self.number_of_atoms):							# BUILD HESSIAN MATRIX
					for k in range(3):								#
						self.hessian.itemset((j*3+k,i*3 + self.x_ind), res.item((j,k)))		#

			except FileNotFoundError:							#
				sys.exit(1)								# END X

			try:										# READ FDM Y 
				with open(file_y,'r') as f:						#
					cm = f.readline()						# READ COMMENT LINE
					for j in range(self.number_of_atoms):				# END READ FORWARD  + dy
						rl  = f.readline()					#
						spl = rl.split()					#
						self.ws_ford.itemset((j,0),spl[1])			#
						self.ws_ford.itemset((j,1),spl[2])			#
						self.ws_ford.itemset((j,2),spl[3])			#
					cm = f.readline()						# READ COMMENT LINE
					for j in range(self.number_of_atoms):				# END READ BACKWARD - dy
						rl  = f.readline()					#
						spl = rl.split()					#
						self.ws_back.itemset((j,0),spl[1])			#
						self.ws_back.itemset((j,1),spl[2])			#
						self.ws_back.itemset((j,2),spl[3])			#
				
				res = np.subtract(self.ws_ford,self.ws_back)				# CALCULATE NUMERICAL 2ND DERIVATIVES
				res = res/self.delta							#
				
				for j in range(self.number_of_atoms):							# BUILD HESSIAN MATRIX
					for k in range(3):								#
						self.hessian.itemset((j*3+k,i*3 + self.y_ind), res.item((j,k)))		#

			except FileNotFoundError:							#
				sys.exit(1)								# END Y

			try:										# READ FDM Z 
				with open(file_z,'r') as f:						#
					cm = f.readline()						# READ COMMENT LINE
					for j in range(self.number_of_atoms):				# END READ FORWARD  + dz
						rl  = f.readline()					#
						spl = rl.split()					#
						self.ws_ford.itemset((j,0),spl[1])			#
						self.ws_ford.itemset((j,1),spl[2])			#
						self.ws_ford.itemset((j,2),spl[3])			#
					cm = f.readline()						# READ COMMENT LINE
					for j in range(self.number_of_atoms):				# END READ BACKWARD - dz
						rl  = f.readline()					#
						spl = rl.split()					#
						self.ws_back.itemset((j,0),spl[1])			#
						self.ws_back.itemset((j,1),spl[2])			#
						self.ws_back.itemset((j,2),spl[3])			#
				
				res = np.subtract(self.ws_ford,self.ws_back)				# CALCULATE NUMERICAL 2ND DERIVATIVES
				res = res/self.delta							#
				
				for j in range(self.number_of_atoms):							# BUILD HESSIAN MATRIX
					for k in range(3):								#
						self.hessian.itemset((j*3+k,i*3 + self.z_ind), res.item((j,k)))		#

			except FileNotFoundError:							#
				sys.exit(1)								# END Z

		
		with open("actual_hessian.dat",'w') as f:									# WRITE ACTUAL HESSIAN 
			for i in range(self.stride):										#
				for j in range(self.stride):									#
					symm_checker = (self.hessian.item((i,j)) - self.hessian.item((j,i)))			# CHECK IF HESSIAN SYMMETRIC
					if abs(symm_checker) >= self.symm_checker_tol :
						print("\n")
						print(" @WARNING!!! : Hessian element (%2d , %2d) is not strictly symmetric !" % (i+1,j+1))
						print(" D/TOL       : %12.9f / %12.9f" % (symm_checker,self.symm_checker_tol) )
						print(" Elem1/Elem2 : %12.9f / %12.9f" % (self.hessian.item((i,j)),self.hessian.item((j,i))))
						print(" Err	    : %12.9f" % (abs(symm_checker)/abs((self.hessian.item((j,i))+self.hessian.item((i,j)))/2.)*100.))
						print("\n")
					f.write("%14.6e" % self.hessian.item((i,j)))						#
				f.write("\n")											#

		self.trans = np.transpose(self.hessian)
		self.ave_hessian = np.add(self.trans,self.hessian)
		self.ave_hessian = self.ave_hessian/2.

		with open("ave_hessian.dat",'w') as f:										# WRITE AVE HESSIAN 
			for i in range(self.stride):										#
				for j in range(self.stride):									#
					f.write("%14.6e" % self.ave_hessian.item((i,j)))					#
				f.write("\n")											#

	def build_mwc_hessian(self):												# MASS WEIGHTED HESSIAN
	
		# @ WARNING !!, THIS METHOD HAS TO BE CALLED AFTER 'read_main_config() / build_hessian()' AT LEAST CARRIED OUT	
		# OTHERWISE, NO ATOMIC MASS DATA AVAILBALE !!

		self.mwc_hessian = self.ave_hessian.copy()
		stride = self.number_of_atoms*3

		for i in range(stride):
			for j in range(stride):
				index_1 = i%3
				index_2 = j%3
				numerator_1 = np.sqrt( self.atomic_mass_dic[ self.main_config_type[index_1] ] )
				numerator_2 = np.sqrt( self.atomic_mass_dic[ self.main_config_type[index_2] ] )
				self.mwc_hessian.itemset((i,j), self.mwc_hessian.item((i,j))/numerator_1/numerator_2)

		with open("mwc_hessian.dat",'w') as f:										# WRITE MWC HESSIAN 
			for i in range(stride):											#
				for j in range(stride):										#
					f.write("%14.6e" % self.mwc_hessian.item((i,j)))					#
				f.write("\n")											#


		self.mode, self.vec = np.linalg.eig(self.mwc_hessian)
		neg=[]								# RECORDING NEG FREQ (IM FREQ)
		for i in range(len(self.mode)):					#
			neg.append(self.mode.item(i) > 0)			# TRUE/FALSE TABLE
			if self.mode.item(i) < 0:
				self.mode.itemset(i,self.mode.item(i)*-1.)
		
		self.mode_wn = self.mode.copy()

		'''
		print("%3s" % ("result by simple mwc hessian"))
		for i in range(len(self.mode)):
			if neg[i] == False :
				self.mode_wn.itemset(i, -1.*np.sqrt( self.mode.item(i)/(self.wavenumber_conv) ) )		
			else:
				self.mode_wn.itemset(i, +1.*np.sqrt( self.mode.item(i)/(self.wavenumber_conv) ) )		

			print("%12.6lf%6s" % (self.mode_wn.item(i),"cm-1"))
		'''
	
	def build_int_hessian(self):
	
		stride = self.number_of_atoms*3

		self.D_T = np.transpose(self.D)					# INTERNAL TRANSFORMATION
		int_hessian_ws = np.dot(self.D_T,self.mwc_hessian)		#
		self.int_hessian = np.dot(int_hessian_ws,self.D)		# H_INT = D_T . H_MWC . D

		with open("int_hessian.dat",'w') as f:										# WRITE MWC HESSIAN 
			for i in range(stride):											#
				for j in range(stride):										#
					f.write("%14.6e" % self.int_hessian.item((i,j)))					#
				f.write("\n")											#
		
		self.mode_int, self.vec_int = np.linalg.eig(self.int_hessian)
		self.mode_wn_int = self.mode_int.copy()				# KEEP THE ORG DIAGONALISATION RESULT
		self.neg=[]							# RECORDING NEG FREQ (IM FREQ)
		for i in range(len(self.mode)):					#
			self.neg.append(self.mode_int.item(i) > 0)		# TRUE/FALSE TABLE
			if self.mode_int.item(i) < 0:
				self.mode_int.itemset(i,self.mode_int.item(i)*-1.)
		

		'''
		print("%3s" % ("result by simple int hessian"))
		for i in range(len(self.mode_int)):
			if neg[i] == False :
				self.mode_wn_int.itemset(i, -1.*np.sqrt( self.mode_int.item(i)/(self.wavenumber_conv) ) )		
			else:
				self.mode_wn_int.itemset(i, +1.*np.sqrt( self.mode_int.item(i)/(self.wavenumber_conv) ) )		

			print("%12.6lf%6s" % (self.mode_wn_int.item(i),"cm-1"))
		'''
	def get_vib_direction_mat(self):
		
		stride = self.number_of_atoms*3

		m_matrix = np.array( [[0. for i in range(stride)] for j in range(stride)] )

		for i in range(self.number_of_atoms):
			m_matrix.itemset((i*3+0,i*3+0),1./np.sqrt( self.atomic_mass_dic[ self.main_config_type[i] ] ))
			m_matrix.itemset((i*3+1,i*3+1),1./np.sqrt( self.atomic_mass_dic[ self.main_config_type[i] ] ))
			m_matrix.itemset((i*3+2,i*3+2),1./np.sqrt( self.atomic_mass_dic[ self.main_config_type[i] ] ))
		
		# calculate M * D * L : M (to mwc H), D (to int H), L (to find NMs)

		tmp = np.dot(m_matrix,self.D)
		self.vib_dir_mat = np.dot(tmp,self.vec_int)

		# normalise the raw displacements in cartesian coordinate

		for i in range(self.number_of_atoms*3):
			tmp_sum = 0.
			for j in range(self.number_of_atoms*3):
				tmp_sum += self.vib_dir_mat.item((j,i))**2.		# for 'i'th mode, displacements of atoms in 'j' row elements
			normalisation_factor = 1./np.sqrt(tmp_sum)			# get normalisation factor

			for j in range(self.number_of_atoms*3):
				self.vib_dir_mat.itemset((j,i), self.vib_dir_mat.item((j,i))*normalisation_factor)	# Normalise

	def res_write(self):

		print(" Normal Mode Frequencies ( Translational / Rotational Motions are projected out )")
		print(" Mode     Freq ( cm**-1 ) ... results computed by interal hessian matrix")
		print("--------------------------------------------------------------------------------------")
		for i in range(len(self.mode_int)):
			if self.neg[i] == False :
				self.mode_wn_int.itemset(i, -1.*np.sqrt( self.mode_int.item(i)/(self.wavenumber_conv) ) )		
			else:
				self.mode_wn_int.itemset(i, +1.*np.sqrt( self.mode_int.item(i)/(self.wavenumber_conv) ) )		
			print("%3.2d%17.6lf" % (i+1,self.mode_wn_int.item(i)))
		print("--------------------------------------------------------------------------------------")
		print("")
		print(" Moment of inertia tensor ( u.Angs^2 )")
		print("")
		print("--------------------------------------------------------------------------------------")
		print("	     x       y      z")
		print("--------------------------------------------------------------------------------------")
		print(" %3s%10.3f%10.3f%10.3f" % ("x",self.inertia_tensor_eval.item((0)),0,0))
		print(" %3s%10.3f%10.3f%10.3f" % ("y",0,self.inertia_tensor_eval.item((1)),0))
		print(" %3s%10.3f%10.3f%10.3f" % ("z",0,0,self.inertia_tensor_eval.item((2))))
		print("--------------------------------------------------------------------------------------")
		print(" Centre of mass (Ang) = %10.4f%10.4f%10.4f" % ( self.com.item((0)), self.com.item((1)), self.com.item((2)) ) )
		print("--------------------------------------------------------------------------------------")

		

		try: 
			with open("vib_displ.out","w") as f:
			
				stride = self.number_of_atoms*3

				for i in range(stride):					# FOR THIS NUMBER OF MODES
					f.write("\t%d\n" % (self.number_of_atoms))
					f.write("normal_mode_# %d - Freq(cm^-1): %8.4f\n" % (i+1,self.mode_wn_int.item(i)))

					for j in range(self.number_of_atoms):
						f.write("%3s%10.6f%10.6f%10.6f\n" % (self.main_config_type[j], self.vib_dir_mat.item((j*3+0,i)),self.vib_dir_mat.item((j*3+1,i)),self.vib_dir_mat.item((j*3+2,i))))

					f.write("\n#--\n\n")
					
			'''
			with open("vib_evec.out","w") as f:
				f.write(" # normal mode corresponding eigenvectors (same order with frequencies)\n")
				for i in range(len(self.mode_int)):
					for j in range(len(self.mode_int)):
						f.write("%18.6e" % (self.vec_int.item((i,j))) )
					f.write("\n")
			'''
		except FileNotFoundError:
			print(FileNotFoundError);	sys.exit(1)
	



	def write(self):			# WRITE GENERAL OUTPUT

	# write '.vesta' dipole vector file

		with open("slam_dip.vesta","w") as f:

			f.write("#VESTA_FORMAT_VERSION 3.3.0\n")
			f.write("MOLECULE\n")
			f.write("TITLE\n")
			f.write(" SLAM_DIPOLE\n")
			f.write("STRUC\n")		# STRUCTURE

			cnt=0
			offset_cnt=0
			for i in range(self.mm_cnt):

				if self.mm_config[i][1] == "c":
					cnt += 1
					offset_cnt += 1
					atom_label = self.mm_config[i][0] + str(cnt)
					f.write("%3d%6.3s%10.4s%12s" % (cnt,self.mm_config[i][0],atom_label,"1.0000"))
					f.write("%12.6f%12.6f%12.6f" % (self.mm_config[i][3],self.mm_config[i][4],self.mm_config[i][5]))
					f.write("%12s%12.2s\n" % ("1","-"))
					f.write("%21.6f%12.6f%12.6f%12.2f\n" % (0.,0.,0.,0.))

			cnt=0
			for i in range(self.qm_cnt):
				cnt += 1
				atom_label = self.qm_config[i][0] + str(cnt)
				f.write("%3d%6.3s%10.4s%12s" % (cnt+offset_cnt,self.qm_config[i][0],atom_label,"1.0000"))
				f.write("%12.6f%12.6f%12.6f" % (self.qm_config[i][3],self.qm_config[i][4],self.qm_config[i][5]))
				f.write("%12s%12.2s\n" % ("1","-"))
				f.write("%21.6f%12.6f%12.6f%12.2f\n" % (0.,0.,0.,0.))
		
			f.write("    0 0 0 0 0 0 0\n")	# END FLAG
			# END OF WRITING STRUC		

			# BOND INFO
			f.write("SBOND\n")
			f.write("%3d%6.3s%6.3s%70s\n" % (1,self.mm_config[0][0],self.qm_config[0][0],"0.00000    2.82146  0  1  1  0  1  0.250  2.000 127 127 127"))
			f.write("    0 0 0 0\n")

			# WRITE VECTOR
			f.write("VECTR\n")

			v_cnt=0
			for i in range(self.mm_cnt):

				if self.mm_config[i][1] == "c":
					v_cnt += 1
				if self.mm_config[i][1] == "s":
				#	v_cnt += 1

					rx = self.mm_config[i][3] - self.mm_config[i-1][3]
					ry = self.mm_config[i][4] - self.mm_config[i-1][4]
					rz = self.mm_config[i][5] - self.mm_config[i-1][5]

					ux = rx*self.mm_config[i][2]*self.v_mod
					uy = ry*self.mm_config[i][2]*self.v_mod
					uz = rz*self.mm_config[i][2]*self.v_mod

					ux=ux*-1; uy=uy*-1.; uz=uz*-1.			

					amp = math.sqrt(ux*ux + uy*uy + uz*uz)
					scl = math.exp(-amp/self.v_rho)
					ux = ux*scl
					uy = uy*scl
					uz = uz*scl

					f.write("%3d%15.6f%12.6f%12.6f%12.1d\n" % (v_cnt,ux,uy,uz,0))
					f.write("%3d" % (v_cnt))
					f.write("    0 0 0 0\n")
					f.write(" 0 0 0 0 0\n")

			for i in range(self.qm_cnt):

				v_cnt += 1
				ux = 2.*self.mo_config[i][2]*self.mo_config[i][3]*self.pos_integral*self.qm_config[i][2]*self.v_mod
				uy = 2.*self.mo_config[i][2]*self.mo_config[i][4]*self.pos_integral*self.qm_config[i][2]*self.v_mod
				uz = 2.*self.mo_config[i][2]*self.mo_config[i][5]*self.pos_integral*self.qm_config[i][2]*self.v_mod

				ux=ux*-1; uy=uy*-1.; uz=uz*-1.			

				amp = math.sqrt(ux*ux + uy*uy + uz*uz)
				scl = math.exp(-amp/self.v_rho)
				ux = ux*scl
				uy = uy*scl
				uz = uz*scl

				f.write("%3d%15.6f%12.6f%12.6f%12.1d\n" % (v_cnt,ux,uy,uz,0))
				f.write("%3d" % (v_cnt))
				f.write("    0 0 0 0\n")
				f.write(" 0 0 0 0 0\n")
			f.write("    0 0 0 0 0\n")
			f.write("    0 0 0 0 0\n")
			# END WRITING VECTOR

			f.write("VECTT\n")
			# VECTOR FORMAT
			for i in range(v_cnt):
				f.write("%3d%9.3f%12d%12d%12d%12.3s\n" % (i+1,self.v_radius,self.v_rgb_r,self.v_rgb_g,self.v_rgb_b,self.v_tail_flag))
			f.write(" 0 0 0 0 0\n")

if __name__ == '__main__':

	dig = slam_vib_dig(sys.argv[1],sys.argv[2])
	dig.read_main_config()
	dig.build_inertia_tensor()
	dig.build_tf_matrix()
	dig.build_hessian()
	dig.build_mwc_hessian()
	dig.build_int_hessian()
	dig.get_vib_direction_mat()
	dig.res_write()
	now = datetime.datetime.now()
	#print(" Calculation is finished : ",end="");	print (now.strftime("%Y-%m-%d %H:%M:%S"))
	print(" Calculation is finished : ")
	print(now.strftime(" %Y-%m-%d %H:%M:%S"))
	print(" Terminating Normal Mode Calculator")
	print("--------------------------------------------------------------------------------------")
