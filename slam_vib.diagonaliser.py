#!/bin/python

import numpy as np
import sys
import datetime

kd = lambda i,j : 1. if (i==j) else 0.

class slam_vib_dig:

	def __init__(self,*argv):
		
		self.number_of_atoms = int(argv[0])
		self.delta	     = float(argv[1])
		self.main_config_file= "geo.txt"
		self.main_config     = None
		self.main_config_type= None		

		self.wavenumber_conv=5.8921E-5

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
				self.main_config_type = [ None for i in range(self.number_of_atoms) ]

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

			for j in range(3):
				for k in range(3):
					self.inertia_tensor.itemset((j,k), self.inertia_tensor.item((j,k)) 
						+ m*( rsq*kd(j,k) - self.main_config_shift.item((i,j))*self.main_config_shift.item((i,k)) ))
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
			self.tf_matrix_ws.itemset((i*3 + 0, 3) , ((P_y.item(i)*X_x.item(2) - P_z.item(i)*X_x.item(1)))/np.sqrt( self.atomic_mass_dic[ self.main_config_type[i] ] ))
			self.tf_matrix_ws.itemset((i*3 + 1, 3) , ((P_y.item(i)*X_y.item(2) - P_z.item(i)*X_y.item(1)))/np.sqrt( self.atomic_mass_dic[ self.main_config_type[i] ] ))
			self.tf_matrix_ws.itemset((i*3 + 2, 3) , ((P_y.item(i)*X_z.item(2) - P_z.item(i)*X_z.item(1)))/np.sqrt( self.atomic_mass_dic[ self.main_config_type[i] ] ))
			# D5
			self.tf_matrix_ws.itemset((i*3 + 0, 4) , ((P_z.item(i)*X_x.item(0) - P_x.item(i)*X_x.item(2)))/np.sqrt( self.atomic_mass_dic[ self.main_config_type[i] ] ))
			self.tf_matrix_ws.itemset((i*3 + 1, 4) , ((P_z.item(i)*X_y.item(0) - P_x.item(i)*X_y.item(2)))/np.sqrt( self.atomic_mass_dic[ self.main_config_type[i] ] ))
			self.tf_matrix_ws.itemset((i*3 + 2, 4) , ((P_z.item(i)*X_z.item(0) - P_x.item(i)*X_z.item(2)))/np.sqrt( self.atomic_mass_dic[ self.main_config_type[i] ] ))
			# D6
			self.tf_matrix_ws.itemset((i*3 + 0, 5) , ((P_x.item(i)*X_x.item(1) - P_y.item(i)*X_x.item(0)))/np.sqrt( self.atomic_mass_dic[ self.main_config_type[i] ] ))
			self.tf_matrix_ws.itemset((i*3 + 1, 5) , ((P_x.item(i)*X_y.item(1) - P_y.item(i)*X_y.item(0)))/np.sqrt( self.atomic_mass_dic[ self.main_config_type[i] ] ))
			self.tf_matrix_ws.itemset((i*3 + 2, 5) , ((P_x.item(i)*X_z.item(1) - P_y.item(i)*X_z.item(0)))/np.sqrt( self.atomic_mass_dic[ self.main_config_type[i] ] ))


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
						print(" D/TOL       : %12.9f/%12.9f" % (symm_checker,self.symm_checker_tol) )
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

		try: 
			with open("vib_evec.out","w") as f:
				f.write(" # normal mode corresponding eigenvectors (same order with frequencies)\n")
				for i in range(len(self.mode_int)):
					for j in range(len(self.mode_int)):
						f.write("%18.6e" % (self.vec_int.item((i,j))) )
					f.write("\n")
		except FileNotFoundError:
			print(FileNotFoundError);	sys.exit(1)
	

if __name__ == '__main__':

	dig = slam_vib_dig(sys.argv[1],sys.argv[2])
	dig.read_main_config()
	dig.build_inertia_tensor()
	dig.build_tf_matrix()
	dig.build_hessian()
	dig.build_mwc_hessian()
	dig.build_int_hessian()
	dig.res_write()
	now = datetime.datetime.now()
	print(" Calculation is finished : ",end="");	print (now.strftime("%Y-%m-%d %H:%M:%S"))
	print(" Terminating Normal Mode Calculator")
	print("--------------------------------------------------------------------------------------")
