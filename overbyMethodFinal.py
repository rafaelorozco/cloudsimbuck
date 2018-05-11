from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph.opengl as gl
import pyqtgraph as pg
import numpy as np
import random
import time
import functools
import noise
import sys

app = QtGui.QApplication([])
w = gl.GLViewWidget()
w.show()
g = gl.GLGridItem()
w.addItem(g)


def add_force(grid_0, grid_force, dt):
	"""from stam section 2.2:
	grid = grid + dt*force
    """

	row, col, lay = grid_0.shape
	grid_0[0:row, 0:col, 0:lay] += dt*grid_force[0:row, 0:col, 0:lay] 

def dissipate(grid_1, grid_0, a, dt):
	"""a is disipation rate
	from stam sec 2.4 and 3.2:
	divide each element in the first array and store the result
	in the new array.
	for element in grid_1:
		element = element / (1+dt*a)
    """

	row, col, lay = grid_0.shape
	grid_1[0:row, 0:col, 0:lay] = grid_0[0:row, 0:col, 0:lay] / (1+a*dt)

def transport(grid_1, grid_0, grid_U, dt):
	"""transport scalars through  field
	Use runge-kuta 2nd order for particle trace
	"""
	d = 1 #size of voxels
	origin = 0
	rows, cols, lays = grid_0.shape
	for i in range(0,rows):	
		for j in range(0,cols):
			for k in range(0,lays):
				#point  = ndarray([origin + (i+0.5)*d,(j+0.5)*d,(k+0.5)*d] #get point in middle of grid cell
				point  = np.array([(i)*d,(j)*d,(k)*d]) #get point
				point_before = trace_particle(point, grid_U, -dt) #find where this point came from
				grid_1[i][j][k] = get_value(point_before, grid_0) #put the source there

	set_bound(grid_1,0)

def get_value(point, grid_U):
	"""Use bilinear interpolation to get value 
	of grid at point. (suppose that grid cells have vol = 1) 

	tested separately (testget_value.py)
	"""
	i_hat = int(np.floor(point[0]))
	j_hat = int(np.floor(point[1]))
	k_hat = int(np.floor(point[2]))

	#relative to grid point start. 
	x_max,y_max,z_max = grid_U.shape
	point_x = point[0] - i_hat
	point_y = point[1] - j_hat
	point_z = point[2] - k_hat

	#edge cases
	if(i_hat+1 > x_max-1 or i_hat < 0):
		return 0
	if(j_hat+1 > y_max-1 or j_hat < 0):
		return 0
	if(k_hat+1 > z_max-1 or k_hat < 0):
		return 0

	#find volume of each section
	vol_a = grid_U[i_hat][j_hat][k_hat]      *(1-point_x)*(1-point_y)*(1-point_z)
	vol_b = grid_U[i_hat+1][j_hat][k_hat]    *(point_x)*(1-point_y)*(1-point_z)
	vol_c = grid_U[i_hat][j_hat+1][k_hat]    *(1-point_x)*(point_y)*(1-point_z)
	vol_d = grid_U[i_hat][j_hat][k_hat+1]    *(1-point_x)*(1-point_y)*(point_z)
	vol_e = grid_U[i_hat+1][j_hat+1][k_hat]  *(point_x)*(point_y)*(1-point_z)
	vol_f = grid_U[i_hat+1][j_hat][k_hat+1]  *(point_x)*(1-point_y)*(point_z)
	vol_g = grid_U[i_hat][j_hat+1][k_hat+1]  *(1-point_x)*(point_y)*(point_z)
	vol_h = grid_U[i_hat+1][j_hat+1][k_hat+1]*(point_x)*(point_y)*(point_z)

	return vol_a + vol_b + vol_c + vol_d + vol_e + vol_f + vol_g +  vol_h

def get_value_vel(point, grid_U):
	"""Use bilinear interpolation to get value 
	of grid at point. (suppose that grid cells have vol = 1) 

	tested separately (testget_value.py)

	This version of get value takes in to consideration a 0.5 vol around the point to help propagates sources further
	than the actual velocity field
	"""
	p_a = point + np.array([(-0.5),(-0.5),(-0.5)])
	p_b = point + np.array([(-0.5),(0.5),(0.5)])
	p_c = point + np.array([(0.5),(-0.5),(0.5)])
	p_d = point + np.array([(0.5),(0.5),(-0.5)])
	p_e = point + np.array([(-0.5),(-0.5),(0.5)])
	p_f = point + np.array([(0.5),(-0.5),(-0.5)])
	p_g = point + np.array([(-0.5),(0.5),(-0.5)])
	p_h = point + np.array([(0.5),(0.5),(0.5)])

	return (get_value(p_a, grid_U)+ 
			get_value(p_b, grid_U)+ 
			get_value(p_c, grid_U)+ 
			get_value(p_d, grid_U)+
			get_value(p_e, grid_U)+
			get_value(p_f, grid_U)+
			get_value(p_g, grid_U)+
			get_value(p_h, grid_U))/8

def trace_particle(point, grid_U, dt):
	"""trace particle through given vector field
	Use runge-kuta 2nd order
	grid_U is components of the velocity field
	if dt is negative it will trace backwards in time

	tested separately (testtrace.py)
	"""
	#k1 = dt * np.array([get_value(point, grid_U[0]), get_value(point, grid_U[1]), get_value(point, grid_U[2])])
	#k2 = dt * np.array([get_value(point+k1, grid_U[0]), get_value(point+k1, grid_U[1]), get_value(point+k1, grid_U[2])])

	#return point + 0.5*(k1+k2)


	k1 = point + dt * np.array([get_value_vel(point, grid_U[0]), get_value_vel(point, grid_U[1]),get_value_vel(point, grid_U[2])])

	return point + dt*0.5*(np.array([get_value_vel(point, grid_U[0]), get_value_vel(point, grid_U[1]), get_value_vel(point, grid_U[2])])
		                            +np.array([get_value_vel(k1, grid_U[0]), get_value_vel(k1, grid_U[1]),get_value_vel(k1, grid_U[2])]))

def project(grid_1, grid_0, dt):
	"""Use poisson method to project field to its zero-divergent component
	"""
	rows, cols, lays = grid_0[0].shape

	#First get divergence of vel field
	grid_div = np.zeros((rows, cols, lays))

	matR1 = np.roll(grid_0[0],-1,axis=0) # mat[i+1,j,k]
	matR2 = np.roll(grid_0[0],1,axis=0) # mat[i-1,j,k]
	matR3 = np.roll(grid_0[1],-1,axis=1) # mat[i,j+1,k]
	matR4 = np.roll(grid_0[1],1,axis=1) # mat[i,j-1,k]
	matR5 = np.roll(grid_0[2],-1,axis=2) # mat[i,j,k+1]
	matR6 = np.roll(grid_0[2],1,axis=2) # mat[i,j,k-1]

	h = 1/(rows*1)
	grid_div[0:rows][0:cols][0:lays] = -0.5*( matR1[0:rows][0:cols][0:lays] - matR2[0:rows][0:cols][0:lays] +
													matR3[0:rows][0:cols][0:lays] - matR4[0:rows][0:cols][0:lays] +
													matR5[0:rows][0:cols][0:lays] - matR6[0:rows][0:cols][0:lays])
	set_bound(grid_div,0)

	#use Gauss-Seidel relaxation to solve linear equations
	grid_sol = np.zeros((rows, cols, lays)) #solution grid
	for n in range(20): #iterate many times or until converges
		matR1 = np.roll(grid_sol,-1,axis=0) # mat[i+1,j,k]
		matR2 = np.roll(grid_sol,1,axis=0) # mat[i-1,j,k]
		matR3 = np.roll(grid_sol,-1,axis=1) # mat[i,j+1,k]
		matR4 = np.roll(grid_sol,1,axis=1) # mat[i,j-1,k]
		matR5 = np.roll(grid_sol,-1,axis=2) # mat[i,j,k+1]
		matR6 = np.roll(grid_sol,1,axis=2) # mat[i,j,k-1]

		grid_sol[0:rows][0:cols][0:lays] = (grid_div[0:rows][0:cols][0:lays] + 
													    matR1[0:rows][0:cols][0:lays] + matR2[0:rows][0:cols][0:lays] +
														matR3[0:rows][0:cols][0:lays] + matR4[0:rows][0:cols][0:lays] +
														matR5[0:rows][0:cols][0:lays] + matR6[0:rows][0:cols][0:lays] ) / 6
		set_bound(grid_sol,0)


	#Set new grid to the old grid minus the solution (sol grid)
	#this makes the the new grid divergence free
	matR1 = np.roll(grid_sol,-1,axis=0) # mat[i+1,j,k]
	matR2 = np.roll(grid_sol,1,axis=0) # mat[i-1,j,k]
	matR3 = np.roll(grid_sol,-1,axis=1) # mat[i,j+1,k]
	matR4 = np.roll(grid_sol,1,axis=1) # mat[i,j-1,k]
	matR5 = np.roll(grid_sol,-1,axis=2) # mat[i,j,k+1]
	matR6 = np.roll(grid_sol,1,axis=2) # mat[i,j,k-1]

	grid_1[0][0:rows][0:cols][0:lays] = grid_0[0][0:rows][0:cols][0:lays] - \
											  0.5*(matR1[0:rows][0:cols][0:lays] - matR2[0:rows][0:cols][0:lays])
	grid_1[1][0:rows][0:cols][0:lays] = grid_0[1][0:rows][0:cols][0:lays] - \
											  0.5*(matR3[0:rows][0:cols][0:lays] - matR4[0:rows][0:cols][0:lays])
	grid_1[2][0:rows][0:cols][0:lays] = grid_0[2][0:rows][0:cols][0:lays] - \
											  0.5*(matR5[0:rows][0:cols][0:lays] - matR6[0:rows][0:cols][0:lays])
													


def v_step(grid_u1, grid_u0, grid_forces, dt):
	"""Execute one velocity step
	"""
	NDIM = 3 #three dimensions so three components of velocity

	#add velocity due to forces present
	for i in range(NDIM):
		add_force(grid_u0[i], grid_forces[i], dt)

	#Transport the velocity field on itself
	for i in range(NDIM):
		transport(grid_u1[i], grid_u0[i], grid_u0, dt)

	#for i in range(NDIM): dont have to do because we can ignore viscosity?
	#	diffuse()

	#project the velocity to keep with divergence free
	project(grid_u0, grid_u1, dt)

def s_step(grid_s1, grid_s0, grid_u,a, dt):
	"""Execute one source step. 
    """

    #transport sources through velocity field
	transport(grid_s1, grid_s0, grid_u, dt)

	#diffuse() dont have to do because we can ignore viscosity?

	#dissipate the sources
	dissipate(grid_s0, grid_s1, a, dt)
	dissipate(grid_s1, grid_s0, a, dt)



def set_bound(grid_0, C):
	"""Set bounds of a 3d grid to C
    """
	row, col, lay = grid_0.shape
	
	grid_0[0, 0:col, 0:lay]     = C
	grid_0[row-1, 0:col, 0:lay] = C

	grid_0[0:row, 0, 0:lay]     = C
	grid_0[0:row, col-1, 0:lay] = C

	grid_0[0:row, 0:col, 0]     = C
	grid_0[0:row, 0:col, lay-1] = C


###__EXTERNAL METHODS_###
def perlin3d(x,y,z):
	"""make 3d perlin noise in a space specified by parameters. 
    """
	
	scale = 10
	height = x
	width = y
	z_axis = z

	point_size = []

	for i in range(height-1):
		for j in range(width-1):
			for k in range(z_axis-1):
				icoord = (float(i) / height) * scale 
				jcoord = (float(j) / width) * scale 
				kcoord = (float(k) / z_axis) * scale 

				point_size.append((noise.pnoise3(icoord,jcoord,kcoord)+1)**2)

	return np.reshape(point_size ,(height-1,width-1,z_axis-1))

def normalize(x):
	"""Normalize array from [0,1]
    """
	x = np.asarray(x)
	return (x - x.min()) / (np.ptp(x))

def normalizeT(x,t):
	"""Normalize array on certain number t
    """

	x = np.asarray(x)
	return (x - x.min()) / (t*np.ptp(x))

###-____OVERBY METHODS___###

def make_heat_sourceA(grid_0,num_sources,size_factor):
	""""first we define heat sources at the ground level" pg 38
    """
	x, y, z= grid_0.shape
	z_max = z // 5
	
	##Put them in random places in the bottom 5th of the space. 
	for i in range(num_sources):
		pos_x = 3#random.randint(1,x)
		pos_y = 3#random.randint(1,y)
		pos_z = 3#:random.randint(1,z_max)
		grid_0[pos_x:pos_x+size_factor, pos_y:pos_y+size_factor, pos_z:pos_z+size_factor] = 50

def make_heat_sourceCenter(grid_0,num_sources,size_factor):
	"""Define heat source in center of grid but in the bottom
    """
	row, col, lay= grid_0.shape

	size_factor = size_factor // 2
	
	pos_x = row // 2
	pos_y = col // 2
	pos_z = size_factor+1

	grid_0[pos_x-size_factor:pos_x+size_factor, pos_y-size_factor:pos_y+size_factor, 0:pos_z+size_factor] = 50


def make_hum_sourceA(grid_0):
	""""we use a perlin noise function scaled by an overall concentration value to define both of our humidty sources" pg 44
    """
	row, col, lay = grid_0.shape

	#bottom third
	lay = (lay // 3)*1

	#bound = 1
	grid_0[0:row, 0:col, 0:lay] = perlin3d(row+1, col+1, lay+1)[0:row, 0:col, 0:lay]

	#normalize by factor t
	grid_0 = normalizeT(grid_0,t=100)
	#print(np.amin(grid_0))
	#print(np.amax(grid_0))



def get_temp(grid_0):
	"""Generate temperature of space.
    """
	row, col, max_z = grid_0.shape

	#We are able to compute the temperature at an point in the field by 
	#multiplying the amont of heat energy present by the pressure 
	#T_local = Heat_energy * Pressure
	for z_i in range(1,max_z):
		grid_0[0:row, 0:col, z_i] = grid_0[0:row, 0:col, z_i] * get_pressure(z_i, max_z)

def get_pressure(curr_z, max_z):
	"""Get pressure at current altitude (curr_z)
	Pressure is defined as a gradient based on altitude (pg 39)
	The presure is 1.0 bar at ground level and descreases exponentially with height to 0.1 Bar 
	at the tropopause. 
    """

    #Pressure depends on total height.
	a = np.log(0.1) / max_z
	return np.exp(curr_z*a*4)


def gen_bouyant_forces(grid_force, grid_temp):
	"""We use temperatre use that to calculate bouyant forces. 
	At each point, if the surrounding temperature is less than 
	the current temperature at that point then we add a possitive bouyant force
	AND VICE VERSA.
    """

	row, col, lay = grid_force.shape

	#constant k 
	K_bouyancy = 1.0
	
	#Find local temperature by averaging the temperature values of all neighbor voxels in the same horizontal plane.
	#all neighbors (8 of them):
	matR1 = np.roll(grid_temp,-1,axis=0) # mat[i+1,j,k]
	matR2 = np.roll(grid_temp,1,axis=0) # mat[i-1,j,k]
	matR3 = np.roll(grid_temp,-1,axis=1) # mat[i,j+1,k]
	matR4 = np.roll(grid_temp,1,axis=1) # mat[i,j-1,k]

	matR5 = np.roll(np.roll(grid_temp,1,axis=0),1,axis=1) # mat[i-1,j-1,k]
	matR6 = np.roll(np.roll(grid_temp,-1,axis=0),-1,axis=1) # mat[i+1,j+1,k] 
	matR7 = np.roll(np.roll(grid_temp,1,axis=0),-1,axis=1) # mat[i-1,j+1,k]
	matR8 = np.roll(np.roll(grid_temp,-1,axis=0),1,axis=1) # mmat[i+1,j-1,k] 

	#F_bouyancy = K_bouyancy(T_local - T_surrounding)
	grid_force[0:row, 0:col, 0:lay] = K_bouyancy * (grid_temp[0:row, 0:col, 0:lay] 
		                                                                      	 - (   matR1[0:row, 0:col, 0:lay] 
																					+  matR2[0:row, 0:col, 0:lay] 
															     					+  matR3[0:row, 0:col, 0:lay] 
																					+  matR4[0:row, 0:col, 0:lay] 
																					+  matR5[0:row, 0:col, 0:lay] 
																					+  matR6[0:row, 0:col, 0:lay] 
																					+  matR7[0:row, 0:col, 0:lay] 
																					+  matR8[0:row, 0:col, 0:lay])/8)



def gen_rel_humidity(grid_curr_hum, grid_temp):
	"""Generate relative humidity at each voxel.
    RelativeHumidity = Current_wv / Total_wv (wv = water vapor)
	To get total wv multiply current temp and press and scale [0,1]
    """

	#grid_source refers to current humidity source
	row, col, lay = grid_curr_hum.shape
	grid_rel_hum = np.zeros((row, col, lay))
	
	grid_norm_temp = np.copy(grid_temp)
	grid_norm_temp = normalize(grid_norm_temp)

	grid_norm_pres = np.zeros((row, col, lay))
	for z_i in range(0,lay):
		grid_norm_pres[0:row, 0:col, z_i] = get_pressure(z_i, lay)
	grid_norm_pres = normalize(grid_norm_pres)

	#print(grid_norm_pres)
	#print(grid_norm_temp)

	#print(np.amax(grid_norm_pres))
	#print(np.amax(grid_norm_temp))

	grid_total_hum = np.zeros((row, col, lay))
	grid_total_hum[0:row, 0:col, 0:lay] =  \
									grid_norm_temp[0:row, 0:col, 0:lay] * grid_norm_pres[0:row, 0:col, 0:lay]

	#print(np.amax(grid_norm_pres))
	#print(np.amax(grid_norm_temp))
	#print("ew:" +str(np.amax(grid_total_hum)))

	grid_total_hum = (grid_total_hum +0.0001)

	#print(np.amax(grid_curr_hum))
	#print("ew:" +str(np.amax(grid_total_hum)))
	
	grid_rel_hum[0:row, 0:col, 0:lay] = \
	grid_curr_hum[0:row, 0:col, 0:lay] / (grid_total_hum[0:row, 0:col, 0:lay])


	return grid_rel_hum

def gen_condensation(grid_cond,grid_rel_hum):
	"""Generate condensation.
	The relative humidity at a given voxel must be at or near saturation for condensation to 
    occur. If concentration of hygroscopic nuclei is high, condensation can occur at lower relative 
    humidity levels. 
    Water_condensation = (RH / 100) * HN * K_c
    """

	#grid_source refers to relative humidity source
	K_c = 1
	Hydro_Nuc = 1
	row, col, lay = grid_cond.shape
	grid_cond[0:row, 0:col, 0:lay] =  K_c*Hydro_Nuc*(1/100)*grid_rel_hum[0:row, 0:col, 0:lay] 
	
	return grid_cond

def gen_lat_heat(grid_source):
	"""Generate latent heat due to condensation.
	This is heat released from the phase transition of water vapor. 
	H_latent = W_condensed * (2.5 * 10**6) * (K_l)
    """

	K_l = 1
	TEMP_K = 2.5*(10^6)
	row, col, lay = grid_source.shape
	grid_lat_heat = np.zeros((row, col, max_z))
	grid_lat_heat[0:row, 0:col, 0:lay] =  K_l*TEMP_K*grid_source[0:row, 0:col, 0:lay] 
	
	return grid_lat_heat

####----SIM CODE

def sim_step():
	"""Perform one step of simulation
    """
	global grid_u0, grid_u1, grid_heat0, grid_heat1, grid_hum0, grid_hum1, grid_cond0, grid_cond1, grid_forces, grid_sources_heat, grid_sources_hum
	a = 0.01 
	dt = 0.1

	##SWAP SOURCE GRIDS
	grid_temp = np.copy(grid_hum1)
	grid_hum1 = np.copy(grid_hum0)
	grid_hum0 = np.copy(grid_temp)

	grid_temp = np.copy(grid_heat1)
	grid_heat1 = np.copy(grid_heat0)
	grid_heat0 = np.copy(grid_temp)

	grid_temp = np.copy(grid_cond1)
	grid_cond1 = np.copy(grid_cond0)
	grid_cond0 = np.copy(grid_temp)

	##SWAP VELOCITY GRIDS
	grid_temp = np.copy(grid_u1)
	grid_u1 = np.copy(grid_u0)
	grid_u0 = np.copy(grid_temp)

	#add density sources to grid
	add_force(grid_hum0, grid_sources_hum, dt*0.1)
	add_force(grid_heat0, grid_sources_heat, dt)

	#add forces due to these sources
	gen_bouyant_forces(grid_forces[2], grid_heat0)

	#Step for velocity field
	v_step(grid_u1, grid_u0, grid_forces, dt)

	#Steps for scalar sources
	s_step(grid_hum1, grid_hum0, grid_u0, a, dt)
	s_step(grid_heat1, grid_heat0, grid_u0, a, dt)
	s_step(grid_cond1, grid_cond0, grid_u0, a, dt)

	#cap humidity at 1. Humidity must be a range from [0,1]
	grid_hum1[grid_hum1 > 1.0] = 1.0

	#generate relative humidity and condensation 
	grid_rel_hum = gen_rel_humidity(grid_hum1, grid_heat1)
	gen_condensation(grid_cond1,grid_rel_hum)

	
####_____INITIALIZATION___###
#total space size
lenI = 12
lenJ = 12
lenK = 32

#Initialize forces matrix
grid_forces = np.zeros((3,lenI, lenJ, lenK))

#Initialize velocity matrices
grid_u0 = np.zeros((3,lenI, lenJ, lenK))
grid_u1 = np.zeros((3,lenI, lenJ, lenK))

#Initialize sources Three: humidity concensation and heat
grid_sources_empty = np.zeros((lenI, lenJ, lenK))
grid_sources_heat = np.zeros((lenI, lenJ, lenK))
grid_sources_hum = np.zeros((lenI, lenJ, lenK))
#grid_sources_cond = np.zeros((lenI, lenJ, lenK))

#grid_h_test = np.ones((lenI, lenJ, lenK))

#define a battery of tests.
def tests():
	

	##TEST HUMID SOURCE
	#make_hum_sourceA(grid_sources)

	##TEST HEAT SOURCE
	#make_heat_sourceA(grid_sources,3,5)

	##TEST MAKE PRESSURE
	#get_temp(grid_sources)

	gen_rel_humidity(grid_curr_hum, grid_temp)

	##TEST MAKE BOUYANCY FORCE
	# make_heat_sourceA(grid_sources,3,5)

	# grid_forces[0].fill(0.1)
	# gen_bouyant_forces(grid_forces[0], grid_sources)

	##TEST RELATIVE HUMIDITY
	# make_heat_sourceA(grid_sources_temp,3,10)
	# make_hum_sourceA(grid_sources_hum)
	# grid_rel_hum = gen_rel_humidity(grid_sources_hum, grid_sources_temp)

	##TEST CONDENSATION
	# make_heat_sourceA(grid_sources_temp,3,10)
	# make_hum_sourceA(grid_sources_hum)
	# grid_rel_hum = gen_rel_humidity(grid_sources_hum, grid_sources_temp)

	# grid_cond = gen_condensation(grid_rel_hum)

	##TEST LATENT HEAT
	# make_heat_sourceA(grid_sources_temp,3,10)
	# make_hum_sourceA(grid_sources_hum)
	# grid_rel_hum = gen_rel_humidity(grid_sources_hum, grid_sources_temp)
	# grid_cond = gen_condensation(grid_rel_hum)
	# grid_lat_heat = gen_lat_heat(grid_cond)

	## test grid swap:
	# grid_a = np.ones((2, 2, 2))
	# grid_b = np.ones((2, 2, 2))
	# grid_a.fill(0.1)
	# print(grid_a)
	# print(grid_b)
	# swap_grid(grid_a, grid_b)
	# print(grid_a)
	# print(grid_b)


#define a humidity source
make_hum_sourceA(grid_sources_hum)

#generate a heat source
make_heat_sourceCenter(grid_sources_heat,1,6)

#Initialize source grids
grid_hum0 = np.copy(grid_sources_empty)
grid_hum1 = np.copy(grid_sources_empty)
grid_heat0 = np.copy(grid_sources_empty)
grid_heat1 = np.copy(grid_sources_empty)
grid_cond0 = np.copy(grid_sources_empty)
grid_cond1 = np.copy(grid_sources_empty)

#___VISUALIZATION__#

#initialize points. Just product a matrix with coords of each point in 3d space
gridAll = np.ones((lenI, lenJ, lenK))
indexes = np.where(gridAll == 1)
indexesAll = np.array([[indexes[0][i],indexes[1][i],indexes[2][i]] for i in range(len(indexes[0]))])
		
grid_s1 = grid_cond1
sp2 = gl.GLScatterPlotItem(pos=indexesAll,pxMode=False)
w.addItem(sp2)
sp2.setData(pos=indexesAll, size=grid_s1.flatten())

# sp5 = gl.GLScatterPlotItem(pos=indexesAll,pxMode=False)
# w.addItem(sp5)
# sp5.setData(pos=indexesAll, size=grid_s1.flatten())

#plot positive bouyancy
sp3 = gl.GLScatterPlotItem(pos=indexesAll,pxMode=False,color=((255, 0, 0,1)))
w.addItem(sp3)
sp3.setData(pos=indexesAll, size=grid_u1[2].flatten())

#plot negative bouyancy
# sp4 = gl.GLScatterPlotItem(pos=indexesAll,pxMode=False,color=((0, 0, 255,1)))
# w.addItem(sp4)
# sp4.setData(pos=indexesAll, size=grid_u1[2].flatten())

#manually set position in 3d space.
w.setCameraPosition(0.12,100.12,6)


#set total number of iterations to run simulation.
totalIterations = 20
numIteration    = 0
def update():
	"""Create a new user.
    Line 2 of comment...
    And so on... 
    """

	
	global numIteration, indexesFinal
	if(numIteration < totalIterations) :

		grid_source = grid_cond1

		grid_source[grid_source > 0.9] = 0.9
		flat_cond = grid_source.flatten()
		color_flat = np.array([(255,255,255,i) for i in flat_cond])
		sp2.setData(pos=indexesAll, size=(grid_source).flatten())


		grid_hum1[grid_hum1 > 1.5] = 1.5
		flat_cond = grid_hum1.flatten()
		#color_flat = np.array([(0,0,255,i) for i in flat_cond])
		#print(color_flat.shape)
		sp3.setData(pos=indexesAll, size=(grid_hum1).flatten())

		# grid_heat1[grid_heat1 > 1.5] = 1.5
		# sp5.setData(pos=indexesAll, size=(grid_heat1).flatten())

		# grid_u1[2][grid_u1[2] > 1.5] = 1.5
		# sp3.setData(pos=indexesAll, size=(grid_forces[2]).flatten())

		# grid_cond1[grid_cond1 > 1.5] = 1.5
		# sp3.setData(pos=indexesAll, size=(grid_cond1).flatten())
		#print(w.cameraPosition())

		#grid_forces[2][grid_forces[2] > 1.5] = 1.5
		#grid_forces[2] *= 10
		#sp3.setData(pos=indexesAll, size=(grid_forces[2]).flatten())

		# grid_neg = np.copy(grid_u1[2])
		# grid_neg = grid_neg*(-1)
		# grid_neg[grid_neg > 1.5] = 1.5
		# sp4.setData(pos=indexesAll, size=(grid_neg).flatten())

		#indexes = np.where(grid_sources > 0)
		#indexesFinal = np.array([[indexes[0][i],indexes[1][i],indexes[2][i]] for i in range(len(indexes[0]))])
		print("iter: " +str(numIteration))
		sim_step()
		numIteration+=1
		

		#After each simulation make a window image grab of the current view.
		#pause program to allow for window grab.
		#time.sleep(5.5)
		w.grabFrameBuffer().save("vid/cond_test_"+str(numIteration).zfill(4)+".png")
	else:
		#End simulation after n simulation steps.
		pass
		

t = QtCore.QTimer()
t.timeout.connect(update)
t.start(1000)  

## Start Qt event loop unless running in interactive mode.
if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, PYQT_VERSION):
        QtGui.QApplication.instance().exec_()