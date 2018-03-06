def add_force(grid_i, force_i, dt):
	#from stam:
	new_grid_i = grid_i + dt*force_i
	return new_grid_i

#a is disipation rate
def dissipate(grid_0, grid_1, a, dt)
	#from stam sec 2.4:
