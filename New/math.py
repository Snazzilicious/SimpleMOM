

def si_prefix_factor( unit_str ):
	prefixes = {"G":1e9, "M":1e6, "K":1e3, "k":1e3, "m":1e-3, "u":1e-6, "n":1e-9}
	
	return prefixes.get( unit_str, default=1.0 )


"""Spherical coordinates convenience functions
"""
def x_angle_z_angle_to_dir( x_angle, z_angle ):
	return np.array([ np.cos(x_angle)*np.sin(z_angle), np.sin(x_angle)*np.sin(z_angle), np.cos(z_angle) ])

def z_angle_hat( x_angle, z_angle ):
	return -x_angle_z_angle_to_dir( x_angle+np.pi, (np.pi/2) - z_angle )

def x_angle_hat( x_angle, z_angle ):
	return x_angle_z_angle_to_dir( x_angle+np.pi/2, np.pi/2 )
