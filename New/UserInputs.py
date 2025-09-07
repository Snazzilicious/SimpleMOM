

"""Excitations
"""
class PlanewaveDef:
	def __init__(self, dir_x_ang, dir_z_ang, pol_x_ang, magnitude):
		self.dir = Util.math.x_angle_z_angle_to_dir( x_ang, z_ang )
		self.pol = magnitude * ( np.sin(pol_x_ang)*z_angle_hat(dir_x_ang,dir_z_ang) + np.cos(pol_x_ang)*x_angle_hat(dir_x_ang,dir_z_ang) )

class SphericalWaveDef:
	def __init__(self, x,y,z):
		self.center = np.array([x,y,z])

"""Observations
"""
class FarfieldDef:
	def __init__(self, x_ang, z_ang ):
		self.dir = Util.math.x_angle_z_angle_to_dir( x_ang, z_ang )

class MonostaticFarfieldDef:
	def __init__(self, x_ang, z_ang ):
		self.dir = Util.math.x_angle_z_angle_to_dir( x_ang, z_ang )



from argparse import ArgumentParser

def get_standard_arg_parser():
	p = ArgumentParser()
	
	p.add_argument('--version', action='store_true', help='print version and exit')
	
	p.add_argument('--input', help='name of input file')
	
	p.add_argument('--output', help='name for output files')
	
	return p


import xmltodict

def read_input_file( fname ):
	f = open( fname, "rb" )
	content = xmltodict.parse( f )
	f.close()
	
	# get global excitations and observations
	
	# for each scene
		# get mesh, frequency, casename
		# for each excitation
			# get
			# for each observation
				# get

	return content
