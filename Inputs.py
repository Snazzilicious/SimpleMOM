"""Module for acquiring inputs from the user.
"""

from argparse import ArgumentParser

"""Interprets metric/si prefixes specified on units. e.g. '1.2G' -> 1.2e9
"""
def parse_si_prefix(s):
	"""Parses a float from a string that may contain a magnitude prefix, e.g. '5.5G' parses to 5.5e9.
	"""
	prefixes = {"G":1e9, "M":1e6, "K":1e3, "k":1e3, "m":1e-3, "u":1e-6, "n":1e-9}
	if s[-1] in prefixes:
		return float(s[:-1]) * prefixes[s[-1]]
	else:
		return float(s)
		

def get_standard_arg_parser():
	"""Returns an argparse.ArgumentParser with arguments common to all applications.
	"""
	p = ArgumentParser()
	
	p.add_argument('--version', action='store_true', help='print version and exit')
	
	p.add_argument('--mesh', help='name of mesh file')
	
	p.add_argument('--frequency', type=parse_si_prefix, help='excitation frequency in Hz')
	
	p.add_argument('--output', help='name for output files')
	
	return p



import xmltodict

def read_input_file( fname ):
	"""Parses input file.
	
	Returns
		content : dict
			Content of the input file organized by xmltodict.
			Keywords are keys and input values are the mapped values.
	"""
	f = open( fname, "rb" )
	content = xmltodict.parse( f )
	f.close()
	
	return content
