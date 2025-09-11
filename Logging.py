"""Routines for printing to screen or log files.
"""

from time import ctime

def stdout(msg):
	"""Prints a message to console with a time stamp. Immediately flushes the standard out buffer.
	
	Arguments
		msg : string
			Message to print to console.
	"""
	print(ctime(), msg, flush=True)


