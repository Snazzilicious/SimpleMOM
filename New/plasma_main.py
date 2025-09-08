
import numpy as np

# load mesh

# get connectivity graph; Interface & vertex
# get cell volumes [ nCells ]
# get interface areas [ nInterfacesPerCell ]

# set initial values [ nCells, nQuants ]

for it in range(numIts):
	
	# compute min mod slopes* [ nCells, nQuants, nDim ]
	# compute all inner face values [ nCells, nQuants, nInterfacesPerCell ]
	# apply boundary conditions; in/outflow, solid boundary
	
	# compute wave speeds* [ nCells, nQuants, nInterfacesPerCell ]
	
	# compute fluxes* [ nCells, nQuants, nInterfacesPerCell ]
	
	# sum up fluxes [ nCells, nQuants ]
	
	# step

# *involves scatter/gather (interprocess communication)
"""
Is there an ESP module that does *simple* cfd?
"""
