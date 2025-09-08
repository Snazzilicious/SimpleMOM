

* Input file
	* xml format like aurora
		* scene w/ mesh, units, frequency
		* excitation
			* plane
				* polarization
			* sphere
				* polarization and amplitude
			* multi
		* Observation angles
	* Need SceneDefinitions module
		* need excitation types
		* specify observations
		* frequency parsing
		* angle sweeps
* Revise input system
	* cmdline > Input file
	* scenes > contained excitations > global excitations
* Remove 'Util' from module names
* C ray trace
	* use ctypes
	* Tracer type for python
	* trace function
		* all combinations of org/dir layouts
	* V&V
	* build system
* Properly organize modules for import
	* set up `__init__.py`
* Docstrings
	* Mesh - done
	* Results - done
	* EM
	* General
	* IPO
* PyTests
* Standardized Output format
	* With multiple frequencies
* MPI support


#### Modules
* Inputs
	* cmdline
	* input file
	* scene element definitions
	* job recipes
* Mesh - done
	* read mesh files -> vertices and faces
	* connectivity
	* basis, centroids
	* write visualization files (vtk)
* EM
	* omega
	* G, \grad G
	* get farfield
	* plane waves
	* spherical waves
* Results
	* write observation arrays
	* read various results files, including Stars and Aurora
	* Handle multiple frequencies in output files
* Geometry
	* spherical to/from cartesian
* Util
	* logging
* MOM
	* Integral kernels
	* Matrix construction
	* RHS construction
* IPO
	* all pairs visibility
	* plane wave visibility
	* sparse solver
* POSBR
	* construct initial rays
	* `const` iterators
	* `const` methods in Ray interfaces
	* Revise or remove concrete Ray types


#### Verification
* mie series expressions on wikipedia
* mie python

