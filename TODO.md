

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
* Docstrings
	* EM
	* POSBR
	* IPO
* PyTests
	* Mesh readers & writers
	* farfields
	* spherical to/from cartesian
* V&V Tests
* build system
	* `g++ -fPIC -shared nanortlib.cpp -o nanortlib.so`
* Standardized Output format
	* With multiple frequencies
* MPI support
* Get mutual import paths sorted out
	* nanortlib
	* pytests


#### Modules
* Inputs
	* cmdline
	* input file
	* scene element definitions
	* job recipes
* Mesh - done
* EM
	* omega
	* G, \grad G
	* get farfield
	* plane waves
	* spherical waves
	* doc strings not all complete
* Results
	* write observation arrays
	* read various results files, including Stars and Aurora
	* Handle multiple frequencies in output files
* Geometry - needs complete docstrings
* Logging - done
* MOM
	* Integral kernels
	* Matrix construction
	* RHS construction
* IPO
	* all pairs visibility
	* plane wave visibility
	* sparse solver
	* Iterative solver
	* IPO initial currents
		* remove from EM
* POSBR
	* construct initial rays
	* `const` iterators
	* `const` methods in Ray interfaces
	* PO currents
		* remove from EM


#### Verification
* mie series expressions on wikipedia
* mie python

