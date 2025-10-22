

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
	* `g -fPIC -shared nanortlib.cpp -o nanortlib.so`
* Standardized Output format
	* With multiple frequencies
* MPI support
	* Distributed Objects and RPCs
		* Pyro and Ray
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
* ReWrite Stars
	* Matrix classes
		* Dense & Low Rank & Zero
			* need to be careful about which operators get called
			* need to map out result type of each binary operation depending on types of inputs
				* add, mult, matmult
				* scalar, Dense, LR, Zero, H
		* Operations
			* LU factor
			* Tri solve
				* L and U
				* Left and Right
			* GEMM
			* Transpose
			* ACA fill
		* Use sqlite for meta data
			* give this a tree-like interface
			* ID, nrows, ncols, rank, start_row, start_col
		* Block matrix
		* OOC Matrix
			* replacements for DenseMatrix and LowRankMatrix
		* Distributed Matrix
	* Error throwing
	* Parallelism mechanism
		* Keep in mind, no branching in algorithms - maybe
		* Note pickling functions & classes just saves the name, not the whole definition
		* Matmul: parallelize over entries in C
			* Note: Result isn't guaranteed to fit in original allocation
			* Note: General add to LR matrix is the hang up
				* Particularly when is a subdivision of a LR matrix 
		* Add: parallelize over entries in C - may not need this
		* LU Factor
		* LU Solve
	* IO mechanism (mainly for OOC matrix)
		* Assurance of honorance of memory limits
		* Assume have a PFS e.g. Lustre
	* Checkpointing system
		* might be easier if have a list of tasks
		* Ideally tracks diffs between check points
		* Ideally transparent, e.g. a decorator
	* Clustering
	* Fill
		* ACA
			* OOC ACA
		* Integrals


#### Verification
* mie series expressions on wikipedia
* mie python

