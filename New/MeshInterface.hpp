
#ifndef MESH_INTERFACE_HPP
#define MESH_INTERFACE_HPP

#include "../Common/Iterators.hpp"

/*
Minimial definition of a mesh required to support the tracer engine.
Can be instantiated as-is to wrap mesh data, or can be used as a base class for a more feature-rich object.
*/
class MeshInterface {
	private:
		// Data format is driven by ray tracers.
		// Currently assuming element-major (as opposed to component-major).
		// This interface only points to the data, it does not store or manage it.
		// There may be a better abstraction than this based on iterators, but that may also be overkill.
		// Iterators proivded here help access elements and their components.
		unsigned _n_facets, _n_vertices;
		unsigned* facets;
		double* vertices;
		
	public:
		MeshInterface( unsigned num_facets, unsigned num_vertices, unsigned* facets_data, double* vertices_data );
		
		// Minimal set of functionality needed by Engine to access the mesh
		unsigned n_facets() const ;
		unsigned n_vertices() const ;
		
		typedef Iterator2D<double*> VertexIterator ;
		typedef Iterator2D<unsigned*> FacetIterator ;
		
		VertexIterator vertices_begin() const ;
		VertexIterator vertices_end() const ;
		
		FacetIterator facets_begin() const ;
		FacetIterator facets_end() const ;
};

#endif /* MESH_INTERFACE_HPP */

