
#include "MeshInterface.hpp"

// Basic constructor
MeshInterface::MeshInterface( unsigned num_facets, unsigned num_vertices, unsigned* facets_data, double* vertices_data )
	: _n_facets(num_facets), _n_vertices(num_vertices), facets(facets_data), vertices(vertices_data) {}

unsigned MeshInterface::n_vertices() const { return _n_vertices; }
unsigned MeshInterface::n_facets() const { return _n_facets; }

MeshInterface::VertexIterator MeshInterface::vertices_begin() const { return VertexIterator( vertices, 3 ); }
MeshInterface::VertexIterator MeshInterface::vertices_end() const { return vertices_begin() + n_vertices(); }

MeshInterface::FacetIterator MeshInterface::facets_begin() const { return FacetIterator( facets, 3 ); }
MeshInterface::FacetIterator MeshInterface::facets_end() const { return facets_begin() + n_facets(); }
