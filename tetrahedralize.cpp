#include "tetrahedralize.h"
#define BOOST_PARAMETER_MAX_ARITY 12

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <boost/function_output_iterator.hpp>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/make_mesh_3.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;
typedef Mesh::Vertex_index Vertex_index;

typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Mesh>::edge_descriptor     edge_descriptor;

typedef CGAL::Polyhedral_mesh_domain_with_features_3<K> Mesh_domain;
// Polyhedron type
typedef CGAL::Mesh_polyhedron_3<K>::type Polyhedron;
// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<
	Tr, Mesh_domain::Corner_index, Mesh_domain::Curve_segment_index> C3t3;
// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
// To avoid verbose function and named parameters call
using namespace CGAL::parameters;
using namespace Eigen;
void tetrahedralize(Eigen::MatrixXf & V, Eigen::MatrixXi &F, Eigen::MatrixXd & Vt, Eigen::MatrixXi &T, float maxTriangle) {
	using namespace std;
	// Load a polyhedron
	Polyhedron poly;
	CGAL::Polyhedron_incremental_builder_3<Polyhedron::HalfedgeDS> B(poly.hds(), true);
	B.begin_surface(V.rows(), F.rows());
	for (int i = 0; i < V.rows();++i) {
		B.add_vertex(Polyhedron::HalfedgeDS::Vertex::Point(V(i,0),V(i,1),V(i,2)));
	}
	for (int i = 0; i < F.rows(); i++) {
		B.begin_facet();
		B.add_vertex_to_facet(F(i,0));
		B.add_vertex_to_facet(F(i, 1));
		B.add_vertex_to_facet(F(i, 2));
		B.end_facet();
	}
	B.end_surface();
	// Create a vector with only one element: the pointer to the polyhedron.
	std::vector<Polyhedron*> poly_ptrs_vector(1, &poly);
	// Create a polyhedral domain, with only one polyhedron,
	// and no "bounding polyhedron", so the volumetric part of the domain will be
	// empty.
	Mesh_domain domain(poly);

	// Get sharp features
	//domain.detect_borders();
	domain.detect_features(90); //includes detection of borders
	// Mesh criteria
	Mesh_criteria criteria(facet_angle = 23,
		facet_size = maxTriangle,
		facet_distance = 0.04, cell_radius_edge_ratio = 2, cell_size = 2);

	// Mesh generation
	C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_perturb(), no_exude());
	//turn back on
	c3t3.rescan_after_load_of_triangulation();

	// Output the facets of the c3t3 to an OFF file. The facets will not be
	// oriented.
	typedef C3t3::Triangulation::Vertex_handle Vertex_handle;
	map<Vertex_handle, unsigned int> vMap;
	int vIndex = 0;


	int numTets = 0;

	for (auto it = c3t3.cells_in_complex_begin(); it != c3t3.cells_in_complex_end(); ++it) {
		numTets++;
	}
	T.resize(numTets, 4);

	Vt.resize(numTets * 4, 3);
	int tI = 0;
	for (auto it = c3t3.cells_in_complex_begin(); it != c3t3.cells_in_complex_end(); ++it) {
		for (int i = 0; i < 4; ++i) {
			auto vh = it->vertex(i);
			auto findex = vMap.find(vh);
			int index = vIndex;
			if (findex == vMap.end()) {
				auto pt = vh->point();
				Vt.row(vIndex) = Vector3d(pt.x(), pt.y(), pt.z());
				vMap.emplace(vh, vIndex);
				vIndex++;
			}
			else {
				index = findex->second;
			}
			T(tI, i) = index;
		}
		tI++;
	}

	int numVerts = 0;
	numVerts = vIndex;
	Vt.conservativeResize(vIndex, 3);

	cout << "tets " << numTets << endl;


	std::size_t inum = 0;
	std::size_t nfacets = 0;
	array<std::size_t, 3> indices = { { 0,0,0 } };
	//vMap.clear();
	/*
	mesh.clear();
	for (auto fit = c3t3.facets_in_complex_begin(),
		end = c3t3.facets_in_complex_end();
		fit != end; ++fit)
	{
		C3t3::Subdomain_index cell_sd = c3t3.subdomain_index(fit->first);
		C3t3::Subdomain_index opp_sd = c3t3.subdomain_index(fit->first->neighbor(fit->second));

		if (cell_sd != 0 && opp_sd != 0) continue;

		++nfacets;
		int j = -1;


		for (int i = 0; i < 4; ++i) {
			if (i != fit->second) {
				auto vh = fit->first->vertex(i);
				auto findex = vMap.find(vh);
				int index;
				if (findex == vMap.end()) {
					index = mesh.getNumVertices();
					mesh.addVertex(ofVec3f(vh->point().x(), vh->point().y(), vh->point().z()));
					vMap.emplace(vh, index);
				}
				else {
					index = findex->second;
				}
				indices[++j] = index;
			}
		}
		if (((cell_sd == 0) == (fit->second % 2 == 1)) == 1) {
			std::swap(indices[0], indices[1]);
		}
		//facet_buffer << "3" << " " << indices[0] << " " << indices[1] << " " << indices[2] << "\n";
		mesh.addIndex(indices[0]);
		mesh.addIndex(indices[2]);
		mesh.addIndex(indices[1]);
	}
	*/
	
}

void tetrahedralize(ofMesh & mesh, Eigen::MatrixXd & Vt, Eigen::MatrixXi &T, float maxTriangle) {
	using namespace std;
	// Load a polyhedron
	Polyhedron poly;
	CGAL::Polyhedron_incremental_builder_3<Polyhedron::HalfedgeDS> B(poly.hds(), true);
	B.begin_surface(mesh.getNumVertices(), mesh.getNumIndices()/3);
	for (ofVec3f & v : mesh.getVertices()) {
		B.add_vertex(Polyhedron::HalfedgeDS::Vertex::Point(v.x,v.y,v.z));
	}
	for (int i = 0; i < mesh.getNumIndices(); i+=3) {
		B.begin_facet();
		B.add_vertex_to_facet(mesh.getIndex(i));
		B.add_vertex_to_facet(mesh.getIndex(i+1));
		B.add_vertex_to_facet(mesh.getIndex(i+2));
		B.end_facet();
	}
	B.end_surface();
	// Create a vector with only one element: the pointer to the polyhedron.
	std::vector<Polyhedron*> poly_ptrs_vector(1, &poly);
	// Create a polyhedral domain, with only one polyhedron,
	// and no "bounding polyhedron", so the volumetric part of the domain will be
	// empty.
	Mesh_domain domain(poly);

	// Get sharp features
	//domain.detect_borders();
	domain.detect_features(90); //includes detection of borders
								// Mesh criteria
	Mesh_criteria criteria(facet_angle = 23,
		facet_size = maxTriangle,
		facet_distance = 0.06, cell_radius_edge_ratio = 2, cell_size = maxTriangle);

	// Mesh generation
	C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_perturb(), no_exude());
	//turn back on
	c3t3.rescan_after_load_of_triangulation();

	// Output the facets of the c3t3 to an OFF file. The facets will not be
	// oriented.
	typedef C3t3::Triangulation::Vertex_handle Vertex_handle;
	map<Vertex_handle, unsigned int> vMap;
	int vIndex = 0;


	int numTets = 0;

	for (auto it = c3t3.cells_in_complex_begin(); it != c3t3.cells_in_complex_end(); ++it) {
		numTets++;
	}
	T.resize(numTets, 4);

	Vt.resize(numTets * 4, 3);
	int tI = 0;
	mesh.clear();
	for (auto it = c3t3.cells_in_complex_begin(); it != c3t3.cells_in_complex_end(); ++it) {
		for (int i = 0; i < 4; ++i) {
			auto vh = it->vertex(i);
			auto findex = vMap.find(vh);
			int index = vIndex;
			if (findex == vMap.end()) {
				auto pt = vh->point();
				Vt.row(vIndex) = Vector3d(pt.x(), pt.y(), pt.z());
				mesh.addVertex(ofVec3f(pt.x(), pt.y(), pt.z()));
				vMap.emplace(vh, vIndex);
				vIndex++;
			}
			else {
				index = findex->second;
			}
			T(tI, i) = index;
		}
		tI++;
	}

	int numVerts = 0;
	numVerts = vIndex;
	Vt.conservativeResize(vIndex, 3);

	cout << "tets " << numTets << endl;


	std::size_t inum = 0;
	std::size_t nfacets = 0;
	array<std::size_t, 3> indices = { { 0,0,0 } };
	//vMap.clear();
	
	
	for (auto fit = c3t3.facets_in_complex_begin(),
	end = c3t3.facets_in_complex_end();
	fit != end; ++fit)
	{
		C3t3::Subdomain_index cell_sd = c3t3.subdomain_index(fit->first);
		C3t3::Subdomain_index opp_sd = c3t3.subdomain_index(fit->first->neighbor(fit->second));

		if (cell_sd != 0 && opp_sd != 0) continue;

		++nfacets;
		int j = -1;

		bool goodFace = true;
		for (int i = 0; i < 4; ++i) {
			if (i != fit->second) {
				auto vh = fit->first->vertex(i);
				auto findex = vMap.find(vh);
				int index;
				if (findex == vMap.end()) {
					index = mesh.getNumVertices();
					//mesh.addVertex(ofVec3f(vh->point().x(), vh->point().y(), vh->point().z()));
					//vMap.emplace(vh, index);
					goodFace = false;
					cout << "bad face. how does this happen?" << endl;
				} else {
					index = findex->second;
				}
				indices[++j] = index;
			}
		}
		if (goodFace) {
			if (((cell_sd == 0) == (fit->second % 2 == 1)) == 1) {
				std::swap(indices[0], indices[1]);
			}
			//facet_buffer << "3" << " " << indices[0] << " " << indices[1] << " " << indices[2] << "\n";
			mesh.addIndex(indices[0]);
			mesh.addIndex(indices[2]);
			mesh.addIndex(indices[1]);
		}
	}
	

}