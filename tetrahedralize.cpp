#include "tetrahedralize.h"
#define BOOST_PARAMETER_MAX_ARITY 12

#include <tetwild/tetwild.h>
#include <floattetwild/FloatTetwild.h>
#include <igl/boundary_facets.h>


using namespace Eigen;
void tetrahedralize(Eigen::MatrixXf & V, Eigen::MatrixXi &F, Eigen::MatrixXd & Vt, Eigen::MatrixXi &T, float maxTriangle) {
	
	
}

void tetrahedralize(ofMesh & mesh, Eigen::MatrixXd & Vt, Eigen::MatrixXi &T, float maxTriangle) {
	

}


void tetrahedralizeWild(ofMesh & mesh, Eigen::MatrixXd & Vt, Eigen::MatrixXi &T, float maxTriangle) {
	GEO::initialize();
	using namespace std;
	// Load a polyhedron
	GEO::Mesh gmesh;
	
	Eigen::MatrixXd V(mesh.getNumVertices(),3);
	Eigen::MatrixXi F(mesh.getNumIndices()/3,3);
	double p[3];
	for (int i = 0; i < mesh.getNumVertices();++i) {
		ofVec3f & v = mesh.getVertex(i);
		p[0] = v.x;
		p[1] = v.y;
		p[2] = v.z;
		gmesh.vertices.create_vertex(p);
		V(i, 0) = v.x;
		V(i, 1) = v.y;
		V(i, 2) = v.z;
	}
	for (int i = 0; i < mesh.getNumIndices(); i += 3) {
		gmesh.facets.create_triangle(mesh.getIndex(i), mesh.getIndex(i + 1), mesh.getIndex(i + 2));
		F(i / 3, 0) = mesh.getIndex(i);
		F(i / 3, 1) = mesh.getIndex(i+1);
		F(i / 3, 2) = mesh.getIndex(i+2);
	}
	gmesh.facets.connect();
	gmesh.facets.compute_borders();
	gmesh.show_stats();

	tetwild::Args args;
	args.write_csv_file = false;
	Eigen::VectorXd A;
	floatTetWild::Parameters params;
	params.is_quiet = false;
	floatTetWild::tetrahedralization(gmesh, params, Vt, T);
	//tetwild::tetrahedralization(V, F, Vt, T, A, args);
	//tetwild::extractSurfaceMesh(Vt, T, V, F);
	igl::boundary_facets(T, F);
	mesh.clear();
	cout << "num tets " << T.rows() << " num pts " << Vt.rows() << endl;
	for (int i = 0; i < Vt.rows(); ++i) {
		ofVec3f pt(Vt(i, 0), Vt(i, 1), Vt(i, 2));
		mesh.addVertex(pt);
	}
	for (int i = 0; i < F.rows(); ++i) {
		mesh.addTriangle(F(i, 0), F(i, 1), F(i, 2));
	}

}