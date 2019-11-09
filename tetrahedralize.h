#pragma once
#include <Eigen/Dense>
#include "ofMesh.h"
void tetrahedralize(Eigen::MatrixXf & V, Eigen::MatrixXi &F, Eigen::MatrixXd & Vt, Eigen::MatrixXi &T, float maxTriangle = 1.0);
void tetrahedralize(ofMesh & mesh, Eigen::MatrixXd & Vt, Eigen::MatrixXi &T, float maxTriangle = 1.0);
void tetrahedralizeWild(ofMesh & mesh, Eigen::MatrixXd & Vt, Eigen::MatrixXi &T, float maxTriangle = 1.0);