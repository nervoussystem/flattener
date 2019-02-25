#pragma once
#include <Eigen/Dense>
#include "ofMesh.h"
#include "ofThread.h"

class TransformThread : public ofThread {
public:
	TransformThread(ofMesh & _patternMeshFlat,
		ofMesh & _srfMesh,
		Eigen::MatrixXd & _Vt,
		Eigen::MatrixXi & _T,
		Eigen::MatrixXd & _Vf) :patternMesh(_patternMeshFlat), srfMesh(_srfMesh), Vt(_Vt), T(_T), Vf(_Vf) {

	}
	void TransformThread::threadedFunction();
	float complete;
	ofMesh & patternMesh;
	ofMesh  srfMesh;
	Eigen::MatrixXd  Vt;
	Eigen::MatrixXi  T;
	Eigen::MatrixXd  Vf;
};
void transformSrf(ofMesh & srf, ofMesh & boundary, Eigen::MatrixXd & V, Eigen::MatrixXi & T, Eigen::MatrixXd & V_deformed);