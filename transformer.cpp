#include "transformer.h"
#include "AABB1.h"
#include "AABB.h"
#include <vector>
using namespace Eigen;
using namespace std;

double determinant(Vector3d & p1, Vector3d & p2, Vector3d & p3) {
	return (p1.y()*p2.z() - p1.z()*p2.y())*p3.x() + (p1.z()*p2.x() - p1.x()*p2.z())*p3.y() + (p1.x()*p2.y() - p1.y()*p2.x())*p3.z();
}

bool isInsideTet(int tIndex, Vector3d & p, Vector3d & bary, MatrixXd & V, MatrixXi & T) {
	Vector3d t1 = V.row(T(tIndex, 0));
	Vector3d t2 = V.row(T(tIndex, 1));
	Vector3d t3 = V.row(T(tIndex, 2));
	Vector3d t4 = V.row(T(tIndex, 3));
	Vector3d v = p - t1;
	t2 -= t1;
	t3 -= t1;
	t4 -= t1;

	double vol = determinant(t2, t3, t4);
	if (abs(vol) > 1e-7) {
		double invVol = 1.0 / vol;
		bary.x() = determinant(v, t3, t4)*invVol;
		if (bary.x() < 0) return false;
		bary.y() = determinant(t2, v, t4)*invVol;
		if (bary.y() < 0) return false;
		bary.z() = determinant(t2, t3, v)*invVol;
		if (bary.z() < 0) return false;
		if (bary.x() + bary.y() + bary.z() > 1.0) return false;
		return true;
	}
	else {
		return false;
	}
}

Vector3d transformDistortPt(Vector3d & p, MatrixXd & V, MatrixXi & T, MatrixXd & V_deformed, aabb::Tree & tetTree, BTree * tree) {
	int Tindex = -1;
	double b2, b3, b4;
	vector<double> minPt = { p.x() - .01, p.y() - .01, p.z() - .01 }, maxPt = { p.x() + .01, p.y() + .01, p.z() + .01 };
	aabb::AABB q(minPt, maxPt);
	Vector3d bary;
	vector<unsigned int> nearest = tetTree.query(q);
	for (int i = 0; i<nearest.size(); ++i) {
		//for (int i = 0; i < T.rows(); ++i) {
		int tIndex = nearest[i];
		if (isInsideTet(tIndex, p, bary, V,T)) {
			Vector3d s1 = V_deformed.row(T(tIndex, 0));
			Vector3d s2 = V_deformed.row(T(tIndex, 1));
			Vector3d s3 = V_deformed.row(T(tIndex, 2));
			Vector3d s4 = V_deformed.row(T(tIndex, 3));

			double b1 = 1 - bary.x() - bary.y() - bary.z();
			Vector3d newPt = s1*b1 + s2* bary.x() + s3* bary.y() + s4* bary.z();
			return newPt;
		}
	}
	
	auto closest = tree->closest_point_and_primitive(BPoint(p.x(), p.y(), p.z()));
	auto info = (closest.second);
	auto base = info.base();
	int i1 = *base;
	++base;
	int i2 = *base;
	++base;
	int i3 = *base;

	Vector3d v1 = V.row(i1);
	Vector3d v2 = V.row(i2);
	Vector3d v3 = V.row(i3);

	v2 -= v1;
	v3 -= v1;

	Vector3d norm = v2.cross(v3);
	float area = norm.norm();
	norm /= area;

	float nDist = norm.dot(p - v1);
	Vector3d c = p - nDist*norm;
	//Vector3d c(closest.first.x(), closest.first.y(), closest.first.z());
	c -= v1;

	float a3 = v2.cross(c).dot(norm);
	float a2 = c.cross(v3).dot(norm);

	a2 /= area;
	a3 /= area;
	float a1 = 1.0 - a2 - a3;


	Vector3d vb1 = V_deformed.row(i1);
	Vector3d vb2 = V_deformed.row(i2);
	Vector3d vb3 = V_deformed.row(i3);

	Vector3d normD = (vb2-vb1).cross((vb3-vb1));
	float areaD = normD.norm();

	normD /= areaD;
	c = vb1*a1 + vb2*a2 + vb3*a3;
	c += normD*nDist;
	//cout << "project" << endl;
	return c;
}

void transformSrf(ofMesh & srf, ofMesh & boundary, MatrixXd & V, MatrixXi & T, MatrixXd & V_deformed) {
	aabb::Tree tetTree;

	vector<double> minV(3), maxV(3);
	vector<aabb::AABB> boxes;
	for (int i = 0; i < T.rows(); ++i) {
		Vector3d t1 = V.row(T(i, 0));
		Vector3d t2 = V.row(T(i, 1));
		Vector3d t3 = V.row(T(i, 2));
		Vector3d t4 = V.row(T(i, 3));
		Vector3d minPt = t1.cwiseMin(t2);
		minPt = minPt.cwiseMin(t3);
		minPt = minPt.cwiseMin(t4);
		Vector3d maxPt = t1.cwiseMax(t2);
		maxPt = maxPt.cwiseMax(t3);
		maxPt = maxPt.cwiseMax(t4);
		minV.assign(minPt.data(), minPt.data() + 3);
		maxV.assign(maxPt.data(), maxPt.data() + 3);
		//tetTree.insertParticle(i, minV,maxV);
		boxes.push_back(aabb::AABB(minV, maxV));
	}
	tetTree.build(boxes);

	BTree * tree = new BTree(Triangle_iterator(boundary.getIndices().begin(), &(boundary.getVertices())), Triangle_iterator(boundary.getIndices().end(), &(boundary.getVertices())));
	tree->accelerate_distance_queries();

	for (int i = 0; i < srf.getNumVertices();++i) {
		ofVec3f v = srf.getVertex(i);
		Vector3d vd(v.x, v.y, v.z);
		vd = transformDistortPt(vd, V, T, V_deformed, tetTree, tree);
		//v.set(vd.x(), vd.y(), vd.z());
		srf.setVertex(i, ofVec3f(vd.x(), vd.y(), vd.z()));
	}
}

struct ofMeshClosestPt{

	double operator() (const array<double, 3> & pt, int tri, array<double,3> & closest) {
		int t3 = tri * 3;
		int i1 = mesh->getIndex(t3);
		int i2 = mesh->getIndex(t3+1);
		int i3 = mesh->getIndex(t3+2);

		ofVec3f & v1 = mesh->getVertex(i1);
		ofVec3f & v2 = mesh->getVertex(i2);
		ofVec3f & v3 = mesh->getVertex(i3);

		ofVec3f p(pt[0], pt[1], pt[2]);
		v2 -= v1;
		v3 -= v1;
		p = v1-p;
		
		float a = v2.dot(v2);
		float b = v2.dot(v3);
		float c = v3.dot(v3);
		float d = v2.dot(p);
		float e = v3.dot(p);

		float det = a*c - b*b;
		float s = b*e - c*d;
		float t = b*d - a*e;
		if (s + t < det)
		{
			if (s < 0.f)
			{
				if (t < 0.f)
				{
					if (d < 0.f)
					{
						s = min(max(-d / a, 0.f), 1.f);
						t = 0.f;
					}
					else
					{
						s = 0.f;
						t = min(max(-e / c, 0.f), 1.f);
					}
				}
				else
				{
					s = 0.f;
					t = min(max(-e / c, 0.f), 1.f);
				}
			}
			else if (t < 0.f)
			{
				s = min(max(-d / a, 0.f), 1.f);
				t = 0.f;
			}
			else
			{
				float invDet = 1.f / det;
				s *= invDet;
				t *= invDet;
			}
		}
		else
		{
			if (s < 0.f)
			{
				float tmp0 = b + d;
				float tmp1 = c + e;
				if (tmp1 > tmp0)
				{
					float numer = tmp1 - tmp0;
					float denom = a - 2 * b + c;
					s = min(max(numer / denom, 0.f), 1.f);
					t = 1 - s;
				}
				else
				{
					t = min(max(-e / c, 0.f), 1.f);
					s = 0.f;
				}
			}
			else if (t < 0.f)
			{
				if (a + d > b + e)
				{
					float numer = c + e - b - d;
					float denom = a - 2 * b + c;
					s = min(max(numer / denom, 0.f), 1.f);
					t = 1 - s;
				}
				else
				{
					s = min(max(-e / c, 0.f), 1.f);
					t = 0.f;
				}
			}
			else
			{
				float numer = c + e - b - d;
				float denom = a - 2 * b + c;
				s = min(max(numer / denom, 0.f), 1.f);
				t = 1.f - s;
			}
		}
		
		ofVec3f close = s * v2 + t * v3;
		double dist = close.distanceSquared(-p);
		close += v1;
		closest[0] = close.x;
		closest[1] = close.y;
		closest[2] = close.z;
		return dist;
	}
	ofMesh * mesh;
};

void TransformThread::threadedFunction() {
	aabb::Tree tetTree;

	vector<double> minV(3), maxV(3);
	vector<aabb::AABB> boxes;
	for (int i = 0; i < T.rows(); ++i) {
		Vector3d t1 = Vt.row(T(i, 0));
		Vector3d t2 = Vt.row(T(i, 1));
		Vector3d t3 = Vt.row(T(i, 2));
		Vector3d t4 = Vt.row(T(i, 3));
		Vector3d minPt = t1.cwiseMin(t2);
		minPt = minPt.cwiseMin(t3);
		minPt = minPt.cwiseMin(t4);
		Vector3d maxPt = t1.cwiseMax(t2);
		maxPt = maxPt.cwiseMax(t3);
		maxPt = maxPt.cwiseMax(t4);
		minV.assign(minPt.data(), minPt.data() + 3);
		maxV.assign(maxPt.data(), maxPt.data() + 3);
		//tetTree.insertParticle(i, minV,maxV);
		boxes.push_back(aabb::AABB(minV, maxV));
	}
	tetTree.build(boxes);

	BTree * tree = new BTree(Triangle_iterator(srfMesh.getIndices().begin(), &(srfMesh.getVertices())), Triangle_iterator(srfMesh.getIndices().end(), &(srfMesh.getVertices())));
	tree->accelerate_distance_queries();

	aabb::Tree srfTree;

	boxes.clear();
	for (int i = 0; i < srfMesh.getNumIndices();) {
		int index = i / 3;
		int i1 = srfMesh.getIndex(i++);
		ofVec3f v = srfMesh.getVertex(i1);
		minV[0] = maxV[0] = v[0];
		minV[1] = maxV[1] = v[1];
		minV[2] = maxV[2] = v[2];
		for (int j = 1; j < 3; ++j) {
			i1 = srfMesh.getIndex(i++);
			ofVec3f & v1 = srfMesh.getVertex(i1);
			minV[0] = min<double>(minV[0], v1.x);
			minV[1] = min<double>(minV[1], v1.y);
			minV[2] = min<double>(minV[2], v1.z);
			maxV[0] = max<double>(maxV[0], v1.x);
			maxV[1] = max<double>(maxV[1], v1.y);
			maxV[2] = max<double>(maxV[2], v1.z);
		}
		//tetTree.insertParticle(i, minV,maxV);
		boxes.push_back(aabb::AABB(minV, maxV));
	}
	srfTree.build(boxes);

	ofMeshClosestPt closestPter;
	closestPter.mesh = &srfMesh;
	int numSrfPts = patternMesh.getNumVertices();
	auto & pts = patternMesh.getVertices();

	uint64_t totalMe = 0;
	uint64_t totalYou = 0;
	for (int i = 0; i < numSrfPts && isThreadRunning(); ++i) {
		ofVec3f & v = pts[i];
		Vector3d vd(v.x, v.y, v.z);
		vd = transformDistortPt(vd, Vt, T, Vf, tetTree, tree);
		auto t1 = ofGetElapsedTimeMicros();
		array<double, 3> aPt = { v.x,v.y,v.z };
		auto closest1 = srfTree.closestPointAndPrimitive(aPt, 2.0, closestPter);
		auto t2 = ofGetElapsedTimeMicros();
		auto closest2 = tree->closest_point_and_primitive(BPoint(v.x,v.y,v.z));
		auto t3 = ofGetElapsedTimeMicros();

		totalMe += t2 - t1;
		totalYou += t3 - t2;
		v.set(vd.x(), vd.y(), vd.z());
		//srf.setVertex(i, ofVec3f(vd.x(), vd.y(), vd.z()));
		complete = i*1.0 / numSrfPts;
	}
	cout << totalMe << " " << totalYou << endl;
	patternMesh.smoothNormals(0);
	patternMesh.getVertices();
}

