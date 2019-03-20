#include "ofApp.h"
#include "tetrahedralize.h"
#include "flatten.h"
#include "transformer.h"
#include "AABB.h"
ofVboMesh srfMesh;
ofVboMesh deformedMesh;
nsHEMesh hemesh;
ofVboMesh patternMesh, patternMeshFlat;

vector<int> flattenIndices;

bool drawFlat = false;
bool drawPattern = false;
bool fromFlatTo3D = false;

bool hasFlattened = true;
bool hasPattern = false;
bool hasDeformed = false;

float angleTolerance = 5;
Eigen::MatrixXd Vt, Vf;
Eigen::MatrixXi T;

TransformThread transThread(patternMeshFlat, srfMesh, Vt,T,Vf);
HTree tree;

ofVec3f minOfV(const ofVec3f & a, const ofVec3f & b) {
	return ofVec3f(min(a.x, b.x), min(a.y, b.y), min(a.z, b.z));
}
ofVec3f maxOfV(const ofVec3f & a, const ofVec3f & b) {
	return ofVec3f(max(a.x, b.x), max(a.y, b.y), max(a.z, b.z));
}

void ofToEigen(ofMesh & mesh, Eigen::MatrixXd & V, Eigen::MatrixXi &F) {
	V.resize(mesh.getNumVertices(), 3);
	F.resize(mesh.getNumIndices() / 3, 3);
	for (int i = 0; i < mesh.getNumVertices(); ++i) {
		ofVec3f v = mesh.getVertex(i);
		V(i, 0) = v.x;
		V(i, 1) = v.y;
		V(i, 2) = v.z;
	}

	for (int i = 0; i < mesh.getNumIndices(); i+=3) {
		F(i / 3, 0) = mesh.getIndex(i);
		F(i / 3, 1) = mesh.getIndex(i+1);
		F(i / 3, 2) = mesh.getIndex(i+2);
	}
}

//--------------------------------------------------------------
void ofApp::setup(){
	loadSrf("test.obj");
	gui.setup();
	srfMesh.smoothNormals(0);

	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	ofToEigen(srfMesh, V, F);
	Vf = Vt;
	//get negative normals
	flattenIndices.clear();
	for (int i = 0; i < srfMesh.getNumVertices(); ++i) {
		ofVec3f n = srfMesh.getNormal(i);
		if (n.z > .5) {
			flattenIndices.push_back(i);
		}
	}
	flatten(Vt,T,Vf,flattenIndices);

	//deformedMesh = srfMesh;
	deformedMesh.getVertices() = srfMesh.getVertices();
	deformedMesh.getIndices() = srfMesh.getIndices();
	for (int i = 0; i < Vf.rows(); ++i) {
		deformedMesh.setVertex(i, ofVec3f(Vf(i, 0), Vf(i, 1), Vf(i, 2)));
	}

	deformedMesh.smoothNormals(0);
}

void ofApp::loadSrf(string filename) {
	srfMesh.load(filename);
	zoom(srfMesh);
	if (srfMesh.getNumVertices() > 0) {
		tetrahedralize(srfMesh, Vt, T,2.5);
		hemesh.loadFromOfMesh(srfMesh);
		hemesh.getFaceNormals();
		tree.clear();
		tree.insert(hemesh.faces.begin(), hemesh.faces.end());
		srfMesh.enableColors();
		srfMesh.getColors().resize(srfMesh.getNumVertices(), ofColor::white);
		hasFlattened = false;
		hasPattern = false;
		hasDeformed = false;
	}
	else {
		cout << "no mesh found: " << filename << endl;
	}
}
//--------------------------------------------------------------
void ofApp::update(){

}

void ofApp::zoom(const ofMesh & mesh) {
	ofVec3f minV, maxV;
	minV = maxV = mesh.getVertex(0);
	for (auto & v : mesh.getVertices()) {
		minV = minOfV(minV, v);
		maxV = maxOfV(maxV, v);
	}
	zoom(minV, maxV);
}

void ofApp::zoom(ofVec3f bboxMin, ofVec3f bboxMax) {
	ofVec3f dir = cam.getLookAtDir();
	ofVec3f up = cam.getUpDir();
	ofVec3f center = (bboxMin + bboxMax)*0.5;
	float dist = bboxMin.distance(bboxMax);
	dir.normalize();
	cam.setTarget(center);
	cam.setPosition(center - dir*dist);
	cam.lookAt(center, up);
	cam.setNearClip(dist*.1);
}

void ofApp::guiFunc() {
	gui.begin();
	ImFont * font = ImGui::GetFont();
	font->Scale = 1;
	ImGui::PushFont(font);
	if (ImGui::BeginMainMenuBar())
	{
		if (ImGui::BeginMenu("File"))
		{
			if (ImGui::MenuItem("Open Template")) {
				auto result = ofSystemLoadDialog();
				if (result.bSuccess) {
					string filename = result.filePath;
					loadSrf(filename);
					flattenIndices.clear();
				}
			}
			if (ImGui::MenuItem("Save Flat Template")) {
				//ofSystemSaveDialog("volume.vdb", "Save VDB file");
				auto result = ofSystemSaveDialog("flat.obj", "Save flattened template");
				if (result.bSuccess) {
					deformedMesh.save(result.filePath);
				}
			}
			if (ImGui::MenuItem("Open Pattern")) {
				auto result = ofSystemLoadDialog();
				if (result.bSuccess) {
					patternMesh.load(result.filePath);
					if (patternMesh.getNumNormals() == 0) {
						patternMesh.smoothNormals(0);
					}
					hasPattern = true;
				}
			}
			if (ImGui::MenuItem("Save Deformed Pattern")) {
				auto result = ofSystemSaveDialog("deformed.obj", "Save deformed pattern");
				if (result.bSuccess) {
					patternMeshFlat.save(result.filePath);
				}
			}
			if (ImGui::MenuItem("Save Configuration")) {
				auto result = ofSystemSaveDialog("thing.flat", "Save flat/unflat configuration");
				if (result.bSuccess) {
					saveConfiguration(result.filePath);
				}
			}
			if (ImGui::MenuItem("Load Configuration")) {
				auto result = ofSystemLoadDialog("Load flat/unflat configuration");
				if (result.bSuccess) {
					loadConfiguration(result.filePath);
				}
			}
			ImGui::EndMenu();
		}
		ImGui::EndMainMenuBar();
	}

	ImGui::Begin("tools");

	if (ImGui::Button("flatten")) {
		if (flattenIndices.size() > 0) {
			Vf = Vt;
			flatten(Vt, T, Vf, flattenIndices);
			hasFlattened = true;

			deformedMesh.getIndices() = srfMesh.getIndices();
			deformedMesh.getVertices() = srfMesh.getVertices();
			for (int i = 0; i < Vf.rows(); ++i) {
				deformedMesh.setVertex(i, ofVec3f(Vf(i, 0), Vf(i, 1), Vf(i, 2)));
			}
			deformedMesh.smoothNormals(0);
		}
		else {

		}
	}
	if (ImGui::Button("deform")) {
		if (patternMesh.getNumVertices() > 0) {
			if (fromFlatTo3D) {
				patternMeshFlat.getIndices() = patternMesh.getIndices();
				patternMeshFlat.getVertices() = patternMesh.getVertices();
				transThread.srfMesh = deformedMesh;
				transThread.Vt = Vf;
				transThread.T = T;
				transThread.Vf = Vt;
				transThread.startThread();
			}
			else {
				patternMeshFlat.getIndices() = patternMesh.getIndices();
				patternMeshFlat.getVertices() = patternMesh.getVertices();
				transThread.srfMesh = srfMesh;
				transThread.Vt = Vt;
				transThread.T = T;
				transThread.Vf = Vf;
				transThread.startThread();
			}
			hasDeformed = true;
		}
	}
	if (ImGui::RadioButton("From 3D to Flat", !fromFlatTo3D)) {
		fromFlatTo3D = false;
	}
	if (ImGui::RadioButton("From Flat to 3D", fromFlatTo3D)) {
		fromFlatTo3D = true;
	}
	ImGui::Separator();
	if(ImGui::Button("clear selection")) {
		flattenIndices.clear();

	}
	ImGui::PushItemWidth(100);
	ImGui::InputFloat("angle tolerance", &angleTolerance);
	ImGui::PopItemWidth();
	ImGui::TextWrapped("Click to select area to flatten. SHIFT to add. CTRL to remove.");

	ImGui::Checkbox("draw deformed", &drawFlat);
	ImGui::Checkbox("draw pattern", &drawPattern);
	ImGui::End();
	if (transThread.isThreadRunning()) {
		ImGui::Begin("flattening");
		ImGui::ProgressBar(transThread.complete);
		if (ImGui::Button("cancel")) {
			transThread.stopThread();
		}
		ImGui::End();
	}
}

//--------------------------------------------------------------
void ofApp::draw(){
	ofBackground(0);
	ofSetColor(255);
	guiFunc();
	gui.end();
	ofEnableDepthTest();
	cam.begin();

	ofEnableLighting();
	ofLight light0;
	light0.enable();
	ofSetColor(200);
	if (drawPattern) {
		if (drawFlat) {
			patternMeshFlat.draw();
		}
		else {
			patternMesh.draw();
		}
	}
	else {
		if (drawFlat) {
			deformedMesh.draw();
		}
		else {
			srfMesh.draw();
		}
	}
	cam.end();
	ofDisableLighting();
}

void ofApp::saveConfiguration(string filename) {
	ofstream out(filename, ios::binary | ios::out);
	unsigned int numVertices = Vt.rows();
	unsigned int numTet = T.rows();
	unsigned int numTri = srfMesh.getNumIndices() / 3;

	out.write((const char *)&numVertices, sizeof(unsigned int));
	out.write((const char *)&numTet, sizeof(unsigned int));
	out.write((const char *)&numTri, sizeof(unsigned int));
	out.write((const char *)Vt.data(), sizeof(double)*numVertices*3);
	out.write((const char *)Vf.data(), sizeof(double)*numVertices * 3);
	out.write((const char *)T.data(), sizeof(int)*numTet * 4);
	out.write((const char*)srfMesh.getIndexPointer(), sizeof(ofIndexType)*numTri * 3);

}

void ofApp::loadConfiguration(string filename) {
	ifstream in(filename, ios::binary | ios::in);
	unsigned int numVertices, numTet, numTri;
	in.read((char *)&numVertices, sizeof(unsigned int));
	in.read((char *)&numTet, sizeof(unsigned int));
	in.read((char *)&numTri, sizeof(unsigned int));
	Vt.resize(numVertices, 3);
	Vf.resize(numVertices, 3);
	T.resize(numTet, 4);

	in.read((char *)Vt.data(), sizeof(double)*numVertices * 3);
	in.read((char *)Vf.data(), sizeof(double)*numVertices * 3);
	in.read((char *)T.data(), sizeof(int)*numTet * 4);

	srfMesh.clear();
	srfMesh.getIndices().resize(numTri * 3);
	in.read((char *)srfMesh.getIndexPointer(), sizeof(ofIndexType)*numTri * 3);

	in.close();
	srfMesh.getVertices().resize(numVertices);
	deformedMesh.getVertices().resize(numVertices);
	deformedMesh.getIndices() = srfMesh.getIndices();
	for (int i = 0; i < Vf.rows(); ++i) {
		srfMesh.setVertex(i, ofVec3f(Vt(i, 0), Vt(i, 1), Vt(i, 2)));
		deformedMesh.setVertex(i, ofVec3f(Vf(i, 0), Vf(i, 1), Vf(i, 2)));
	}

	deformedMesh.smoothNormals(0);
	srfMesh.smoothNormals(0);

	hemesh.loadFromOfMesh(srfMesh);
	hemesh.getFaceNormals();
	tree.clear();
	tree.insert(hemesh.faces.begin(), hemesh.faces.end());
	srfMesh.enableColors();
	srfMesh.getColors().resize(srfMesh.getNumVertices(), ofColor::white);
	

	hasFlattened = true;
	hasDeformed = false;

}
//--------------------------------------------------------------
void ofApp::keyPressed(int key){
	if (key == 'f') {
		drawFlat = !drawFlat;
	}
	else if (key == 'i') {
		drawPattern = !drawPattern;
	}
}

//--------------------------------------------------------------
void ofApp::keyReleased(int key){

}

//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y ){

}

//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button){
	if (button == OF_MOUSE_BUTTON_LEFT && !drawFlat && !drawPattern) {
		ofVec3f p1 = cam.screenToWorld(ofVec3f(x,y,-1));
		ofVec3f p2 = cam.screenToWorld(ofVec3f(x, y, 1));
		ofVec3f dir = p2 - p1;
		auto intersect = tree.first_intersected_primitive(BRay(BPoint(p1.x,p1.y,p1.z),BPoint(p2.x,p2.y,p2.z)));
		if (intersect) {
			nsHEFace * start = *(intersect.get());
			ofVec3f norm = start->normal;
			float maxDot = cos(ofDegToRad(angleTolerance));
			for (auto f : hemesh.faces) {
				f->flag = false;
			}
			start->flag = true;
			if (!(ofGetKeyPressed(OF_KEY_SHIFT) || ofGetKeyPressed(OF_KEY_CONTROL))) {
				for (auto v : hemesh.vertices) {
					v->val = 0;
				}
			}
			bool addPts = !ofGetKeyPressed(OF_KEY_CONTROL);
			list<nsHEFace * > faceStack;
			faceStack.push_back(start);
			while (faceStack.size() > 0) {
				nsHEFace * f = faceStack.front();
				faceStack.pop_front();

				nsHEdge * he = f->edge;
				for (int i = 0; i < 3; ++i) {
					if (addPts) {
						he->vertex->val = 1;
					}
					else {
						he->vertex->val = 0;
					}
					if (he->pair->face != NULL) {
						nsHEFace * neighbor = he->pair->face;
						if (!neighbor->flag) {
							neighbor->flag = true;
							if (neighbor->normal.dot(norm) > maxDot) {
								faceStack.push_back(neighbor);
							}
						}
					}
					he = he->next;
				}

			}
			flattenIndices.clear();
			for (int i = 0; i < hemesh.vertices.size(); ++i) {
				if (hemesh.vertices[i]->val == 1) {
					flattenIndices.push_back(i);
					srfMesh.setColor(i, ofColor::red);
				}
				else {
					srfMesh.setColor(i, ofColor::white);
				}
			}
			
			
		}
	}
}

//--------------------------------------------------------------
void ofApp::mouseEntered(int x, int y){

}

//--------------------------------------------------------------
void ofApp::mouseExited(int x, int y){

}

//--------------------------------------------------------------
void ofApp::windowResized(int w, int h){

}

//--------------------------------------------------------------
void ofApp::gotMessage(ofMessage msg){

}

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo){ 

		if (dragInfo.files.size() > 0) {
			
		}

}
