#pragma once

#include "ofMain.h"
#include "ofxImGui.h"
class ofApp : public ofBaseApp{

	public:
		ofxImGui::Gui gui;
		ofEasyCam cam;
		void loadSrf(string filename);
		void setup();
		void update();
		void draw();
		void guiFunc();
		void zoom(ofVec3f bboxMin, ofVec3f bboxMax);
		void zoom(const ofMesh & mesh);
		void keyPressed(int key);
		void keyReleased(int key);
		void mouseMoved(int x, int y );
		void mouseDragged(int x, int y, int button);
		void mousePressed(int x, int y, int button);
		void mouseReleased(int x, int y, int button);
		void mouseEntered(int x, int y);
		void mouseExited(int x, int y);
		void windowResized(int w, int h);
		void dragEvent(ofDragInfo dragInfo);
		void gotMessage(ofMessage msg);
		
};
