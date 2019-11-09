#include "ofMain.h"
#include "ofApp.h"

//========================================================================
int main(int argc, char * argv[]){
	ofSetupOpenGL(1024,768,OF_WINDOW);			// <-------- setup the GL context

												// vector for storing args
	vector<string> myArgs;

	if (argc > 1) {
		for (int i = 0; i < argc; i++) {
			cout << "arg num " << i << " : " << argv[i] << endl;
			myArgs.push_back(argv[i]);
		}
	}

	// create a window
	//ofAppGlutWindow window;

	// set width, height, mode (OF_WINDOW or OF_FULLSCREEN)
	//ofSetupOpenGL(&window, 1024, 768, OF_WINDOW);

	ofApp *app = new ofApp();
	// add args to app instance
	app->arguments = myArgs;
	ofRunApp(app);

}
