//
//  fluid.cpp
//  
//
//  Created by Matthew Visco/Thomas Yopes/Justin Kay on 4/10/13.
//
//

#include "fluid.h"

inline float sqr(float x) { return x*x; }

using namespace std;
using namespace glm;

//****************************************************
// Some Classes
//****************************************************

class Viewport;

class Viewport {
public:
    int w, h; // width and height
    
};

//****************************************************
// Global Variables
//****************************************************
Viewport viewport;

//****************************************************
// Simple init function
//****************************************************
void initScene(){
    glClearColor (0.0, 0.0, 0.0, 0.0);
    glClearDepth(1.0);
}

//****************************************************
// reshape viewport if the window is resized
//****************************************************
void myReshape(int w, int h) {
    viewport.w = w;
    viewport.h = h;

    float ratio =  w * 1.0 / h;
    
    // Use the Projection Matrix
	glMatrixMode(GL_PROJECTION);
    
    // Reset Matrix
	glLoadIdentity();
    
	// Set the viewport to be the entire window
	glViewport(0, 0, viewport.w, viewport.h);
    
	// Set the correct perspective.
	gluPerspective(45,ratio,1,100);
    
    
	// Get Back to the Modelview
	glMatrixMode(GL_MODELVIEW);
}

//****************************************************
// function that does the actual drawing of stuff
//***************************************************
void myDisplay(void) {
    
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // Clear Color and Depth Buffers
    glLoadIdentity ();     // make sure transformation is "zero'd"
    
    
    glFlush();
    glutSwapBuffers();					// swap buffers (we earlier set double buffer)
}

int main(int argc, char *argv[]) {

    //This initializes glut
    glutInit(&argc, argv);
    
    //This tells glut to use a double-buffered window with red, green, and blue channels
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    
    // Initalize theviewport size
    viewport.w = 400;
    viewport.h = 400;
    
    //The size and position of the window
    glutInitWindowSize(viewport.w, viewport.h);
    glutInitWindowPosition(0,0);
    glutCreateWindow(argv[0]);
    
    initScene();							// quick function to set up scene
    
    glutDisplayFunc(myDisplay);				// function to run when its time to draw something
    glutReshapeFunc(myReshape);				// function to run when the window gets resized
    glutIdleFunc(myDisplay);
    
    glutMainLoop();	// infinite loop that will keep drawing and resizing
    
    return 0;
}