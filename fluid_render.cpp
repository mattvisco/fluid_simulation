    //
//  fluid.cpp
//  
//
//  Created by Matthew Visco/Thomas Yopes/Justin Kay on 4/10/13.
//
//

#include "fluid_render.h"


inline float sqr(float x) { return x*x; }

using namespace std;
using namespace glm;

//
// Some Classes
//

class Viewport;

class Viewport {
public:
    int w, h; // width and height
    
};

//
// Global Variables
//
Viewport viewport;
float gridX, gridY, gridZ, cellSize;
Grid grid;
vector<Particle> particles;
Simulator simulator;

//
// Simple init function
//
void initScene(){
    glClearColor (0.0, 0.0, 0.0, 0.0);
    glClearDepth(1.0);
    
    Grid grid(gridX, gridY, gridZ, cellSize);
    for (int i=0;i<grid.xcells;i++){
        for (int j=0; j<grid.ycells;j++){
            for (int h=0; h<grid.zcells; h++){
                grid.pressures[i][j][h]=10.0;
                grid.xvelocityOld[i][j][h]=2.0;
                grid.yvelocityOld[i][j][h]=2.0;
                grid.zvelocityOld[i][j][h]=2.0;
            }
        }
    }
    grid.xvelocityOld[grid.xcells][grid.ycells-1][grid.zcells-1]=3.0;
    grid.yvelocityOld[grid.xcells-1][grid.ycells][grid.zcells-1]=3.0;
    grid.zvelocityOld[grid.xcells-1][grid.ycells-1][grid.zcells]=3.0;
    grid.computePressure();
    //Simulator simulator(particles);
    //Particle h;
}

//
// reshape viewport if the window is resized
//
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

//
// function that does the actual drawing of stuff
//
void myDisplay(void) {
    
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // Clear Color and Depth Buffers
    glLoadIdentity ();     // make sure transformation is "zero'd"
    
    // Render all particles
    glBegin(GL_POINTS);
    for (vector<Particle>::iterator particle = particles.begin(); particle != particles.end(); ++particle) {
        glColor3f(0,0,1);
        glVertex3f(((Particle)*particle).pos.x, ((Particle)*particle).pos.y, ((Particle)*particle).pos.z);
    }
    glEnd();
    
    // Simulates the movement of all particles
    // at current time step
    simulator.simulate();
    
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
    

    // Currently grid size hard coded in,
    // later parsed from command line
    gridX = 1000, gridY = 1000, gridZ = 1000;
    cellSize = 100;
    // Cutty scene set up, later make robust function
    Particle particle1(vec3(1,0,0),vec3(1,0,0),vec3(1,0,0),vec3(1,0,0),3.0,3.0);
    particles.push_back(particle1);
    Particle particle2(vec3(100,100,100),vec3(100,100,100),vec3(100,100,100),vec3(100,100,100),4.0,4.0);
    particles.push_back(particle2);
    Particle particle3(vec3(500,500,0),vec3(500,500,0),vec3(500,500,0),vec3(500,500,0),7.0,7.0);
    particles.push_back(particle3);
    Particle particle4(vec3(0,300,700),vec3(0,300,700),vec3(0,300,700),vec3(0,300,700),10.0,10.0);
    particles.push_back(particle4);
    
    
    initScene();							// quick function to set up scene
    
    
    //glutDisplayFunc(myDisplay);				// function to run when its time to draw something
    //glutReshapeFunc(myReshape);				// function to run when the window gets resized
    //glutIdleFunc(myDisplay);
    
    //glutMainLoop();	// infinite loop that will keep drawing and resizing
    
    return 0;
}

