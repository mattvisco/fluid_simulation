    //
//  fluid.cpp
//  
//
//  Created by Matthew Visco/Thomas Yopes/Justin Kay on 4/10/13.
//
//

#include "fluid_render.h"

#define KB_SPACE 32

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
vector<Particle> particles;
Simulator simulator;
float xrot,yrot,zrot,xtrans,ytrans,zoom;

//
// Simple init function
//
void initScene(){
    glClearColor (0.0, 0.0, 0.0, 0.0);
    glClearDepth(1.0);
    
    simulator = Simulator(&particles, gridX, gridY, gridZ, cellSize);
}

//orient the camera
void camera() {
    glTranslatef(0.0+xtrans,0.0+ytrans,zoom-30.0);
    glRotatef(xrot+25,1.0,0.0,0.0);
    glRotatef(yrot-35,0.0,1.0,0.0);
    glRotatef(zrot,0.0,0.0,1.0);
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
    glMatrixMode(GL_MODELVIEW);		// indicate we are specifying camera transformations
    glLoadIdentity ();     // make sure transformation is "zero'd"
    camera();
    
    /* Allow particles to blend with each other. */
    glDepthMask(GL_TRUE);
    
    glPointSize(8.0);
    
    // Render all particles
    glBegin(GL_POINTS);
    for (int i = 0; i < particles.size(); i++) {
        glColor3f(0,0,1);
        glVertex3f(particles[i].pos.x, particles[i].pos.y, particles[i].pos.z);
    }
    glEnd();
    
    simulator.simulate();
    
    glFlush();
    glutSwapBuffers();					// swap buffers (we earlier set double buffer)
    
    //cout << "zoom " << zoom << " xrot " << xrot << " yrot " << yrot << " zrot " << zrot << "\n";
}

void setupParticles() {
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 10; j++) {
            for (int k = 0; k < 10; k++) {
                Particle p(vec3(i+(rand() % 99) * 0.01,j+(rand() % 99) * 0.01,k+(rand() % 99) * 0.01),vec3(0,0,0),vec3(1,0,0),vec3(1,0,0),1,1);
                Particle p2(vec3(i + (rand() % 99) * 0.01,j + (rand() % 99) * 0.01,k+(rand() % 99) * 0.01),vec3(0,0,0),vec3(1,0,0),vec3(1,0,0),1,1);
                particles.push_back(p);
                particles.push_back(p2);
            }
        }
    }
}

// functions to do stuff with key presses
void keyboard(unsigned char key, int x, int y) {
    switch (key) {
        case '=':
            zoom+=1;
            break;
        case '-':
            zoom-=1;
            break;
    }
}

void special(int key, int x, int y)
{
    int mod = glutGetModifiers();
    
    switch (key) {
        case KB_SPACE:
            exit(0);
            break; 
        case GLUT_KEY_RIGHT:
            if (mod == GLUT_ACTIVE_SHIFT) {
                xtrans+=.5;
            } else {
                yrot+=5;
            }
            break;
        case GLUT_KEY_LEFT:
            if (mod == GLUT_ACTIVE_SHIFT) {
                xtrans-=.5;
            } else {
                yrot-=5;
            }
            break;
        case GLUT_KEY_UP:
            if (mod == GLUT_ACTIVE_SHIFT) {
                ytrans+=.5;
            } else {
                xrot+=5;
            }
            break;
        case GLUT_KEY_DOWN:
            if (mod == GLUT_ACTIVE_SHIFT) {
                ytrans-=.5;
            } else {
                xrot-=5;
            }
            break;
    }
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
    gridX = 10, gridY = 10, gridZ = 10;
    cellSize = 1;
    
    setupParticles();
    // Cutty scene set up, later make robust function
//    Particle particle1(vec3(1,0,0),vec3(0,0,0),vec3(1,0,0),vec3(1,0,0),3.0,3.0);
//    particles.push_back(particle1);
//    Particle particle2(vec3(0,0,1),vec3(0,0,0),vec3(100,100,100),vec3(100,100,100),4.0,4.0);
//    particles.push_back(particle2);
//    Particle particle3(vec3(0,1,0),vec3(0,0,0),vec3(500,500,0),vec3(500,500,0),7.0,7.0);
//    particles.push_back(particle3);
//    Particle particle4(vec3(1,1,7),vec3(0,0,0),vec3(0,300,700),vec3(0,300,700),10.0,10.0);
//    particles.push_back(particle4);
    
    initScene();							// quick function to set up scene
    
    // Simulates the movement of all particles
    // at current time step
//    simulator.simulate();
//    cout << "particle1 " << particles[0].pos.x << " " << particles[0].pos.y << " " << particles[0].pos.z << " " << particles[0].vel.x << " " << particles[0].vel.y << " " << particles[0].vel.z << "\n";
//    cout << "particle2 " << particles[1].pos.x << " " << particles[1].pos.y << " " << particles[1].pos.z << " " << particles[1].vel.x << " " << particles[1].vel.y << " " << particles[1].vel.z << "\n"; 
//    cout << "particle3 " << particles[2].pos.x << " " << particles[2].pos.y << " " << particles[2].pos.z << " " << particles[2].vel.x << " " << particles[2].vel.y << " " << particles[2].vel.z << "\n"; 
//    cout << "particle4 " << particles[3].pos.x << " " << particles[3].pos.y << " " << particles[3].pos.z << " " << particles[3].vel.x << " " << particles[3].vel.y << " " << particles[3].vel.z << "\n"; 
    
    glutDisplayFunc(myDisplay);				// function to run when its time to draw something
    glutReshapeFunc(myReshape);				// function to run when the window gets resized
    glutIdleFunc(myDisplay);
    glutKeyboardFunc(keyboard);               // to do key things
    glutSpecialFunc(special);
    
    glutMainLoop();	// infinite loop that will keep drawing and resizing
    
    return 0;
}

