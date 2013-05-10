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
ofstream outputfile;

int frame_number = 0;
unsigned char* g_video_memory_start = NULL;
unsigned char* g_video_memory_ptr = NULL;
int g_video_seconds_total = 10;
int g_video_fps = 25;

// write into a file
char name[1024];



// Simple init function
//
void initScene(){
    glClearColor (0.0, 0.0, 0.0, 0.0);
    glClearDepth(1.0);
    //Simulator simulator(particles);
    //Particle h;
    simulator = Simulator(&particles, gridX, gridY, gridZ, cellSize);
}

// Reserves memeory space for files
void reserve_video_memory () {
    // 480 MB at 800x800 resolution 230.4 MB at 640x480 resolution
    g_video_memory_ptr = (unsigned char*)malloc (viewport.w * viewport.h * 3 * g_video_fps * g_video_seconds_total);
    g_video_memory_start = g_video_memory_ptr;
}

// Dumps a frame into a png file
void dumpFrames() {
    // save name will have number
    sprintf (name, "video_frame_%03ld.raw", frame_number);
    std::ofstream file;
    file.open (name, std::ios::out | std::ios::binary);
    // natural order is upside down so flip y
    int bytes_in_row = viewport.w * 3;
    int bytes_left = viewport.w * viewport.h * 3;
    while (bytes_left > 0) {
        int start_of_row = bytes_left - bytes_in_row;
        // write the row
        for (int i = 0; i < bytes_in_row; i++) {
            file << g_video_memory_ptr[start_of_row + i];
        }
        bytes_left -= bytes_in_row;
    }
    file.close ();
    // invoke ImageMagick to convert from .raw to .png
    char command[2056];
    sprintf (command, "convert -depth 8 -size %ix%i rgb:video_frame_%03ld.raw video_frame_%03ld.png", viewport.w, viewport.h, frame_number, frame_number);
    printf ("%s\n", command);
    system (command);
    // delete the .raw
    sprintf (command, "rm -f %s", name);
    system (command);
    
}

//orient the camera
void camera() {
    glTranslatef(-11.5+xtrans,-17.5+ytrans,zoom-70.0);
    glRotatef(xrot - 10,1.0,0.0,0.0);
    glRotatef(yrot,0.0,1.0,0.0);
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
    frame_number++;
    
    /* Allow particles to blend with each other. */
    glDepthMask(GL_TRUE);
    
    glPointSize(8.0);
    
    // Render all particles
    glBegin(GL_POINTS);
    for (int i = 0; i < particles.size(); i++) {
        glColor3f(particles[i].color.x,particles[i].color.y,particles[i].color.z);
        //glColor3f(0,0,1);
        glVertex3f(particles[i].pos.x, particles[i].pos.y, particles[i].pos.z);
//        outputfile << "sphere" << i << " = cmds.polySphere(ax=(0,0,0),r=3)\n";
//        //outputfile << particles[i].pos.x << " " << particles[i].pos.y << " " << particles[i].pos.z << "\n";
//        outputfile << "cmds.polyMoveVertex(sphere" << i << ", t=(" << particles[i].pos.x*10 << "," << particles[i].pos.y*10 << "," << particles[i].pos.z*10 << "))\n";
    }
//    outputfile << "mm.eval('renderWindowRender redoPreviousRender renderView')\n";
//    outputfile << "cmds.renderWindowEditor(editor, e=True, writeImage='image" << frame_number << "')\n";
//    for (int s = 0; s < particles.size(); s++) {
//        outputfile << "cmds.delete(sphere" << s << ")\n";
//    }
    //outputfile << "FRAME\n";
    glEnd();
    
    
    simulator.simulate();
    
    // Stores Frames in a pointer to convert to Image later
//    glReadPixels (0, 0, viewport.w, viewport.h, GL_RGB, GL_UNSIGNED_BYTE, g_video_memory_ptr);
//    dumpFrames();

    
    glFlush();
    glutSwapBuffers();					// swap buffers (we earlier set double buffer)
    }

void setupParticles() {
    for (int i = 15; i < 30; i++) {
        for (int j = 15; j < 30; j++) {
            for (int k = 15; k < 30; k++) {
                Particle p(vec3(i+(rand() % 99) * 0.01,j+(rand() % 99) * 0.01,k+(rand() % 99) * 0.01),vec3(0,0,0),vec3(0,0,1),vec3(0,0,1),1,1);
                Particle p2(vec3(i + (rand() % 99) * 0.01,j + (rand() % 99) * 0.01,k+(rand() % 99) * 0.01),vec3(0,0,0),vec3(0,0,1),vec3(1,0,0),1,1);
                Particle p3(vec3(i+(rand() % 99) * 0.01,j+(rand() % 99) * 0.01,k+(rand() % 99) * 0.01),vec3(0,0,0),vec3(0,0,1),vec3(0,0,1),1,1);
                Particle p4(vec3(i + (rand() % 99) * 0.01,j + (rand() % 99) * 0.01,k+(rand() % 99) * 0.01),vec3(0,0,0),vec3(0,0,1),vec3(0,0,1),1,1);
                Particle p5(vec3(i+(rand() % 99) * 0.01,j+(rand() % 99) * 0.01,k+(rand() % 99) * 0.01),vec3(0,0,0),vec3(0,0,1),vec3(0,0,1),1,1);
                Particle p6(vec3(i + (rand() % 99) * 0.01,j + (rand() % 99) * 0.01,k+(rand() % 99) * 0.01),vec3(0,0,0),vec3(0,0,1),vec3(0,0,1),1,1);
                Particle p7(vec3(i+(rand() % 99) * 0.01,j+(rand() % 99) * 0.01,k+(rand() % 99) * 0.01),vec3(0,0,0),vec3(0,0,1),vec3(0,0,1),1,1);
                Particle p8(vec3(i + (rand() % 99) * 0.01,j + (rand() % 99) * 0.01,k+(rand() % 99) * 0.01),vec3(0,0,0),vec3(0,0,1),vec3(0,0,1),1,1);
                particles.push_back(p);
                particles.push_back(p2);
//                particles.push_back(p3);
//                particles.push_back(p4);
//                particles.push_back(p5);
//                particles.push_back(p6);
//                particles.push_back(p7);
//                particles.push_back(p8);
            }
        }
    }
}

// functions to do stuff with key presses
void keyboard(unsigned char key, int x, int y) {
    switch (key) {
        case '=':
            zoom+=1;
            cout << "zoom" << zoom << "\n";
            break;
        case '-':
            zoom-=1;
            cout << "zoom" << zoom << "\n";
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
                cout << "x" << xtrans << "\n";
            } else {
                yrot+=5;
            }
            break;
        case GLUT_KEY_LEFT:
            if (mod == GLUT_ACTIVE_SHIFT) {
                xtrans-=.5;
                cout << "x" << xtrans << "\n";

            } else {
                yrot-=5;
            }
            break;
        case GLUT_KEY_UP:
            if (mod == GLUT_ACTIVE_SHIFT) {
                ytrans+=.5;
                cout << "y" << ytrans << "\n";

            } else {
                xrot+=5;
            }
            break;
        case GLUT_KEY_DOWN:
            if (mod == GLUT_ACTIVE_SHIFT) {
                ytrans-=.5;
                cout << "y" << ytrans << "\n";
            } else {
                xrot-=5;
            }
            break;
    }
}
 
int main(int argc, char *argv[]) {
    
//    outputfile.open ("/Users/mattvisco/Documents/School/CS184/fluid_simulation/PIC_Column.py");
//    outputfile << "import maya.cmds as cmds\n";
//    outputfile << "import maya.mel as mm\n";
//    outputfile << "import os\n";
//    outputfile << "os.chdir('/Users/mattvisco/Documents/School/CS184/fluid_simulation/renderings/PIC_Column')\n";

    //This initializes glut
    glutInit(&argc, argv);
    
    //This tells glut to use a double-buffered window with red, green, and blue channels
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    
    // Initalize theviewport size
    viewport.w = 1000;
    viewport.h = 1000;
    
    //The size and position of the window
    glutInitWindowSize(viewport.w, viewport.h);
    glutInitWindowPosition(0,0);
    glutCreateWindow(argv[0]);
    
    // Currently grid size hard coded in,
    // later parsed from command line
    gridX = 30, gridY = 30, gridZ = 30;
    cellSize = 1;
    
    setupParticles();
    
    //reserve_video_memory ();
    
    initScene();							// quick function to set up scene
    


    glutDisplayFunc(myDisplay);				// function to run when its time to draw something
    glutReshapeFunc(myReshape);				// function to run when the window gets resized
    glutIdleFunc(myDisplay);

    glutKeyboardFunc(keyboard);               // to do key things
    glutSpecialFunc(special);



    
    glutMainLoop();	// infinite loop that will keep drawing and resizing
    
    return 0;
}

