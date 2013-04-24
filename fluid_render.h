//
//  fluid.h
//  
//
//  Created by Matthew Visco/Thomas Yopes/Justin Kay on 4/10/13.
//
//

#ifndef ____fluid_render__
#define ____fluid_render__

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>

#include "fluid_render.h"
#include "particle.h"
#include "fluid_simulator.h"
#include "grid.h"



#ifdef _WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif

#ifdef OSX
#include <GLUT/glut.h>
#include <OpenGL/glu.h>
#else
#include <GL/glut.h>
#include <GL/glu.h>
#endif

#include "glm/glm.hpp"
#include <time.h>
#include <math.h>

#define PI 3.14159265
#define epsilon .0001


#endif 
//defined(____fluid__)
