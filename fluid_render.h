//
//  fluid.h
//  
//
//  Created by Matthew Visco/Thomas Yopes/Justin Kay on 4/10/13.
//
//

#ifndef ____fluid__
#define ____fluid__

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>

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


#endif /* defined(____fluid__) */
