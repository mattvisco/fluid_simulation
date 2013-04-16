//
//  fluid_simulator.h
//  
//
//  Created by Matthew Visco on 4/15/13.
//
//

#ifndef ____fluid_simulator__
#define ____fluid_simulator__

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>

#include "particle.h"
#include "grid.h"

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif

#include "glm/glm.hpp"
#include <time.h>
#include <math.h>

#define PI 3.14159265
#define epsilon .0001

#endif /* defined(____fluid_simulator__) */
