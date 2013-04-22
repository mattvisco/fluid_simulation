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

#include "fluid_simulator.h"
#include "particle.h"
#include "grid.h"
#include "GridCell.h"


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

using namespace std;
using namespace glm;

class Simulator {
public:
    vector<Particle>* particles;
    Grid curr_grid;
    Grid new_grid;
    Simulator () {};
    Simulator (vector<Particle>*, Grid);
    void simulate();
};

#endif 
//defined(____fluid_simulator__)
