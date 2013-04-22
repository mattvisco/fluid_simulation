//
//  particle.h
//  
//
//  not created by Thomas Rodman Yopes on 4/15/13.
//
//

#ifndef ____particle__
#define ____particle__

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include "glm/glm.hpp"

using namespace std;
using namespace glm;

#include "particle.h"
//#include "fluid_simulator.h"
//#include "fluid_render.h"
//#include "grid.h"

class Particle;

class Particle {
public:
    vec3 pos, vel, color, force;
    float mass, den;
    Particle* copy;
    Particle(vec3 pos, vec3 vel, vec3 color, vec3 force, float mass, float den);
    Particle(void) {};
};

#endif 
//defined(____particle__)
