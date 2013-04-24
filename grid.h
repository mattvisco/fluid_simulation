//
//  grid.h
//
//

#ifndef ____grid__
#define ____grid__
//#define COMPLEX std::complex<double>

#include <vector>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <cmath>
#include <stdio.h>

#include <sys/time.h>
#include "glm/glm.hpp"
#include <time.h>

#include "UMFPACK/Include/umfpack.h"
#include <Accelerate/Accelerate.h>


#define PI 3.14159265
#define epsilon .0001

#include <iostream>
#include <math.h>


#include "grid.h"
#include "particle.h"


using namespace std;
using namespace glm;

class Grid {
public:
    float h; // size of cell
    float xdim; // in real world size
    float ydim;
    float zdim;
    float boundaryC;
    int xcells; // number of cells
    int ycells;
    int zcells;
    vector<Particle> particles;
    vector<vector<vector<float> > > pressures; 
    vector<vector<vector<float> > > xvelocityOld; 
    vector<vector<vector<float> > > yvelocityOld; 
    vector<vector<vector<float> > > zvelocityOld; 
    vector<vector<vector<float> > > xvelocityNew;
    vector<vector<vector<float> > > yvelocityNew;
    vector<vector<vector<float> > > zvelocityNew;
    vector<vector<vector<vector<Particle> > > > particleCopies;
    void computePressure();
    Grid (float, float, float, float);
    Grid (void) {};
    void setParticles(vector<Particle>);
    void setupVector(vector<vector<vector<float> > >&, int, int, int);
    void clearParticleCopies();
    void setupParticleGrid();
    vector<Particle> getNeighbors(Particle);
    vec3 getCell(Particle&);
    void clearGrid();
        
};

#endif 
//defined(____grid__)
