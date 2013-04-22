//
//  grid.h
//
//

#ifndef ____grid__
#define ____grid__

#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>

#include <sys/time.h>
#include "SparseLib++/1.7/include/iohb.h"
#include "glm/glm.hpp"
//#include "gmare.h"
//#include "ilupres.h"
//#include "compcol1.h"
//#include "iohb.h"
//#include "vector.h"
//#include "blas1.h"
#include <time.h>

#define PI 3.14159265
#define epsilon .0001

#include <iostream>
#include <math.h>


#include "grid.h"
#include "particle.h"
#include "GridCell.h"

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
    void computePressue();
    Grid (float, float, float, float);
    Grid (void) {};
    void setParticles(vector<Particle>);
    vector<Particle> getNeighbors(Particle);
protected: //?? maybe private dgaf
    vec3 getCell(float, float, float);
    void clearGrid();
        
};

#endif 
//defined(____grid__)
