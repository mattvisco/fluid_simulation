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

#include "glm/glm.hpp"
#include <time.h>

#define PI 3.14159265
#define epsilon .0001

#include <iostream>
#include <math.h>
#include "fluid_simulator.h"
#include "particle.h"

using namespace std;
using namespace glm;

class Grid {
public:
    float h; // size of cell
    float xdim; // in real world size
    float ydim;
    float zdim;
    int xcells;
    int ycells;
    int zcells;
    vector<Particle> particles;
    vector<vector<vector<vector<GridCell> > > > grid; // 3d vector of grid cells
    Grid (float, float, float, float);
    Grid (void);
    void setParticles(vector<Particle>);
    vector<Particle> getNeighbors(Particle);
protected: //?? maybe private dgaf
    vec3 getCell(float, float, float);
    void clearGrid();
};

#endif /* defined(____grid__) */
