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

#include "particle.h"

#define X_AXIS 0
#define Y_AXIS 1
#define Z_AXIS 2
#define GRAVITY -9.81f // meters/second^2
#define DENSITY 1.0f
#define KCFL 1 // constant for CFL condition for timestep
#define DAMPENING 0.4f

using namespace std;
using namespace glm;

class Grid {
public:
    Grid (float, float, float, float);
    Grid (void) {};
    float h; // size of cell
    float xdim; // in real world size
    float ydim;
    float zdim;
    int xcells; // number of cells
    int ycells;
    int zcells;
    float timeStep; // in milliseconds ?
    float maxVelocity;
    vector<Particle>* particles;
    vector<vector<vector<float> > > pressures; 
    vector<vector<vector<float> > > xvelocityOld; 
    vector<vector<vector<float> > > yvelocityOld; 
    vector<vector<vector<float> > > zvelocityOld; 
    vector<vector<vector<float> > > xvelocityNew; // new velocities store differences
    vector<vector<vector<float> > > yvelocityNew;
    vector<vector<vector<float> > > zvelocityNew;
    vector<vector<vector<vector<Particle> > > > particleCopies;
    void setParticles(vector<Particle>*);
    void setupVector(vector<vector<vector<float> > >&, int, int, int);
    void clearParticleCopies();
    void setupParticleGrid();
    vector<Particle> getNeighbors(Particle);
    vec3 getCell(Particle&);
    void clearGrid();
    vector<Particle> getNeighbors(float,float,float,float);
    vector<vector<Particle> > getCellNeighbors(float,float,float);
    float distance(vec3,vec3);
    void storeOldVelocities();
    float weightedAverage(vector<Particle>,vec3,int);
    void computeNonAdvection();
    float computeGravityToAdd();
    void computeTimeStep();
    float computePressureToAdd(int,int,int,int);
    vec3 getInterpolatedVelocityDifference(vec3);
    vec3 getInterpolatedVelocity(vec3);
    float getInterpolatedValue(float,float,float,vector<vector<vector<float> > >);
    void updateParticleVels();
};

#endif 
//defined(____grid__)
