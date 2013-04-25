//
//  fluid_simulator.cpp
//  
//
//
//

#include "fluid_simulator.h"


inline float sqr(float x) { return x*x; }

//using namespace std;
//using namespace glm;

Simulator::Simulator(vector<Particle>* particles, float gridX, float gridY, float gridZ, float gridH) {
    grid = Grid(gridX, gridY, gridZ, gridH);
    Simulator::particles = particles;
    grid.setParticles(particles);    
}

void Simulator::simulate() {
    grid.setupParticleGrid(); // puts particles in correct grid cells
    grid.storeOldVelocities(); // stores the weighted averages of particles at grid points
    grid.computeTimeStep();
    grid.computeNonAdvection();
    grid.updateParticleVels();
    moveParticles();
}

// move the particles using rk2 on the velocity field
void Simulator::moveParticles() {
    for (int i = 0; i < (*particles).size(); i++) {
//        (*particles)[i].pos += grid.timeStep*(grid.getInterpolatedVelocity((*particles)[i].pos + grid.timeStep/2*(grid.getInterpolatedVelocity((*particles)[i].pos)+grid.getInterpolatedVelocityDifference((*particles)[i].pos)))+grid.getInterpolatedVelocityDifference((*particles)[i].pos + grid.timeStep/2*(grid.getInterpolatedVelocity((*particles)[i].pos)+grid.getInterpolatedVelocityDifference((*particles)[i].pos))));
        
        // the sketchy way
        (*particles)[i].pos += grid.timeStep*(*particles)[i].vel;
        
        // check for boundary collisions
        // move the particles in the normal direction outside the solid
        // reverse and dampen normal velocity
        if ((*particles)[i].pos.x > grid.xdim) {
            (*particles)[i].pos.x = grid.xdim;
            (*particles)[i].vel.x *= -1;
            (*particles)[i].vel *= DAMPENING;
        } else if ((*particles)[i].pos.x < 0) {
            (*particles)[i].pos.x = 0;
            (*particles)[i].vel.x *= -1;
            (*particles)[i].vel *= DAMPENING;
        }
        if ((*particles)[i].pos.y > grid.ydim) {
            (*particles)[i].pos.y = grid.ydim;
            (*particles)[i].vel.y *= -1;
            (*particles)[i].vel *= DAMPENING;
        } else if ((*particles)[i].pos.y < 0) {
            (*particles)[i].pos.y = 0;
            (*particles)[i].vel.y *= -1;
            (*particles)[i].vel *= DAMPENING;
        }
        if ((*particles)[i].pos.z > grid.zdim) {
            (*particles)[i].pos.z = grid.zdim;
            (*particles)[i].vel.z *= -1;
            (*particles)[i].vel *= DAMPENING;
        } else if ((*particles)[i].pos.z < 0) {
            (*particles)[i].pos.z = 0;
            (*particles)[i].vel.z *= -1;
            (*particles)[i].vel *= DAMPENING;
        }
    }
}
