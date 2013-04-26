//
//  fluid_simulator.cpp
//  
//
//
//

#include "fluid_simulator.h"


inline float sqr(float x) { return x*x; }


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
       //(*particles)[i].pos += grid.timeStep*(grid.getInterpolatedVelocity((*particles)[i].pos + grid.timeStep/2*(grid.getInterpolatedVelocity((*particles)[i].pos)+grid.getInterpolatedVelocityDifference((*particles)[i].pos)))+grid.getInterpolatedVelocityDifference((*particles)[i].pos + grid.timeStep/2*(grid.getInterpolatedVelocity((*particles)[i].pos)+grid.getInterpolatedVelocityDifference((*particles)[i].pos))));
        
        //the kind of sketchy way?
        vec3 oldvel = (*particles)[i].vel - grid.getInterpolatedVelocityDifference((*particles)[i].pos);
        (*particles)[i].pos += grid.timeStep*(oldvel + grid.getInterpolatedVelocityDifference((*particles)[i].pos + (grid.timeStep/2*(*particles)[i].vel)));
        
        //the sketchy way
        //(*particles)[i].pos += grid.timeStep*(*particles)[i].vel;
        
        // check for boundary collisions
        // move the particles in the normal direction outside the solid
        // reverse and dampen normal velocity
        bool dampen = false;
        if ((*particles)[i].pos.x > grid.xdim) {
            (*particles)[i].pos.x = grid.xdim;
            (*particles)[i].vel.x *= -1;
            dampen = true;
        } else if ((*particles)[i].pos.x < 0) {
            (*particles)[i].pos.x = 0;
            (*particles)[i].vel.x *= -1;
            dampen = true;
        }
        if ((*particles)[i].pos.y > grid.ydim) {
            (*particles)[i].pos.y = grid.ydim;
            (*particles)[i].vel.y *= -1;
            dampen = true;
        } else if ((*particles)[i].pos.y < 0) {
            (*particles)[i].pos.y = 0;
            (*particles)[i].vel.y *= -1;
            dampen = true;
        }
        if ((*particles)[i].pos.z > grid.zdim) {
            (*particles)[i].pos.z = grid.zdim;
            (*particles)[i].vel.z *= -1;
            dampen = true;
        } else if ((*particles)[i].pos.z < 0) {
            (*particles)[i].pos.z = 0;
            (*particles)[i].vel.z *= -1;
            dampen = true;
        }
        
        if (dampen) {
            (*particles)[i].vel *= DAMPENING;
        }
    }
}
