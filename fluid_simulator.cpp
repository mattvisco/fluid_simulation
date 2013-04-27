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
    grid.extrapolateVelocities();
    grid.updateParticleVels();
    moveParticles();
    
    
    // color particles based on cell pressure
    for (int i = 0; i < grid.xcells; i++) {
        for (int j = 0; j < grid.ycells; j++) {
            for (int k = 0; k < grid.zcells; k++) {
                float pressure = grid.pressures[i][j][k];
                float r = 0,g = 0,b = 0;
                if (pressure < .1) {
                    r = 1;
                    g = 1;
                    b = 1;
                } else if (pressure < 31190) {
                    b = pressure / 31190;
                } else if (pressure < 62380) {
                    g = pressure / 62380;
                } else if (pressure < 93570) {
                    r = pressure / 93570;
                }
                vector<Particle> parts = grid.particleCopies[i][j][k];
                for (int s = 0; s < parts.size(); s++) {
                    (*(parts[s].copy)).color = vec3(r,g,b);
                }
            }
        }
    }
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
