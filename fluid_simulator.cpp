//
//  fluid_simulator.cpp
//  
//
//
//

#include "fluid_simulator.h"


inline float sqr(float x) { return x*x; }
float xMax = 0;

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
    checkDivergence();
    grid.updateParticleVels();
    moveParticles();
    //cout << "\n";

    //pressureColorMap(); // Testing purposes only
}

// Testing purposes only
void Simulator::checkDivergence() {

    vector<vec3> fluids = grid.getFluids();
    int size=fluids.size();
    int i,j,k;
    for (int s = 0; s < size; s++) {
        i = fluids[s].x;
        j = fluids[s].y;
        k = fluids[s].z;
        
        if (i == 9 && j == 9 && k == 0) {
            //cout << "New divergence" << "\n";
            //cout << i << " " << j << " " << k << " " << grid.divergence(vec3(i,j,k)) << "\n";
        }
            //cout << "x: " << grid.xvelocityNew[i+1][j][k]-grid.xvelocityNew[i][j][k]<< "\n";
        //cout << "y: " << grid.yvelocityNew[i][j+1][k]-grid.yvelocityNew[i][j][k]<< "\n";
        //cout << "z: " << grid.zvelocityNew[i][j][k+1]-grid.zvelocityNew[i][j][k]<< "\n";
    }
}

// Colors particles based on cell pressure
void Simulator::pressureColorMap() {
    for (int i = 0; i < grid.xcells; i++) {
        for (int j = 0; j < grid.ycells; j++) {
            for (int k = 0; k < grid.zcells; k++) {
                
                float pressure = grid.pressures[i][j][k];

                float r = 0,g = 0,b = 0;
                if (pressure < 32) {
                    r = 1;
                    g = 0;
                    b = 1;
                }else if (pressure < 64) {
                    r = 0;
                    g = 1;
                    b = 1;
                } else if (pressure < 96){
                    r = 1;
                    g = 1;
                    b = 0;
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
        
        // move particles through updated velocity field
        (*particles)[i].pos += grid.timeStep*grid.getInterpolatedVelocityNew((*particles)[i].pos + grid.timeStep/2.0f*grid.getInterpolatedVelocityNew((*particles)[i].pos));
        
        //the sketchy way
        //(*particles)[i].pos += grid.timeStep*(*particles)[i].vel;
        
        // check for boundary collisions
        // move the particles in the normal direction outside the solid
        // reverse and dampen normal velocity
        bool dampen = false;
        if ((*particles)[i].pos.x >= grid.xdim) {
            (*particles)[i].pos.x = grid.xdim-0.01;
            (*particles)[i].vel.x *= -1;
            //(*particles)[i].vel.x *= DAMPENING;
            dampen = true;
        } else if ((*particles)[i].pos.x <= 0) {
            (*particles)[i].pos.x = 0.01;
            (*particles)[i].vel.x *= -1;
            //(*particles)[i].vel.x *= DAMPENING;
            dampen = true;
        }
        if ((*particles)[i].pos.y >= grid.ydim) {
            (*particles)[i].pos.y = grid.ydim-0.01;
            (*particles)[i].vel.y *= -1;
            //(*particles)[i].vel.y *= DAMPENING;
            dampen = true;
        } else if ((*particles)[i].pos.y <= 0) {
            (*particles)[i].pos.y = 0.01;
            (*particles)[i].vel.y *= -1;
            //(*particles)[i].vel.y *= DAMPENING;
            dampen = true;
        }
        if ((*particles)[i].pos.z >= grid.zdim) {
            (*particles)[i].pos.z = grid.zdim-0.01;
            (*particles)[i].vel.z *= -1;
            //(*particles)[i].vel.z *= DAMPENING;
            dampen = true;
        } else if ((*particles)[i].pos.z <= 0) {
            (*particles)[i].pos.z = 0.01;
            (*particles)[i].vel.z *= -1;
            //(*particles)[i].vel.z *= DAMPENING;
            dampen = true;
        }
        
//        if (dampen) {
//            (*particles)[i].vel *= DAMPENING;
//        }
    }
}
