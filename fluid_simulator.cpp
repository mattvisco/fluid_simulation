//
//  fluid_simulator.cpp
//  
//
//  Created by Matthew Visco on 4/15/13.
//
//

#include "fluid_simulator.h"


inline float sqr(float x) { return x*x; }

//using namespace std;
//using namespace glm;

Simulator::Simulator(vector<Particle>& particles) {
    Simulator::particles = particles;
    grid.setupParticleGrid();
    
}

void Simulator::simulate() {
//    for (vector<Particle>::iterator particle = particles->begin(); particle != particles->end(); ++particle) {
//        // may have to call *particle -- not sure yet
//        vector<Particle> neighbors = grid.getNeighbors(*particle);
//        
//    }
}


