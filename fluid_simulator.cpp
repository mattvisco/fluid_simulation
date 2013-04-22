//
//  fluid_simulator.cpp
//  
//
//  Created by Matthew Visco on 4/15/13.
//
//

#include "fluid_simulator.h"
//#include "particle.h"
//#incldue "grid.h"

inline float sqr(float x) { return x*x; }

//using namespace std;
//using namespace glm;

Simulator::Simulator(vector<Particle>* particles, Grid grid) {
    Simulator::particles = particles;
    Simulator::curr_grid = grid;
    Simulator::new_grid = grid;
    grid.setParticles(*particles);
    
}
void Simulator::simulate() {
    for (vector<Particle>::iterator particle = particles->begin(); particle != particles->end(); ++particle) {
        // may have to call *particle -- not sure yet
        vector<Particle> neighbors = curr_grid.getNeighbors(*particle);
        
    }
}


