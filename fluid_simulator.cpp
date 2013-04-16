//
//  fluid_simulator.cpp
//  
//
//  Created by Matthew Visco on 4/15/13.
//
//

#include "fluid_simulator.h"

inline float sqr(float x) { return x*x; }

using namespace std;
using namespace glm;

class Simulator {
public:
    Vector<Particle>* particles;
    Grid curr_grid;
    Grid new_grid;
    Simulator (Vector<Particle>*, Grid);
    void simulate();
};
Simulator::Simulator(Vector<Particle>* particles, Grid grid) {
    Simulator::particles = particles;
    Simulator::curr_grid = grid;
    Simulator::new_grid = grid;
    grid.setParticles(particles);
    
}
void Simulator::simulate() {
    for (vector<Particle> particle = particles.begin(); particle != particles.end(); ++particle) {
        // may have to call *particle -- not sure yet
        vector<Particle> neighbors = curr_grid.neighbors(particle);
        
    }
}