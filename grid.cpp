//
//  grid.cpp
//  
//
//  Created by Justin Kay on 4/15/13.
//  
//

#include "grid.h"

Grid::Grid(float xdim, float ydim, float zdim, float h) {
    Grid::xdim = xdim;
    Grid::ydim = ydim;
    Grid::zdim = zdim;
    Grid::h = h;
    
    xcells = (int)xdim/h;
    ycells = (int)ydim/h;
    zcells = (int)zdim/h;
    
    // set up size of grid and initialize the linked lists for each cell
    grid.resize(xcells);
    for (int i = 0; i < xcells; i++) {
        grid[i].resize(ycells);
        for (int j = 0; j < ycells; j++) {
            grid[i][j].resize(zcells);
        }
    }
}

// put all the particles in the grid
void Grid::setParticles(vector<Particle> particles) {
    Grid::particles = particles;
    
    // clear the grid of old particles - free the allocated memory
    clearGrid();
    
    for (int i = 0; i < particles.size(); i++) {
        vec3 cell = getCell(particles[i].pos[0], particles[i].pos[1], particles[i].pos[2]);
        // store copies for spacial locality - I DONT KNOW IF THIS IS HAPPENING? ??!!?1
        grid[(int)cell[0]][(int)cell[1]][(int)cell[2]].push_back(particles[i]);
    }
}

// get the cell of a particle at position (x,y,z) in space
vec3 Grid::getCell(float x, float y, float z) {
    float xcell = (float)floor(x/h);
    float ycell = (float)floor(y/h);
    float zcell = (float)floor(z/h);
    vec3 cell(xcell, ycell, zcell);
    return cell;
}

void clearGrid() {
    for (int i = 0; i < xcells; i++) {
        for (int j = 0; j < ycells; j++) {
            for (int k = 0; k < zcells; k++) {
                grid[i][j][k].clear();
            }
        }
    }
}

vector<Particle> getNeighbors(Particle p) {
    vec3 cell = getCell(particles[i].pos[0], particles[i].pos[1], particles[i].pos[2]);
}