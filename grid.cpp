//
//  grid.cpp
//  staggered MAC grid that stores GridCells containing pressures, velocities, and particle locations

//  pressures[i][j][k] = pressure[i[j][k]
//  xvelocityOld[i][j][k] = xvel[i-1/2][j][k]
//  yvelocityOld[i][j][k] = yvel[i][j-1/2][k]
//  zvelocityOld[i][j][k] = zvel[i][j][k-1/2]

#include "grid.h"


Grid::Grid(float xdim, float ydim, float zdim, float h) {
    Grid::xdim = xdim;
    Grid::ydim = ydim;
    Grid::zdim = zdim;
    Grid::h = h;
    
    Grid::xcells = (int)xdim/h;
    Grid::ycells = (int)ydim/h;
    Grid::zcells = (int)zdim/h;
    
    // set up each 3d vector and initialize the entries of each cell
    setupVector(pressures, xcells, ycells, zcells);
    setupVector(xvelocityOld, xcells+1, ycells, zcells);
    setupVector(yvelocityOld, xcells, ycells+1, zcells);
    setupVector(zvelocityOld, xcells, ycells, zcells+1);
    
    // set up each 3d vector and initialize the entries of each cell
    setupVector(xvelocityNew, xcells+1, ycells, zcells);
    setupVector(yvelocityNew, xcells, ycells+1, zcells);
    setupVector(zvelocityNew, xcells, ycells, zcells+1);
    
    // set up 3d vector that stores copies of particles
    particleCopies.resize(xcells);
    for (int i = 0; i < xcells; i++) {
        particleCopies[i].resize(ycells);
        for (int j = 0; j < ycells; j++) {
            particleCopies[i][j].resize(zcells);
        }
    }
}

void Grid::setupVector(vector<vector<vector<float> > >& vec, int xsize, int ysize, int zsize) {
    vec.resize(xsize);
    for (int i = 0; i < xsize; i++) {
        vec[i].resize(ysize);
        for (int j = 0; j < ysize; j++) {
            vec[i][j].resize(zsize);
            for (int k = 0; k < zsize; k++) {
                vec[i][j][k] = 0.0f;
            }
        }
    }
}

// clear particleCopies and copy all new particles into it
// store a pointer in each particle copy to its original particle (done in Particle constructor)
void Grid::setupParticleGrid() {
    // clear the grid of old copies of particles
    clearParticleCopies();
    
    for (int i = 0; i < particles.size(); i++) {
        vec3 cell = getCell(particles[i]);
        Particle copy(&particles[i]);
        particleCopies[(int)cell.x][(int)cell.y][(int)cell.z].push_back(copy);
    }
}

// get the cell of a particle at position (x,y,z) in space
// if outside grid, puts in grid
vec3 Grid::getCell(Particle& particle) {
    float xcell = std::min(std::max(0.0f, (float)floor(particle.pos.x/h)), (float)xdim-1.0f);
    float ycell = std::min(std::max(0.0f, (float)floor(particle.pos.y/h)), (float)ydim-1.0f);
    float zcell = std::min(std::max(0.0f, (float)floor(particle.pos.z/h)), (float)zdim-1.0f);
    vec3 cell(xcell, ycell, zcell);
    return cell;
}

// clear the particleCopies grid
void Grid::clearParticleCopies() {
    for (int i = 0; i < xcells; i++) {
        for (int j = 0; j < ycells; j++) {
            for (int k = 0; k < zcells; k++) {
                particleCopies[i][j][k].clear();
            }
        }
    }
}

//vec3 getXPos(int i, int j, int k) {
//    
//}

//vector<Particle> Grid::getNeighbors(int x, int y, int z, float radius) {
//    
//}



