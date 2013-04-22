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
    Grid:boundaryC=1.0;
    
    Grid::xcells = (int)xdim/h;
    Grid::ycells = (int)ydim/h;
    Grid::zcells = (int)zdim/h;
    
    // set up each 3d vector and initialize the entries of each cell
    setupVector(pressures&, xcells, ycells, zcells);
    setupVector(xvelocityOld&, xcells+1, ycells, zcells);
    setupVector(yvelocityOld&, xcells, ycells+1, zcells);
    setupVector(zvelocityOld&, xcells, ycells, zcells+1);
    
    // set up each 3d vector and initialize the entries of each cell
    setupVector(xvelocityNew&, xcells+1, ycells, zcells);
    setupVector(yvelocityNew&, xcells, ycells+1, zcells);
    setupVector(zvelocityNew&, xcells, ycells, zcells+1);
    
    // set up 3d vector that stores copies of particles
    particleCopies.resize(xcells);
    for (int i = 0; i < xcells; i++) {
        particleCopies[i].resize(ycells);
        for (int j = 0; j < ycells; j++) {
            particleCopies[i][j].resize(zcells);
        }
    }
}

void Grid::setupVector(&vector<vector<vector<float>>> vec, int xsize, int ysize, int zsize) {
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
        vec3 cell = getCell(particles[i].pos.x, particles[i].pos.y, particles[i].pos.z);
        Particle copy(particles[i]);
        particleCopies[(int)cell.x][(int)cell.y][(int)cell.z].push_back(copy);
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


void computePressure(){
    float deltaT=1.0;
    float density=1.0;
    int size=xcells*ycells*zcells;
    Vector_double x(size,0.0);
    Vector_double b(size,0,0);
    vector<float> val;
    vector<int> ival;
    vector<int> jval;
    for (int i=0;i<xcells;i++){
        for (int j=0; j<ycells; j++){
            for (int k=0; k<zcells; k++){
                int n=0;
                if (i==0){
                    if (particleCopies[i+1][j][k].size()!=0){
                        val.push_back(-1/(h*h));
                        jval.push_back(k+zcells*j+zcells*ycells*(i+1));
                        ival.push_back(k+zcells*j+zcells*ycells*i);
                        n++;
                    }
                } else if (i==xcells-1){
                    if (particleCopies[i-1][j][k].size()!=0){
                        val.push_back(-1/(h*h));
                        jval.push_back(k+zcells*j+zcells*ycells*(i-1));
                        ival.push_back(k+zcells*j+zcells*ycells*i);
                        n++;
                    }
                } else {
                    if (particleCopies[i-1][j][k].size()!=0){
                        val.push_back(-1/(h*h));
                        jval.push_back(k+zcells*j+zcells*ycells*(i-1));
                        ival.push_back(k+zcells*j+zcells*ycells*i);
                        n++;
                    }
                    if (particleCopies[i+1][j][k].size()!=0){
                        val.push_back(-1/(h*h));
                        jval.push_back(k+zcells*j+zcells*ycells*(i+1));
                        ival.push_back(k+zcells*j+zcells*ycells*i);
                        n++;
                    }
                }
                if (j==0){
                    if (particleCopies[i][j+1][k].size()!=0){
                        val.push_back(-1/(h*h));
                        jval.push_back(k+zcells*(j+1)+zcells*ycells*i);
                        ival.push_back(k+zcells*j+zcells*ycells*i);
                        n++;
                    }
                } else if (j==ycells-1){
                    if (particleCopies[i][j-1][k].size()!=0){
                        val.push_back(-1/(h*h));
                        jval.push_back(k+zcells*(j-1)+zcells*ycells*i);
                        ival.push_back(k+zcells*j+zcells*ycells*i);
                        n++;
                    }
                } else {
                    if (particleCopies[i][j-1][k].size()!=0){
                        val.push_back(-1/(h*h));
                        jval.push_back(k+zcells*(j-1)+zcells*ycells*i);
                        ival.push_back(k+zcells*j+zcells*ycells*i);
                        n++;
                    }
                    if (particleCopies[i][j+1][k].size()!=0){
                        val.push_back(-1/(h*h));
                        jval.push_back(k+zcells*(j+1)+zcells*ycells*i);
                        ival.push_back(k+zcells*j+zcells*ycells*i);
                        n++;
                    }
                }
                if (k==0){
                    if (particleCopies[i][j][k+1].size()!=0){
                        val.push_back(-1/(h*h));
                        jval.push_back(k+1+zcells*j+zcells*ycells*i);
                        ival.push_back(k+zcells*j+zcells*ycells*i);
                        n++;
                    }
                } else if (k==zcells-1){
                    if (particleCopies[i][j][k-1].size()!=0){
                        val.push_back(-1/(h*h));
                        jval.push_back(k-1+zcells*j+zcells*ycells*i);
                        ival.push_back(k+zcells*j+zcells*ycells*i);
                        n++;
                    }
                } else {
                    if (particleCopies[i][j][k-1].size()!=0){
                        val.push_back(-1/(h*h));
                        jval.push_back(k-1+zcells*j+zcells*ycells*i);
                        ival.push_back(k+zcells*j+zcells*ycells*i);
                        n++;
                    }
                    if (particleCopies[i][j][k+1].size()!=0){
                        val.push_back(-1/(h*h));
                        jval.push_back(k+1+zcells*j+zcells*ycells*i);
                        ival.push_back(k+zcells*j+zcells*ycells*i);
                        n++;
                    }
                }
                if (particleCopies[i][j][k].size()!=0){
                    val.push_back(n/(h*h));
                    ival.push_back(k+zcells*j+zcells*ycells*i);
                    jval.push_back(k+zcells*j+zcells*ycells*i);
                }
                b[k+zcells*j+zcells*ycells*i]=(-1)*(xvelocityOld[i+1][j][k]-xvelocityOld[i][j][k]+yvelocityOld[i][j+1][k]-yvelocityOld[i][j][k]+zvelocityOld[i][j][k+1]-zvelocityOld[i][j][k])/h;
            }
        }
    }
    Coord_Mat_double A(size, size, val.size(), val, ival, jval);
    ICPreconditioner D(A);
    int result = CG(A,x,b,D,150,1e-6);
    for (int ii=0;ii<xcells;ii++){
        for (int jj=0;jj<ycells;jj++){
            for (int kk=0;kk<zells;kk++){
                pressures[ii][jj][kk]=(1/density)*(deltaT)*x[kk+zcells*jj+zcells*ycells*ii];
            }
        }
    }
}

vector<Particle> Grid::getNeighbors(Particle p) {
    int i=0;
    vec3 cell = getCell(particles[i].pos[0], particles[i].pos[1], particles[i].pos[2]);
vec3 getXPos(int i, int j, int k) {
    
}
    





