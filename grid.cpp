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
    setupVector(pressures, xcells, ycells, zcells);
    setupVector(xvelocityOld, xcells+1, ycells, zcells);
    setupVector(yvelocityOld, xcells, ycells+1, zcells);
    setupVector(zvelocityOld, xcells, ycells, zcells+1);
    
    // set up each 3d vector and initialize the entries of each cell
    setupVector(xvelocityNew, xcells+1, ycells, zcells);
    setupVector(yvelocityNew, xcells, ycells+1, zcells);
    setupVector(zvelocityNew, xcells, ycells, zcells+1);
    
    // sets up two 3d vectors that stores respectively
    // copies of particles and what component is in grid
    particleCopies.resize(xcells);
    gridComponents.resize(xcells);
    for (int i = 0; i < xcells; i++) {
        particleCopies[i].resize(ycells);
        gridComponents[i].resize(ycells);
        for (int j = 0; j < ycells; j++) {
            particleCopies[i][j].resize(zcells);
            gridComponents[i][j].resize(zcells);
        }
    }
}

void Grid::setParticles(vector<Particle>* particles) {
    Grid::particles = particles;
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
// also store the velocity of the particle with the maximum velocity for timestep calculations
void Grid::setupParticleGrid() {
    // clear the grid of old copies of particles
    clearParticleCopies();
    maxVelocity = 0.0f;
    
    for (int i = 0; i < (*particles).size(); i++) {
        vec3 cell = getCell((*particles)[i]);
        Particle copy((*particles)[i]);
        particleCopies[(int)cell.x][(int)cell.y][(int)cell.z].push_back(copy);
        
        // TODO -- Have this also iterate over some solid array
        gridComponents[(int)cell.x][(int)cell.y][(int)cell.z] = FLUID;
        if (length((*particles)[i].vel) > maxVelocity) {
            maxVelocity = length((*particles)[i].vel);
        }
    }
}

// get the cell of a particle at position (x,y,z) in space
// if outside grid, puts in grid
vec3 Grid::getCell(Particle& particle) {
    float xcell = std::min(std::max(0.0f, (float)floor((particle.pos.x)/h)), std::max((float)xcells-1.0f, 0.0f));
    float ycell = std::min(std::max(0.0f, (float)floor((particle.pos.y)/h)), std::max((float)ycells-1.0f, 0.0f));
    float zcell = std::min(std::max(0.0f, (float)floor((particle.pos.z)/h)), std::max((float)zcells-1.0f, 0.0f));
    vec3 cell(xcell, ycell, zcell);
    return cell;
}

// clear the particleCopies and gridComponents grid
void Grid::clearParticleCopies() {
    for (int i = 0; i < xcells; i++) {
        for (int j = 0; j < ycells; j++) {
            for (int k = 0; k < zcells; k++) {
                particleCopies[i][j][k].clear();
                gridComponents[i][j][k] = NONE;
            }
        }
    }
}


void Grid::computePressure(){
    int size=xcells*ycells*zcells;
    double x[size];
    double b[size];
    vector<float> val;
    vector<float> tempval;
    vector<int> ival;
    vector<int> tempival;
    vector<int> rval;
    int colCounter=0;
    rval.push_back(colCounter);
    for (int i=0;i<xcells;i++){
        for (int j=0; j<ycells; j++){
            for (int k=0; k<zcells; k++){
                int n=0;
                int ns=6;
                bool fx, fy, fz;
                if (i>0){
                    if (gridComponents[i-1][j][k]==FLUID){
                        //printf("#1\n");
                        val.push_back((1/DENSITY)*(timeStep)*(-1)/(h*h));
                        ival.push_back(k+zcells*j+zcells*ycells*(i-1));
                        n++;
                    }
                }
                if (j>0){
                    if (gridComponents[i][j-1][k]==FLUID){
                        val.push_back((1/DENSITY)*(timeStep)*(-1)/(h*h));
                        ival.push_back(k+zcells*(j-1)+zcells*ycells*i);
                        n++;
                    }
                    
                }
                if (k>0) {
                    if (gridComponents[i][j][k-1]==FLUID){
                        val.push_back((1/DENSITY)*(timeStep)*(-1)/(h*h));
                        ival.push_back(k-1+zcells*j+zcells*ycells*i);
                        n++;
                    }
                    
                }
                if (k<zcells-1){
                    if (gridComponents[i][j][k+1]==FLUID){
                        tempval.push_back((1/DENSITY)*(timeStep)*(-1)/(h*h));
                        tempival.push_back(k+1+zcells*j+zcells*ycells*i);
                        n++;
                    }
                }
                if (j<ycells-1){
                    if (gridComponents[i][j+1][k]==FLUID){
                        tempval.push_back((1/DENSITY)*(timeStep)*(-1)/(h*h));
                        tempival.push_back(k+zcells*(j+1)+zcells*ycells*i);
                        n++;
                    }
                    
                } if (i<xcells-1){
                    if (gridComponents[i+1][j][k]==FLUID){
                        tempval.push_back((1/DENSITY)*(timeStep)*(-1)/(h*h));
                        tempival.push_back(k+zcells*j+zcells*ycells*(i+1));
                        n++;
                    }
                }
                if (i==0){
                    ns--;
                    if (j==0){
                        ns--;
                        if (k==0){
                            ns--;
                        } else if (k==zcells-1){
                            ns--;
                        }
                    } else if (j==ycells-1){
                        ns--;
                    }
                } else if (i==xcells-1){
                    ns--;
                }
                if (gridComponents[i][j][k]==FLUID){
                    val.push_back((1/DENSITY)*(timeStep)*ns/(h*h));
                    ival.push_back(k+zcells*j+zcells*ycells*i);
                    n++;
                }
                for (int r =0;r<tempval.size();r++){
                    val.push_back(tempval[r]);
                    ival.push_back(tempival[r]);
                }
                tempval.clear();
                tempival.clear();
                colCounter+=n;
                rval.push_back(colCounter);
                b[k+zcells*j+zcells*ycells*i] = (-1)*(xvelocityOld[i+1][j][k]-xvelocityOld[i][j][k]
                                                     +yvelocityOld[i][j+1][k]-yvelocityOld[i][j][k]
                                                     +zvelocityOld[i][j][k+1]-zvelocityOld[i][j][k])/h;
            }
        }
    }
    int c=val.size();
    double value[c];
    int ivalue[c];
    int pvalue[rval.size()];
    for (int g=0;g<ival.size();g++){
        ivalue[g]=ival[g];
        value[g]=1;
    }
    for (int m=0;m<rval.size();m++){
        pvalue[m]=rval[m];
    }
    int status;
    double *null = (double *) NULL ;
    int i,j,k;
    void *Symbolic, *Numeric ;
    status = umfpack_di_symbolic (size, size, pvalue, ivalue, value, &Symbolic, null, null) ;
    status = umfpack_di_numeric (pvalue, ivalue, value, Symbolic, &Numeric, null, null);
    umfpack_di_free_symbolic (&Symbolic) ;
    status = umfpack_di_solve (UMFPACK_At, pvalue, ivalue, value, x, b, Numeric, null, null) ;
    umfpack_di_free_numeric (&Numeric) ;
    for (i = 0; i<xcells; i++){
        for (j = 0; j<ycells; j++){
            for (k = 0; k<zcells; k++){
                if (isfinite(x[k+j*zcells+i*ycells*zcells])){
                    pressures[i][j][k]=x[k+j*zcells+i*ycells*zcells];
                } else {
                    pressures[i][j][k]=0.0;
                }
            }
        }
    }
    for (i=0;i<size;i++) {printf("x[%i]= %f\n",i,x[i]);}
    
}


// get all particles within a radius of a point x,y,z
// neighbors include the cell we're in
// radius = cell width
vector<Particle> Grid::getNeighbors(float x, float y, float z, float radius) {
    vector<vector<Particle> > cellNeighbors = getCellNeighbors(x,y,z);
    vector<Particle> neighbors;
    for (int i = 0; i < cellNeighbors.size(); i++) {
        for (int j = 0; j < cellNeighbors[i].size(); j++) {
            if (distance(cellNeighbors[i][j].pos, vec3(x,y,z)) <= radius) { // < vs. <= ? --- this is a sphere.....
                neighbors.push_back(cellNeighbors[i][j]);
            }
        }
    }
    return neighbors;
}



float Grid::distance(vec3 p1, vec3 p2) {
    return length(p1 - p2);
}


vector<vector<Particle> > Grid::getCellNeighbors(float x, float y, float z) {
    float xcell = std::min(std::max(0.0f, (float)floor((x)/h)), std::max((float)xcells-1.0f, 0.0f));
    float ycell = std::min(std::max(0.0f, (float)floor((y)/h)), std::max((float)ycells-1.0f, 0.0f));
    float zcell = std::min(std::max(0.0f, (float)floor((z)/h)), std::max((float)zcells-1.0f, 0.0f));
    
    float xplus = std::min(xcell + 1.0f, xcells - 1.0f);
    float yplus = std::min(ycell + 1.0f, ycells - 1.0f);
    float zplus = std::min(zcell + 1.0f, zcells - 1.0f);
    
    float xminus = std::max(xcell - 1.0f, 0.0f);
    float yminus = std::max(ycell - 1.0f, 0.0f);
    float zminus = std::max(zcell - 1.0f, 0.0f);
    
    vector<vector<Particle> > cellParticles;
    for (int i = (int)xminus; i <= (int)xplus; i++) {
        for (int j = (int)yminus; j <= (int)yplus; j++) {
            for (int k = (int)zminus; k <= (int)zplus; k++) {
                cellParticles.push_back(particleCopies[i][j][k]);
            }
        }
    }
    return cellParticles;
}

// at each velocity grid point (in center of cell faces), compute and store a weighted average of particle velocites within a **sphere** of radius h
void Grid::storeOldVelocities() {
    // x
    for (int i=0; i < xcells+1; i++) {
        for (int j=0; j < ycells; j++) {
            for (int k=0; k < zcells; k++) {
                vec3 xpt(i-0.5f, (float)j, (float)k);
                vector<Particle> xneighbors = getNeighbors(xpt.x, xpt.y, xpt.z, h);
                xvelocityOld[i][j][k] = weightedAverage(xneighbors, xpt, X_AXIS);
            }
        }
    }
    // y
    for (int i=0; i < xcells; i++) {
        for (int j=0; j < ycells+1; j++) {
            for (int k=0; k < zcells; k++) {
                vec3 ypt((float)i, j-0.5f, (float)k);
                vector<Particle> yneighbors = getNeighbors(ypt.x, ypt.y, ypt.z, h);
                yvelocityOld[i][j][k] = weightedAverage(yneighbors, ypt, Y_AXIS);
            }
        }
    }
    // z
    for (int i=0; i < xcells; i++) {
        for (int j=0; j < ycells; j++) {
            for (int k=0; k < zcells+1; k++) {
                vec3 zpt((float)i, (float)j, k-0.5f);
                vector<Particle> zneighbors = getNeighbors(zpt.x, zpt.y, zpt.z, h);
                zvelocityOld[i][j][k] = weightedAverage(zneighbors, zpt, Z_AXIS);
            }
        }
    }
}

float Grid::weightedAverage(vector<Particle> particles, vec3 pt, int AXIS) {
    float avg = 0.0f;
    float totalDst = 0.0f;
    for (int i = 0; i < particles.size(); i++) {
        float dst = distance(particles[i].pos, pt);
        switch (AXIS) {
            case X_AXIS:
                avg += (float)particles[i].vel.x*dst;
                break;
            case Y_AXIS:
                avg += (float)particles[i].vel.y*dst;
                break;
            case Z_AXIS:
                avg += (float)particles[i].vel.z*dst;
                break;
            default:
                break;
        }
        totalDst += dst;
    }
    if (totalDst > 0.0f) {
        avg /= totalDst;
    }
    return avg;
}

//Do all the non-advection steps of a standard water simulator on the grid.
void Grid::computeNonAdvection() {
    computePressure();
    // x
    for (int i=0; i < xcells+1; i++) {
        for (int j=0; j < ycells; j++) {
            for (int k=0; k < zcells; k++) {
                xvelocityNew[i][j][k] = 0.0f;
                xvelocityNew[i][j][k] += computePressureToAdd(i,j,k,X_AXIS);
            }
        }
    }
    // y
    for (int i=0; i < xcells; i++) {
        for (int j=0; j < ycells+1; j++) {
            for (int k=0; k < zcells; k++) {
                yvelocityNew[i][j][k] = 0.0f;
                yvelocityNew[i][j][k] += computeGravityToAdd();
                yvelocityNew[i][j][k] += computePressureToAdd(i,j,k,Y_AXIS);
            }
        }
    }
    // z
    for (int i=0; i < xcells; i++) {
        for (int j=0; j < ycells; j++) {
            for (int k=0; k < zcells+1; k++) {
                zvelocityNew[i][j][k] = 0.0f;
                zvelocityNew[i][j][k] += computePressureToAdd(i,j,k,Z_AXIS);
            }
        }
    }
}

float Grid::computeGravityToAdd() {
    return timeStep*GRAVITY;
}

// compute the change in velocity due to pressure for a cell face
float Grid::computePressureToAdd(int i, int j, int k, int AXIS) {
    
    switch (AXIS) {
        case X_AXIS:
            if (i-1 < 0 || i >= xcells) {
                return 0.0f;
            } else {
                return -timeStep*1.0f/DENSITY*(pressures[i][j][k]-pressures[i-1][j][k])/h;
            }
        case Y_AXIS:
            if (j-1 < 0 || j >= ycells) {
                return 0.0f;
            } else {
                return -timeStep*1.0f/DENSITY*(pressures[i][j][k]-pressures[i][j-1][k])/h;
            }
        case Z_AXIS:
            if (k-1 < 0 || k >= zcells) {
                return 0.0f;
            } else {
                return -timeStep*1.0f/DENSITY*(pressures[i][j][k]-pressures[i][j][k-1])/h;
            }
        default:
            return 0.0f;
    }
}

void Grid::computeTimeStep() {
    if (maxVelocity > 0) {
        timeStep = std::min(KCFL*h/maxVelocity,0.1f);
    } else {
        timeStep = 0.1;
    }
}

vec3 Grid::getInterpolatedVelocityDifference(vec3 pt) {
    vec3 velDif;
    velDif.x = getInterpolatedValue(pt.x/h, pt.y/h-0.5f, pt.z/h-0.5f, xvelocityNew);
    velDif.y = getInterpolatedValue(pt.x/h-0.5f, pt.y/h, pt.z/h-0.5f, yvelocityNew);
    velDif.z = getInterpolatedValue(pt.x/h-0.5f, pt.y/h-0.5f, pt.z/h, zvelocityNew);
    return velDif;
}

vec3 Grid::getInterpolatedVelocity(vec3 pt) {
    vec3 vel;
    vel.x = getInterpolatedValue(pt.x/h, pt.y/h-0.5f, pt.z/h-0.5f, xvelocityOld);
    vel.y = getInterpolatedValue(pt.x/h-0.5f, pt.y/h, pt.z/h-0.5f, yvelocityOld);
    vel.z = getInterpolatedValue(pt.x/h-0.5f, pt.y/h-0.5f, pt.z/h, zvelocityOld);
    return vel;
}

float Grid::getInterpolatedValue(float x, float y, float z, vector<vector<vector<float> > > values) {
    int i = floor(x);
    int j = floor(y);
    int k = floor(z);
    float value = 0.0f;
    float totalWeight = 0.0f;
    if (i >= 0 && j >= 0 && k >= 0 && i < values.size() && j < values[i].size() && k < values[i][j].size()) {
        value += (i+1-x) * (j+1-y) * (k+1-z) * values[i][j][k];
        totalWeight += (i+1-x) * (j+1-y) * (k+1-z);
    }
    if (i+1 >= 0 && j >= 0 && k >= 0 && i+1 < values.size() && j < values[i+1].size() && k < values[i+1][j].size()) {
        value += (x-i) * (j+1-y) * (k+1-z) * values[i+1][j][k];
        totalWeight += (x-i) * (j+1-y) * (k+1-z);
    }
    if (i >= 0 && j+1 >= 0 && k >= 0 && i < values.size() && j+1 < values[i].size() && k < values[i][j+1].size()) {
        value += (i+1-x) * (y-j) * (k+1-z) * values[i][j+1][k];
        totalWeight += (i+1-x) * (y-j) * (k+1-z);
    }
    if (i+1 >= 0 && j+1 >= 0 && k >= 0 && i+1 < values.size() && j+1 < values[i+1].size() && k < values[i+1][j+1].size()) {
        value += (x-i) * (y-j) * (k+1-z) * values[i+1][j+1][k];
        totalWeight += (x-i) * (y-j) * (k+1-z);
    }
    if (i >= 0 && j >= 0 && k+1 >= 0 && i < values.size() && j < values[i].size() && k+1 < values[i][j].size()) {
        value += (i+1-x) * (j+1-y) * (z-k) * values[i][j][k+1];
        totalWeight += (i+1-x) * (j+1-y) * (z-k);
    }
    if (i+1 >= 0 && j >= 0 && k+1 >= 0 && i+1 < values.size() && j < values[i+1].size() && k+1 < values[i+1][j].size()) {
        value += (x-i) * (j+1-y) * (z-k) * values[i+1][j][k+1];
        totalWeight += (x-i) * (j+1-y) * (z-k);
    }
    if (i >= 0 && j+1 >= 0 && k+1 >= 0 && i < values.size() && j+1 < values[i].size() && k+1 < values[i][j+1].size()) {
        value += (i+1-x) * (y-j) * (z-k) * values[i][j+1][k+1];
        totalWeight += (i+1-x) * (y-j) * (z-k);
    }
    if (i+1 >= 0 && j+1 >= 0 && k+1 >= 0 && i+1 < values.size() && j+1 < values[i+1].size() && k+1 < values[i+1][j+1].size()) {
        value += (x-i) * (y-j) * (z-k) * values[i+1][j+1][k+1];
        totalWeight += (x-i) * (y-j) * (z-k);
    }
    
    if (totalWeight > 0) {
        value = value / totalWeight;
    }
    return value;
}

void Grid::updateParticleVels() {
    for (int i = 0; i < (*particles).size(); i++) {
        (*particles)[i].vel += getInterpolatedVelocityDifference((*particles)[i].pos);
    }
}
