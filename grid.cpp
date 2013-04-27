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
    
    // sets up two 3d vectors that stores respectively
    // copies of particles and what component is in grid
    // also initializes layers vector for velocity extrapolation
    particleCopies.resize(xcells);
    gridComponents.resize(xcells);
    layers.resize(xcells);
    for (int i = 0; i < xcells; i++) {
        particleCopies[i].resize(ycells);
        gridComponents[i].resize(ycells);
        layers[i].resize(ycells);
        for (int j = 0; j < ycells; j++) {
            particleCopies[i][j].resize(zcells);
            gridComponents[i][j].resize(zcells);
            layers[i][j].resize(zcells);
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
// default is AIR
void Grid::clearParticleCopies() {
    for (int i = 0; i < xcells; i++) {
        for (int j = 0; j < ycells; j++) {
            for (int k = 0; k < zcells; k++) {
                particleCopies[i][j][k].clear();
                gridComponents[i][j][k] = AIR;
            }
        }
    }
}

// get a vector of coordinates of fluid cells
vector<vec3> Grid::getFluids() {
    vector<vec3> fluids;
    for (int i = 0; i < xcells; i++) {
        for (int j = 0; j < ycells; j++) {
            for (int k = 0; k < zcells; k++) {
                if (gridComponents[i][j][k] == FLUID)
                    fluids.push_back(vec3((float)i,(float)j,(float)k));
            }
        }
    }
    return fluids;
}

bool Grid::isNeighbor(vec3 p1, vec3 p2) {
    int i1 = (int)p1.x;
    int j1 = (int)p1.y;
    int k1 = (int)p1.z;
    int i2 = (int)p2.x;
    int j2 = (int)p2.y;
    int k2 = (int)p2.z;
    
    bool left = (i1 == i2-1) && (j1 == j2) && (k1 == k2);
    bool right = (i1 == i2+1) && (j1 == j2) && (k1 == k2);
    bool top = (i1 == i2) && (j1 == j2-1) && (k1 == k2);
    bool bottom = (i1 == i2) && (j1 == j2+1) && (k1 == k2);
    bool front = (i1 == i2) && (j1 == j2) && (k1 == k2-1);
    bool back = (i1 == i2) && (j1 == j2) && (k1 == k2+1);    
    return left || right || top || bottom || back || front;
}

int Grid::getNonSolidNeighbors(vec3 p) {
    int i = (int)p.x;
    int j = (int)p.y;
    int k = (int)p.z;
    int nonSolids = 0;
    if (i+1 < xcells && gridComponents[i+1][j][k] != SOLID)
        nonSolids++;
    if (i-1 >= 0 && gridComponents[i-1][j][k] != SOLID)
        nonSolids++;
    if (j+1 < ycells && gridComponents[i][j+1][k] != SOLID)
        nonSolids++;
    if (j-1 >= 0 && gridComponents[i][j-1][k] != SOLID)
        nonSolids++;
    if (k+1 < zcells && gridComponents[i][j][k+1] != SOLID)
        nonSolids++;
    if (k-1 >= 0 && gridComponents[i][j][k-1] != SOLID)
        nonSolids++;
    return nonSolids;
}

int Grid::getAirNeighbors(vec3 p) {
    int i = (int)p.x;
    int j = (int)p.y;
    int k = (int)p.z;
    int airs = 0;
    if (i+1 < xcells && gridComponents[i+1][j][k] == AIR)
        airs++;
    if (i-1 >= 0 && gridComponents[i-1][j][k] == AIR)
        airs++;
    if (j+1 < ycells && gridComponents[i][j+1][k] == AIR)
        airs++;
    if (j-1 >= 0 && gridComponents[i][j-1][k] == AIR)
        airs++;
    if (k+1 < zcells && gridComponents[i][j][k+1] == AIR)
        airs++;
    if (k-1 >= 0 && gridComponents[i][j][k-1] == AIR)
        airs++;
    return airs;
}

vector<vec3> Grid::getvec3Neighbors(vec3 p) {
    int i = (int)p.x;
    int j = (int)p.y;
    int k = (int)p.z;
    vector<vec3> flds;
    if (i+1 < xcells)
        flds.push_back(vec3((float)i+1,(float)j,(float)k));
    if (i-1 >= 0)
        flds.push_back(vec3((float)i-1,(float)j,(float)k));
    if (j+1 < ycells)
        flds.push_back(vec3((float)i,(float)j+1,(float)k));
    if (j-1 >= 0)
        flds.push_back(vec3((float)i,(float)j-1,(float)k));
    if (k+1 < zcells)
        flds.push_back(vec3((float)i,(float)j,(float)k+1));
    if (k-1 >= 0)
        flds.push_back(vec3((float)i,(float)j,(float)k-1));
    return flds;
}

// get divergence of cell u
// not sure if forward difference is right
float Grid::divergence(vec3 u) {
    int i = (int)u.x;
    int j = (int)u.y;
    int k = (int)u.z;
    return xvelocityOld[i+1][j][k]-xvelocityOld[i][j][k] + yvelocityOld[i][j+1][k]-yvelocityOld[i][j][k] + zvelocityOld[i][j][k+1]-zvelocityOld[i][j][k];
}

void Grid::computePressure(){
    vector<vec3> fluids = getFluids();
    int size=fluids.size();
    double x[size];
    double b[size];
    double a[size][size];
    
    // zero out x
    for (int i = 0; i < size; i++) {
        x[i] = 0.0;
    }
        
    // for each row in a, store negative how many neighbors are fluid and 1s in columns of neighbors
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (isNeighbor(fluids[i], fluids[j])) {
                a[i][j] = 1.0;
            } else {
                a[i][j] = 0.0;
            }
        }
        a[i][i] = -1.0*getNonSolidNeighbors(fluids[i]);
    }
    
    // for each entry in b
    for (int i = 0; i < size; i++) {
        b[i] = (DENSITY*h/timeStep)*divergence(fluids[i]) - getAirNeighbors(fluids[i])*PATM;
    }
    
    // set up arrays for sparse matrix solve
    //pvalue - total nonzero column entries so far - initial 0
    //ivalue - row number of values in value
    vector<int> pvaluevect;
    int count = 0;
    pvaluevect.push_back(count);
    vector<double> valuevect;
    vector<int> ivaluevect;
    // go by columns
    for (int j = 0; j < size; j++) {
        int numInCol = 0;
        for (int i = 0; i < size; i++) {
            if (a[i][j] != 0) {
                numInCol++;
                valuevect.push_back(a[i][j]);
                ivaluevect.push_back(i);
            }
        }
        count += numInCol;
        pvaluevect.push_back(count);
    }
    
    // store these in arrays
    int* pvalue = &pvaluevect[0];
    int* ivalue = &ivaluevect[0];
    double* value = &valuevect[0];

    int status;
    double *null = (double *) NULL ;
    void *Symbolic, *Numeric ;
    status = umfpack_di_symbolic (size, size, pvalue, ivalue, value, &Symbolic, null, null) ;
    status = umfpack_di_numeric (pvalue, ivalue, value, Symbolic, &Numeric, null, null);
    umfpack_di_free_symbolic (&Symbolic) ;
    status = umfpack_di_solve (UMFPACK_A, pvalue, ivalue, value, x, b, Numeric, null, null) ;
    umfpack_di_free_numeric (&Numeric) ;
    
    // store new pressure values
    setupVector(pressures, xcells, ycells, zcells);
    for (int s = 0; s < size; s++) {
        int i = (int)fluids[s].x;
        int j = (int)fluids[s].y;
        int k = (int)fluids[s].z;
        pressures[i][j][k] = x[s];
    }
    
    //for (i=0;i<size;i++) {printf("x[%i]= %f\n",i,x[i]);}
}


// get all particles within a cube of size radius of a point x,y,z
// neighbors include the cell we're in
// radius = cell width
vector<Particle> Grid::getNeighbors(float x, float y, float z, float radius) {
    vector<vector<Particle> > cellNeighbors = getCellNeighbors(x,y,z);
    vector<Particle> neighbors;
    for (int i = 0; i < cellNeighbors.size(); i++) {
        for (int j = 0; j < cellNeighbors[i].size(); j++) {
            if (cellNeighbors[i][j].pos.x >= x-radius && cellNeighbors[i][j].pos.x < x+radius && cellNeighbors[i][j].pos.y >= y-radius && cellNeighbors[i][j].pos.y < y+radius && cellNeighbors[i][j].pos.z >= z-radius && cellNeighbors[i][j].pos.z < z+radius) {
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

// at each velocity grid point (in center of cell faces), compute and store a weighted average of particle velocites within a cube of radius h
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
        if (dst != 0) {
            dst = 1.0/dst;
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
                xvelocityNew[i][j][k] += KPRES*computePressureToAdd(i,j,k,X_AXIS);
            }
        }
    }
    // y
    for (int i=0; i < xcells; i++) {
        for (int j=0; j < ycells+1; j++) {
            for (int k=0; k < zcells; k++) {
                yvelocityNew[i][j][k] = 0.0f;
                yvelocityNew[i][j][k] += computeGravityToAdd();
                yvelocityNew[i][j][k] += KPRES*computePressureToAdd(i,j,k,Y_AXIS);
            }
        }
    }
    // z
    for (int i=0; i < xcells; i++) {
        for (int j=0; j < ycells; j++) {
            for (int k=0; k < zcells+1; k++) {
                zvelocityNew[i][j][k] = 0.0f;
                zvelocityNew[i][j][k] += KPRES*computePressureToAdd(i,j,k,Z_AXIS);
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
            if (i-1 < 0) {
                return 0.0f;
                //return -timeStep*1.0f/DENSITY*(pressures[i][j][k])/h;
            } else if (i >= xcells) {
                return 0.0f;
                //return -timeStep*1.0f/DENSITY*(-pressures[i-1][j][k])/h;
            } else {
                return -timeStep*1.0f/DENSITY*(pressures[i][j][k]-pressures[i-1][j][k])/h;
            }
        case Y_AXIS:
            if (j-1 < 0) {
                return 0.0f;
                //return -timeStep*1.0f/DENSITY*(pressures[i][j][k])/h;
            } else if (j >= ycells) {
                return 0.0f;
                //return -timeStep*1.0f/DENSITY*(-pressures[i][j-1][k])/h;
            } else {
                return -timeStep*1.0f/DENSITY*(pressures[i][j][k]-pressures[i][j-1][k])/h;
            }
        case Z_AXIS:
            if (k-1 < 0) {
                return 0.0f;
                //return -timeStep*1.0f/DENSITY*(pressures[i][j][k])/h;
            } else if (k >= zcells) {
                return 0.0f;
                //return -timeStep*1.0f/DENSITY*(-pressures[i][j][k-1])/h;
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

// set layers to 0 for fluid cells, -1 for air cells and solid cells
void Grid::updateLayers() {
    for (int i = 0; i < xcells; i++) {
        for (int j = 0; j < ycells; j++) {
            for (int k = 0; k < zcells; k++) {
                if (gridComponents[i][j][k] == FLUID)
                    layers[i][j][k] = 0;
                else {
                    layers[i][j][k] = -1;
                }
            }
        }
    }
}

float Grid::avgNeighbLayers(vector<vec3> neighbs, int i, int AXIS) {
    float avg = 0.0f;
    int tot = 0;
    for (int n = 0; n < neighbs.size(); n++) {
        if (layers[(int)neighbs[n].x][(int)neighbs[n].y][(int)neighbs[n].z] == i-1) {
            switch (AXIS) {
                case X_AXIS:
                    avg += xvelocityNew[(int)neighbs[n].x][(int)neighbs[n].y][(int)neighbs[n].z];
                    tot++;
                    break;
                case Y_AXIS:
                    avg += yvelocityNew[(int)neighbs[n].x][(int)neighbs[n].y][(int)neighbs[n].z];
                    tot++;
                    break;
                case Z_AXIS:
                    avg += zvelocityNew[(int)neighbs[n].x][(int)neighbs[n].y][(int)neighbs[n].z];
                    tot++;
                    break;
            }
        }
    }
    if (tot > 0) {
        avg = avg / tot;
    }
    return avg;
}

//check if neighbs contains a cell such that n.layer == l-1
bool Grid::hasl1Neighbor(vector<vec3> neighbs, int l) {
    for (int n = 0; n < neighbs.size(); n++) {
        if (layers[(int)neighbs[n].x][(int)neighbs[n].y][(int)neighbs[n].z] == l-1) {
            return true;
        }
    }
    return false;
}

// propogate known fluid velocities into air cells surrounding the fluid
// define a buffer zone of at least 2 cells
void Grid::extrapolateVelocities() {
    updateLayers();
    for (int l = 1; l <= std::max(2, KCFL); l++) {
        for (int i = 0; i < xcells; i++) {
            for (int j = 0; j < ycells; j++) {
                for (int k = 0; k < zcells; k++) {
                    if (layers[i][j][k] == -1) {
                        vector<vec3> neighbs = getvec3Neighbors(vec3((float)i,(float)j,(float)k));
                            if (hasl1Neighbor(neighbs, l)) {
                                // x
                                if (i-1 >= 0 && gridComponents[i-1][j][k] != FLUID) {
                                    xvelocityNew[i][j][k] = avgNeighbLayers(neighbs, l, X_AXIS);
                                }
                                // y
                                if (j-1 >= 0 && gridComponents[i][j-1][k] != FLUID) {
                                    yvelocityNew[i][j][k] = avgNeighbLayers(neighbs, l, Y_AXIS);
                                }
                                // z
                                if (k-1 >= 0 && gridComponents[i][j][k-1] != FLUID) {
                                    zvelocityNew[i][j][k] = avgNeighbLayers(neighbs, l, Z_AXIS);
                                }
                                layers[i][j][k] = l;
                            }
                        }
                    }
                }
            }
        }
}
