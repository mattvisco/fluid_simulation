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

void Grid::computePressure(){
    float deltaT=1.0;
    float density=1.0;
    int size=xcells*ycells*zcells;
    double x[size];
    double b[size];
    vector<float> val;
    vector<int> ival;
    vector<int> rval;
    rval.push_back(0);
    rval.push_back(2);
    ival.push_back(0);
    ival.push_back(1);
    val.push_back(1);
    val.push_back(1);
    for (int i=0;i<xcells;i++){
        for (int j=0; j<ycells; j++){
            for (int k=0; k<zcells; k++){
                int n=0;
                if (i==0){
                    if (particleCopies[i+1][j][k].size()!=0){
                        printf("#1\n");
                        //val.push_back((1/density)*(deltaT)*(-1)/(h*h));
                        val.push_back(1);
                        ival.push_back(k+zcells*j+zcells*ycells*(i+1));
                        n++;
                    }
                } else if (i==xcells-1){
                    if (particleCopies[i-1][j][k].size()!=0){
                        printf("#2\n");
                        //val.push_back((1/density)*(deltaT)*(-1)/(h*h));
                        val.push_back(1);
                        ival.push_back(k+zcells*j+zcells*ycells*(i-1));
                        n++;
                    }
                } else {
                    if (particleCopies[i-1][j][k].size()!=0){
                        printf("#3\n");
                        //val.push_back((1/density)*(deltaT)*(-1)/(h*h));
                        val.push_back(1);
                        ival.push_back(k+zcells*j+zcells*ycells*(i-1));
                        n++;
                    }
                    if (particleCopies[i+1][j][k].size()!=0){
                        printf("#4\n");
                        //val.push_back((1/density)*(deltaT)*(-1)/(h*h));
                        val.push_back(1);
                        ival.push_back(k+zcells*j+zcells*ycells*(i+1));
                        n++;
                    }
                }
                if (j==0){
                    if (particleCopies[i][j+1][k].size()!=0){
                        printf("#5\n");
                        //val.push_back((1/density)*(deltaT)*(-1)/(h*h));
                        val.push_back(1);
                        ival.push_back(k+zcells*(j+1)+zcells*ycells*i);
                        n++;
                    }
                } else if (j==ycells-1){
                    if (particleCopies[i][j-1][k].size()!=0){
                        printf("#6\n");
                        //val.push_back((1/density)*(deltaT)*(-1)/(h*h));
                        val.push_back(1);
                        ival.push_back(k+zcells*(j-1)+zcells*ycells*i);
                        n++;
                    }
                } else {
                    if (particleCopies[i][j-1][k].size()!=0){
                        printf("#7\n");
                        //val.push_back((1/density)*(deltaT)*(-1)/(h*h));
                        val.push_back(1);
                        ival.push_back(k+zcells*(j-1)+zcells*ycells*i);
                        n++;
                    }
                    if (particleCopies[i][j+1][k].size()!=0){
                        printf("#8\n");
                        //val.push_back((1/density)*(deltaT)*(-1)/(h*h));
                        val.push_back(1);
                        ival.push_back(k+zcells*(j+1)+zcells*ycells*i);
                        n++;
                    }
                }
                if (k==0){
                    if (particleCopies[i][j][k+1].size()!=0){
                        printf("#9\n");
                        //val.push_back((1/density)*(deltaT)*(-1)/(h*h));
                        val.push_back(1);
                        ival.push_back(k+1+zcells*j+zcells*ycells*i);
                        n++;
                    }
                } else if (k==zcells-1){
                    if (particleCopies[i][j][k-1].size()!=0){
                        printf("#10\n");
                        //val.push_back((1/density)*(deltaT)*(-1)/(h*h));
                        val.push_back(1);
                        ival.push_back(k-1+zcells*j+zcells*ycells*i);
                        n++;
                    }
                } else {
                    if (particleCopies[i][j][k-1].size()!=0){
                        printf("#11\n");
                        //val.push_back((1/density)*(deltaT)*(-1)/(h*h));
                        val.push_back(1);
                        ival.push_back(k-1+zcells*j+zcells*ycells*i);
                        n++;
                    }
                    if (particleCopies[i][j][k+1].size()!=0){
                        printf("#12\n");
                        //val.push_back((1/density)*(deltaT)*(-1)/(h*h));
                        val.push_back(1);
                        ival.push_back(k+1+zcells*j+zcells*ycells*i);
                        n++;
                    }
                }
                if (particleCopies[i][j][k].size()!=0){
                    //val.push_back(n/(h*h));
                    val.push_back(1);
                    ival.push_back(k+zcells*j+zcells*ycells*i);
                    rval.push_back(n+1);
                } else {
                    rval.push_back(n);
                }
                //b[k+zcells*j+zcells*ycells*i]=(-1)*(xvelocityOld[i+1][j][k]-xvelocityOld[i][j][k]+yvelocityOld[i][j+1][k]-yvelocityOld[i][j][k]+zvelocityOld[i][j][k+1]-zvelocityOld[i][j][k])/h;
                if (k+zcells*j+zcells*ycells*i==0 || k+zcells*j+zcells*ycells*i==1){
                    b[k+zcells*j+zcells*ycells*i]=1;
                } else {
                    b[k+zcells*j+zcells*ycells*i]=0;
                }
            }
        }
    }
    
    int it=150;
    int tol=1.e-6;
    int c=val.size();
    double value[c];
    int ivalue[c];
    int pvalue[rval.size()];
    for (int g=0;g<c;g++){
        value[g]=val[g];
        ivalue[g]=ival[g];
    }
    for (int h=0;h<rval.size();h++){
        pvalue[h]=rval[h];
    }
    
    int status;
    //int n=5;
    //int Ap[]={0,2,5,9, 10, 12} ;
    //int Ai[]={0, 1, 0, 2, 4, 1, 2, 3, 4, 2, 1, 4}; double Ax [ ] = {2., 3., 3., -1., 4., 4., -3., 1., 2., 2., 6., 1.} ; double b [ ] = {8., 45., -3., 3., 19.} ;
    //double x [5] ;
    double *null = (double *) NULL ;
    int i ;
    int *Nr, *Ni, *P, *Q;
    double *Nv;
    void *Symbolic, *Numeric ;
    (void) umfpack_di_symbolic (size, size, pvalue, ivalue, value, &Symbolic, null, null) ;
    printf("Checkpoint1\n");
    (void) umfpack_di_numeric (pvalue, ivalue, value, Symbolic, &Numeric, null, null) ;
    printf("Checkpoint2\n");
    umfpack_di_free_symbolic (&Symbolic) ;
    printf("Checkpoint3\n");
    status = umfpack_di_solve (UMFPACK_At, pvalue, ivalue, value, x, b, Numeric, null, null) ;
    if (status==UMFPACK_OK){
        printf("yes\n");
    } else {
        printf("%i\n",status);
    }
    printf("Checkpoint4\n");
    //umfpack_di_free_numeric (&Numeric) ;
    for (i = 0 ; i < size ; i++) printf ("x [%d] = %g\n", i, x [i]) ;
    /*
    int value={
    double *null = (double *) NULL;
    int i;
    void *Symbolic, *Numeric;
    int *Nr, *Ni, *P, *Q;
    double *Nv;
    int status;
     */
    //int status = umfpack_di_transpose(size,size,rvalue,ivalue,value,P,Q,Nr,Ni,Nv);
    //if (status==UMFPACK_OK){
    //    printf("YES");
    //} else {
    //    printf("NO\n %i", status);
    //    exit(0);
    //}
    //status = umfpack_di_symbolic(size, size, Nr, Ni, Nv, &Symbolic, null, null);
    //status = umfpack_di_numeric(rvalue, ivalue, value, Symbolic, &Numeric, null, null);
    //umfpack_di_free_symbolic (&Symbolic);
    //(void) umfpack_di_solve(UMFPACK_A, rvalue, ivalue, value, x, b, Numeric, null, null);
    //umfpack_di_free_numeric(&Numeric);
    /*
    for (int ii=0;ii<xcells;ii++){
        for (int jj=0;jj<ycells;jj++){
            for (int kk=0;kk<zcells;kk++){
                pressures[ii][jj][kk]=x[kk+zcells*jj+zcells*ycells*ii];
            }
        }
    }
     */
     
     
}
/*
vector<Particle> Grid::getNeighbors(Particle p) {
    int i=0;
    vec3 cell = getCell(particles[i].pos[0], particles[i].pos[1], particles[i].pos[2]);
vec3 getXPos(int i, int j, int k) {
    
}
    
*/





