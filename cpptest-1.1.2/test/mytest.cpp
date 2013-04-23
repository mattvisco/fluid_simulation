// ---
//
// $Id: mytest.cpp,v 1.5 2008/07/11 16:49:43 hartwork Exp $
//
// CppTest - A C++ Unit Testing Framework
// Copyright (c) 2003 Niklas Lundell
//
// ---
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the
// Free Software Foundation, Inc., 59 Temple Place - Suite 330,
// Boston, MA 02111-1307, USA.
//
// ---
//
// Test program demonstrating all assert types and output handlers.
//
// ---

#include <cstdlib>
#include <cstring>
#include <iostream>

#ifdef _MSC_VER
#pragma warning (disable: 4290)
#endif

#include "../src/cpptest.h"
#include "../../particle.h"
#include "../../grid.h"
#include "../../grid.cpp"
#include "../../Particle.cpp"


using namespace std;

// Tests unconditional fail asserts
//
class GridTestSuite : public Test::Suite
{
public:
    Grid grid;
    vector<Particle> particles, particles2, particles3;
    float xdim, ydim, zdim, h;
	GridTestSuite()
	{
        TEST_ADD(GridTestSuite::testSetup)
		TEST_ADD(GridTestSuite::testSetupVector)
        TEST_ADD(GridTestSuite::testGetCell)
        TEST_ADD(GridTestSuite::testClearParticleCopies)
        TEST_ADD(GridTestSuite::testSetParticles)
        TEST_ADD(GridTestSuite::testSetupParticleGrid)
        TEST_ADD(GridTestSuite::testSetupParticleGrid2)
        TEST_ADD(GridTestSuite::testDistance)
        TEST_ADD(GridTestSuite::testGetNeighbors)
        TEST_ADD(GridTestSuite::testWeightedAverage)
        TEST_ADD(GridTestSuite::testStoreOldVelocities)  
	}
    
protected:
    virtual void setup() {
        grid = Grid(10.0f,10.0f,10.0f,1.0f);
        Particle particle1(vec3(0,0,0),vec3(1,0,0),vec3(1,0,0),vec3(1,0,0),3.0,3.0);
        Particle particle2(vec3(10,10,10),vec3(100,100,100),vec3(100,100,100),vec3(100,100,100),4.0,4.0);
        Particle particle3(vec3(5,5,5),vec3(500,500,0),vec3(500,500,0),vec3(500,500,0),7.0,7.0);
        Particle particle4(vec3(0.5,3.2,7.7),vec3(0,300,700),vec3(0,300,700),vec3(0,300,700),10.0,10.0);
        particles.push_back(particle1);
        particles.push_back(particle2);
        particles.push_back(particle3);
        particles.push_back(particle4);
        
        Particle particle5(vec3(5,5,5),vec3(1,1,1),vec3(1,0,0),vec3(1,0,0),3.0,3.0);
        Particle particle6(vec3(5.4,5,5),vec3(1,0,0),vec3(1,0,0),vec3(1,0,0),3.0,3.0);
        Particle particle7(vec3(5.3,5.4,5.1),vec3(1,0,0),vec3(1,0,0),vec3(1,0,0),3.0,3.0);
        Particle particle8(vec3(5.01,4.5,5.49),vec3(1,0,0),vec3(1,0,0),vec3(1,0,0),3.0,3.0);
        particles2.push_back(particle5);
        particles2.push_back(particle6);
        particles2.push_back(particle7);
        particles2.push_back(particle8);
        
        for(int i = 0; i < grid.xcells; i++) {
            for(int j = 0; j < grid.ycells; j++) {
                for(int k = 0; k < grid.zcells; k++) {
                    Particle justinsucksparticle(vec3(i,j,k),vec3(1,1,1),vec3(1,0,0),vec3(1,0,0),3.0,3.0);
                    particles3.push_back(justinsucksparticle);
                }
            }
        }
        
    }
    
    virtual void tear_down() {
        grid.clearParticleCopies();
        particles.clear();
        particles2.clear();
    }
    
private:
    void testSetup() {
        TEST_ASSERT_MSG(grid.h == 1.0f, "grid h");
        TEST_ASSERT_MSG(grid.xdim == 10.0f, "grid xdim");
        TEST_ASSERT_MSG(grid.ydim == 10.0f, "grid ydim");
        TEST_ASSERT_MSG(grid.zdim == 10.0f, "grid zdim");
        TEST_ASSERT_MSG(grid.xcells == 10, "grid xcells");
        TEST_ASSERT_MSG(grid.ycells == 10, "grid ycells");
        TEST_ASSERT_MSG(grid.zcells == 10, "grid grid.zcells");
    }
    
	void testSetupVector() {
        grid.setupVector(grid.pressures, 10, 10, 10);
        TEST_ASSERT_MSG(grid.pressures.size() == 10, "pressures x10");
        for (int i = 0; i < 10; i++) {
            TEST_ASSERT_MSG(grid.pressures[i].size() == 10, "pressures y10");
            for (int j = 0; j < 10; j++) {
                TEST_ASSERT_MSG(grid.pressures[i][j].size() == 10, "pressures z10");
            }
        }
        
        grid.setupVector(grid.xvelocityOld, 11, 10, 10);
        TEST_ASSERT_MSG(grid.xvelocityOld.size() == 11, "xvelold x11");
        for (int i = 0; i < 10; i++) {
            TEST_ASSERT_MSG(grid.xvelocityOld[i].size() == 10, "xvelold y10");
            for (int j = 0; j < 10; j++) {
                TEST_ASSERT_MSG(grid.xvelocityOld[i][j].size() == 10, "xvelold z10");
            }
        }
        
        grid.setupVector(grid.xvelocityNew, 11, 10, 10);
        TEST_ASSERT_MSG(grid.xvelocityNew.size() == 11, "xvelnew x11");
        for (int i = 0; i < 10; i++) {
            TEST_ASSERT_MSG(grid.xvelocityNew[i].size() == 10, "xvelnew y10");
            for (int j = 0; j < 10; j++) {
                TEST_ASSERT_MSG(grid.xvelocityNew[i][j].size() == 10, "xvelnew z10");
            }
        }
        
        grid.setupVector(grid.yvelocityOld, 10, 11, 10);
        TEST_ASSERT_MSG(grid.yvelocityOld.size() == 10, "yvelold x10");
        for (int i = 0; i < 10; i++) {
            TEST_ASSERT_MSG(grid.yvelocityOld[i].size() == 11, "yvelold y11");
            for (int j = 0; j < 10; j++) {
                TEST_ASSERT_MSG(grid.yvelocityOld[i][j].size() == 10, "yvelold z10");
            }
        }
        
        grid.setupVector(grid.yvelocityNew, 10, 11, 10);
        TEST_ASSERT_MSG(grid.yvelocityNew.size() == 10, "yvelnew x10");
        for (int i = 0; i < 10; i++) {
            TEST_ASSERT_MSG(grid.yvelocityNew[i].size() == 11, "yvelnew y11");
            for (int j = 0; j < 10; j++) {
                TEST_ASSERT_MSG(grid.yvelocityNew[i][j].size() == 10, "yvelnew z10");
            }
        }
        
        grid.setupVector(grid.zvelocityOld, 10, 10, 11);
        TEST_ASSERT_MSG(grid.zvelocityOld.size() == 10, "zvelold x10");
        for (int i = 0; i < 10; i++) {
            TEST_ASSERT_MSG(grid.zvelocityOld[i].size() == 10, "zvelold y10");
            for (int j = 0; j < 10; j++) {
                TEST_ASSERT_MSG(grid.zvelocityOld[i][j].size() == 11, "zvelold z11");
            }
        }
        
        grid.setupVector(grid.zvelocityNew, 10, 10, 11);
        TEST_ASSERT_MSG(grid.zvelocityNew.size() == 10, "zvelnew x10");
        for (int i = 0; i < 10; i++) {
            TEST_ASSERT_MSG(grid.zvelocityNew[i].size() == 10, "zvelnew y10");
            for (int j = 0; j < 10; j++) {
                TEST_ASSERT_MSG(grid.zvelocityNew[i][j].size() == 11, "zvelnew z11");
            }
        }
    }
    
    void testGetCell() {
        Particle particle1(vec3(0,0,0),vec3(1,0,0),vec3(1,0,0),vec3(1,0,0),3.0,3.0);
        vec3 cell1(0.0,0.0,0.0);
        Particle particle2(vec3(10,10,10),vec3(100,100,100),vec3(100,100,100),vec3(100,100,100),4.0,4.0);
        vec3 cell2(9.0,9.0,9.0);
        Particle particle3(vec3(5,5,5),vec3(500,500,0),vec3(500,500,0),vec3(500,500,0),7.0,7.0);
        vec3 cell3(5.0,5.0,5.0);
        Particle particle4(vec3(0.5,3.2,7.7),vec3(0,300,700),vec3(0,300,700),vec3(0,300,700),10.0,10.0);
        vec3 cell4(1.0,3.0,8.0);
        
        vec3 c = grid.getCell(particle1);        
        TEST_ASSERT_MSG(cellsEqual(c, cell1), "getcell1");
        c = grid.getCell(particle2);        
        TEST_ASSERT_MSG(cellsEqual(c, cell2), "getcell2");
        c = grid.getCell(particle3);        
        TEST_ASSERT_MSG(cellsEqual(c, cell3), "getcell3");
        c = grid.getCell(particle4);        
        TEST_ASSERT_MSG(cellsEqual(c, cell4), "getcell4");
    }
    
    bool cellsEqual(const glm::vec3 &vecA, const glm::vec3 &vecB) {
        const double epsilion = 0.0001;
        return    vecA.x == vecB.x && vecA.y == vecB.y && vecA.z == vecB.z;
    }
    
    void testClearParticleCopies() {
        grid.clearParticleCopies();
        for (int i = 0; i < grid.xcells; i++) {
            for (int j = 0; j < grid.ycells; j++) {
                for (int k = 0; k < grid.zcells; k++) {
                    TEST_ASSERT_MSG(grid.particleCopies[i][j][k].size()==0, "clear not 0");
                }
            }
        }
    }
    
    void testSetParticles() {
        grid.setParticles(&particles);
        TEST_ASSERT_MSG((*(grid.particles)).size() == 4, "particles size");
        for (int i = 0; i < (*(grid.particles)).size(); i++) {
            (*(grid.particles))[i].pos.x += 1;
            TEST_ASSERT_MSG(particles[i].pos.x == (*(grid.particles))[i].pos.x, "reset pointer x");
        }
        
    }
    
    //only tests 1 particle per cell
    void testSetupParticleGrid() {
        grid.setParticles(&particles);
        grid.setupParticleGrid();
        for (int i = 0; i < particles.size(); i++) {
            vec3 cell = grid.getCell(particles[i]);
            TEST_ASSERT_MSG(&(grid.particleCopies[(int)cell.x][(int)cell.y][(int)cell.z][0]) != &(particles[i]), "particleCopies contains copies");
            TEST_ASSERT_MSG(grid.particleCopies[(int)cell.x][(int)cell.y][(int)cell.z][0].pos.x == particles[i].pos.x, "particleCopies posx");
            TEST_ASSERT_MSG(grid.particleCopies[(int)cell.x][(int)cell.y][(int)cell.z][0].pos.y == particles[i].pos.y, "particleCopies posy");
            TEST_ASSERT_MSG(grid.particleCopies[(int)cell.x][(int)cell.y][(int)cell.z][0].pos.z == particles[i].pos.z, "particleCopies posz");
        }
    }
    
    //test multiple particles per cell
    void testSetupParticleGrid2() {
        grid.setParticles(&particles2);
        grid.setupParticleGrid();
        TEST_ASSERT_MSG(grid.particleCopies[5][5][5].size() == 4, "4 particles same cell");
    }
    
    void testDistance() {
        vec3 p1(1.0f, 1.0f, 1.0f);
        vec3 p2(2.0f, 1.0f, 1.0f);
        float dist = grid.distance(p1,p2);
        TEST_ASSERT_MSG(dist == 1.0f, "dist fucked");
    }
    
    void testGetNeighbors() {
        grid.setParticles(&particles2);
        grid.setupParticleGrid();
        TEST_ASSERT_MSG(grid.getNeighbors(4.5, 5.0, 5.0, 1.0).size() == 4, "x.5 fucked");
        TEST_ASSERT_MSG(grid.getNeighbors(5.0, 4.5, 5.0, 1.0).size() == 4, "y.5 fucked");
        TEST_ASSERT_MSG(grid.getNeighbors(5.0, 5.0, 4.5, 1.0).size() == 3, "z.5 fucked");
    }
    
    // TODO - dont be baased?
    void testWeightedAverage() {
        vec3 xpt(4.5, 5.0, 5.0);
        vec3 ypt(5.0, 4.5, 5.0);
        vec3 zpt(5.0, 5.0, 4.5);
        float xresult = grid.weightedAverage(particles2, xpt, X_AXIS);
        float yresult = grid.weightedAverage(particles2, ypt, Y_AXIS);
        float zresult = grid.weightedAverage(particles2, zpt, Z_AXIS);
        vector<Particle> nurb;
        float nurbresult = grid.weightedAverage(nurb, xpt, X_AXIS);
    }
    
    void testStoreOldVelocities() {
        grid.setParticles(&particles3);
        grid.setupParticleGrid();
        grid.storeOldVelocities();
        
        // x
        for (int i=0; i < grid.xcells+1; i++) {
            for (int j=0; j < grid.ycells; j++) {
                for (int k=0; k < grid.zcells; k++) {
                    TEST_ASSERT_MSG(grid.xvelocityOld[i][j][k] == 1, "x is 1");
                }
            }
        }
        // y
        for (int i=0; i < grid.xcells; i++) {
            for (int j=0; j < grid.ycells+1; j++) {
                for (int k=0; k < grid.zcells; k++) {
                    TEST_ASSERT_MSG(grid.yvelocityOld[i][j][k] == 1, "y is 1");
                }
            }
        }
        // z
        for (int i=0; i < grid.xcells; i++) {
            for (int j=0; j < grid.ycells; j++) {
                for (int k=0; k < grid.zcells+1; k++) {
                    TEST_ASSERT_MSG(grid.zvelocityOld[i][j][k] == 1, "z is 1");
                }
            }
        }
    }
    
};

enum OutputType
{
	Compiler,
	Html,
	TextTerse,
	TextVerbose
};

static void
usage()
{
	cout << "usage: mytest [MODE]\n"
    << "where MODE may be one of:\n"
    << "  --compiler\n"
    << "  --html\n"
    << "  --text-terse (default)\n"
    << "  --text-verbose\n";
	exit(0);
}

static auto_ptr<Test::Output>
cmdline(int argc, char* argv[])
{
	if (argc > 2)
		usage(); // will not return
	
	Test::Output* output = 0;
	
	if (argc == 1)
		output = new Test::TextOutput(Test::TextOutput::Verbose);
	else
	{
		const char* arg = argv[1];
		if (strcmp(arg, "--compiler") == 0)
			output = new Test::CompilerOutput;
		else if (strcmp(arg, "--html") == 0)
			output =  new Test::HtmlOutput;
		else if (strcmp(arg, "--text-terse") == 0)
			output = new Test::TextOutput(Test::TextOutput::Terse);
		else if (strcmp(arg, "--text-verbose") == 0)
			output = new Test::TextOutput(Test::TextOutput::Verbose);
		else
		{
			cout << "invalid commandline argument: " << arg << endl;
			usage(); // will not return
		}
	}
	
	return auto_ptr<Test::Output>(output);
}

// Main test program
//
int
main(int argc, char* argv[])
{
	try
	{
		// Demonstrates the ability to use multiple test suites
		//
		Test::Suite ts;
		ts.add(auto_ptr<Test::Suite>(new GridTestSuite));
        
		// Run the tests
		//
		auto_ptr<Test::Output> output(cmdline(argc, argv));
		ts.run(*output, true);
        
		Test::HtmlOutput* const html = dynamic_cast<Test::HtmlOutput*>(output.get());
		if (html)
			html->generate(cout, true, "MyTest");
	}
	catch (...)
	{
		cout << "unexpected exception encountered\n";
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

