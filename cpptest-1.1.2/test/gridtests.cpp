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
#include "../../grid.h"

using namespace std;

// Tests unconditional fail asserts
//
class GridTestSuite : public Test::Suite
{
public:
    Grid grid;
    vector<Particle> particles;
    float xdim, ydim, zdim, h;
	GridTestSuite()
	{
		TEST_ADD(GridTestSuite::testSetupVector)
	}
    
protected:
    void setup() {
        h = 1;
        xdim = 10;
        ydim = 10;
        zdim = 10;
        Grid grid(xdim, ydim, zdim, h);
    }
    
private:
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
		ts.add(auto_ptr<Test::Suite>(new FailTestSuite));
		ts.add(auto_ptr<Test::Suite>(new CompareTestSuite));
		ts.add(auto_ptr<Test::Suite>(new ThrowTestSuite));
        
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

