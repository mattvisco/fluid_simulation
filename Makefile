qCC = g++
ifeq ($(shell sw_vers 2>/dev/null | grep Mac | awk '{ print $$2}'),Mac)
	CFLAGS = -g -DGL_GLEXT_PROTOTYPES -I./include/ -I/usr/X11/include -DOSX
	LDFLAGS = -framework GLUT -framework OpenGL \
    	-L"/System/Library/Frameworks/OpenGL.framework/Libraries" \
    	-lGL -lGLU -lm -lstdc++
else
	CFLAGS = -g -DGL_GLEXT_PROTOTYPES -Iglut-3.7.6-bin
	LDFLAGS = -lglut -lGLU
endif
	
RM = /bin/rm -f 
all: clean fluid_render
fluid_render: fluid_render.o fluid_simulator.o grid.o particle.o GridCell.o
	$(CC) $(CFLAGS) -o fluid_render fluid_render.o fluid_simulator.o grid.o GridCell.o particle.o $(LDFLAGS)
fluid_render.o: fluid_render.cpp fluid_render.h fluid_simulator.cpp fluid_simulator.h grid.cpp grid.h GridCell.cpp GridCell.h particle.cpp particle.h 
	$(CC) $(CFLAGS) -c fluid_render.cpp -o fluid_render.o
fluid_simulator.o: fluid_simulator.cpp fluid_simulator.h particle.cpp particle.h grid.cpp grid.h GridCell.cpp GridCell.h
	$(CC) $(CFLAGS) -c fluid_simulator.cpp -o fluid_simulator.o
grid.o: grid.cpp grid.h particle.cpp particle.h GridCell.cpp GridCell.h
	$(CC) $(CFLAGS) -c grid.cpp -o grid.o 
GridCell.o: GridCell.cpp GridCell.h particle.cpp particle.h
	$(CC) $(CFLAGS) -c GridCell.cpp -o GridCell.o
particle.o: particle.cpp particle.h
	$(CC) $(CFLAGS) -c particle.cpp -o particle.o

clean:
	$(RM) *.o fluid_render particle grid fluid_simulator
	