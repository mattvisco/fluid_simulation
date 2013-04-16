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
all: main 
main: fluid_render.o 
	$(CC) $(CFLAGS) -o fluid_render fluid_render.o $(LDFLAGS) 
fluid_render.o: fluid_render.cpp
	$(CC) $(CFLAGS) -c fluid_render.cpp -o fluid_render.o
clean:
	$(RM) *.o fluid_render