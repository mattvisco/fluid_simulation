//
//  Particle.cpp
//  
//
//  Created by Thomas Rodman Yopes on 4/15/13.
//
//

#include "particle.h"
#include <vector>
#include "glm/glm.hpp"


class Particle;

class Particle {
public:
    vec3 pos, vel, color, force;
    float mass, den, pres;
    Particle(vec3 pos, vec3 vel, vec3, color, vec3 force, float mass, float den, float pres) {Particle::pos=pos; Particle::vel=vel; Particle::color=color; Particle::force=force; Particle::mass=mass; Particle::den=den; Particle::pres=pres;}
    void setPos(vec3);
    void setVel(vec3);
    void setColor(vec3);
    void setForce(vec3);
    void setMass(float);
    void setDen(float)
    void setPres(float);
}

void setPos(vec3 npos){
    pos=npos;
}

void setVel(vec3 nvel){
    vel=nvel;
}

void setColor(vec3 ncolor){
    color=ncolor;
}

void setForce(vec3 nforce){
    force=nforce;
}

void setMass(float nmass){
    mass=nmass;
}

void setDen(float nden){
    den=nden;
}




