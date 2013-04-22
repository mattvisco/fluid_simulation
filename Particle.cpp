//
//  particle.cpp
//  
//
//  Created by Thomas Rodman Yopes on 4/15/13.
//
//

#include "particle.h"

Particle::Particle(vec3 pos, vec3 vel, vec3 color, vec3 force, float mass, float den) {
    Particle::pos=pos;
    Particle::vel=vel;
    Particle::color=color;
    Particle::force=force;
    Particle::mass=mass;
    Particle::den=den;
}

void Particle::setPos(vec3 npos){
    Particle::pos=npos;
}

void Particle::setVel(vec3 nvel){
    Particle::vel=nvel;
}

void Particle::setColor(vec3 ncolor){
    Particle::color=ncolor;
}

void Particle::setForce(vec3 nforce){
    Particle::force=nforce;
}

void Particle::setMass(float nmass){
    Particle::mass=nmass;
}

void Particle::setDen(float nden){
    Particle::den=nden;
}




