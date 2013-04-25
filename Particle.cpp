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

Particle::Particle(Particle* particle) {
    pos = (*particle).pos;
    vel = (*particle).vel;
    color = (*particle).color;
    force = (*particle).force;
    mass = (*particle).mass;
    den = (*particle).den;
    copy = particle;
}