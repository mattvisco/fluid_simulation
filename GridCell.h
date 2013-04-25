//
//  GridCell.h
//  
//
//  Created by Justin Kay on 4/18/13.
//  Copyright (c) 2013 UC Berkeley. All rights reserved.
//

#ifndef _GridCell_h
#define _GridCell_h
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iostream>
#include "glm/glm.hpp"

#include "particle.h"
#include "GridCell.h"


class GridCell {
public:
    vector<Particle> particles;
    vec3 velocity;
    vec3 velocityDifferences;
    float pressure;
    float h; //size of cell dimension
};


#endif
