//
//  GridCell.cpp
//  
//
//

#include <iostream>

using namespace std;
using namespace glm;

class GridCell {
public:
    vector<Particle> particles;
    vec3 velocity;
    vec3 velocityDifferences;
    float pressure;
    
};