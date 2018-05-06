/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#ifndef CRVEDGERESHAPE_H
#define CRVEDGERESHAPE_H

#include "crvAdapt.h"

namespace crv {

class Adapt;

class EdgeReshape
{
  public:
    void Init(Adapt* a, int t);
    void setSimplex(ma::Entity* e);
    bool reshape();
    void rePosition(ma::Entity* edge);
    bool isValid(ma::Entity* edge);
    void cancel();
    Adapt* adapter;
    int tag;
    ma::Mesh* mesh;
    ma::Entity* simplex;
    ma::Entity* edges[6];
    ma::Entity* candidateEdge;
    int ne;
    std::vector<ma::Vector> oldPositions;
};

}

#endif
