/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#ifndef CRVVERTEXMOVE_H
#define CRVVERTEXMOVE_H

#include "crvAdapt.h"

namespace crv {

class Adapt;

class VertexMove
{
  public:
    void Init(Adapt* a, int t, ma::Vector d);
    bool setSimplex(ma::Entity* s, ma::Entity* e);
    ma::Entity* findTargetVertex();
    bool move();
    bool isValid(ma::Entity* edge);
    void cancel();
    Adapt* adapter;
    int tag;
    ma::Mesh* mesh;
    ma::Entity* simplex;
    ma::Entity* refEdge;
    ma::Vector oldPosition;
    ma::Vector dir;
    ma::Entity* target;
    std::vector<ma::Entity*> upEdges;
    std::vector<ma::Vector> upEdgePositions;
};

}

#endif
