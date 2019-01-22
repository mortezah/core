/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#ifndef CRVEDGERESHAPEHO_H
#define CRVEDGERESHAPEHO_H

#include "crvAdapt.h"

namespace crv {

class EdgeReshapeHO2D
{
};

class EdgeReshapeHO3D
{
  public:
    EdgeReshapeHO3D(ma::Mesh* m, ma::Entity* e, int o) :
      mesh(m), edge(e), p(o)
    {
      ma::Entity* vs[2];
      mesh->getDownward(edge, 0, vs);
      vertA = vs[0];
      vertB = vs[1];
      faceListA.clear();
      faceListB.clear();
    }
    ~EdgeReshapeHO3D() {}

  protected:
    void setFaceLists();
  private:
    ma::Mesh* mesh;
    ma::Entity* edge;
    int p;
    ma::Entity* vertA;
    ma::Entity* vertB;
    std::vector<ma::Entity*> faceListA;
    std::vector<ma::Entity*> faceListB;
};

}

#endif
