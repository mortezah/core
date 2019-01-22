/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#include "crvEdgeReshapeHO.h"
#include "crvShape.h"
#include "crvTables.h"
#include <pcu_util.h>

namespace crv {

/* static ma::Entity* getFirstFaceInTetSharingEdge( */
/*   ma::Mesh* m, ma::Entity* t, ma::Entity* e) */
/* { */
/*   ma::Entity* fs[4]; */
/*   m->getDownward(t, 2, fs); */
/*   for (int i = 0; i < 4; i++) { */
/*     ma::Entity* es[3]; */
/*     m->getDownward(fs[i], 1, es); */
/*     if (apf::findIn(es, 3, e) > -1) */
/*       return fs[i]; */
/*   } */
/*   return 0; */
/* } */

static ma::Entity* getOtherFaceInTetSharingEdge(
  ma::Mesh* m, ma::Entity* t, ma::Entity* f, ma::Entity* e)
{
  ma::Entity* fs[4];
  m->getDownward(t, 2, fs);
  for (int i = 0; i < 4; i++) {
    if (fs[i] == f) continue;
    ma::Entity* es[3];
    m->getDownward(fs[i], 1, es);
    if (apf::findIn(es, 3, e) > -1)
      return fs[i];
  }
  return 0;
}

static std::vector<ma::Entity*> getOrderedTets(
    ma::Mesh* m, ma::Entity* e)
{
  std::vector<ma::Entity*> res;
  res.clear();

  // assert edge is an internal edge
  // TODO: complete this

  // get the first face (name it faceStar)
  ma::Entity* faceStar = m->getUp(e, 0);
  ma::Entity* firstTet = m->getUp(faceStar, 0);
  res.push_back(firstTet);

  // get the tet opposing faceStar
  faceStar = getOtherFaceInTetSharingEdge(m, firstTet, faceStar, e);
  ma::Entity* nextTet = 0;
  ma::Entity* tet0 = m->getUp(faceStar, 0);
  if (tet0 != firstTet)
    nextTet = tet0;
  else
    nextTet = m->getUp(faceStar, 1);

  while (nextTet != firstTet)
  {
    res.push_back(nextTet);
    faceStar = getOtherFaceInTetSharingEdge(m, nextTet, faceStar, e);
    PCU_ALWAYS_ASSERT(faceStar);
    ma::Entity* tet0 = m->getUp(faceStar, 0);
    if (tet0 != nextTet)
      nextTet = tet0;
    else
      nextTet = m->getUp(faceStar, 1);
  }

  PCU_ALWAYS_ASSERT(res.size() == m->countUpward(e));
  return res;
}

void EdgeReshapeHO3D::setFaceLists()
{

}

}
