/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#include "crvVertexMove.h"
#include "crvShape.h"
#include "crvTables.h"
#include <pcu_util.h>

namespace crv {

static bool shareVert(ma::Mesh* m, ma::Entity* e1, ma::Entity* e2)
{
  ma::Entity* vs1[2];
  ma::Entity* vs2[2];
  m->getDownward(e1, 0, vs1);
  m->getDownward(e2, 0, vs2);

  return (vs1[0] == vs2[0]) ||
	 (vs1[0] == vs2[1]) ||
	 (vs1[1] == vs2[0]) ||
	 (vs1[1] == vs2[1]);
}

static void getJacobianAtEdgeVerts(ma::Mesh* m, ma::Entity* s, ma::Entity* e,
    double& J1, double& J2)
{
  ma::Entity* ev[2];
  m->getDownward(e, 0, ev);

  ma::Entity* sv[4];
  m->getDownward(s, 0, sv);
  apf::MeshElement* me = apf::createMeshElement(m, s);

  {
    ma::Entity* vert = ev[0];
    int idx = apf::findIn(sv, 4, vert);
    ma::Vector xi(0,0,0);
    switch (idx) {
      case 0:
	xi = ma::Vector(0.0, 0.0, 0.0);
	break;
      case 1:
	xi = ma::Vector(1.0, 0.0, 0.0);
	break;
      case 2:
	xi = ma::Vector(0.0, 1.0, 0.0);
	break;
      case 3:
	xi = ma::Vector(0.0, 0.0, 1.0);
	break;
    }
    ma::Matrix J;
    apf::getJacobian(me, xi, J);
    J1 = apf::getJacobianDeterminant(J, 3);
  }
  {
    ma::Entity* vert = ev[1];
    int idx = apf::findIn(sv, 4, vert);
    ma::Vector xi(0,0,0);
    switch (idx) {
      case 0:
	xi = ma::Vector(0.0, 0.0, 0.0);
	break;
      case 1:
	xi = ma::Vector(1.0, 0.0, 0.0);
	break;
      case 2:
	xi = ma::Vector(0.0, 1.0, 0.0);
	break;
      case 3:
	xi = ma::Vector(0.0, 0.0, 1.0);
	break;
    }
    ma::Matrix J;
    apf::getJacobian(me, xi, J);
    J2 = apf::getJacobianDeterminant(J, 3);
  }

  if (J1 < 0 && J2 < 0) {
    for (int i = 0; i < 4; i++) {
      if (isBoundaryEntity(m, sv[i]))
      	printf("vertex %d is a boundary entity \n", i);
      else
      	printf("vertex %d is NOT a boundary entity \n", i);
    }
  }


  apf::destroyMeshElement(me);
}


void VertexMove::Init(Adapt* a, int t, ma::Vector d)
{
  adapter = a;
  mesh = a->mesh;
  tag = t;
  simplex = 0;
  refEdge = 0;
  dir = d;
  target = 0;
  upEdges.clear();
  upEdgePositions.clear();
}

bool VertexMove::setSimplex(ma::Entity* s, ma::Entity* e)
{
  simplex = s;
  refEdge = e;
  double J1, J2;
  getJacobianAtEdgeVerts(mesh, simplex, refEdge, J1, J2);
  if (J1>0 || J2>0)
    return false;
  return true;
}


ma::Entity* VertexMove::findTargetVertex()
{

  double J1, J2;
  getJacobianAtEdgeVerts(mesh, simplex, refEdge, J1, J2);
  PCU_ALWAYS_ASSERT(J1<=0 && J2<=0);

  // find the vert
  ma::Entity* sv[4];
  mesh->getDownward(simplex, 0, sv);
  int n = 0;
  for (int i = 0; i < 4; i++) {
    if (!isBoundaryEntity(mesh, sv[i])) {
      n++;
      target = sv[i];
    }
  }
  PCU_ALWAYS_ASSERT(n == 1);
  PCU_ALWAYS_ASSERT(target);

  oldPosition = ma::getPosition(mesh, target);
}

static ma::Entity* getTetFaceOppositeVert(ma::Mesh* m, ma::Entity* e, ma::Entity* v)
{
  ma::Downward faces;
  ma::Entity* oppositeFace = 0;
  int nDownFaces = m->getDownward(e, 2, faces);
  for (int i = 0; i < nDownFaces; i++) {
    ma::Downward verts;
    int nDownVerts = m->getDownward(faces[i], 0, verts);
    int j;
    for (j = 0; j < nDownVerts; j++) {
      if (v == verts[j])
      	break;
    }
    if (j == nDownVerts)
      oppositeFace = faces[i];
    else
      continue;
  }

  // make sure that oppositeFace is what it is meant to be!
  ma::Downward verts;
  int numDownVerts = m->getDownward(oppositeFace, 0, verts);
  bool flag = true;
  for (int i = 0; i < numDownVerts; i++) {
    if (v == verts[i]){
      flag = false;
      break;
    }
  }

  PCU_ALWAYS_ASSERT(flag);

  return oppositeFace;
}

bool VertexMove::move()
{
  findTargetVertex();
  for (int i = 0; i < mesh->countUpward(target); i++) {
    ma::Entity* upEdge = mesh->getUpward(target, i);
    upEdges.push_back(upEdge);
    upEdgePositions.push_back(ma::getPosition(mesh, upEdge));
  }

  ma::Entity* face = getTetFaceOppositeVert(mesh, simplex, target);
  ma::Entity* fe[3];
  mesh->getDownward(face, 1, fe);
  ma::Vector normal = computeFaceNormalAtEdgeInTet(mesh, simplex, face, fe[0], ma::Matrix(1,0,0,0,1,0,0,0,1))
                     +computeFaceNormalAtEdgeInTet(mesh, simplex, face, fe[1], ma::Matrix(1,0,0,0,1,0,0,0,1))
                     +computeFaceNormalAtEdgeInTet(mesh, simplex, face, fe[2], ma::Matrix(1,0,0,0,1,0,0,0,1));
  normal = normal / normal.getLength();

  // get the physical coord of the face center
  ma::Vector curvedCentroid;
  ma::Vector linearCentroid;

  apf::MeshElement* me = apf::createMeshElement(mesh, face);
  apf::mapLocalToGlobal(me, ma::Vector(1./3., 1./3., 1./3.), curvedCentroid);
  apf::destroyMeshElement(me);

  apf::getLinearCentroid(mesh, face);

  ma::Vector diff = linearCentroid - curvedCentroid;
  PCU_ALWAYS_ASSERT(diff*normal < 0);

  diff = diff * (-1.);

  /* double delta = diff.getLength() / 10; */
  double delta = (diff*normal) / 20;
  ma::Entity* se[6];
  mesh->getDownward(simplex, 1, se);
  for (int i = 0; i < 10; i++) {
    /* ma::Vector vChange = diff * i * delta / diff.getLength(); */
    ma::Vector vChange = normal * i * delta;
    ma::Vector newPosition = oldPosition + vChange;
    mesh->setPoint(target, 0, newPosition);
    ma::Vector eChange1 = vChange * 1.1;
    std::vector<ma::Entity*> es;
    std::vector<double> ls;
    double lMax = 0.;
    for (std::size_t j = 0; j < upEdges.size(); j++) {
      if (findIn(se, 6, upEdges[j]) > -1)
	mesh->setPoint(upEdges[j], 0, upEdgePositions[j] + eChange1);
      else {
      	es.push_back(upEdges[j]);
	ma::Entity* vs[2];
	mesh->getDownward(upEdges[j], 0, vs);
      	double l = (ma::getPosition(mesh, vs[0]) - ma::getPosition(mesh, vs[1])).getLength();
      	if (l > lMax)
      	  lMax = l;
      	ls.push_back(l);
      }
    }

    for (std::size_t j = 0; j < ls.size(); j++) {
      double lp = ls[j]/lMax;
      ma::Vector eChange2 = vChange / 2.25 / (2.0 - lp*lp);
      mesh->setPoint(es[j], 0, ma::getPosition(mesh, es[j]) + eChange2);
    }




    /* double J1, J2; */
    /* getJacobianAtEdgeVerts(mesh, simplex, refEdge, J1, J2); */
    /* if (J1 > 0 || J2 > 0) { */

    if (isValid(refEdge)) {
      printf("Vertex moving was successful\n");
      return true;
    }
  }

  /* printf("Vertex moving was not successful\n"); */
  /* return false; */
  return true;
}



bool VertexMove::isValid(ma::Entity* edge)
{
  apf::Adjacent adjacent;
  mesh->getAdjacent(edge,3,adjacent);
  Quality* qual = makeQuality(mesh, 2);
  for (std::size_t i = 0; i < adjacent.getSize(); ++i){
    if (qual->getQuality(adjacent[i]) < 0){
      delete qual;
      return false;
    }
  }
  delete qual;
  return true;
}


void VertexMove::cancel()
{
  mesh->setPoint(target, 0, oldPosition);
  for (std::size_t i = 0; i < upEdges.size(); i++) {
    mesh->setPoint(upEdges[i], 0, upEdgePositions[i]);
  }
}

}
