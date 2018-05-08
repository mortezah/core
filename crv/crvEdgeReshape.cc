/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#include "crvEdgeReshape.h"
#include "crvShape.h"
#include "crvTables.h"
#include <pcu_util.h>

namespace crv {

static int markEdges(ma::Mesh* m, ma::Entity* e, int tag,
    ma::Entity* edges[6])
{
  if ( tag <= 1 ) // if its valid, or not checked, don't worry about it
    return 0;
  int dim = (tag-2)/6;
  int index = (tag-2) % 6;
  int n = 0;
  int md = m->getDimension();

  switch (dim) {
    case 0:
    {
      // if we have an invalid vertex, operate on its edges
      ma::Downward ed;
      m->getDownward(e,1,ed);
      n = md;
      if(md == 2){
        edges[0] = ed[index];
        edges[1] = ed[(index+2) % 3];
      } else {
        PCU_ALWAYS_ASSERT(index < 4);
        edges[0] = ed[vertEdges[index][0]];
        edges[1] = ed[vertEdges[index][1]];
        edges[2] = ed[vertEdges[index][2]];
      }
    }
    break;
    case 1:
    {
      // if we have a single invalid edge, operate on it
      ma::Downward ed;
      m->getDownward(e,1,ed);
      edges[0] = ed[index];
      n = 1;
    }
    break;
    case 2:
    {
      // if we have an invalid face, operate on its edges
      ma::Downward ed, faces;
      m->getDownward(e,2,faces);
      m->getDownward(faces[index],1,ed);
      n = 3;
      edges[0] = ed[0];
      edges[1] = ed[1];
      edges[2] = ed[2];
    }
    break;
    case 3:
      m->getDownward(e,1,edges);
      n = 6;
      break;
    default:
      fail("invalid quality tag in markEdges\n");
      break;
  }

  return n;
}



void EdgeReshape::Init(Adapt* a, int t, ma::Vector d)
{
  adapter = a;
  mesh = a->mesh;
  tag = t;
  simplex = 0;
  edges[0] = edges[1] = edges[2] = edges[3] = edges[4] = edges[5] = 0;
  refEdge = 0;
  ne = 0;
  oldPositions.clear();
  dir = d;
}

void EdgeReshape::setSimplex(ma::Entity* s, ma::Entity* e)
{
  ne = markEdges(mesh,s,tag,edges);
  simplex = s;
  refEdge = e;
  for (int i = 0; i < ne; i++) {
    ma::Vector pos;
    mesh->getPoint(edges[i],0,pos);
    oldPositions.push_back(pos);
  }
}

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

ma::Entity* EdgeReshape::findCandidateEdge()
{
  ma::Entity* edge = 0;
  double maxCosineAngle = -1.0;
  printf("ne is %d tag is %d\n", ne, tag);
  for (int i = 0; i < ne; i++) {
    if (isBoundaryEntity(mesh, edges[i])) {
      printf("boundary edge\n");
      continue;
    }
    if (!shareVert(mesh, edges[i], refEdge)) {
      printf("edges don't share a vertex\n");
      continue;
    }
    /* if (refEdge == edges[i]) continue; */
    double cosineAngle = apf::computeCosAngleInTet(mesh, simplex, refEdge, edges[i],
    	ma::Matrix(1.,0.,0.,0.,1.,0.,0.,0.,1.));
    printf("cos angle is %f \n", cosineAngle);
    if (cosineAngle > maxCosineAngle) {
      maxCosineAngle = cosineAngle;
      edge = edges[i];
    }
  }
  printf("=========\n");
  /* PCU_ALWAYS_ASSERT(edge); */
  return edge;
}

bool EdgeReshape::reshape()
{
  ma::Entity* edgeToTry = findCandidateEdge();
  if (!edgeToTry)
    return false;
  rePosition(edgeToTry);
  if (isValid(edgeToTry)) {
    printf("++++ edge reshape was successful ++++\n");
    return true;
  }
  else {
    printf("++++ edge reshape wasn't successful ++++\n");
    return false;
  }
}

void EdgeReshape::rePosition(ma::Entity* edge)
{
  ma::Entity* vs[4];
  ma::Entity* es[6];

  ma::Vector pivotPoint;
  ma::Vector edgeVectors[3];
  mesh->getDownward(simplex,0,vs);
  mesh->getDownward(simplex,1,es);

  // pick a pivotVert, the vertex with the worse jacobian determinant
  ma::Entity* pivotVert;
  int pivotIndex;
  {
    apf::MeshElement* me = apf::createMeshElement(mesh,simplex);

    ma::Entity* edgeVerts[2];
    mesh->getDownward(edge,0,edgeVerts);
    apf::Matrix3x3 J;
    pivotIndex = apf::findIn(vs,4,edgeVerts[0]);
    PCU_ALWAYS_ASSERT(pivotIndex >= 0);

    ma::Vector xi = crv::elem_vert_xi[apf::Mesh::TET][pivotIndex];
    apf::getJacobian(me,xi,J);

    double j = apf::getJacobianDeterminant(J,3);
    pivotVert = edgeVerts[0];

    int index = apf::findIn(vs,4,edgeVerts[1]);
    PCU_ALWAYS_ASSERT(index >= 0);
    xi = crv::elem_vert_xi[apf::Mesh::TET][index];
    apf::getJacobian(me,xi,J);
    if (apf::getJacobianDeterminant(J,3) < j){
      pivotVert = edgeVerts[1];
      pivotIndex = index;
    }
    apf::destroyMeshElement(me);
  }

  mesh->getPoint(pivotVert,0,pivotPoint);

  // local, of edges around vert, [0,2]
  int edgeIndex = 0;

  for (int i = 0; i < 3; ++i){
    // theres only one point, so reuse this...
    edgeVectors[i] = ma::getPosition(mesh,es[vertEdges[pivotIndex][i]])
		    - pivotPoint;
    if (es[vertEdges[pivotIndex][i]] == edge)
      edgeIndex = i;
  }

  PCU_ALWAYS_ASSERT(edgeIndex >= 0);

  ma::Entity* edge1 = es[vertEdges[pivotIndex][(edgeIndex+1)%3]];
  ma::Entity* edge2 = es[vertEdges[pivotIndex][(edgeIndex+2)%3]];

  ma::Vector t1 = computeEdgeTangentAtVertex(mesh, edge1, pivotVert, ma::Matrix(1,0,0,0,1,0,0,0,1));
  ma::Vector t2 = computeEdgeTangentAtVertex(mesh, edge2, pivotVert, ma::Matrix(1,0,0,0,1,0,0,0,1));

  ma::Vector normal = apf::cross(t1, t2);
  double length = normal.getLength();
  double validity = edgeVectors[edgeIndex]*normal;

  if(validity > 1e-10)
    return;

  normal = normal/length;

  /* mirror the vector edgeVectors[edgeIndex] with respect to the plane
    * perpendicular to the normal. The parameter alpha scales the normal
    * (to the plane) component of the mirrored vector.
    */
  /* normal = ma::Vector(0,0,0) - dir; */
  /* normal = dir; */
  double alpha = 1.5;
  double delta = dir.getLength()/100;
  int n = 65;

  for (int i = 0; i < 150; i++) {
    ma::Vector newPoint = pivotPoint + edgeVectors[edgeIndex] +
      dir * i * delta / dir.getLength();
    mesh->setPoint(edge,0,newPoint);
    if (isValid(edge)) {
      printf("++++ found a valid reshape ++++\n");
      break;
    }
  }
  /* ma::Vector newPoint = pivotPoint + edgeVectors[edgeIndex] + */
  /*   dir * (1+alpha) * std::abs(edgeVectors[edgeIndex]*normal) * std::abs(dir*normal); */

  /* ma::Vector newPoint = pivotPoint + edgeVectors[edgeIndex] - */
  /*   normal * (normal * edgeVectors[edgeIndex]) * (1 + alpha) / length / length; */

  /* mesh->setPoint(edge,0,newPoint); */

}

bool EdgeReshape::isValid(ma::Entity* edge)
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


void EdgeReshape::cancel()
{
  for (int i = 0; i < ne; i++) {
    mesh->setPoint(edges[i],0,oldPositions[i]);
  }

}

}
