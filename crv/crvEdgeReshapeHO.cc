/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "crvEdgeReshapeHO.h"
#include "crvMath.h"
#include "crvTables.h"
#include "crvBezier.h"
#include <iostream>

namespace crv {

// this can be re-written using "edgeFaces" in crvTables.h
static ma::Entity* getTetOppFaceSharingEdge(
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

static ma::Entity* getSharedEdge(
    ma::Mesh* m, ma::Entity* f0, ma::Entity* f1)
{
  ma::Entity* f0es[3];
  ma::Entity* f1es[3];
  m->getDownward(f0, 1, f0es);
  m->getDownward(f1, 1, f1es);
  for (int i = 0; i < 3; i++) {
    if (apf::findIn(f1es, 3, f0es[i]) > -1)
      return f0es[i];
  }
  return 0;
}

static ma::Entity* getSharedVert(
    ma::Mesh* m, ma::Entity* e0, ma::Entity* e1)
{
  ma::Entity* e0vs[2];
  ma::Entity* e1vs[2];
  m->getDownward(e0, 0, e0vs);
  m->getDownward(e1, 0, e1vs);
  for (int i = 0; i < 2; i++) {
    if (apf::findIn(e1vs, 2, e0vs[i]) > -1)
      return e0vs[i];
  }
  return 0;
}

static ma::Entity* getEdgeOppVert(
    ma::Mesh* m, ma::Entity* e, ma::Entity* v)
{
  ma::Entity* evs[2];
  m->getDownward(e, 0, evs);
  if (v == evs[0])
    return evs[1];
  else
    return evs[0];
  return 0;
}

static int getFaceDir(
    ma::Mesh* m, ma::Entity* f, ma::Entity* v0, ma::Entity* v1, ma::Entity* v2)
{
  ma::Entity* vs[3];
  m->getDownward(f, 0, vs);
  PCU_ALWAYS_ASSERT(apf::findIn(vs, 3, v0) > -1);
  PCU_ALWAYS_ASSERT(apf::findIn(vs, 3, v1) > -1);
  PCU_ALWAYS_ASSERT(apf::findIn(vs, 3, v2) > -1);

  int dir = -1;
  if ((vs[0]==v0) && (vs[1]==v1) && (vs[2]==v2))
    dir = 0;
  else if ((vs[0]==v1) && (vs[1]==v2) && (vs[2]==v0))
    dir = 1;
  else if ((vs[0]==v2) && (vs[1]==v0) && (vs[2]==v1))
    dir = 2;
  else if ((vs[0]==v2) && (vs[1]==v1) && (vs[2]==v0))
    dir = 3;
  else if ((vs[0]==v1) && (vs[1]==v0) && (vs[2]==v2))
    dir = 4;
  else if ((vs[0]==v0) && (vs[1]==v2) && (vs[2]==v1))
    dir = 5;
  else
    dir = 6;
  PCU_ALWAYS_ASSERT((dir >= 0) && (dir <= 5));
  return dir;
}

static apf::Vector3 map0(double t, double s)
{
  return apf::Vector3((1-s)*(1-t), (1-s)*t, 0);
}

static apf::Vector3 map1(double t, double s)
{
  return apf::Vector3((1-s)*t, s, 0);
}

static apf::Vector3 map2(double t, double s)
{
  return apf::Vector3(s, (1-s)*(1-t), 0);
}

static apf::Vector3 map3(double t, double s)
{
  return apf::Vector3((1-s)*(1-t), s, 0);
}

static apf::Vector3 map4(double t, double s)
{
  return apf::Vector3(s, (1-s)*t, 0);
}

static apf::Vector3 map5(double t, double s)
{
  return apf::Vector3((1-s)*t, (1-s)*(1-t), 0);
}

static MapToParent all_permutations[6] =
{map0,
 map1,
 map2,
 map3,
 map4,
 map5
};

static ma::Vector getCenterPoint(std::vector<ma::Vector> poly)
{
  ma::Vector center(0,0,0);
  for (size_t i = 0; i < poly.size(); i++)
    center += poly[i];
  return center/poly.size();

}

static ma::Vector getCrossSection(std::vector<ma::Vector> strand, ma::Vector x)
{
  /* printf("strand is\n"); */
  /* for (size_t i = 0; i < strand.size(); i++) { */
  /*   printf("%f, %f, %f\n", strand[i].x(), */
			   /* strand[i].y(), */
			   /* strand[i].z()); */
  /* } */

  // get the end points and compute the normal
  ma::Vector a = strand[0];
  ma::Vector b = strand[strand.size() - 1];
  ma::Vector normal = b - a;
  normal = normal / normal.getLength();

  PCU_ALWAYS_ASSERT(x*normal > a*normal);
  PCU_ALWAYS_ASSERT(b*normal > x*normal);

  int idx = 0; // holds index of the 1st point in strand that is beyond x
  for (size_t i = 0; i < strand.size(); i++) {
    ma::Vector s = strand[i];
    if (s*normal < x*normal)
      continue;
    else {
      idx = i;
      break;
    }
  }
  PCU_ALWAYS_ASSERT(idx > 0);
  ma::Vector before = strand[idx-1];
  ma::Vector after  = strand[idx];

  ma::Vector p = before + (after - before) *
    (((x - before) * normal) / ((after - before) * normal));

  return p;
}


void CrvEdgeReshapeHO::orderTets()
{
  ma::Entity* currentFace = mesh->getUpward(edge, 0);
  apf::Up up;
  mesh->getUp(currentFace, up);
  PCU_ALWAYS_ASSERT(up.n == 2);
  ma::Entity* firstTet = up.e[0];
  ma::Entity* nextTet  = up.e[1];
  tets.push_back(firstTet);

  while (nextTet != firstTet) {
    tets.push_back(nextTet);
    currentFace = getTetOppFaceSharingEdge(mesh, nextTet, currentFace, edge);
    PCU_ALWAYS_ASSERT(currentFace);
    apf::Up up;
    mesh->getUp(currentFace, up);
    PCU_ALWAYS_ASSERT(up.n == 2);
    if (nextTet != up.e[0])
      nextTet = up.e[0];
    else
      nextTet = up.e[1];
  }
}

void CrvEdgeReshapeHO::computeFaceLists()
{
  // first create the faceLists
  for (size_t i = 0; i < tets.size(); i++) {
    ma::Entity* es[6];
    ma::Entity* fs[4];
    mesh->getDownward(tets[i], 1, es);
    mesh->getDownward(tets[i], 2, fs);
    int idx = apf::findIn(es, 6, edge);
    ma::Entity* f0 = fs[edgeFaces[oppEdges[idx]][0]];
    ma::Entity* f1 = fs[edgeFaces[oppEdges[idx]][1]];
    ma::Entity* f0vs[3];
    ma::Entity* f1vs[3];
    mesh->getDownward(f0, 0, f0vs);
    mesh->getDownward(f1, 0, f1vs);
    if (apf::findIn(f0vs, 3, vertA) > -1) {
      faceListA.push_back(f0);
      faceListB.push_back(f1);
    }
    else {
      faceListA.push_back(f1);
      faceListB.push_back(f0);
    }
  }
  PCU_ALWAYS_ASSERT(faceListA.size() == tets.size());
  PCU_ALWAYS_ASSERT(faceListB.size() == tets.size());

  // then create the direction lists
  ma::Entity* v0 = vertA;
  ma::Entity* v1 = vertB;
  int size = tets.size();
  for (int i = 0; i < size; i++) {
    ma::Entity* currentE =
      getSharedEdge(mesh, faceListA[i], faceListB[i]);
    ma::Entity* nextE =
      getSharedEdge(mesh, faceListA[(i+1)%size], faceListB[(i+1)%size]);
    PCU_ALWAYS_ASSERT(currentE);
    PCU_ALWAYS_ASSERT(nextE);
    ma::Entity* v3 = getSharedVert(mesh, currentE, nextE);
    PCU_ALWAYS_ASSERT(v3);
    ma::Entity* v2 = getEdgeOppVert(mesh, currentE, v3);
    PCU_ALWAYS_ASSERT(v2);
    dirListA.push_back(getFaceDir(mesh, faceListA[i], v0, v2, v3));
    dirListB.push_back(getFaceDir(mesh, faceListB[i], v1, v2, v3));
  }
}

std::vector<ma::Vector> CrvEdgeReshapeHO::getPolygon(ma::Vector x)
{
  // array of edge nodes in parent [0 to 1] coordinates including verts
  std::vector<double> eNodes;
  eNodes.clear();
  eNodes.push_back(0.0);
  for (int i = 0; i < p-1; i++) {
    ma::Vector xi;
    getBezierNodeXi(apf::Mesh::EDGE, p, i, xi);
    eNodes.push_back((xi[0]+1)/2);
  }
  eNodes.push_back(1.0);

  std::vector<ma::Vector> poly;
  for (size_t i = 0; i < faceListA.size(); i++) {
    apf::MeshElement* meA = apf::createMeshElement(mesh, faceListA[i]);
    apf::MeshElement* meB = apf::createMeshElement(mesh, faceListB[i]);
    for (size_t j = 0; j < eNodes.size()-1; j++) {
      // add the coord of vertA
      std::vector<ma::Vector> strand;
      strand.push_back(ma::getPosition(mesh, vertA));
      // add the coords for internal nodes
      double t = eNodes[j];
      // first for faceA
      for (size_t k = eNodes.size()-2; k > 0; k--) {
      	double s = eNodes[k];
	ma::Vector xi = all_permutations[dirListA[i]](t, s);
	ma::Vector pt;
	apf::mapLocalToGlobal(meA, xi, pt);
	strand.push_back(pt);
      }
      // then for faceB
      for (size_t l = 0; l < eNodes.size()-1; l++) {
      	double s = eNodes[l];
	ma::Vector xi = all_permutations[dirListB[i]](t, s);
	ma::Vector pt;
	apf::mapLocalToGlobal(meB, xi, pt);
	strand.push_back(pt);
      }
      // add the coord of vertB
      strand.push_back(ma::getPosition(mesh, vertB));
      PCU_ALWAYS_ASSERT(int(strand.size()) == 2*p+1);
      poly.push_back(getCrossSection(strand, x));

    }
    apf::destroyMeshElement(meA);
    apf::destroyMeshElement(meB);
  }
  return poly;
}

void CrvEdgeReshapeHO::computeInterpolatingPoints()
{
  // setup vertA and vertB
  ma::Entity* vs[2];
  mesh->getDownward(edge, 0, vs);
  vertA = vs[0];
  vertB = vs[1];

  ma::Vector xA = ma::getPosition(mesh, vertA);
  ma::Vector xB = ma::getPosition(mesh, vertB);
  // first do the operations that are necessary for all the nodes
  // like computing the set "faceListA" and "faceListB"
  orderTets();
  computeFaceLists();
  // there will be p-1 (internal) interpolating points on the edge
  for (int i = 0; i < p-1; i++) {
    ma::Vector nodeXi;
    getBezierNodeXi(apf::Mesh::EDGE, p, i, nodeXi);
    double xi = (nodeXi[0] + 1.) / 2.;
    std::vector<ma::Vector> polygon = getPolygon(xA + (xB - xA) * xi);

    polygons.push_back(polygon);

    /* printf("Polygon is \n"); */
    /* for (size_t j = 0; j < polygon.size(); j++) { */
    /*   printf("%f, %f, %f\n", polygon[j].x(), */
			     /* polygon[j].y(), */
			     /* polygon[j].z()); */
    /* } */


    ma::Vector centerPoint = getCenterPoint(polygon);

    interpolatingPoints.push_back(centerPoint);
  }

  // set the nodal values to be the interpolating points for now
  apf::FieldShape* fs = mesh->getShape();
  int ne = fs->countNodesOn(apf::Mesh::EDGE);
  for (int i = 0; i < ne; i++) {
    mesh->setPoint(edge, i, interpolatingPoints[i]);
  }


}

void CrvEdgeReshapeHO::computeControlPoints()
{
  apf::FieldShape* fs = mesh->getShape();
  int n = fs->getEntityShape(apf::Mesh::EDGE)->countNodes();
  int ne = fs->countNodesOn(apf::Mesh::EDGE);
  PCU_ALWAYS_ASSERT(int(interpolatingPoints.size()) == ne);


  // bezier coeffs
  apf::NewArray<double> c;
  getBezierTransformationCoefficients(p, apf::Mesh::EDGE, c);

  convertInterpolationPoints(mesh, edge, n, ne, c);

  /* // construct "nodes" */
  /* apf::NewArray<ma::Vector> nodes(n); */
  /* nodes[0] = ma::getPosition(mesh, vertA); */
  /* for (int i = 1; i < n-1; i++) { */
  /*   nodes[i] = interpolatingPoints[i-1]; */
  /* } */
  /* nodes[n-1] = ma::getPosition(mesh, vertB); */


  /* apf::NewArray<ma::Vector> newNodes(ne); */

  /* convertInterpolationPoints(n, ne, nodes, c, newNodes); */

  /* for (int i = 0; i < ne; i++) { */
  /*   controlPoints.push_back(newNodes[i]); */
  /* } */

  /* PCU_ALWAYS_ASSERT(int(controlPoints.size()) == ne); */

  /* for (int i = 0; i < ne; i++) { */
  /*   mesh->setPoint(edge, i, controlPoints[i]); */
  /* } */

}


// crvEdgeReshapeHOS member functions //
void CrvEdgeReshapeHOS::getStrands()
{
  std::vector<ma::Entity*> fs;
  apf::Up up;
  mesh->getUp(edge, up);
  for (int i = 0; i < up.n; i++) {
    int type = mesh->getModelType(mesh->toModel(up.e[i]));
    if (type == 2) fs.push_back(up.e[i]);
  }
  PCU_ALWAYS_ASSERT(int(fs.size()) == 2);

  // edges of fs[0] corresponds to upper strand
  // edges of fs[1] corresponds to lower strand

  ma::Entity* edges[2][2] = {{0, 0}, {0, 0}};

  int dirs[2][2];

  // compute the upper ones
  for (int i = 0; i < 2; i++) {
    ma::Entity* face = fs[i];
    ma::Entity* es[3];
    mesh->getDownward(face, 1, es);
    for (int k = 0; k < 3; k++) {
      if (es[k] == edge) continue;
      if (getSharedVert(mesh, es[k], edge) == vertA) {
	edges[i][0] = es[k];
	ma::Entity* vs[2];
	mesh->getDownward(es[k], 0, vs);
	if (vs[0] == vertA)
	  dirs[i][0] = 1;
	else
	  dirs[i][0] = -1;
      }
      if (getSharedVert(mesh, es[k], edge) == vertB) {
	edges[i][1] = es[k];
	ma::Entity* vs[2];
	mesh->getDownward(es[k], 0, vs);
	if (vs[0] == vertB)
	  dirs[i][1] = -1;
	else
	  dirs[i][1] = 1;
      }
    }
    PCU_ALWAYS_ASSERT(edges[i][0]);
    PCU_ALWAYS_ASSERT(edges[i][1]);
  }

  // array of edge nodes including verts
  std::vector<ma::Vector> eNodes;
  eNodes.clear();
  eNodes.push_back(ma::Vector(-1.0, 0.0, 0.0));
  for (int i = 0; i < p-1; i++) {
    ma::Vector xi;
    getBezierNodeXi(apf::Mesh::EDGE, p, i, xi);
    eNodes.push_back(xi);
  }
  eNodes.push_back(ma::Vector(1.0, 0.0, 0.0));

  for (int i = 0; i < 2; i++) { // strands 0 and 1
    apf::MeshElement* m0 = apf::createMeshElement(mesh, edges[i][0]);
    apf::MeshElement* m1 = apf::createMeshElement(mesh, edges[i][1]);
    // add the points from the edge 0 in the strand
    for (size_t j = 0; j < eNodes.size(); j++) {
      ma::Vector xi = eNodes[j] * dirs[i][0];
      ma::Vector x;
      apf::mapLocalToGlobal(m0, xi, x);
      strands[i].push_back(x);
    }
    // add the points from the edge 1 in the strand
    for (size_t j = 1; j < eNodes.size(); j++) {
      ma::Vector xi = eNodes[j] * dirs[i][1];
      ma::Vector x;
      apf::mapLocalToGlobal(m1, xi, x);
      strands[i].push_back(x);
    }
    apf::destroyMeshElement(m0);
    apf::destroyMeshElement(m1);
  }
}

std::vector<ma::Vector> CrvEdgeReshapeHOS::getEndPoints(ma::Vector x)
{
  std::vector<ma::Vector> endPoints;
  for (int i = 0; i < 2; i++) {
    endPoints.push_back(getCrossSection(strands[i], x));
  }
  return endPoints;
}

void CrvEdgeReshapeHOS::computeInterpolatingPoints()
{
  // setup vertA and vertB
  ma::Entity* vs[2];
  mesh->getDownward(edge, 0, vs);
  vertA = vs[0];
  vertB = vs[1];

  ma::Vector xA = ma::getPosition(mesh, vertA);
  ma::Vector xB = ma::getPosition(mesh, vertB);

  // do things that are used by other member functions
  // like setting up the strands (upper and lower)
  getStrands();

  // there will be p-1 (internal) interpolating points on the edge
  for (int i = 0; i < p-1; i++) {
    ma::Vector nodeXi;
    getBezierNodeXi(apf::Mesh::EDGE, p, i, nodeXi);
    double xi = (nodeXi[0] + 1.) / 2.;
    std::vector<ma::Vector> endPoints = getEndPoints(xA + (xB - xA) * xi);

    /* printf("Polygon is \n"); */
    /* for (size_t j = 0; j < polygon.size(); j++) { */
    /*   printf("%f, %f, %f\n", polygon[j].x(), */
			     /* polygon[j].y(), */
			     /* polygon[j].z()); */
    /* } */


    ma::Vector centerPoint = (endPoints[0] + endPoints[1]) / 2.;

    // make sure that interpolating points are on the model boundary
    //

    interpolatingPoints.push_back(centerPoint);
  }

  // set the nodal values to be the interpolating points for now
  apf::FieldShape* fs = mesh->getShape();
  int ne = fs->countNodesOn(apf::Mesh::EDGE);
  for (int i = 0; i < ne; i++) {
    mesh->setPoint(edge, i, interpolatingPoints[i]);
  }

}

void CrvEdgeReshapeHOS::computeControlPoints()
{
  apf::FieldShape* fs = mesh->getShape();
  int n = fs->getEntityShape(apf::Mesh::EDGE)->countNodes();
  int ne = fs->countNodesOn(apf::Mesh::EDGE);
  PCU_ALWAYS_ASSERT(int(interpolatingPoints.size()) == ne);


  // bezier coeffs
  apf::NewArray<double> c;
  getBezierTransformationCoefficients(p, apf::Mesh::EDGE, c);


  convertInterpolationPoints(mesh, edge, n, ne, c);

  if (mesh->getModelTag(mesh->toModel(edge)) == 42) {
    apf::Up up;
    mesh->getUp(edge, up);
    for (int i = 0; i < up.n; i++) {
      if (mesh->getModelType(mesh->toModel(edge)) == 3) continue;
      apf::NewArray<apf::Vector3> l;
      apf::Element* elem =
	  apf::createElement(mesh->getCoordinateField(),up.e[i]);
      apf::getVectorNodes(elem,l);
      std::cout << "Control Points Are:" << std::endl;

      for (int j = 0; j < (p+1)*(p+2)/2; j++) {
        printf("%f, %f, %f\n", l[j].x(), l[j].y(), l[j].z());
      }

      apf::destroyElement(elem);
    }
  }

  // blend up
  /// upward adjacent faces first
  if (fs->hasNodesIn(2)) {
    std::cout << std::endl;
    std::cout << "HERE" << std::endl;
    std::cout << std::endl;
    n = fs->getEntityShape(apf::Mesh::TRIANGLE)->countNodes();
    ne = fs->countNodesOn(apf::Mesh::TRIANGLE);
    c.deallocate();
    getInternalBezierTransformationCoefficients(mesh,p,1,
	apf::Mesh::TRIANGLE,c);

/*     printf("cs are\n"); */
/*     for (size_t j = 0; j < 300; j++) { */
/*       std::cout << c[j] << " "; */
/*     } */
/*     std::cout << std::endl; */


    apf::Up up;
    mesh->getUp(edge, up);
    for (int i = 0; i < up.n; i++) {
      printf("blending %d upward face\n", i);
      printf("before points are\n");
      for (int j = 0; j < ne; j++) {
	ma::Vector pt;
	mesh->getPoint(up.e[i], j, pt);
	std::cout << pt << std::endl;
      }
      convertInterpolationPoints(mesh, up.e[i], n-ne, ne, c);
      printf("after points are\n");
      for (int j = 0; j < ne; j++) {
	ma::Vector pt;
	mesh->getPoint(up.e[i], j, pt);
	std::cout << pt << std::endl;
      }
    }
  }
  /// upward adjacent tets next
  if (fs->hasNodesIn(3)) {
    std::cout << std::endl;
    std::cout << "HERE" << std::endl;
    std::cout << std::endl;
    n = fs->getEntityShape(apf::Mesh::TET)->countNodes();
    ne = fs->countNodesOn(apf::Mesh::TET);
    c.deallocate();
    getInternalBezierTransformationCoefficients(mesh,p,1,
	apf::Mesh::TET,c);

    apf::Adjacent adj;
    mesh->getAdjacent(edge, 3, adj);
    for (size_t i = 0; i < adj.getSize(); i++) {
      convertInterpolationPoints(mesh, adj[i], n-ne, ne, c);
    }
  }


  /* // construct "nodes" */
  /* apf::NewArray<ma::Vector> nodes(n); */
  /* nodes[0] = ma::getPosition(mesh, vertA); */
  /* for (int i = 1; i < n-1; i++) { */
  /*   nodes[i] = interpolatingPoints[i-1]; */
  /* } */
  /* nodes[n-1] = ma::getPosition(mesh, vertB); */


  /* apf::NewArray<ma::Vector> newNodes(ne); */

  /* convertInterpolationPoints(n, ne, nodes, c, newNodes); */

  /* for (int i = 0; i < ne; i++) { */
  /*   controlPoints.push_back(newNodes[i]); */
  /* } */

  /* PCU_ALWAYS_ASSERT(int(controlPoints.size()) == ne); */

  /* for (int i = 0; i < ne; i++) { */
  /*   mesh->setPoint(edge, i, controlPoints[i]); */
  /*   /1* mesh->setPoint(edge, i, nodes[i+1]); *1/ */
  /* } */

}


} // namespace crv
