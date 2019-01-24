/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef CRVEDGERESHAPEHO_H
#define CRVEDGERESHAPEHO_H

#include "crv.h"
#include <apf.h>
#include <pcu_util.h>
#include <vector>

/** \file crvEdgeReshapeHO.h
    \brief main class for edge reshaping for higher order meshes */

namespace crv {

typedef ma::Vector (*MapToParent)(double t, double s);

// HO stands for Higher Order
// this class is used for internal edges
class CrvEdgeReshapeHO
{
public:
  CrvEdgeReshapeHO(ma::Mesh* m, ma::Entity* e, int order) :
    mesh(m), edge(e), p(order)
  {
    vertA = 0;
    vertB = 0;
    tets.clear();
    faceListA.clear();
    faceListB.clear();
    dirListA.clear();
    dirListB.clear();
    interpolatingPoints.clear();
    controlPoints.clear();
  }
  ~CrvEdgeReshapeHO()
  {}
public:
  void orderTets();
  void computeFaceLists();
  std::vector<ma::Vector> getPolygon(ma::Vector x);
  void computeInterpolatingPoints();
  void computeControlPoints();
public:
  ma::Mesh* mesh;
  ma::Entity* edge;
  int p;
  ma::Entity* vertA;
  ma::Entity* vertB;
  std::vector<ma::Entity*> tets; // ordered list of upward adj tets to edge
  std::vector<ma::Entity*> faceListA;
  std::vector<ma::Entity*> faceListB;
  std::vector<int> dirListA;
  std::vector<int> dirListB;
  std::vector<std::vector<ma::Vector> > polygons;
  std::vector<ma::Vector> interpolatingPoints;
  std::vector<ma::Vector> controlPoints;
};

// HO stands for Higher Order and S stands for Surface
// this class is used for edges on the model surface
class CrvEdgeReshapeHOS
{
public:
  CrvEdgeReshapeHOS(ma::Mesh* m, ma::Entity* e, int order) :
    mesh(m), edge(e), p(order)
  {
    vertA = 0;
    vertB = 0;
    strands[0].clear();
    strands[1].clear();
    controlPoints.clear();
  }
  ~CrvEdgeReshapeHOS()
  {}
public:
  void getStrands();
  std::vector<ma::Vector> getEndPoints(ma::Vector x);
  void computeInterpolatingPoints();
  void computeControlPoints();
public:
  ma::Mesh* mesh;
  ma::Entity* edge;
  int p;
  ma::Entity* vertA;
  ma::Entity* vertB;
  std::vector<ma::Vector> strands[2];
  std::vector<ma::Vector> interpolatingPoints;
  std::vector<ma::Vector> controlPoints;
};

}

#endif
