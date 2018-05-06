/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef CRVCURVEENTITIES_H
#define CRVCURVEENTITIES_H

/** \file crvCurveEntities.h
  * \brief main file for curving edges and faces */

#include<maOperator.h>
#include"crvSnap.h"
#include"crvShape.h"
#include"crvBezier.h"
#include"crvEdgeReshape.h"

namespace crv {

static void convertInterpolationPoints2(int n, int ne,
    apf::NewArray<apf::Vector3>& nodes,
    apf::NewArray<double>& c,
    apf::NewArray<apf::Vector3>& newNodes){

  for(int i = 0; i < ne; ++i)
    newNodes[i].zero();

  for( int i = 0; i < ne; ++i)
    for( int j = 0; j < n; ++j)
      newNodes[i] += nodes[j]*c[i*n+j];

}

static void convertInterpolationPoints2(apf::Mesh2* m, apf::MeshEntity* e,
    int n, int ne, apf::NewArray<double>& c){

  apf::NewArray<apf::Vector3> l, b(ne);
  apf::Element* elem =
      apf::createElement(m->getCoordinateField(),e);
  apf::getVectorNodes(elem,l);

  crv::convertInterpolationPoints2(n,ne,l,c,b);

  for(int i = 0; i < ne; ++i)
    m->setPoint(e,i,b[i]);

  apf::destroyElement(elem);
}


static void snapToInterpolate2(apf::Mesh2* m, apf::MeshEntity* e)
{
  PCU_ALWAYS_ASSERT(m->canSnap());
  int type = m->getType(e);
  if(type == apf::Mesh::VERTEX){
    apf::Vector3 p, pt(0,0,0);
    apf::ModelEntity* g = m->toModel(e);
    m->getParamOn(g,e,p);
    m->snapToModel(g,p,pt);
    m->setPoint(e,0,pt);
    return;
  }
  apf::FieldShape * fs = m->getShape();
  int non = fs->countNodesOn(type);
  apf::Vector3 p, xi, pt(0,0,0);
  for(int i = 0; i < non; ++i){
    apf::ModelEntity* g = m->toModel(e);
    fs->getNodeXi(type,i,xi);
    if(type == apf::Mesh::EDGE)
      transferParametricOnEdgeSplit(m,e,0.5*(xi[0]+1.),p);
    else
      transferParametricOnTriSplit(m,e,xi,p);
    m->snapToModel(g,p,pt);
    m->setPoint(e,i,pt);
  }
}

/* static int markEdges(ma::Mesh* m, ma::Entity* e, int tag, */
/*     ma::Entity* edges[6]) */
/* { */
/*   if ( tag <= 1 ) // if its valid, or not checked, don't worry about it */
/*     return 0; */
/*   int dim = (tag-2)/6; */
/*   int index = (tag-2) % 6; */
/*   int n = 0; */
/*   int md = m->getDimension(); */

/*   switch (dim) { */
/*     case 0: */
/*     { */
/*       // if we have an invalid vertex, operate on its edges */
/*       ma::Downward ed; */
/*       m->getDownward(e,1,ed); */
/*       n = md; */
/*       if(md == 2){ */
/*         edges[0] = ed[index]; */
/*         edges[1] = ed[(index+2) % 3]; */
/*       } else { */
/*         PCU_ALWAYS_ASSERT(index < 4); */
/*         edges[0] = ed[vertEdges[index][0]]; */
/*         edges[1] = ed[vertEdges[index][1]]; */
/*         edges[2] = ed[vertEdges[index][2]]; */
/*       } */
/*     } */
/*     break; */
/*     case 1: */
/*     { */
/*       // if we have a single invalid edge, operate on it */
/*       ma::Downward ed; */
/*       m->getDownward(e,1,ed); */
/*       edges[0] = ed[index]; */
/*       n = 1; */
/*     } */
/*     break; */
/*     case 2: */
/*     { */
/*       // if we have an invalid face, operate on its edges */
/*       ma::Downward ed, faces; */
/*       m->getDownward(e,2,faces); */
/*       m->getDownward(faces[index],1,ed); */
/*       n = 3; */
/*       edges[0] = ed[0]; */
/*       edges[1] = ed[1]; */
/*       edges[2] = ed[2]; */
/*     } */
/*     break; */
/*     case 3: */
/*       m->getDownward(e,1,edges); */
/*       n = 6; */
/*       break; */
/*     default: */
/*       fail("invalid quality tag in markEdges\n"); */
/*       break; */
/*   } */

/*   return n; */
/* } */

class CurveEdge : public ma::Operator
{
public:
  CurveEdge(Adapt* a)
  {
    adapter = a;
    mesh = a->mesh;
    edge = 0;
    fs = mesh->getShape();
    qual = makeQuality(mesh, 2);
    ns = nf = 0;
  }
  virtual ~CurveEdge()
  {
    delete qual;
  }
  virtual int getTargetDimension() {return 1;}
  virtual bool shouldApply(ma::Entity* e)
  {
    int md = mesh->getModelType(mesh->toModel(e));
    if (md == 0 || md == 3) return false;
    if (!fs->hasNodesIn(1)) return false;
    edge = e;
    return true;
  }
  virtual bool requestLocality(apf::CavityOp* o)
  {
    return o->requestLocality(&edge,1);
  }
  virtual void apply()
  {
    // get the old node locations
    int nn = fs->countNodesOn(apf::Mesh::simplexTypes[1]);
    std::vector<ma::Vector> oldPositions;
    for (int i = 0; i < nn; i++) {
      ma::Vector pos;
      mesh->getPoint(edge, i, pos);
      oldPositions.push_back(pos);
    }

    snapToInterpolate2(mesh, edge);

    int order = fs->getOrder();
    int n = fs->getEntityShape(apf::Mesh::simplexTypes[1])->countNodes();
    int ne = fs->countNodesOn(apf::Mesh::simplexTypes[1]);
    apf::NewArray<double> c;
    getBezierTransformationCoefficients(order,
    	apf::Mesh::simplexTypes[1], c);
    convertInterpolationPoints2(mesh,edge,n,ne,c);

    // check the validity here
    apf::Adjacent adj;
    mesh->getAdjacent(edge, 3, adj);

    bool isCavityValid = true;
    int invalidCount = 0;
    for (std::size_t i = 0; i < adj.getSize(); ++i){
      if (qual->getQuality(adj[i]) < 0){
      	isCavityValid = false;
      	invalidCount++;
      	/* break; */
      }
    }


    // if cavity is ok then return
    if (isCavityValid) {
      ns++;
      return;
    }

    // at this point there is an invalid element in the cavity
    // find it and operate on it if possible
    for (std::size_t i = 0; i < adj.getSize(); ++i){
      int qualityTag = qual->checkValidity(adj[i]);
      /* int nc = markEdges(mesh, adj[i], qualityTag, edges); */
      if (qualityTag < 2) continue;
      EdgeReshape er;
      er.Init(adapter, qualityTag);
      er.setSimplex(adj[i]);
      if (er.reshape()) {
      	// check the outer cavity
      	bool isCavityValid = true;
      	int invalidCount = 0;
      	for (std::size_t j = 0; j < adj.getSize(); ++j){
      	  /* int v = qual->checkValidity(adj[j]); */
      	  /* if (v >= 2) { */
      	  double v = qual->getQuality(adj[j]);
      	  if (v < 0) {
      	    isCavityValid = false;
      	    invalidCount++;
      	    /* break; */
	  }
	}
	if (isCavityValid) {
	  ns++;
	  return;
	}
	else {
	  er.cancel();
	}
      }
      else {
      	er.cancel();
      }
    }


    nf++;
    for (int i = 0; i < nn; i++) {
      ma::Vector pos = oldPositions[i];
      mesh->setPoint(edge, i, pos);
    }
  }
private:
  Adapt* adapter;
  ma::Mesh* mesh;
  ma::Entity* edge;
  apf::FieldShape* fs;
  Quality* qual;
public:
  int ns;
  int nf;
};


}

#endif
