/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#include <PCU.h>
#include "maSize.h"
#include "apfMatrix.h"
#include <apfShape.h>
#include <cstdlib>
#include <pcu_util.h>

namespace ma {

SizeField::~SizeField()
{
}

IdentitySizeField::IdentitySizeField(Mesh* m):
  mesh(m)
{
}

double IdentitySizeField::measure(Entity* e)
{
  apf::MeshElement* me = apf::createMeshElement(mesh,e);
  double x = apf::measure(me);
  apf::destroyMeshElement(me);
  return x;
}

bool IdentitySizeField::shouldSplit(Entity*)
{
  return false;
}

bool IdentitySizeField::shouldCollapse(Entity*)
{
  return false;
}

void IdentitySizeField::interpolate(
    apf::MeshElement*,
    Vector const&,
    Entity*, int)
{
}

void IdentitySizeField::getTransform(
        apf::MeshElement*,
        Vector const&,
        Matrix& t)
{
  t = Matrix(1,0,0,
             0,1,0,
             0,0,1);
}

double IdentitySizeField::getWeight(Entity*)
{
  return 1.0;
}

static void orthogonalizeR(Matrix& R)
{
  /* by the way, the principal direction vectors
     are in the columns, so lets do our cleanup
     work on the transpose */
  Matrix RT = transpose(R);
  /* we have to make sure this remains a rotation matrix.
     lets assume that component-wise interpolation gave
     a somewhat decent answer.
     (if not, the truly proper solution is to interpolate
      quaternion components, see also Lie groups)
     R[0] should be the interpolated main principal direction
     vector. lets first normalize it. */
  RT[0] = RT[0].normalize();
  /* now the next principal direction should roughly align
     with R[1], but first it must be orthogonal to R[0].
     so, remove its component along R[0] (remember R[0] is unit length)*/
  RT[1] = RT[1] - RT[0]*(RT[0]*RT[1]);
  /* and normalize it too */
  RT[1] = RT[1].normalize();
  /* finally, forget what was in R[2] and just make it the cross
     product of R[0] and R[1] to make a nice frame */
  RT[2] = apf::cross(RT[0],RT[1]);
  /* and back out of the transpose */
  /* if you're wondering, this is why R and h are stored separately:
     so we can do this interpolation dance on the rotation matrix */
  R = transpose(RT);
}

static void orthogonalEigenDecompForSymmetricMatrix(Matrix const& A, Vector& v, Matrix& R)
{
  /* here we assume A to be real symmetric 3x3 matrix, 
   * we should be able to get 3 orthogonal eigen vectors
   * we also normalize the eigen vectors */
  double eigenValues[3];
  Vector eigenVectors[3];
  
  apf::eigen(A, eigenVectors, eigenValues);
  
  Matrix RT(eigenVectors); // eigen vectors are stored in the rows of RT

  RT[0] = RT[0].normalize();
  RT[1] = RT[1] - RT[0]*(RT[0]*RT[1]);
  RT[1] = RT[1].normalize();
  RT[2] = apf::cross(RT[0],RT[1]);

  v = Vector(eigenValues);
  R = transpose(RT);
}

/* the length, area, or volume of
   the parent element for this
   entity type */
static double parentMeasure[apf::Mesh::TYPES] =
{0.0     //vert
,2.0     //edge
,1.0/2.0 //tri
,4.0     //quad
,1.0/6.0 //tet
,-42.0   //hex - definitely not sure
,-42.0   //prism - definitely not sure
,8.0/3.0 //pyramid
};

class SizeFieldIntegrator : public apf::Integrator
{
  public:
    SizeFieldIntegrator(SizeField* sF):
      Integrator(2),
      measurement(0),
      sizeField(sF),
      meshElement(0),
      dimension(0)
    {
    }
    void process(apf::MeshElement* me)
    {
      this->inElement(me);
      int np = countIntPoints(me,this->order);
      for (int p=0; p < np; ++p)
      {
        ipnode = p;
        Vector point;
        getIntPoint(me,this->order,p,point);
        double w = getIntWeight(me,this->order,p);
        this->atPoint(point,w,0);
      }
      this->outElement();
    }
    void inElement(apf::MeshElement* me)
    {
      meshElement = me;
      dimension = apf::getDimension(me);
    }
    void atPoint(Vector const& p , double w, double )
    {
      Matrix Q;
      sizeField->getTransform(meshElement,p,Q);
      Matrix J;
      apf::getJacobian(meshElement,p,J);
/* transforms the rows of J, the differential tangent vectors,
   into the metric space, then uses the generalized determinant */
      double dV2 = apf::getJacobianDeterminant(J*Q,dimension);
      measurement += w*dV2;
    }
    double measurement;
    SizeField* sizeField;
  private:
    apf::MeshElement* meshElement;
    int dimension;
};

struct MetricSizeField : public SizeField
{
  double measure(Entity* e)
  {
    SizeFieldIntegrator sFI(this); 
    apf::MeshElement* me = apf::createMeshElement(mesh, e);
    sFI.process(me);
    apf::destroyMeshElement(me);
    return sFI.measurement;
  }
  bool shouldSplit(Entity* edge)
  {
    return this->measure(edge) > 1.5;
  }
  bool shouldCollapse(Entity* edge)
  {
    return this->measure(edge) < 0.5;
  }
  double getWeight(Entity* e)
  {
    /* parentMeasure is used to normalize */
    return measure(e) / parentMeasure[mesh->getType(e)];
  }
  Mesh* mesh;
};

AnisotropicFunction::~AnisotropicFunction()
{
}

IsotropicFunction::~IsotropicFunction()
{
}

/* struct IsoWrapper : public AnisotropicFunction */
/* { */
/*   IsoWrapper() */
/*   { */
/*   } */
/*   IsoWrapper(IsotropicFunction* f) */
/*   { */
/*     function = f; */
/*   } */
/*   void getValue(const Vector& x, Matrix& r, Vector& h) */
/*   { */
/*     r = Matrix(1,0,0, */
/*                0,1,0, */
/*                0,0,1); */
/*     double s = function->getValue(x); */
/*     h = Vector(s,s,s); */
/*   } */
/*   IsotropicFunction* function; */
/* }; */

struct BothEval
{
  BothEval()
  {
  }
  BothEval(AnisotropicFunction* f, Mesh* m) :
    function(f), mesh(m)
  {
    // if shape of mesh is interpolating (i.e., Lagrange)
    // and order is <= 3 use coordinateField
    if (mesh->getShape()->isInterpolating() &&
	mesh->getShape()->getOrder() <= 3) {
      useCoordinates = true;
      shape = mesh->getShape();
    }
    else { // else use a Lagrange shape of order of the mesh up to order 3
      useCoordinates = false;
      int order = mesh->getShape()->getOrder();
      order = order <= 3 ? order : 3;
      shape = apf::getLagrange(order);
    }
    cachedEnt = 0;
  }
  void updateCache(Entity* e)
  {
    if (e == cachedEnt)
      return;
    int non = shape->countNodesOn(mesh->getType(e));
    cachedSizeVecs.setSize(non);
    cachedFrameMats.setSize(non);
    Vector coords;
    if (useCoordinates)
      for (int i = 0; i < non; i++) {
	mesh->getPoint(e, i, coords);
	function->getValue(coords, cachedFrameMats[i], cachedSizeVecs[i]);
      }
    else {
      Vector xi;
      apf::MeshElement* me = apf::createMeshElement(mesh, e);
      for (int i = 0; i < non; i++) {
      	shape->getNodeXi(mesh->getType(e), i, xi);
	apf::mapLocalToGlobal(me, xi, coords);
	function->getValue(coords, cachedFrameMats[i], cachedSizeVecs[i]);
      }
      apf::destroyMeshElement(me);
    }
    cachedEnt = e;
  }
  void getSizeAtNode(Entity* e, int n, Vector& s)
  {
    updateCache(e);
    s = cachedSizeVecs[n];
  }
  void getFrameAtNode(Entity* e, int n, Matrix& f)
  {
    updateCache(e);
    f = cachedFrameMats[n];
  }
  void getAllSizes(Entity* e, apf::DynamicArray<Vector>& s)
  {
     updateCache(e);
     s = cachedSizeVecs;
  }
  void getAllFrames(Entity* e, apf::DynamicArray<Matrix>& s)
  {
     updateCache(e);
     s = cachedFrameMats;
  }
  AnisotropicFunction* function;
  Mesh* mesh;
  bool useCoordinates;
  apf::FieldShape* shape;
  Entity* cachedEnt;
  apf::DynamicArray<Vector> cachedSizeVecs;
  apf::DynamicArray<Matrix> cachedFrameMats;
};

struct SizesEval : public apf::Function
{
  SizesEval()
  {
  }
  SizesEval(BothEval* b)
  {
    both = b;
  }
  void eval(Entity* e, double* result)
  {
    apf::DynamicArray<Vector> s;
    both->getAllSizes(e, s);
    for (std::size_t i = 0; i < s.getSize(); i++) {
      s[i].toArray(result+3*i);
    }
  }
  BothEval* both;
};

struct FrameEval : public apf::Function
{
  FrameEval()
  {
  }
  FrameEval(BothEval* b)
  {
    both = b;
  }
  void eval(Entity* e, double* result)
  {
    apf::DynamicArray<Matrix> f;
    both->getAllFrames(e, f);
    for (std::size_t i = 0; i < f.getSize(); i++) { // ith matrix in the list
      for (int j = 0; j < 3; j++) { // jth row of the matrix
	((Vector)f[i][j]). //cast jth row of the matrix to Vector so toArray can be used
	  toArray(result+9*i+3*j);
      }
    }
  }
  BothEval* both;
};

struct LogMEval : public apf::Function
{
  LogMEval()
  {
  }
  LogMEval(AnisotropicFunction* f, Mesh* m) :
    function(f), mesh(m)
  {
    // if shape of mesh is Lagrange and order is <= 3 use coordinateField
    if (mesh->getShape()->isInterpolating() &&
	mesh->getShape()->getOrder() <= 3) {
      useCoordinates = true;
      shape = mesh->getShape();
    }
    else {
      useCoordinates = false;
      int order = mesh->getShape()->getOrder();
      order = order <= 3 ? order : 3;
      shape = apf::getLagrange(order);
    }
    cachedEnt = 0;
  }
  void updateCache(Entity* e)
  {
    if (e == cachedEnt)
      return;

    int non = shape->countNodesOn(mesh->getType(e));
    cachedLogMs.setSize(non);
    Vector coords; // physical coordinate of the node
    Vector h;      // vector of sizes
    Matrix R;      // principal directions
    if (useCoordinates)
      for (int i = 0; i < non; i++) {
	mesh->getPoint(e, i, coords);
	function->getValue(coords, R, h);
        Matrix S( -2*log(h[0]), 0          , 0,
	           0          ,-2*log(h[1]), 0,
	           0          , 0          ,-2*log(h[2]));
	cachedLogMs[i] = R*S*transpose(R);
      }
    else {
      Vector xi;     // parent coordinate of the node
      apf::MeshElement* me = apf::createMeshElement(mesh, e);
      for (int i = 0; i < non; i++) {
      	shape->getNodeXi(mesh->getType(e), i, xi);
	apf::mapLocalToGlobal(me, xi, coords);
	function->getValue(coords, R, h);
        Matrix S( -2*log(h[0]), 0          , 0,
	           0          ,-2*log(h[1]), 0,
	           0          , 0          ,-2*log(h[2]));
	cachedLogMs[i] = R*S*transpose(R);
      }
      apf::destroyMeshElement(me);
    }
    cachedEnt = e;
  }
  void getLogMAtNode(Entity* e, int n, Matrix& f)
  {
    updateCache(e);
    f = cachedLogMs[n];
  }
  void getAllLogMs(Entity* e, apf::DynamicArray<Matrix>& f)
  {
    updateCache(e);
    f = cachedLogMs;
  }
  void eval(Entity* e, double* result)
  {
    apf::DynamicArray<Matrix> f;
    getAllLogMs(e, f);
    for (std::size_t i = 0; i < f.getSize(); i++) { // ith matrix in the list
      for (int j = 0; j < 3; j++) { // jth row of the matrix
	((Vector)f[i][j]). //cast jth row of the matrix to Vector so toArray can be used
	  toArray(result+9*i+3*j);
      }
    }
  }
  AnisotropicFunction* function;
  Mesh* mesh;
  bool useCoordinates;
  Entity* cachedEnt;
  apf::FieldShape* shape;
  apf::DynamicArray<Matrix> cachedLogMs;
};

struct IsoEval : public apf::Function
{
  IsoEval()
  {
  }
  IsoEval(IsotropicFunction* f, Mesh* m) :
    function(f), mesh(m)
  {
    // if shape of mesh is Lagrange and order is <= 3 use coordinateField
    if (mesh->getShape()->isInterpolating() &&
	mesh->getShape()->getOrder() <= 3) {
      useCoordinates = true;
      shape = mesh->getShape();
    }
    else {
      useCoordinates = false;
      int order = mesh->getShape()->getOrder();
      order = order <= 3 ? order : 3;
      shape = apf::getLagrange(order);
    }
    cachedEnt = 0;
  }
  void updateCache(Entity* e)
  {
    if (e == cachedEnt)
      return;

    int non = shape->countNodesOn(mesh->getType(e));
    cachedSizes.setSize(non);
    Vector coords; // physical coordinate of the node
    if (useCoordinates)
      for (int i = 0; i < non; i++) {
	mesh->getPoint(e, i, coords);
	cachedSizes[i] = function->getValue(coords);
      }
    else {
      Vector xi;     // parent coordinate of the node
      apf::MeshElement* me = apf::createMeshElement(mesh, e);
      for (int i = 0; i < non; i++) {
      	shape->getNodeXi(mesh->getType(e), i, xi);
	apf::mapLocalToGlobal(me, xi, coords);
	cachedSizes[i] = function->getValue(coords);
      }
      apf::destroyMeshElement(me);
    }
    cachedEnt = e;
  }
  void getSizeAtNode(Entity* e, int n, double& f)
  {
    updateCache(e);
    f = cachedSizes[n];
  }
  void getAllSizes(Entity* e, apf::DynamicArray<double>& f)
  {
    updateCache(e);
    f = cachedSizes;
  }
  void eval(Entity* e, double* result)
  {
    apf::DynamicArray<double> f;
    getAllSizes(e, f);
    for (std::size_t i = 0; i < f.getSize(); i++) {
      result[i] = f[i];
    }
  }
  IsotropicFunction* function;
  Mesh* mesh;
  bool useCoordinates;
  Entity* cachedEnt;
  apf::FieldShape* shape;
  apf::DynamicArray<double> cachedSizes;
};


struct AnisoSizeField : public MetricSizeField
{
  AnisoSizeField()
  {
  }
  AnisoSizeField(Mesh* m, AnisotropicFunction* f):
    bothEval(f, m),
    sizesEval(&bothEval),
    frameEval(&bothEval)
  {
    mesh = m;
    hField = apf::createUserField(m, "ma_sizes", apf::VECTOR,
        apf::getLagrange(bothEval.shape->getOrder()), &sizesEval);
    rField = apf::createUserField(m, "ma_frame", apf::MATRIX,
        apf::getLagrange(bothEval.shape->getOrder()), &frameEval);
  }
  ~AnisoSizeField()
  {
    apf::destroyField(hField);
    apf::destroyField(rField);
  }
  void init(Mesh* m, apf::Field* sizes, apf::Field* frames)
  {
    mesh = m;
    hField = sizes;
    rField = frames;
  }
  void getTransform(
      apf::MeshElement* me,
      Vector const& xi,
      Matrix& Q)
  {
    apf::Element* hElement = apf::createElement(hField,me);
    apf::Element* rElement = apf::createElement(rField,me);
    Vector h;
    Matrix R;
    apf::getVector(hElement,xi,h);
    apf::getMatrix(rElement,xi,R);
    apf::destroyElement(hElement);
    apf::destroyElement(rElement);
    orthogonalizeR(R);
    Matrix S(1/h[0],0,0,
             0,1/h[1],0,
             0,0,1/h[2]);
    Q = R*S;
  }
  void interpolate(
      apf::MeshElement* parent,
      Vector const& xi,
      Entity* ent, int node)
  {
    PCU_ALWAYS_ASSERT_VERBOSE(
    	node < apf::getShape(rField)->countNodesOn(mesh->getType(ent)),
    	"node out of range for the Fields!");
    apf::Element* rElement = apf::createElement(rField,parent);
    apf::Element* hElement = apf::createElement(hField,parent);
    Vector h;
    apf::getVector(hElement,xi,h);
    Matrix R;
    apf::getMatrix(rElement,xi,R);
    orthogonalizeR(R);
    this->setValue(ent,node,R,h);
    apf::destroyElement(hElement);
    apf::destroyElement(rElement);
  }
  void setValue(
      Entity* e,
      int node,
      Matrix const& r,
      Vector const& h)
  {
    apf::setMatrix(rField,e,node,r);
    apf::setVector(hField,e,node,h);
  }
  void setIsotropicValue(
      Entity* e,
      int node,
      double value)
  {
    this->setValue(e,node,
                   Matrix(1,0,0,
                          0,1,0,
                          0,0,1),
                   Vector(value,value,value));
  }
  apf::Field* hField;
  apf::Field* rField;
  BothEval bothEval;
  SizesEval sizesEval;
  FrameEval frameEval;
};

struct LogAnisoSizeField : public MetricSizeField
{
  LogAnisoSizeField()
  {
  }
  LogAnisoSizeField(Mesh* m, AnisotropicFunction* f):
    logMEval(f, m)
  {
    mesh = m;
    logMField = apf::createUserField(m, "ma_logM", apf::MATRIX,
        apf::getLagrange(logMEval.shape->getOrder()), &logMEval);
  }
  ~LogAnisoSizeField()
  {
    apf::destroyField(logMField);
  }
  void init(Mesh* m, apf::Field* sizes, apf::Field* frames)
  {
    mesh = m;
    logMField = apf::createField(m, "ma_logM", apf::MATRIX, apf::getLagrange(1));
    Entity* v;
    Iterator* it = m->begin(0);
    while ( (v = m->iterate(it)) ) {
      Vector h;
      Matrix f;
      apf::getVector(sizes, v, 0, h);
      apf::getMatrix(frames, v, 0, f);
      Vector s(log(1/h[0]/h[0]), log(1/h[1]/h[1]), log(1/h[2]/h[2]));
      Matrix S(s[0], 0   , 0,
              0    , s[1], 0,
              0    , 0   , s[2]);
      apf::setMatrix(logMField, v, 0, f * S * transpose(f));
    }
  }
  void getTransform(
      apf::MeshElement* me,
      Vector const& xi,
      Matrix& Q)
  {
    apf::Element* logMElement = apf::createElement(logMField,me);
    Matrix logM;
    apf::getMatrix(logMElement,xi,logM);
    apf::destroyElement(logMElement);
    Vector v;
    Matrix R;
    orthogonalEigenDecompForSymmetricMatrix(logM, v, R);
    Matrix S( sqrt(exp(v[0])), 0, 0,
              0, sqrt(exp(v[1])), 0,
              0, 0, sqrt(exp(v[2])));
    Q = R*S;
  }
  void interpolate(
      apf::MeshElement* parent,
      Vector const& xi,
      Entity* ent, int node)
  {
    PCU_ALWAYS_ASSERT_VERBOSE(
    	node < apf::getShape(logMField)->countNodesOn(mesh->getType(ent)),
    	"node out of range for the Fields!");
    apf::Element* logMElement = apf::createElement(logMField,parent);
    Matrix logM;
    apf::getMatrix(logMElement,xi,logM);
    this->setValue(ent,node,logM);
    apf::destroyElement(logMElement);
  }
  void setValue(
      Entity* e,
      int node,
      Matrix const& logM)
  {
    apf::setMatrix(logMField,e,node,logM);
  }
  void setIsotropicValue(
      Entity* e,
      int node,
      double value)
  {
    this->setValue(e, node,
                   Matrix(value,0,0,
                          0,value,0,
                          0,0,value));
  }
  apf::Field* logMField;
  LogMEval logMEval;
};

struct IsoSizeField : public MetricSizeField
{
  IsoSizeField()
  {
  }
  IsoSizeField(Mesh* m, IsotropicFunction* f) :
    isoEval(f, m)
  {
    mesh = m;
    sField = apf::createUserField(m, "ma_sizes", apf::SCALAR,
    	apf::getLagrange(isoEval.shape->getOrder()), &isoEval);
  }
  ~IsoSizeField()
  {
    apf::destroyField(sField);
  }
  void init(Mesh* m, apf::Field* sizes)
  {
    mesh = m;
    sField = sizes;
  }
  void getTransform(
      apf::MeshElement* me,
      Vector const& xi,
      Matrix& Q)
  {
    apf::Element* el = apf::createElement(sField,me);
    double h = apf::getScalar(el,xi);
    Q = Matrix(1./h, 0.  , 0.,
               0.,   1./h, 0.,
               0.,   0.  , 1./h);
    apf::destroyElement(el);
  }
  void interpolate(
      apf::MeshElement* parent,
      Vector const& xi,
      Entity* ent, int node)
  {
    PCU_ALWAYS_ASSERT_VERBOSE(
    	node < apf::getShape(sField)->countNodesOn(mesh->getType(ent)),
    	"node out of range for the Fields!");
    apf::Element* el = apf::createElement(sField,parent);
    double h = apf::getScalar(el,xi);
    this->setValue(ent,node,h);
    apf::destroyElement(el);
  }
  void setValue(
      Entity* e,
      int node,
      double h)
  {
    apf::setScalar(sField,e,node,h);
  }
  IsoEval isoEval;
  apf::Field* sField;
};

SizeField* makeSizeField(Mesh* m, apf::Field* sizes, apf::Field* frames,
    bool logInterpolation)
{
  // logInterpolation is "false" by default
  if (! logInterpolation) {
    AnisoSizeField* anisoF = new AnisoSizeField();
    anisoF->init(m, sizes, frames);
    return anisoF;
  }
  else {
    LogAnisoSizeField* logAnisoF = new LogAnisoSizeField();
    logAnisoF->init(m, sizes, frames);
    return logAnisoF;
  }
}

SizeField* makeSizeField(Mesh* m, AnisotropicFunction* f,
    bool logInterpolation)
{
  // logInterpolation is "false" by default
  if(! logInterpolation)
    return new AnisoSizeField(m, f);
  else
    return new LogAnisoSizeField(m,f);
}

SizeField* makeSizeField(Mesh* m, IsotropicFunction* f)
{
  return new IsoSizeField(m, f);
}

SizeField* makeSizeField(Mesh* m, apf::Field* size)
{
  IsoSizeField* isoSize = new IsoSizeField();
  isoSize->init(m, size);
  return isoSize;
}

double getAverageEdgeLength(Mesh* m)
{
  IdentitySizeField sizeField(m);
  double sums[2];
  double& length_sum = sums[0];
  double& edge_count = sums[1];
  length_sum = edge_count = 0;
  apf::MeshIterator* it = m->begin(1);
  Entity* e;
  while ((e = m->iterate(it)))
  {
    length_sum += sizeField.measure(e);
    edge_count += 1.0;
  }
  m->end(it);
  PCU_Add_Doubles(sums,2);
  return length_sum / edge_count;
}

double getMaximumEdgeLength(Mesh* m, SizeField* sf)
{
  if (!sf)
    sf = new IdentitySizeField(m);
  double maxLength = 0.0;
  apf::MeshIterator* it = m->begin(1);
  Entity* e;
  while ((e = m->iterate(it)))
  {
    if (!m->isOwned(e))
      continue;
    double length = sf->measure(e);
    if (length > maxLength)
      maxLength = length;
  }
  m->end(it);
  PCU_Max_Doubles(&maxLength,1);
  return maxLength;
}


}
