#include <ma.h>
#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <PCU.h>
#include <lionPrint.h>
#include <apfNumbering.h>
#include <apfShape.h>
#ifdef HAVE_SIMMETRIX
#include <gmi_sim.h>
#include <SimUtil.h>
#include <MeshSim.h>
#include <SimModel.h>
#endif
#include <pcu_util.h>
#include <cstdlib>
#include <iostream>

static void getEdgeLenAndCnt(
    ma::Mesh* m,
    ma::Entity* v,
    ma::SizeField* isf,
    double& len, int& cnt){
  len = 0;
  cnt = 0;
  apf::Up edges;
  m->getUp(v, edges);
  for(int eIdx=0; eIdx < edges.n; eIdx++) {
    if( m->isOwned(edges.e[eIdx])) {
      cnt++;
      len += isf->measure(edges.e[eIdx]);
    }
  }
}

static void getAnisoSizeX(
    ma::Mesh* m,
    int dir,
    apf::Field* sizes,
    apf::Field* frames)
{
  PCU_ALWAYS_ASSERT(dir>=0 && dir<3);
  PCU_ALWAYS_ASSERT(apf::getValueType(sizes) == apf::VECTOR);
  PCU_ALWAYS_ASSERT(apf::getValueType(frames) == apf::MATRIX);
  ma::SizeField* isf = new ma::IdentitySizeField(m);
  apf::Field* fLen = createFieldOn(m, "incidentEdgeLength", apf::SCALAR);
  apf::Field* fCnt = createFieldOn(m, "incidentEdgeCount", apf::SCALAR);
  double len;
  int cnt;
  ma::Entity* vtx;
  apf::MeshIterator* itr = m->begin(0);
  while( (vtx = m->iterate(itr)) ) {
    getEdgeLenAndCnt(m, vtx, isf, len, cnt);
    apf::setScalar(fLen, vtx, 0, len);
    apf::setScalar(fCnt, vtx, 0, cnt);
  }
  m->end(itr);
  apf::accumulate(fLen);
  apf::accumulate(fCnt);
  apf::synchronize(fLen);
  apf::synchronize(fCnt);

  itr = m->begin(0);
  while( (vtx = m->iterate(itr)) ) {
    int cnt = static_cast<int>(apf::getScalar(fCnt, vtx, 0));
    double l = apf::getScalar(fLen, vtx, 0) / cnt;
    const double f = 1.8;
    double sz[3] = {l,l,l};
    for(int i=0; i<3; i++)
      if( i == dir )
	sz[i] /= f;

    apf::setVector(sizes, vtx, 0, ma::Vector(sz));
    apf::setMatrix(frames, vtx, 0,
    	ma::Matrix(1.,0.,0.,0.,1.,0.,0.,0.,1.));
  }
  apf::destroyField(fLen);
  apf::destroyField(fCnt);
  delete isf;
}


int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  lion_set_verbosity(1);
  if (argc != 5) {
    if(0==PCU_Comm_Self())
      std::cerr << "usage: " << argv[0] 
        << " <model file> <in mesh> <split direction=[0-2]> <out mesh> \n";
    return EXIT_FAILURE;
  }
#ifdef HAVE_SIMMETRIX
  MS_init();
  SimModel_start();
  Sim_readLicenseFile(NULL);
  gmi_sim_start();
  gmi_register_sim();
#endif
  gmi_register_mesh();
  ma::Mesh* m = apf::loadMdsMesh(argv[1],argv[2]);
  apf::Field* sizes = createFieldOn(m, "sizes", apf::VECTOR);
  apf::Field* frames = createFieldOn(m, "frames", apf::MATRIX);
  getAnisoSizeX(m, atoi(argv[3]) ,sizes, frames);
  ma::Input* in = ma::configure(m, sizes, frames);
  in->shouldRunPreZoltanRib = true;
  in->shouldRunPreParma = true;
  in->shouldRunMidParma = true;
  in->shouldRunPostParma = true;
  in->maximumIterations = 10;
  if (in->shouldSnap) {
    in->shouldSnap = false;
    PCU_ALWAYS_ASSERT(in->shouldTransferParametric);
  }
  in->shouldFixShape = false;
  ma::adapt(in);
  m->verify();
  /* delete ansx; */
  apf::destroyField(sizes);
  apf::destroyField(frames);
  m->writeNative(argv[4]);
  m->destroyNative();
  apf::destroyMesh(m);
#ifdef HAVE_SIMMETRIX
  gmi_sim_stop();
  Sim_unregisterAllKeys();
  SimModel_stop();
  MS_exit();
#endif
  PCU_Comm_Free();
  MPI_Finalize();
}
