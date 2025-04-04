#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// concrete class encapsulating surface fluxes as presented 
// by an Earth System Model (e.g. UKESM)

#include "ESMSurfaceFlux.H"
#include "ReadLevelData.H"
#include "FillFromReference.H"
#include "ReflectGhostCells.H"
#include "AmrIceBase.H"
#include "IceConstants.H"
#include "computeSum.H"
#include "NamespaceHeader.H"

SurfaceFlux* ESMSurfaceFlux::new_surfaceFlux() 
{
  // relies on the default copy constructor
  ESMSurfaceFlux* ptr = new ESMSurfaceFlux(*this);
  return static_cast<SurfaceFlux*>(ptr);
}

void ESMSurfaceFlux::surfaceThicknessFlux
(LevelData<FArrayBox>& a_flux,
 const AmrIceBase& a_amrIce, 
 int a_level, Real a_dt)
{
  CH_TIME("ESMSurfaceFlux::surfaceThicknessFlux");
  pout() << "ESMSurfaceFlux::surfaceThicknessFlux " << std::endl;
  
  const RealVect& levelDx = a_amrIce.dx(a_level);
  FillFromReference(a_flux, *m_uniform_source, levelDx ,m_dx, true); // need to replace *m_uniform_source with the correct data
  
  const ProblemDomain& domain = a_flux.disjointBoxLayout().physDomain();
  for (int dir = 0; dir < SpaceDim; ++dir)
    {
      if (!(domain.isPeriodic(dir))){
	ReflectGhostCells(a_flux, domain, dir, Side::Lo);
	ReflectGhostCells(a_flux, domain, dir, Side::Hi);
      }
    }
  a_flux.exchange();

}

#include "NamespaceFooter.H"