#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// Class encapsulating surface fluxes determined  
// from an Earth System Model (e.g. UKESM)

// Features include taking ESM grid points from lat/lon
// and converting to polar stereographic projections
// and converting surface fluxes defined in elevation tiles
// to fluxes defined on the ice sheet grid

#ifndef _ESMSURFACEFLUX_H_
#define _ESMSURFACEFLUX_H_

#include "SurfaceFlux.H"
#include "ParmParse.H"
#include "NamespaceHeader.H"

class ESMSurfaceFlux : public SurfaceFlux
{
public:
  
  ESMSurfaceFlux(ParmParse& a_pp);

  /// destructor
  virtual ~ESMSurfaceFlux() {};

  /// factory method
  /** return a pointer to a new SurfaceFlux object
   */
  virtual SurfaceFlux* new_surfaceFlux();

  /// define source term for thickness evolution and place it in flux
  /** 
   * Is this something that I need to implement?
  */
  virtual void surfaceThicknessFlux(LevelData<FArrayBox>& a_flux,
				                        const AmrIceBase& a_amrIce, 
				                        int a_level, Real a_dt);

protected:

  /// member variables
  std::string m_file_format;
  std::string m_name;

  /// function to read in ESM data
  /// do I need that here? Or should this go into the constructor?

  /// function to convert lat/lon data to polar stereographic
  void regridESMtoBISICLES(const LevelData<FArrayBox>& a_input,
                                       LevelData<FArrayBox>& a_output,
                                       const RealVect& a_dx,
                                       const ProblemDomain& a_domain);

  void regridBISICLEStoESM(const LevelData<FArrayBox>& a_input,
                                       LevelData<FArrayBox>& a_output,
                                       const RealVect& a_dx,
                                       const ProblemDomain& a_domain)

  };

#include "NamespaceFooter.H"

#endif