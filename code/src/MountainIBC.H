#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

//
//  MountainIBC.H
// ============
//
// IceThicknessIBC-derived class implementing the "dino" mountain geometry
// used in the SUHMO coupled subglacial hydrology test case.
//
// Bed topography is a tilted plane with superimposed Gaussian bumps
// (the "dino" shape). Ice thickness is a linear function of position.
// Intended for use in coupled BISICLES-SUHMO test problems.
//

#ifndef _MOUNTAINIBC_H_
#define _MOUNTAINIBC_H_

#include "IceThicknessIBC.H"
#include "ReflectGhostCells.H"


#include "NamespaceHeader.H"

/// Mountain ("dino") geometry IBC for coupled BISICLES-SUHMO test cases
/**
   Implements the same bed topography and ice thickness profile as
   SUHMO's MountainIBC, enabling consistent geometry across both models.

   Bed: tilted plane (ax*x + by*y + cst) with Gaussian bump features.
   Ice thickness: linear H = 2*(ax*x + by*y + cst + 100), clamped to >= 0.

   Velocity boundary conditions are read from bc.lo_bc / bc.hi_bc
   (0 = Dirichlet, 1 = Neumann, 2 = slip-wall) per direction,
   with bc.new_bc = true to enable per-direction control.
*/
class MountainIBC : public IceThicknessIBC
{
public:
  /// Constructor
  MountainIBC();

  /// Destructor
  virtual ~MountainIBC();

  /// Define the object (calls base PhysIBC::define)
  virtual void define(const ProblemDomain& a_domain,
                      const Real&          a_dx);

  /// Factory method
  virtual IceThicknessIBC* new_thicknessIBC();

  /// initialize() is not used for this IBC
  virtual void initialize(LevelData<FArrayBox>& a_U);

  /// Set up initial ice geometry (bed topography + thickness)
  virtual void initializeIceGeometry(LevelSigmaCS& a_coords,
                                     const RealVect& a_dx,
                                     const RealVect& a_domainSize,
                                     const Real& a_time,
                                     const LevelSigmaCS* a_crseCoords,
                                     const int a_refRatio);

  /// Set boundary fluxes (no-inflow at boundaries)
  virtual void primBC(FArrayBox&            a_WGdnv,
                      const FArrayBox&      a_Wextrap,
                      const FArrayBox&      a_W,
                      const int&            a_dir,
                      const Side::LoHiSide& a_side,
                      const Real&           a_time);

  /// Set boundary slopes (one-sided differences, do nothing)
  virtual void setBdrySlopes(FArrayBox&       a_dW,
                             const FArrayBox& a_W,
                             const int&       a_dir,
                             const Real&      a_time);

  /// Artificial viscosity BC (not used)
  virtual void artViscBC(FArrayBox&       a_F,
                         const FArrayBox& a_U,
                         const FArrayBox& a_divVel,
                         const int&       a_dir,
                         const Real&      a_time);

  /// Return velocity solve BC
  virtual RefCountedPtr<CompGridVTOBC> velocitySolveBC();

protected:

   // ice sheet thickness at domain boundaries
  Real          m_boundaryThickness;

  RealVect m_domainSize;

  // have bc's been set up?
  bool m_isBCsetUp;

  /// set up bc's based on parmparse inputs
  void setupBCs();

  RefCountedPtr<CompGridVTOBC> m_velBCs;

private:
  void operator=(const MountainIBC&) { MayDay::Error("invalid operator"); }
  MountainIBC(const MountainIBC&)    { MayDay::Error("invalid operator"); }
};

#include "NamespaceFooter.H"
#endif
