#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "MountainIBC.H"
#include "ParmParse.H"
#include "BoxIterator.H"
#include "LevelSigmaCS.H"

#include "NamespaceHeader.H"

MountainIBC::MountainIBC() : m_boundaryThickness(0.0)
{
  m_isBCsetUp = false;
  m_isDefined = false;
}

MountainIBC::~MountainIBC()
{
}

void
MountainIBC::define(const ProblemDomain& a_domain,
                    const Real&          a_dx)
{
  PhysIBC::define(a_domain, a_dx);
}

IceThicknessIBC*
MountainIBC::new_thicknessIBC()
{
  MountainIBC* retval = new MountainIBC();

  retval->m_boundaryThickness = m_boundaryThickness;

  return static_cast<IceThicknessIBC*>(retval);
}

void
MountainIBC::initialize(LevelData<FArrayBox>& a_U)
{
  MayDay::Error("MountainIBC::initialize not implemented");
}

/// Set up initial ice geometry
/**
   Bed topography: tilted plane + Gaussian bumps (the "dino" shape).
   Matches the geometry initialised by SUHMO's MountainIBC so that
   both models start from a consistent state.

   Parameters (hard-coded to match SUHMO):
     ax = 1.5e-3, by = -1.5e-3, cst = 100 m   (tilted base)
     Five Gaussian features ("dino" shape)

   Ice thickness: H = max( 2*(ax*x + by*y + cst + 100), 0 ) m
*/
void
MountainIBC::initializeIceGeometry(LevelSigmaCS& a_coords,
                                   const RealVect& a_dx,
                                   const RealVect& a_domainSize,
                                   const Real& a_time,
                                   const LevelSigmaCS* a_crseCoords,
                                   const int a_refRatio)
{
  const DisjointBoxLayout& grids = a_coords.grids();
  const RealVect& dx = a_coords.dx();

  LevelData<FArrayBox>& levelH  = a_coords.getH();
  const LevelData<FArrayBox>& levelZbRef = a_coords.getTopography();
  LevelData<FArrayBox> levelZb(grids, 1, levelZbRef.ghostVect());

  a_coords.setSeaLevel(0.0);

  // Parameters shared with SUHMO MountainIBC
  const Real ax  =  1.5e-3;
  const Real by  = -1.5e-3;
  const Real cst = 100.0;

  DataIterator dit = grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& zB = levelZb[dit];
      FArrayBox& H  = levelH[dit];

      BoxIterator bit(zB.box());
      for (bit.begin(); bit.ok(); ++bit)
        {
          IntVect iv = bit();
          Real x_loc = (iv[0] + 0.5) * dx[0];
          Real y_loc = (iv[1] + 0.5) * dx[1];

          // ── Bed topography ──────────────────────────────────────────────

          // Step 1: tilted base
          Real step_1 = std::max(ax*x_loc + by*y_loc + cst, 0.0);

          // Step 2: middle finger of the dino
          Real step_2 = 0.0;
          {
            Real angle = 35.0;
            Real x_bd = x_loc - 100000.0;
            Real y_bd = y_loc;
            Real x_t  =  x_bd * std::cos(angle * M_PI / 180.0)
                       + y_bd * std::sin(angle * M_PI / 180.0);
            Real y_t  = -x_bd * std::sin(angle * M_PI / 180.0)
                       + y_bd * std::cos(angle * M_PI / 180.0);
            Real A_gauss = 250.0 * std::exp(-0.5 / (50000.0*50000.0) * y_t*y_t);
            Real sigma   = 12000.0
                         - 3000.0 * std::min(1.0 - (50000.0 - y_t)/50000.0, 1.0);
            step_2 = A_gauss * std::exp(-0.5 / (sigma*sigma) * x_t*x_t);
          }

          // Step 3: side finger of the dino
          Real step_3 = 0.0;
          {
            Real angle = 90.0;
            Real x_bd = x_loc - 85000.0;
            Real y_bd = y_loc;
            Real x_t  =  x_bd * std::cos(angle * M_PI / 180.0)
                       + y_bd * std::sin(angle * M_PI / 180.0);
            Real y_t  = -x_bd * std::sin(angle * M_PI / 180.0)
                       + y_bd * std::cos(angle * M_PI / 180.0);
            Real A_gauss = 250.0 * std::exp(-0.5 / (30000.0*30000.0) * y_t*y_t);
            step_3 = A_gauss * std::exp(-0.5 / (10000.0*10000.0) * x_t*x_t);
          }

          // Step 4: subsidiary finger of the dino
          Real step_4 = 0.0;
          {
            Real angle = 60.0;
            Real x_bd = x_loc - 55000.0;
            Real y_bd = y_loc;
            Real x_t  =  x_bd * std::cos(angle * M_PI / 180.0)
                       + y_bd * std::sin(angle * M_PI / 180.0);
            Real y_t  = -x_bd * std::sin(angle * M_PI / 180.0)
                       + y_bd * std::cos(angle * M_PI / 180.0);
            Real A_gauss = 100.0 * std::exp(-0.5 / (20000.0*20000.0) * y_t*y_t);
            step_4 = A_gauss * std::exp(-0.5 / (5000.0*5000.0) * x_t*x_t);
          }

          // Step 5: other side of the dino
          Real step_5 = 0.0;
          {
            Real angle = 2.0;
            Real x_bd = x_loc - 100000.0;
            Real y_bd = y_loc - 20000.0;
            Real x_t  =  x_bd * std::cos(angle * M_PI / 180.0)
                       + y_bd * std::sin(angle * M_PI / 180.0);
            Real y_t  = -x_bd * std::sin(angle * M_PI / 180.0)
                       + y_bd * std::cos(angle * M_PI / 180.0);
            Real A_gauss = 300.0 * std::exp(-0.5 / (35000.0*35000.0) * y_t*y_t);
            step_5 = A_gauss * std::exp(-0.5 / (6000.0*6000.0) * x_t*x_t);
          }

          zB(iv, 0) = std::max(step_1 + step_2 + step_3 + step_4 + step_5,
                               0.0);

          // ── Ice thickness ─────────────────────────────────────────────
          // Linear profile matching SUHMO initializeData (no random noise)
          H(iv, 0) = std::max(2.0*(ax*x_loc + by*y_loc + cst + 100.0), 0.0);

        } // end cell loop
    } // end box loop

  a_coords.setTopography(levelZb);
}

void
MountainIBC::primBC(FArrayBox&            a_WGdnv,
                    const FArrayBox&      a_Wextrap,
                    const FArrayBox&      a_W,
                    const int&            a_dir,
                    const Side::LoHiSide& a_side,
                    const Real&           a_time)
{
  // Clamp thickness to >= 0 at domain boundaries (no ice inflow).
  if (!m_domain.isPeriodic(a_dir))
    {
      int lohisign = sign(a_side);
      Box tmp = a_WGdnv.box();
      tmp.shiftHalf(a_dir, lohisign);

      // Guard against multi-layer ghost regions (see MarineIBC::primBC comment)
      Box ghostBox = (a_side == Side::Lo)
                   ? adjCellLo(m_domain.domainBox(), a_dir, 1)
                   : adjCellHi(m_domain.domainBox(), a_dir, 1);
      ghostBox &= tmp;

      if (!ghostBox.isEmpty() && !m_domain.contains(tmp))
        {
          tmp &= m_domain;
          Box boundaryBox = (a_side == Side::Lo) ? bdryLo(tmp, a_dir)
                                                 : bdryHi(tmp, a_dir);
          BoxIterator bit(boundaryBox);
          for (bit.begin(); bit.ok(); ++bit)
            {
              const IntVect& iv = bit();
              a_WGdnv(iv, 0) = std::max(0.0, a_Wextrap(iv, 0));
            }
        }
    }
}

void
MountainIBC::setBdrySlopes(FArrayBox&       a_dW,
                           const FArrayBox& a_W,
                           const int&       a_dir,
                           const Real&      a_time)
{
  // One-sided differences are fine; do nothing.
}

void
MountainIBC::artViscBC(FArrayBox&       a_F,
                       const FArrayBox& a_U,
                       const FArrayBox& a_divVel,
                       const int&       a_dir,
                       const Real&      a_time)
{
  MayDay::Error("MountainIBC::artViscBC not implemented");
}

RefCountedPtr<CompGridVTOBC>
MountainIBC::velocitySolveBC()
{
  if (!m_isBCsetUp)
    {
      setupBCs();
    }
  return m_velBCs;
}


void
MountainIBC::setupBCs()
{
  ParmParse ppBC("bc");
  Vector<int> loBCvect(SpaceDim, 1);  // default: Neumann
  Vector<int> hiBCvect(SpaceDim, 1);
  ppBC.queryarr("lo_bc", loBCvect, 0, SpaceDim);
  ppBC.queryarr("hi_bc", hiBCvect, 0, SpaceDim);

  m_velBCs = RefCountedPtr<CompGridVTOBC>(new CompGridVTOBC());

  for (int dir = 0; dir < SpaceDim; ++dir)
    {
      // Low side
      for (int comp = 0; comp < SpaceDim; ++comp)
        {
          // 0 = Dirichlet, 1 = Neumann, 2 = slip-wall (Dirichlet normal only)
          if (loBCvect[dir] == 0)
            {
              m_velBCs->m_bcDiri[dir][0][comp] = true;
            }
          else if (loBCvect[dir] == 1)
            {
              m_velBCs->m_bcDiri[dir][0][comp] = false;
            }
          else if (loBCvect[dir] == 2)
            {
              m_velBCs->m_bcDiri[dir][0][comp] = (comp == dir);
            }
          else
            {
              MayDay::Error("MountainIBC: unknown lo_bc type");
            }
        }

      // High side
      for (int comp = 0; comp < SpaceDim; ++comp)
        {
          if (hiBCvect[dir] == 0)
            {
              m_velBCs->m_bcDiri[dir][1][comp] = true;
            }
          else if (hiBCvect[dir] == 1)
            {
              m_velBCs->m_bcDiri[dir][1][comp] = false;
            }
          else if (hiBCvect[dir] == 2)
            {
              m_velBCs->m_bcDiri[dir][1][comp] = (comp == dir);
            }
          else
            {
              MayDay::Error("MountainIBC: unknown hi_bc type");
            }
        }
    }

  m_isBCsetUp = true;
}

#include "NamespaceFooter.H"
