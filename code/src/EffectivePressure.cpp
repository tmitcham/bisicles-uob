#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "EffectivePressure.H"
#include "BasalFrictionRelationF_F.H"  // for FORT_BFRICTIONLEGUYEFFPRES
#include "IceConstants.H"
#include "ReadLevelData.H"
#include "FillFromReference.H"
#include "BoxIterator.H"
#include "ParmParse.H"
#include "NamespaceHeader.H"

EffectivePressure*
EffectivePressure::parse(const char* a_prefix)
{
  EffectivePressure* ptr = NULL;
  ParmParse pp(a_prefix);

  std::string type;
  if (!pp.query("effectivePressure", type))
    {
      return NULL;
    }

  if (type == "hydrostatic")
    {
      ParmParse ppHydro("HydrostaticEffectivePressure");

      Real p = 1.0;
      ppHydro.query("p", p);

      Real wtMax = -1.0;
      ppHydro.query("maxTillWaterDepth", wtMax);

      Real pwFactor = 0.97;
      ppHydro.query("tillPressureFactor", pwFactor);

      ptr = new HydrostaticEffectivePressure(p, wtMax, pwFactor);
    }
  else if (type == "LevelData")
    {
      ParmParse ppFile("LevelDataEffectivePressure");

      std::string file;
      ppFile.get("file", file);

      std::string variable = "effectivePressure";
      ppFile.query("variable", variable);

      ptr = new LevelDataEffectivePressure(file, variable);
    }
  else
    {
      MayDay::Error("EffectivePressure::parse -- unknown type");
    }

  return ptr;
}

// This is the same code as in the original PressureLimitedBasalFrictionRelation::comuteAlpha(), 
// but eventually it will be removed from that class and only appear here.
// When that happens, may need to move the fortran routine FORT_BFRICTIONLEGUYEFFPRES 
// to a new location.

void
HydrostaticEffectivePressure::computeN(FArrayBox&            a_N,
                                       const LevelSigmaCS&   a_coords,
                                       const DataIterator&    a_dit,
                                       int                    a_level,
                                       const Box&             a_box) const
{
  const FArrayBox& thckOverFlotation = a_coords.getThicknessOverFlotation()[a_dit];
  const FArrayBox& thck              = a_coords.getH()[a_dit];

  const Real eps = 1.0e-6;

  if (Abs(m_p - 1.0) < eps)
    {
      // p == 1:  N = rho_i * g * (h - h_f)
      a_N.copy(thckOverFlotation);
      a_N *= a_coords.iceDensity() * a_coords.gravity();
    }
  else
    {
      // p != 1:  N = rho_i * g * h^(1-p) * (h - h_f)^p   [Leguy 2014]
      Real rhog = a_coords.iceDensity() * a_coords.gravity();
      FORT_BFRICTIONLEGUYEFFPRES(CHF_FRA1(a_N, 0),
                                 CHF_CONST_FRA1(thckOverFlotation, 0),
                                 CHF_CONST_FRA1(thck, 0),
                                 CHF_CONST_REAL(m_p),
                                 CHF_CONST_REAL(rhog),
                                 CHF_BOX(a_box));
    }

  // prevent zero N
  a_N += 1.0e-10;

  // Optionally reduce N by till-water-depth effect
  // (Van Pelt & Oerlemans 2012, Eqn 2)
  if (m_maxTillWaterDepth > 0.0 && m_tillWaterDepth != NULL)
    {
      const Real& wtMax = m_maxTillWaterDepth;
      const FArrayBox& wt = (*(*m_tillWaterDepth)[a_level])[a_dit];
      for (BoxIterator bit(a_N.box()); bit.ok(); ++bit)
        {
          const IntVect& iv = bit();
          Real pw = m_tillPressureFactor * a_N(iv)
                    * std::min(wt(iv), wtMax) / wtMax;
          a_N(iv) = std::max(0.0, std::min(a_N(iv), a_N(iv) - pw));
        }
    }
}


// New class to read effective pressure from a file, e.g. SUHMO plot file.
// Uses the ReadLevelData infrastructure

LevelDataEffectivePressure::LevelDataEffectivePressure(const std::string& a_file,
                                                       const std::string& a_variable)
  : m_file(a_file),
    m_variable(a_variable),
    m_loaded(false),
    m_dx(0.0)
{
}

void
LevelDataEffectivePressure::loadData() const
{
  if (m_loaded) return;

  pout() << "LevelDataEffectivePressure: reading N from " << m_file << std::endl;

  Vector<std::string> names(1, m_variable);
  readLevelData(m_cachedData, m_dx, m_file, names, 1);

  pout() << "LevelDataEffectivePressure: loaded, dx=" << m_dx
         << ", " << m_cachedData.size() << " level(s)" << std::endl;

  m_loaded = true;
}


void
LevelDataEffectivePressure::fillLevel(const LevelSigmaCS& a_coords,
                                      int                 a_level) const
{
  // Already mapped for this level?
  if (m_mappedData.find(a_level) != m_mappedData.end())
    return;

  loadData();

  // Pick the best available source level
  int srcLevel = std::min(a_level, (int)m_cachedData.size() - 1);
  const LevelData<FArrayBox>& srcLD = *m_cachedData[srcLevel];

  // Build destination LevelData on the BISICLES grids for this level
  const DisjointBoxLayout& destGrids = a_coords.grids();
  RefCountedPtr<LevelData<FArrayBox>> destLD(
    new LevelData<FArrayBox>(destGrids, 1, IntVect::Zero));

  // Source and destination dx
  RealVect srcDx  = m_dx * RealVect::Unit;
  RealVect destDx = a_coords.dx();

  pout() << "LevelDataEffectivePressure::fillLevel: level " << a_level
         << ", srcDx=" << srcDx[0] << ", destDx=" << destDx[0] << std::endl;

  FillFromReference(*destLD, srcLD, destDx, srcDx, true /*verbose*/);

  // Make sure there are no negative values of N
  for (DataIterator dit = destLD->dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& N = (*destLD)[dit];
      const Box& box = destGrids[dit];
      for (BoxIterator bit(box); bit.ok(); ++bit)
        {
          if (N(bit()) < 0.0)  N(bit()) = 0.0;
        }
      N += 1.0e-10;
    }

  m_mappedData[a_level] = destLD;
}


void
LevelDataEffectivePressure::computeN(FArrayBox&             a_N,
                                     const LevelSigmaCS&    a_coords,
                                     const DataIterator&    a_dit,
                                     int                    a_level,
                                     const Box&             a_box) const
{
  // Ensure the mapped data exists for this level
  fillLevel(a_coords, a_level);

  // Copy from the pre-filled mapped data for this patch
  const LevelData<FArrayBox>& mappedLD = *m_mappedData.at(a_level);
  a_N.copy(mappedLD[a_dit], a_box, 0, a_box, 0, 1);
}

#include "NamespaceFooter.H"