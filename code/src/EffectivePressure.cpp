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
      ParmParse ppFile("EffectivePressureFromFile");

      std::string file;
      ppFile.get("file", file);

      std::string variable = "effectivePressure";
      ppFile.query("variable", variable);

      ptr = new EffectivePressureFromFile(file, variable);
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

EffectivePressureFromFile::EffectivePressureFromFile(const std::string& a_file,
                                                     const std::string& a_variable)
  : m_file(a_file),
    m_variable(a_variable)
{
}

void
EffectivePressureFromFile::computeN(FArrayBox&            a_N,
                                    const LevelSigmaCS&   a_coords,
                                    const DataIterator&    a_dit,
                                    int                    a_level,
                                    const Box&             a_box) const
{
  
  // Read the data from file
  pout() << "EffectivePressureFromFile: reading N from " << m_file << std::endl;

  // Build the single-variable name list expected by readLevelData.
  Vector<std::string> names(1, m_variable);

  Real dx;

  RefCountedPtr<LevelData<FArrayBox> > TF(new LevelData<FArrayBox>);
  Vector<RefCountedPtr<LevelData<FArrayBox> > > vectData(1, TF);
  readLevelData(vectData,dx,m_file,names,1);

  pout() << "EffectivePressureFromFile: loaded " << dx << std::endl;

  // For a level-0-only file the index is 0.
  // For multi-level files, clamp to the available levels.
  int srcLevel = std::min(a_level, (int)vectData.size() - 1);

  const LevelData<FArrayBox>& srcLD = *vectData[srcLevel];

  // Copy overlapping region.  If the source grids don't match the
  // destination exactly, Chombo's copyTo will do the right thing
  // (intersection-based copy).  Areas of a_box that aren't covered
  // will keep their previous value, so we initialise to zero first.
  a_N.setVal(0.0);

  // Walk the source boxes and copy any intersection.
  for (DataIterator srcDit = srcLD.dataIterator(); srcDit.ok(); ++srcDit)
    {
      const FArrayBox& srcFab = srcLD[srcDit()];
      Box overlap = a_box & srcFab.box();
      if (!overlap.isEmpty())
        {
          a_N.copy(srcFab, overlap, 0, overlap, 0, 1);
        }
    }

  // Ensure non-negative (SUHMO can produce small negative N near margins).
  for (BoxIterator bit(a_box); bit.ok(); ++bit)
    {
      const IntVect& iv = bit();
      if (a_N(iv) < 0.0)
        a_N(iv) = 0.0;
    }

  // Same regularisation as the hydrostatic model.
  a_N += 1.0e-10;
}

#include "NamespaceFooter.H"