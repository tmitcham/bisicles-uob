#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _EFFECTIVEPRESSURE_H_
#define _EFFECTIVEPRESSURE_H_

#include "FArrayBox.H"
#include "LevelData.H"
#include "LevelSigmaCS.H"
#include "RefCountedPtr.H"
#include "AmrIceBase.H"
#include "NamespaceHeader.H"

/// Base class for effective pressure models.
/**
   Provides the effective pressure field N to basal friction
   relations.  
   Concrete classes implement specific ways of computing 
   or obtaining N, either from a hydrostatic approximation,
   or from a file provided by, for example, SUHMO.
*/

class EffectivePressure
{
public:

  virtual ~EffectivePressure() {}

  /// Fill a_N(a_box, 1 component) with effective pressure [Pa].
  /**
     \param a_N      output: cell-centered effective pressure (Pa)
     \param a_coords geometry data for the level (thickness, bed, etc.)
     \param a_dit    DataIterator index for level data access
     \param a_level  AMR level index
     \param a_box    cell-centered box over which to compute
  */
  virtual void computeN(FArrayBox&            a_N,
                        const LevelSigmaCS&   a_coords,
                        const DataIterator&    a_dit,
                        int                    a_level,
                        const Box&             a_box) const = 0;

  virtual EffectivePressure* getNewEffectivePressure() const = 0;

  static EffectivePressure* parse(const char* a_prefix);
};


/// Hydrostatic effective pressure
/**
   Implements the two existing formulae:

   - When p == 1:  \f$ N = \rho_i g (h - h_f) \f$
   - When p != 1:  \f$ N = \rho_i g \, h^{1-p} (h - h_f)^p \f$  (Leguy 2014)

   Also supports the Van Pelt & Oerlemans 2012 till-water-depth
   reduction, as in the original PressureLimitedBasalFrictionRelation.
*/
class HydrostaticEffectivePressure : public EffectivePressure
{
public:

  HydrostaticEffectivePressure(Real a_p                  = 1.0,
                               Real a_maxTillWaterDepth  = -1.0,
                               Real a_tillPressureFactor = 0.97)
    : m_p(a_p),
      m_maxTillWaterDepth(a_maxTillWaterDepth),
      m_tillPressureFactor(a_tillPressureFactor),
      m_tillWaterDepth(NULL)
  {}

  virtual ~HydrostaticEffectivePressure() {}

  /// Allow the thermodynamics module to supply till water depth.
  void setTillWaterDepth(Vector<LevelData<FArrayBox>*>* a_twd)
  {
    m_tillWaterDepth = a_twd;
  }

  virtual void computeN(FArrayBox&            a_N,
                        const LevelSigmaCS&   a_coords,
                        const DataIterator&    a_dit,
                        int                    a_level,
                        const Box&             a_box) const;

  virtual EffectivePressure* getNewEffectivePressure() const
  {
    return new HydrostaticEffectivePressure(m_p, m_maxTillWaterDepth,
                                            m_tillPressureFactor);
  }

private:
  Real m_p;                  // exponent (!= 1 => Leguy 2014)
  Real m_maxTillWaterDepth;  // Van Pelt & Oerlemans till reduction (-1 = off)
  Real m_tillPressureFactor; // alpha in Van Pelt & Oerlemans Eqn 2
  Vector<LevelData<FArrayBox>*>* m_tillWaterDepth;
};


/// Effective pressure read from a SUHMO (or any Chombo HDF5) plot file.
/**
   Reads N from a variable in a Chombo AMR HDF5 file.  SUHMO writes
   effectivePressure directly.
*/
class EffectivePressureFromFile : public EffectivePressure
{
public:

  EffectivePressureFromFile(const std::string& a_file,
                            const std::string& a_variable = "effectivePressure");

  virtual ~EffectivePressureFromFile() {}

  virtual void computeN(FArrayBox&            a_N,
                        const LevelSigmaCS&   a_coords,
                        const DataIterator&    a_dit,
                        int                    a_level,
                        const Box&             a_box) const;

  virtual EffectivePressure* getNewEffectivePressure() const
  {
    return new EffectivePressureFromFile(m_file, m_variable);
  }

private:

  std::string m_file;
  std::string m_variable;

};


#include "NamespaceFooter.H"
#endif