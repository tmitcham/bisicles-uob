#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _SUBHYDRO_H_
#define _SUBHYDRO_H_

#include "AmrIce.H"
#include "NamespaceHeader.H"

/// Subglacial hydrology model (a call to SUHMO)

class SubglacialHydro
{

public:
  ~SubglacialHydro();
  SubglacialHydro();

};
  
/// Contains a SubglacialHydro, responds to changes in the ice sheet state
class SubglacialHydroIceObserver : public AmrIce::Observer 
{
  SubglacialHydro* m_subHydroPtr;

public:

  SubglacialHydroIceObserver();
  
  ~SubglacialHydroIceObserver();
  
  /// General purpose notify method: do something in response to a signal from an AmrIce
  void notify(AmrIce::Observer::Notification a_n, AmrIce& a_amrIce);
  
};

#include "NamespaceFooter.H"

#endif
