#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "SubglacialHydro.H"
#include "AmrIce.H"
#include "NamespaceHeader.H"

SubglacialHydroIceObserver::SubglacialHydroIceObserver()
{
  m_subHydroPtr = new SubglacialHydro();

  pout() << "Creating a SubglacialHydroIceObserver object" << std::endl;

}

SubglacialHydroIceObserver::~SubglacialHydroIceObserver()
{
  if (m_subHydroPtr != NULL)
    {
      delete m_subHydroPtr;
      m_subHydroPtr = NULL;

      pout() << "Deleting the SubglacialHydroIceObserver object" << std::endl;

    }
}

void SubglacialHydroIceObserver::notify(AmrIce::Observer::Notification a_n, AmrIce& a_amrIce)
{

  pout() <<  "SubglacialHydroIceObserver::notify (" << a_n << ")" << std::endl;

  // write some dummy code here to do something in response to a signal from an AmrIce
  // e.g. that could be calling the SUHMO model

}


SubglacialHydro::~SubglacialHydro()
{
  pout() << "Calling SubglacialHydro object destructor" << std::endl;
}

SubglacialHydro::SubglacialHydro()
{
  pout() << "Calling SubglacialHydro object constructor" << std::endl;
}
