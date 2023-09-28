#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "AMRSubglacial.H"
#include "DamageConstitutiveRelation.H"
#include "AmrIce.H"
#include "FineInterp.H"
#include "CoarseAverageFace.H"
#include "AdvectPhysics.H"
#include "PatchGodunov.H"
#include "PiecewiseLinearFillPatch.H"
#include "DivergenceF_F.H"
#include "ParmParse.H"
#include "LevelMappedDerivatives.H"
#include "NyeCrevasseF_F.H"
#include "SigmaCSF_F.H"
#include "CH_HDF5.H"
//Michael added
#include "CellToEdge.H"
#include "BoxIterator.H"
#include "computeSum.H"
//#include "BiCGStabSolver.H"
#include "RelaxSolver.H"
#include "VCAMRPoissonOp2.H"
#include "AMRMultiGrid.H"
#include "MultilevelLinearOp.H"
#include "IceConstants.H"

#include<cstdio>

#include "NamespaceHeader.H"


//TODO: Introduce member variables m_tol (H==0 tolerance)

/*------------------------------------------------------------------------------*/
//		Start: Define Operator functions
extern
AMRLevelOpFactory<LevelData<FArrayBox> >*
defineOperatorFactory_phi(
                      Vector<DisjointBoxLayout>&             a_grids,
                      Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_aCoef,
                      Vector<RefCountedPtr<LevelData<FluxBox> > >& a_bCoef,
                      SubglacialParameters&                  a_params)
{

  VCAMRPoissonOp2Factory* opFactory = new VCAMRPoissonOp2Factory;

  opFactory->define(a_params.coarsestDomain,
                    a_grids,
                    a_params.refRatio,
                    a_params.coarsestDx,
                    &ParseBC_phi,
                    a_params.alpha,
                    a_aCoef,
                    a_params.beta,
                    a_bCoef);

  if (a_params.coefficient_average_type >= 0)
    {
      opFactory->m_coefficient_average_type
        = a_params.coefficient_average_type;
    }

  return (AMRLevelOpFactory<LevelData<FArrayBox> >*) opFactory;
}


extern
AMRLevelOpFactory<LevelData<FArrayBox> >*
defineOperatorFactory_hcav(
                      Vector<DisjointBoxLayout>&               a_grids,
                      Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_aCoef,
                      Vector<RefCountedPtr<LevelData<FluxBox> > >&   a_bCoef,
                      SubglacialParameters&                     a_params)
{

  VCAMRPoissonOp2Factory* opFactory = new VCAMRPoissonOp2Factory;

  opFactory->define(a_params.coarsestDomain,
                    a_grids,
                    a_params.refRatio,
                    a_params.coarsestDx,
                    &ParseBC_hcav,
                    a_params.alpha,
                    a_aCoef,
                    a_params.beta,
                    a_bCoef);

  if (a_params.coefficient_average_type >= 0)
    {
      opFactory->m_coefficient_average_type
        = a_params.coefficient_average_type;
    }

  return (AMRLevelOpFactory<LevelData<FArrayBox> >*) opFactory;
}



void ParseValue(Real* pos,
                int* dir,
                Side::LoHiSide* side,
                Real* a_values)
{
  Real bcVal = 0.0;
  a_values[0]=bcVal;
}





void ParseBC_phi(FArrayBox& a_state,
             const Box& a_valid,
             const ProblemDomain& a_domain,
             Real a_dx,
             bool a_homogeneous)
{
  if (!a_domain.domainBox().contains(a_state.box()))
    {

      // if (!GlobalBCRS::s_areBCsParsed_phi)
      //   {
      //     ParmParse pp;
      //     pp.getarr("bc_lo_phi", GlobalBCRS::s_bcLo_phi, 0, SpaceDim);
      //     pp.getarr("bc_hi_phi", GlobalBCRS::s_bcHi_phi, 0, SpaceDim);
      //     GlobalBCRS::s_areBCsParsed_phi = true;
      //   }

      Box valid = a_valid;
      for (int i=0; i<CH_SPACEDIM; ++i)
        {
          Box ghostBoxLo = adjCellBox(valid, i, Side::Lo, 1);
          Box ghostBoxHi = adjCellBox(valid, i, Side::Hi, 1);
          if (!a_domain.domainBox().contains(ghostBoxLo))
            {
             NeumBC(a_state,
                         valid,
                         a_dx,
                         a_homogeneous,
                         ParseValue,
                         i,
                         Side::Lo);
            }

          if (!a_domain.domainBox().contains(ghostBoxHi))
            {
             NeumBC(a_state,
                         valid,
                         a_dx,
                         a_homogeneous,
                         ParseValue,
                         i,
                         Side::Hi);
            }
        }
    }
}


void ParseBC_hcav(FArrayBox& a_state,
             const Box& a_valid,
             const ProblemDomain& a_domain,
             Real a_dx,
             bool a_homogeneous)
{
  if (!a_domain.domainBox().contains(a_state.box()))
    {

      // if (!GlobalBCRS::s_areBCsParsed_phi)
      //   {
      //     ParmParse pp;
      //     pp.getarr("bc_lo_phi", GlobalBCRS::s_bcLo_phi, 0, SpaceDim);
      //     pp.getarr("bc_hi_phi", GlobalBCRS::s_bcHi_phi, 0, SpaceDim);
      //     GlobalBCRS::s_areBCsParsed_phi = true;
      //   }

      Box valid = a_valid;
      for (int i=0; i<CH_SPACEDIM; ++i)
        {
          Box ghostBoxLo = adjCellBox(valid, i, Side::Lo, 1);
          Box ghostBoxHi = adjCellBox(valid, i, Side::Hi, 1);
          if (!a_domain.domainBox().contains(ghostBoxLo))
            {
             NeumBC(a_state,
                         valid,
                         a_dx,
                         a_homogeneous,
                         ParseValue,
                         i,
                         Side::Lo);
            }

          if (!a_domain.domainBox().contains(ghostBoxHi))
            {
             NeumBC(a_state,
                         valid,
                         a_dx,
                         a_homogeneous,
                         ParseValue,
                         i,
                         Side::Hi);
            }
        }
    }
}


//		End: Operator functions
/*------------------------------------------------------------------------------*/
//		Start: Observer functions


SubglacialIceObserver::SubglacialIceObserver()
{
  m_subglacialPtr = new AMRSubglacial();
}

SubglacialIceObserver::~SubglacialIceObserver()
{
  if (m_subglacialPtr != NULL)
    {
      delete m_subglacialPtr;
      m_subglacialPtr = NULL;
    }
}


AMRSubglacial& SubglacialIceObserver::subglacial() const
{
  return *m_subglacialPtr;
}



void SubglacialIceObserver::notify(AmrIce::Observer::Notification a_n, AmrIce& a_amrIce)
{

  pout() <<  "SubglacialIceObserver::notify" << std::endl;

  if (a_n == AmrIce::Observer::PreVelocitySolve)
    {
      // TODO: the velocity solve will need the effective pressure field
      m_subglacialPtr->define(a_amrIce.grids(), a_amrIce.refRatios(),
			  a_amrIce.finestLevel(), a_amrIce.dx(0), a_amrIce.amrGeometry());
      //TODO: update effective pressure in friction object based on new grids
      if (m_subglacialPtr->m_effectivePressureDefined)
        {
          // Commented out by TM
          //a_amrIce.setEffectivePressure(m_subglacialPtr->getEffectivePressure());
          //a_amrIce.setTimeBasalFriction();
          pout()<<"Called set effective pressure in PreVelocitySolve"<<endl;
        };
    }
  else if (a_n == AmrIce::Observer::PostVelocitySolve)
    {
      // the velocity solve will need the pressurefield (add method to get N form m_hydraulicPotential) field
      const Vector<LevelData<FluxBox>* >& faceVel = a_amrIce.faceVelocities();
      const Vector<LevelData<FArrayBox>* >& cellVel = a_amrIce.amrVelocity();
      const Vector<RefCountedPtr<LevelSigmaCS> >& geometry = a_amrIce.amrGeometry();

      Real m_time_step = a_amrIce.dt();//years
      // First velocity solve happens before the timestep is initialised, in
      // this case don't run subglacial model
      if (m_time_step <10e10)
        {
         //Set initial state, hopefully this doesn't happen multiple times
          if(m_subglacialPtr->m_time_step ==0)
            {
              m_subglacialPtr->setInitialValues(a_amrIce, geometry);
              m_subglacialPtr->m_time_step +=1;
              m_subglacialPtr->m_effectivePressureDefined = true;
              // Commented out by TM
              //a_amrIce.setEffectivePressure(m_subglacialPtr->getEffectivePressure());
              pout()<<"Setting initial values"<<endl;
            }
          else
            {
              m_subglacialPtr->runLargeTimestep(m_time_step, a_amrIce, geometry, cellVel, faceVel);
              // Commented out by TM
              //a_amrIce.setEffectivePressure(m_subglacialPtr->getEffectivePressure());
            }
        } 

    }


}

void SubglacialIceObserver::setEffectivePressure(Vector<LevelData<FArrayBox>*> a_effectivePressure)
{
   m_subglacialPtr->setObserverEffectivePressure(a_effectivePressure);
   m_effectivePressureDefined = true;
}



#ifdef CH_USE_HDF5
void SubglacialIceObserver::addPlotVars(Vector<std::string>& a_vars)
{
  m_subglacialPtr->addPlotVars(a_vars);
}

void SubglacialIceObserver::writePlotData(LevelData<FArrayBox>& a_data, int a_level)
{
  m_subglacialPtr->writePlotData(a_data, a_level);
}

/// fill a_var with the names of variables to add to the checkpoint file
void SubglacialIceObserver::addCheckVars(Vector<std::string>& a_vars)
{
  m_subglacialPtr->addCheckVars(a_vars);
}

/// copy level a_level checkpoint data to  LevelData<FArrayBox>& a_data
void SubglacialIceObserver::writeCheckData(HDF5Handle& a_handle, int a_level)
{
  m_subglacialPtr->writeCheckData(a_handle, a_level);
}

/// read level a_level checkpoint data from  LevelData<FArrayBox>& a_data
void SubglacialIceObserver::readCheckData(HDF5Handle& a_handle, HDF5HeaderData& a_header, int a_level, const DisjointBoxLayout& a_grids)
{
  m_subglacialPtr->readCheckData(a_handle, a_header, a_level, a_grids);
}
#endif


AMRSubglacial::~AMRSubglacial()
{

  for (int lev = 0; lev < m_hydraulicPotential.size(); lev++)
  {
    if (m_hydraulicPotential[lev] != NULL)
	  {
	  delete m_hydraulicPotential[lev]; m_hydraulicPotential[lev] = NULL;
	  }
  }

  for (int lev = 0; lev < m_cavityHeight.size(); lev++)
  {
    if (m_cavityHeight[lev] != NULL)
	  {
	  delete m_cavityHeight[lev]; m_cavityHeight[lev] = NULL;
	  }
  }
  
  for (int lev = 0; lev < m_effectivePressure.size(); lev++)
    {
      if (m_effectivePressure[lev] != NULL)
	{
	  delete m_effectivePressure[lev]; m_effectivePressure[lev] = NULL;
	}
    }
}
AMRSubglacial::AMRSubglacial()
{
  m_time_step = 0;
  m_time = 0.0;
}

const Vector<LevelData<FArrayBox>*> AMRSubglacial::getEffectivePressure()
{
  return m_effectivePressure;
}

const LevelData<FArrayBox>* AMRSubglacial::hydraulicPotential(int a_level) const
{
  if (!(m_hydraulicPotential.size() > a_level))
    {
      std::string msg("AMRSubglacial::hydraulicPotential !(m_hydraulicPotential.size() > a_level)");
      pout() << msg <<endl;
      CH_assert((m_hydraulicPotential.size() > a_level));
      MayDay::Error(msg.c_str());
    }

  LevelData<FArrayBox>* ptr = m_hydraulicPotential[a_level];
  if (ptr == NULL)
    {
      std::string msg("AMRSubglacial::hydraulicPotential m_hydraulicPotential[a_level] == NULL");
      pout() << msg << endl;
      CH_assert(ptr != NULL);
      MayDay::Error(msg.c_str());
    }
  return ptr;
}


const LevelData<FArrayBox>* AMRSubglacial::cavityHeight(int a_level) const
{
  if (!(m_cavityHeight.size() > a_level))
    {
      std::string msg("AMRSubglacial::cavityHeight !(m_cavityHeight.size() > a_level)");
      pout() << msg <<endl;
      CH_assert((m_cavityHeight.size() > a_level));
      MayDay::Error(msg.c_str());
    }

  LevelData<FArrayBox>* ptr = m_cavityHeight[a_level];
  if (ptr == NULL)
    {
      std::string msg("AMRSubglacial::cavityHeight m_cavityHeight[a_level] == NULL");
      pout() << msg << endl;
      CH_assert(ptr != NULL);
      MayDay::Error(msg.c_str());
    }
  return ptr;
}

const LevelData<FArrayBox>* AMRSubglacial::effectivePressure(int a_level) const
{
  if (!(m_effectivePressure.size() > a_level))
    {
      std::string msg("AMRSubglacial::effectivePressure !(m_effectivePressure.size() > a_level)");
      pout() << msg <<endl;
      CH_assert((m_effectivePressure.size() > a_level));
      MayDay::Error(msg.c_str());
    }

  LevelData<FArrayBox>* ptr = m_effectivePressure[a_level];
  if (ptr == NULL)
    {
      std::string msg("AMRSubglacial::effectivePressure m_effectivePressure[a_level] == NULL");
      pout() << msg << endl;
      CH_assert(ptr != NULL);
      MayDay::Error(msg.c_str());
    }
  return ptr;
}

void AMRSubglacial::setInitialValues(
       AmrIce& a_amrIce,
       const Vector<RefCountedPtr<LevelSigmaCS> >& a_geometry)
{
//TODO Michael: For now we set the initial hydraulic potential value, hcav will stay as 0.5
// Later, set from data     
    // fill ghost data
    for (int lev=0; lev<= m_finestLevel; lev++)
    {
      if (lev > 0)
	    {

	     int nGhost = m_hydraulicPotential[lev]->ghostVect()[0];
	     PiecewiseLinearFillPatch pwl(m_grids[lev],  m_grids[lev-1], 1,
				       m_grids[lev].physDomain(),m_ratio[lev-1],nGhost);
	     // since we're not subcycling, don't need to interpolate in time
	     Real time_interp_coeff = 0.0;
	     pwl.fillInterp(*m_hydraulicPotential[lev], *m_hydraulicPotential[lev-1], *m_hydraulicPotential[lev-1],
			 time_interp_coeff,0, 0, 1);

             int nGhost_hcav = m_cavityHeight[lev]->ghostVect()[0];
	     PiecewiseLinearFillPatch pwl_hcav(m_grids[lev],  m_grids[lev-1], 1,
				       m_grids[lev].physDomain(),m_ratio[lev-1],nGhost_hcav);
	     // since we're not subcycling, don't need to interpolate in time
	     Real time_interp_coeff_hcav = 0.0;
	     pwl_hcav.fillInterp(*m_cavityHeight[lev], *m_cavityHeight[lev-1], *m_cavityHeight[lev-1],
			 time_interp_coeff_hcav,0, 0, 1);
             
             int nGhost_N = m_effectivePressure[lev]->ghostVect()[0];
	     PiecewiseLinearFillPatch pwl_N(m_grids[lev],  m_grids[lev-1], 1,
				       m_grids[lev].physDomain(),m_ratio[lev-1],nGhost_N);
	     // since we're not subcycling, don't need to interpolate in time
	     Real time_interp_coeff_N = 0.0;
	     pwl_N.fillInterp(*m_effectivePressure[lev], *m_effectivePressure[lev-1], *m_effectivePressure[lev-1],
               time_interp_coeff_N,0, 0, 1);

	    }
      m_hydraulicPotential[lev]->exchange();
      m_cavityHeight[lev]->exchange();
      m_effectivePressure[lev]->exchange();
    }
    for (int lev = 0; lev <= m_finestLevel; lev++)
       	{
          const LevelSigmaCS& levelCoords = *a_geometry[lev];

          Real iceDensity = levelCoords.iceDensity();
          Real waterDensity = levelCoords.waterDensity();

          DataIterator dit = (*m_hydraulicPotential[lev]).dataIterator();
       	  for (dit.begin(); dit.ok(); ++dit)
       	    {
       	      FArrayBox& hydraulicPotential = (*m_hydraulicPotential[lev])[dit];
              FArrayBox& effectivePressure  = (*m_effectivePressure[lev])[dit];
             
              const FArrayBox& thckOverFlotation = levelCoords.getThicknessOverFlotation()[dit];
              const FArrayBox& thisZb = levelCoords.getTopography()[dit];
              const FArrayBox& thisZs = levelCoords.getSurfaceHeight()[dit];
              const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
              Real tol = 0.001;


              effectivePressure.copy( thckOverFlotation );
              effectivePressure *= iceDensity * m_g;
              //TODO Michael: Not ideal, but how do I avoid the "not in box issue?"
              Box thisBox = thisZs.box();

              BoxIterator bit(thisBox);
              for (bit.begin(); bit.ok(); ++bit)
                {
                  const IntVect& iv = bit();
                  if (mask(iv) == GROUNDEDMASKVAL)
                    {
                      //TODO: at the moment N(0) = 100Pa (small, arbitary), intialise from observations?
                      //effectivePressure(iv)  = iceDensity*m_g*thckOverFlotation(iv);
                      hydraulicPotential(iv) = iceDensity*m_g*thisZs(iv) + m_g *(waterDensity - iceDensity)*thisZb(iv) - effectivePressure(iv)+100;
                    }
                  else 
                    {
                      effectivePressure(iv)  = 1.0;
                      hydraulicPotential(iv) = 1.0;
                    }
                }
	     }
         }
    for (int lev = 0; lev <= m_finestLevel; lev++)
       	{
          const LevelSigmaCS& levelCoords = *a_geometry[lev];
        
          DataIterator dit = (*m_cavityHeight[lev]).dataIterator();
       	  for (dit.begin(); dit.ok(); ++dit)
       	    {
       	      FArrayBox& cavityHeight = (*m_cavityHeight[lev])[dit];
              
              const FArrayBox& thisZs = levelCoords.getSurfaceHeight()[dit];
              Box thisBox = thisZs.box();

              BoxIterator bit(thisBox);
              for (bit.begin(); bit.ok(); ++bit)
                {
                  const IntVect& iv = bit();
                  cavityHeight(iv) =  0.2;                    
                }
	     }
         }

}

void AMRSubglacial::define(const Vector<DisjointBoxLayout>& a_grids,
 const Vector<int>& a_ratio,
 int a_finestLevel,
 const RealVect& a_crseDx,
 const Vector<RefCountedPtr<LevelSigmaCS> >& a_geometry)
{
  
  ParmParse ppSub("subglacial_params");
  m_Ahat = 2e-17; //Pa^(-n)/yr
  ppSub.query("Ahat",m_Ahat );
  m_g = 9.8;//9.7e15 m/yr^2
  ppSub.query("g",m_g );
  m_n = 3.0; //Glen's flow law exponent
  ppSub.query("n",m_n );


  m_sigma = 0.001; // 10e-4
  ppSub.query("sigma",m_sigma );
  m_basalMelt = 10e-4;//basal melt m/yr
  ppSub.query("basalMelt",m_basalMelt );
  m_alpha = 3.0; //hydraulic permeability exponent
  ppSub.query("alpha",m_alpha );
  m_epsilon_d = 1e5;//m^2 /yr
  ppSub.query("epsilon_d",m_epsilon_d );
  m_h0  = 1.0;  // typical cavity thickness m
  ppSub.query("h0",m_h0 );
  m_hr = 1.0; // obstacale height m
  ppSub.query("hr",m_hr );
  m_K = 9.45e6; //9.45e6 typical hydraulic permeability m/yr
  ppSub.query("K",m_K );
  m_lr = 5.0; //obstacle length scale m
  ppSub.query("lr",m_lr );
  //m_N0 = 9e5;//9e5typical hydraulic pressure scale Pa (time comp?)
  m_N0 = 9e5;//9e5typical hydraulic pressure scale Pa (time comp?)
  ppSub.query("N0",m_N0 );

  m_BigaPhiCoef = 1000000000; // Big number phi boundary condition
  ppSub.query("BigaPhiCoef",m_BigaPhiCoef );
  m_BigNotGrounded = 100001; // Big number for outside domain conditions
  ppSub.query("BigNotGrounded",m_BigNotGrounded );
  
  if (m_dx.size() > 0)
    {
      if ( !(m_dx[0] == a_crseDx) )
	     {
	      std::string msg("AMRSubglacial::define, incompatible mesh");
	      pout() << msg << std::endl;
	      CH_assert(m_dx[0] == a_crseDx);
	      MayDay::Error(msg.c_str());
	     }
    }

  //update the mesh hierarchy
  m_finestLevel = a_finestLevel;
  m_grids.resize(m_finestLevel  + 1);
  m_ratio.resize(m_finestLevel + 1, 1);
  for (int lev = 0; lev <= m_finestLevel ; lev++)
   {
    m_grids[lev] = a_grids[lev];
    if (lev < m_finestLevel )
	   {
	    m_ratio[lev] = a_ratio[lev];
  	 }
   }

  m_dx.resize(m_finestLevel + 1);
  m_dx[0] = a_crseDx;
  for (int lev = 1; lev <= m_finestLevel;  lev++)
    {
      m_dx[lev] = m_dx[lev-1] / Real(m_ratio[lev-1]);
    }


  //copy any previous data pointers
  Vector<LevelData<FArrayBox>* > prevPhi;
  Vector<LevelData<FArrayBox>* > prevhcav;
  Vector<LevelData<FArrayBox>* > prevN;

  prevPhi.resize(m_hydraulicPotential.size());
  prevhcav.resize(m_cavityHeight.size());
  prevN.resize(m_effectivePressure.size());
  
  for (int lev = 0; lev < m_hydraulicPotential.size(); lev++)
    {
      prevPhi[lev]  = m_hydraulicPotential[lev];
      prevhcav[lev] = m_cavityHeight[lev];
      prevN[lev]    = m_effectivePressure[lev];
    }

  //initialize the data
  m_hydraulicPotential.resize(m_finestLevel + 1);
  m_cavityHeight.resize(m_finestLevel + 1);
  m_effectivePressure.resize(m_finestLevel + 1);
  
  for (int lev = 0; lev <= m_finestLevel; lev++)
    {
      m_hydraulicPotential[lev] = new LevelData<FArrayBox>(m_grids[lev],  SUBGLACIAL_N_COMP , SUBGLACIAL_N_GHOST * IntVect::Unit);
      m_cavityHeight[lev] = new LevelData<FArrayBox>(m_grids[lev],  SUBGLACIAL_N_COMP , SUBGLACIAL_N_GHOST * IntVect::Unit);
      m_effectivePressure[lev] = new LevelData<FArrayBox>(m_grids[lev],  SUBGLACIAL_N_COMP , SUBGLACIAL_N_GHOST * IntVect::Unit);
    }
  if (prevPhi.size() == 0)
    {
      // for now, set to zero, but need to set initial conditions from input data of some sort,
      // or read checkpoints
      for (int lev = 0; lev <= m_finestLevel; lev++)
       	{
          const LevelSigmaCS& levelCoords = *a_geometry[lev];
         
       	  for (DataIterator dit(m_grids[lev]);dit.ok();++dit)
       	    {
       	      FArrayBox& hydraulicPotential = (*m_hydraulicPotential[lev])[dit];
              FArrayBox& effectivePressure  = (*m_effectivePressure[lev])[dit];
	      
              hydraulicPotential.setVal(100.0); 
              effectivePressure.setVal(100.0);
              
              FArrayBox& cavityHeight = (*m_cavityHeight[lev])[dit];
	       cavityHeight.setVal(0.75);
	     }
	 }
    }

  if (prevPhi.size() > 0)
    {
      // if previous solution data exists, interpolate onto the new mesh
      for (int lev = 0; lev <= m_finestLevel; lev++)
        {
          //FIXME : only doing this to avoid stupid value outside the domain, should be a BC thing....
          for (DataIterator dit(m_grids[lev]);dit.ok();++dit)
           {
             FArrayBox& hydraulicPotential = (*m_hydraulicPotential[lev])[dit];
	     hydraulicPotential.setVal(100.0);
             FArrayBox& cavityHeight = (*m_cavityHeight[lev])[dit];
	     cavityHeight.setVal(1.0);
             FArrayBox& effectivePressure = (*m_effectivePressure[lev])[dit];
	     effectivePressure.setVal(100.0);
	   }
	   if (lev > 0)
	    {
	      FineInterp fi(m_grids[lev], m_hydraulicPotential[lev]->nComp(),
			    m_ratio[lev-1], m_grids[lev].physDomain());
	      fi.interpToFine(*m_hydraulicPotential[lev], *m_hydraulicPotential[lev-1]);
	      //need to fill the ghost cells on CF-interfaces now, because *m_damage[lev] has to be
	      //in good shape at all times if the interface with DamageConstitutiveRelation is to work
	      PiecewiseLinearFillPatch ghostFiller(m_grids[lev],m_grids[lev-1],m_hydraulicPotential[lev]->nComp(),
						   m_grids[lev-1].physDomain(), m_ratio[lev-1], m_hydraulicPotential[lev]->ghostVect()[0]);
	      ghostFiller.fillInterp(*m_hydraulicPotential[lev], *m_hydraulicPotential[lev-1],*m_hydraulicPotential[lev-1],1.0,0,0,m_hydraulicPotential[lev]->nComp());



              FineInterp fi2(m_grids[lev], m_cavityHeight[lev]->nComp(),
			    m_ratio[lev-1], m_grids[lev].physDomain());
	      fi2.interpToFine(*m_cavityHeight[lev], *m_cavityHeight[lev-1]);
	      //need to fill the ghost cells on CF-interfaces now, because *m_damage[lev] has to be
	      //in good shape at all times if the interface with DamageConstitutiveRelation is to work
	      PiecewiseLinearFillPatch ghostFiller2(m_grids[lev],m_grids[lev-1],m_cavityHeight[lev]->nComp(),
						   m_grids[lev-1].physDomain(), m_ratio[lev-1], m_cavityHeight[lev]->ghostVect()[0]);
	      ghostFiller2.fillInterp(*m_cavityHeight[lev], *m_cavityHeight[lev-1],*m_cavityHeight[lev-1],1.0,0,0,m_cavityHeight[lev]->nComp());
	      
             //Effective pressure
              FineInterp fi3(m_grids[lev], m_effectivePressure[lev]->nComp(),
			    m_ratio[lev-1], m_grids[lev].physDomain());
	      fi3.interpToFine(*m_effectivePressure[lev], *m_effectivePressure[lev-1]);
	      //need to fill the ghost cells on CF-interfaces now, because *m_damage[lev] has to be
	      //in good shape at all times if the interface with DamageConstitutiveRelation is to work
	      PiecewiseLinearFillPatch ghostFiller3(m_grids[lev],m_grids[lev-1],m_effectivePressure[lev]->nComp(), m_grids[lev-1].physDomain(), m_ratio[lev-1], m_effectivePressure[lev]->ghostVect()[0]);
	      ghostFiller3.fillInterp(*m_effectivePressure[lev], *m_effectivePressure[lev-1],*m_effectivePressure[lev-1],1.0,0,0,m_effectivePressure[lev]->nComp());
	     
             }
	    Interval ival(0,m_hydraulicPotential[lev]->nComp()-1);
	    Interval ival2(0,m_cavityHeight[lev]->nComp()-1);
	    Interval ival3(0,m_effectivePressure[lev]->nComp()-1);

	    if (prevPhi.size() > lev && prevPhi[lev])
             {
	        prevPhi[lev]->copyTo(ival,  *m_hydraulicPotential[lev], ival);
	        prevhcav[lev]->copyTo(ival2,  *m_cavityHeight[lev], ival2);
	        prevN[lev]->copyTo(ival3,  *m_effectivePressure[lev], ival3);
             }
	    m_hydraulicPotential[lev]->exchange();
	    m_cavityHeight[lev]->exchange();
	    m_effectivePressure[lev]->exchange();
	  }
      //free old data
      for (int lev =0; lev < prevPhi.size(); lev++)
        {
	  if (prevPhi[lev] != NULL)
	    {
	      delete prevPhi[lev]; prevPhi[lev] = NULL;
	    }
	  if (prevhcav[lev] != NULL)
	    {
	     delete prevhcav[lev]; prevhcav[lev] = NULL;
	    }
	  if (prevhcav[lev] != NULL)
	    {
	     delete prevN[lev]; prevN[lev] = NULL;
	    }
	}

     }
}


void AMRSubglacial::runLargeTimestep(Real a_dt,
       AmrIce& a_amrIce,
       const Vector<RefCountedPtr<LevelSigmaCS> >& a_geometry,
       const Vector<LevelData<FArrayBox>* >& a_cellVel,
       const Vector<LevelData<FluxBox>* >& a_faceVel)
{
  //TODO: Adaptively determine number of subiterations based on fastest relative
  //subglacial process
    // fill ghost data
    for (int lev=0; lev<= m_finestLevel; lev++)
    {
      if (lev > 0)
	    {

	     int nGhost = m_hydraulicPotential[lev]->ghostVect()[0];
	     PiecewiseLinearFillPatch pwl(m_grids[lev],  m_grids[lev-1], 1,
				       m_grids[lev].physDomain(),m_ratio[lev-1],nGhost);
	     // since we're not subcycling, don't need to interpolate in time
	     Real time_interp_coeff = 0.0;
	     pwl.fillInterp(*m_hydraulicPotential[lev], *m_hydraulicPotential[lev-1], *m_hydraulicPotential[lev-1],
			 time_interp_coeff,0, 0, 1);

             int nGhost_hcav = m_cavityHeight[lev]->ghostVect()[0];
	     PiecewiseLinearFillPatch pwl_hcav(m_grids[lev],  m_grids[lev-1], 1,
				       m_grids[lev].physDomain(),m_ratio[lev-1],nGhost_hcav);
	     // since we're not subcycling, don't need to interpolate in time
	     Real time_interp_coeff_hcav = 0.0;
	     pwl_hcav.fillInterp(*m_cavityHeight[lev], *m_cavityHeight[lev-1], *m_cavityHeight[lev-1],
			 time_interp_coeff_hcav,0, 0, 1);
	     
             int nGhost_N = m_effectivePressure[lev]->ghostVect()[0];
	     PiecewiseLinearFillPatch pwl3(m_grids[lev],  m_grids[lev-1], 1,
				       m_grids[lev].physDomain(),m_ratio[lev-1],nGhost_N);
	     // since we're not subcycling, don't need to interpolate in time
	     pwl3.fillInterp(*m_effectivePressure[lev], *m_effectivePressure[lev-1], *m_effectivePressure[lev-1],
			 time_interp_coeff,0, 0, 1);


	    }
      m_hydraulicPotential[lev]->exchange();
      m_cavityHeight[lev]->exchange();
      m_effectivePressure[lev]->exchange();
    }

    //Define runoff profile
    Vector<LevelData<FArrayBox>* >  runoff;
    runoff.resize(m_finestLevel + 1);
    for (int lev=0; lev<= m_finestLevel; lev++)
    {
      runoff[lev] = new LevelData<FArrayBox>(m_grids[lev],1,IntVect::Unit);
    }

    // Define phiold and hcavold
    Vector<LevelData<FArrayBox>* > phiOld;
    Vector<LevelData<FArrayBox>* > hcavOld;

    phiOld.resize(m_finestLevel+1);
    hcavOld.resize(m_finestLevel+1);
    for (int lev = 0; lev <= m_finestLevel ; lev++)
    {
      phiOld[lev] = new LevelData<FArrayBox>(m_grids[lev], 1, IntVect::Unit);
      hcavOld[lev] = new LevelData<FArrayBox>(m_grids[lev], 1, IntVect::Unit);

    }

    Vector<RefCountedPtr<LevelData<FArrayBox>>>  aCoef_phi(m_finestLevel+1);
    for (int lev=0; lev<= m_finestLevel; lev++)
    {
      aCoef_phi[lev] = RefCountedPtr<LevelData<FArrayBox> > (new LevelData<FArrayBox>(m_grids[lev],1,IntVect::Unit));
    }

    Vector<RefCountedPtr<LevelData<FluxBox>>>  bCoef_phi(m_finestLevel+1);
    for (int lev=0; lev<= m_finestLevel; lev++)
    {
      bCoef_phi[lev] = RefCountedPtr<LevelData<FluxBox> > (new LevelData<FluxBox>(m_grids[lev], 1, IntVect::Zero));
    }

    Vector<LevelData<FArrayBox>*>  rhs_phi;
    rhs_phi.resize(m_finestLevel + 1);
    for (int lev=0; lev<= m_finestLevel; lev++)
    {
      rhs_phi[lev] = new LevelData<FArrayBox>(m_grids[lev],1,IntVect::Unit);
    }

    //For now prescribe the parameters here (Eventually these will come from inputs)
    SubglacialParameters subParams;

    subParams.coarsestDomain = m_grids[0].physDomain();
    subParams.refRatio = a_amrIce.refRatios();
    subParams.coarsestDx = m_dx[0][0];
    subParams.alpha = 1.0;
    subParams.beta = 1.0;
    subParams.verbosity = 6;
    subParams.coefficient_average_type = 0;


    Vector<ProblemDomain>  vectDomains(m_finestLevel+1);
    for (int lev=0; lev<= m_finestLevel; lev++)
    {
      vectDomains[lev] = m_grids[lev].physDomain();
    }

    Vector<RealVect> vectDx(m_finestLevel+1);
    Vector<Real> pre_vectDx = a_amrIce.amrDx();
    for (int lev=0; lev<= m_finestLevel; lev++)
    {
      vectDx[lev] = RealVect(pre_vectDx[lev]);
    }



    Vector<RefCountedPtr<LevelData<FArrayBox>>>  aCoef_hcav(m_finestLevel+1);
    for (int lev=0; lev<= m_finestLevel; lev++)
    {
      aCoef_hcav[lev] = RefCountedPtr<LevelData<FArrayBox> > (new LevelData<FArrayBox>(m_grids[lev],1,IntVect::Unit));
    }

    Vector<RefCountedPtr<LevelData<FluxBox>>>  bCoef_hcav(m_finestLevel+1);
    for (int lev=0; lev<= m_finestLevel; lev++)
    {
      bCoef_hcav[lev] = RefCountedPtr<LevelData<FluxBox> > (new LevelData<FluxBox>(m_grids[lev], 1, IntVect::Zero));
    }

    Vector<LevelData<FArrayBox>*>  rhs_hcav;
    rhs_hcav.resize(m_finestLevel + 1);
    for (int lev=0; lev<= m_finestLevel; lev++)
    {
      rhs_hcav[lev] = new LevelData<FArrayBox>(m_grids[lev],1,IntVect::Unit);
    }


    //TODO Michael: choose the subglacial timestep adaptively
   // int tmax_multiple= 30;
    // int n_subIts     = 200; 
    int n_subIts = 10; // edited by TM for initial testing
    Real elapsedTime = 0.0;
   // Real small_dt    = tmax_multiple*a_dt/n_subIts;
    

    Real small_dt    = min(a_dt/n_subIts, 0.005);
    
    // n_subIts = floor(a_dt/small_dt); commented out by TM to control iterations

    Real flux_out = 0.0;
    Real total_runoff = 0.0;

    //TODO: Make this an if "save_hydrograph" option
    //FILE *fp;
    //fp = fopen("hydrograph_data.csv", "a");
    //fprintf(fp, "time, runoff, flux out \n"); // commented out by TM

    for (int subIt =0; subIt < n_subIts; subIt++)
    {
    subglacialAdvance(small_dt,
           aCoef_phi,
           bCoef_phi,
           rhs_phi,
           aCoef_hcav,
           bCoef_hcav,
           rhs_hcav,
           phiOld,
           hcavOld,
           vectDomains,
           vectDx,
           subParams,
           runoff,
    	   a_geometry,
    	   a_cellVel,
    	   a_faceVel);

      m_time += small_dt;
      elapsedTime += small_dt;
      setEffectivePressure(a_geometry, elapsedTime, small_dt);

      //calculate total runoff and flux out of boundary
      Vector<LevelData<FArrayBox>* >  source;
      source.resize(m_finestLevel + 1);
      for (int lev=0; lev<= m_finestLevel; lev++)
       {
          source[lev] = new LevelData<FArrayBox>(m_grids[lev],1,IntVect::Unit);
       }  
      computeSinkTerm(source, m_hydraulicPotential, a_amrIce, a_geometry, a_dt);

      flux_out = compute_fluxOutRate(m_hydraulicPotential, a_amrIce, a_geometry, runoff,  small_dt );

      total_runoff = computeSum(runoff,
                               a_amrIce.refRatios(),
                               pre_vectDx[0],
                               Interval(0,0),
                               0);


      pout() <<  "Iteration number " << subIt << "out of" << n_subIts <<std::endl;
      pout() <<  "For this iteration the flux out is " << flux_out<<std::endl;
      pout() <<  "For this iteration the total runoff is " << total_runoff<<std::endl;
      
   //   fprintf(fp, "%f, %e, %e \n", elapsedTime, total_runoff, flux_out);

    }//end subiterations
    
    //fprintf(fp, "%f, %e, %e \n", elapsedTime, total_runoff, flux_out);
    //fclose(fp);


  for (int lev=0; lev<= m_finestLevel; lev++)
  {
    if (runoff[lev] != NULL)
          {
            delete runoff[lev]; runoff[lev] = NULL;
          }

    if (phiOld[lev] != NULL)
          {
            delete phiOld[lev]; phiOld[lev] = NULL;
          }

    if (hcavOld[lev] != NULL)
          {
            delete hcavOld[lev]; hcavOld[lev] = NULL;
          }

    if (rhs_phi[lev] != NULL)
          {
            delete rhs_phi[lev]; rhs_phi[lev] = NULL;
          }
    if (rhs_hcav[lev] != NULL)
          {
            delete rhs_hcav[lev]; rhs_hcav[lev] = NULL;
          }
  }

  m_time_step++;


}

void AMRSubglacial::subglacialAdvance(Real a_dt,
       Vector<RefCountedPtr<LevelData<FArrayBox> > > a_aCoef_phi,
       Vector<RefCountedPtr<LevelData<FluxBox> > > a_bCoef_phi,
       Vector<LevelData<FArrayBox> * > a_rhs_phi,
       Vector<RefCountedPtr<LevelData<FArrayBox> > > a_aCoef_hcav,
       Vector<RefCountedPtr<LevelData<FluxBox> > > a_bCoef_hcav,
       Vector<LevelData<FArrayBox>* > a_rhs_hcav,
       Vector<LevelData<FArrayBox>* > a_phiOld,
       Vector<LevelData<FArrayBox>* > a_hcavOld,
       Vector<ProblemDomain> a_vectDomains,
       Vector<RealVect> a_vectDx,
       SubglacialParameters a_subParams,
       Vector<LevelData<FArrayBox>* > a_runoff,
       const Vector<RefCountedPtr<LevelSigmaCS> >& a_geometry,
       const Vector<LevelData<FArrayBox>* >& a_cellVel,
       const Vector<LevelData<FluxBox>* >& a_faceVel)
{

    setRunoff(a_runoff, a_geometry, m_time);

    for (int lev = 0; lev <= m_finestLevel ; lev++)
    {
      //TODO Michael: Check that this is doing a deep copy
      for(DataIterator dit(m_grids[lev]); dit.ok(); ++dit )
        {
          (*a_phiOld[lev])[dit].copy((*m_hydraulicPotential[lev])[dit]);
          (*a_hcavOld[lev])[dit].copy((*m_cavityHeight[lev])[dit]);
        }
      //TODO Michael: figure out what exchange does
      a_phiOld[lev]->exchange();
      a_hcavOld[lev]->exchange();
    }

    compute_aCoefPhi(a_aCoef_phi, m_hydraulicPotential, m_cavityHeight, a_geometry, a_cellVel, a_runoff, a_dt);
    compute_bCoefPhi(a_bCoef_phi, m_hydraulicPotential, m_cavityHeight, a_geometry, a_cellVel, a_runoff, a_dt);
    compute_RhsPhi(a_rhs_phi, m_hydraulicPotential, m_cavityHeight, a_phiOld, a_geometry, a_cellVel, a_runoff, a_dt);

    // need to get only the grids relevent to the problem

    Vector<DisjointBoxLayout> reduced_grids(m_finestLevel+1);
    for (int lev = 0; lev <= m_finestLevel ; lev++)
      {
       reduced_grids[lev] = m_grids[lev];
      }
    RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > > opFactory_phi
      = RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > >
        (defineOperatorFactory_phi(reduced_grids, a_aCoef_phi, a_bCoef_phi, a_subParams));

    int lBase = 0;
    MultilevelLinearOp<FArrayBox> mlOp_hcav;
    MultilevelLinearOp<FArrayBox> mlOp_phi;
    int numMGIter = 1;
    //pp.query("numMGIterations", numMGIter);

    mlOp_phi.m_num_mg_iterations = numMGIter;
    mlOp_hcav.m_num_mg_iterations = numMGIter;
    int numMGSmooth = 4;
    mlOp_phi.m_num_mg_smooth = numMGSmooth;
    mlOp_hcav.m_num_mg_smooth = numMGSmooth;
    int preCondSolverDepth = -1;
    mlOp_phi.m_preCondSolverDepth = preCondSolverDepth;
    mlOp_hcav.m_preCondSolverDepth = preCondSolverDepth;


    //Real tolerance = 1.0e-7;
    Real tolerance = 1.0e-1; // Updated by TM for initial hydro model testing
    //pp.query("tolerance", tolerance);

    int max_iter = 30;
    //pp.query("max_iterations", max_iter);

    mlOp_phi.define(reduced_grids, a_subParams.refRatio, a_vectDomains,
              a_vectDx, opFactory_phi, lBase);


    RelaxSolver<Vector<LevelData<FArrayBox>* > > solver_phi;

    bool homogeneousBC = false;
    solver_phi.define(&mlOp_phi, homogeneousBC);
    solver_phi.m_verbosity = a_subParams.verbosity;
    solver_phi.m_normType = 2;
    solver_phi.m_eps = tolerance;
    solver_phi.m_imax = max_iter;


    solver_phi.solve(m_hydraulicPotential, a_rhs_phi);

    //Now solve hcav Problem
    compute_aCoefHcav(a_aCoef_hcav, m_hydraulicPotential, m_cavityHeight, a_geometry, a_cellVel, a_runoff, a_dt);
    compute_bCoefHcav(a_bCoef_hcav, m_hydraulicPotential, m_cavityHeight, a_geometry, a_cellVel, a_runoff, a_dt);
    compute_RhsHcav(a_rhs_hcav, m_hydraulicPotential, m_cavityHeight, a_hcavOld, a_geometry,a_cellVel, a_runoff, a_dt);


    RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > > opFactory_hcav
        = RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > >
          (defineOperatorFactory_hcav(reduced_grids, a_aCoef_hcav, a_bCoef_hcav, a_subParams));

    mlOp_hcav.define(reduced_grids, a_subParams.refRatio, a_vectDomains,
                  a_vectDx, opFactory_hcav, lBase);


    RelaxSolver<Vector<LevelData<FArrayBox>* > > solver_hcav;


    solver_hcav.define(&mlOp_hcav, homogeneousBC);
    solver_hcav.m_verbosity = a_subParams.verbosity;
    solver_hcav.m_normType = 2;
    solver_hcav.m_eps = tolerance;
    solver_hcav.m_imax = max_iter;

    solver_hcav.solve(m_cavityHeight, a_rhs_hcav);

}


///Prescribe runoff from surface.
/**
  Eventually the runoff should be prescribed using an external model of some sort (based on routing etc.)
 */
void AMRSubglacial::setRunoff(Vector<LevelData<FArrayBox>* >& a_runoff,
			      const Vector<RefCountedPtr<LevelSigmaCS> >& a_geometry,
			      Real a_t)
{

  for (int lev = 0 ; lev <= m_finestLevel; lev++)
    {
      const LevelSigmaCS& levelCoords = *a_geometry[lev];
      for (DataIterator dit(m_grids[lev]);dit.ok();++dit)
      	{
	       FArrayBox& netRunoff = (*a_runoff[lev])[dit];
	       const FArrayBox& thick = levelCoords.getH()[dit];
       	       const FArrayBox& surf = levelCoords.getSurfaceHeight()[dit];
               const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];

	       for (BoxIterator bit(netRunoff.box());bit.ok();++bit)
	        {
	         const IntVect& iv = bit();
	         Real runoff;
	         if( mask(iv) == GROUNDEDMASKVAL)
                  {
                     runoff = 1.0;
                  }
                 else
                  {
                      runoff = 0.0;
                  }
                 //TODO: load a runoff model?
                 netRunoff(iv) = runoff;
	        }
      	}
    }
}


void AMRSubglacial::setEffectivePressure(
         const Vector<RefCountedPtr<LevelSigmaCS> >& a_geometry,
         const Real a_elapsedTime,
         const Real a_small_dt)
{
    for (int lev = 0; lev<=m_finestLevel; lev++)
      {
        const LevelSigmaCS& levelCoords = *a_geometry[lev];
        const DisjointBoxLayout& grids = levelCoords.grids();

        DataIterator dit = (*m_effectivePressure[lev]).dataIterator();
        for (dit.begin(); dit.ok(); ++dit)
          {
            FArrayBox& thisN = (*m_effectivePressure[lev])[dit];
            //TODO MM: same issue as phi initialisation; how do you avoid the 'not in box' issue?

            FArrayBox& thisPHI = (*m_hydraulicPotential[lev])[dit];
            
            const FArrayBox& thisZb = levelCoords.getTopography()[dit];
            const FArrayBox& thisZs = levelCoords.getSurfaceHeight()[dit];
                
            Box thisBox = thisZs.box();
            
            Real iceDensity = levelCoords.iceDensity();
            Real waterDensity = levelCoords.waterDensity();
            
            const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
            Real tol = 0.001;
      
            BoxIterator bit(thisBox);
            for (bit.begin(); bit.ok(); ++bit)
              {
                const IntVect& iv = bit();
                if (mask(iv) == GROUNDEDMASKVAL)
                  {
                    //Real oldN    = thisN(iv);
                    thisN(iv)    = (iceDensity*m_g*thisZs(iv)+ m_g*(waterDensity - iceDensity)*thisZb(iv) - thisPHI(iv));
                    //thisN(iv)    = oldN*(a_elapsedTime - a_small_dt) + a_small_dt*(iceDensity*m_g*thisZs(iv)+ m_g*(waterDensity - iceDensity)*thisZb(iv) - thisPHI(iv));
                    //thisN(iv)    = thisN(iv)/a_elapsedTime;
                  }
                else
                  {
                    thisN(iv) = 0.0;
                  }
              }
          }
      }
}

void AMRSubglacial::setObserverEffectivePressure(
                        Vector<LevelData<FArrayBox>*>  a_observerEffectivePressure)
{//TODO: A bit unwieldy for now, can be shortened. Just needs to set the argument (the observer
//Effective pressure) to the calculated effective pressure.
    if (m_effectivePressureDefined)
      {
        for (int lev = 0; lev<=m_finestLevel; lev++)
          {
  
            DataIterator dit = (*m_effectivePressure[lev]).dataIterator();
            for (dit.begin(); dit.ok(); ++dit)
              {
                FArrayBox& thisN = (*m_effectivePressure[lev])[dit];
                //TODO MM: same issue as phi initialisation; how do you avoid the 'not in box' issue?
  
                FArrayBox& thisObserverN = (*a_observerEffectivePressure[lev])[dit];
                
                    
                Box thisBox = thisObserverN.box();
                
                Real tol = 0.001;
          
                BoxIterator bit(thisBox);
                for (bit.begin(); bit.ok(); ++bit)
                  {
                    const IntVect& iv = bit();
                    thisObserverN(iv) =  thisN(iv);
                  }
              }
          }

     } 
   else
     {
    
     }
}



/*----------------------------------------------------------------*/
//          Set coefficients for phi
/*----------------------------------------------------------------*/

void AMRSubglacial::compute_aCoefPhi(Vector<RefCountedPtr<LevelData<FArrayBox> > > a_aCoef,
         Vector<LevelData<FArrayBox>* >&   a_phi,
         Vector<LevelData<FArrayBox>* >&   a_hcav,
         const Vector<RefCountedPtr<LevelSigmaCS> >& a_geometry,
         const Vector<LevelData<FArrayBox>* >&   a_cellVel,
         Vector<LevelData<FArrayBox>* >&   a_runoff,
         Real a_dt
        )
{
  for (int lev = 0 ; lev <= m_finestLevel; lev++){
    const  LevelSigmaCS& levelCoords = *a_geometry[lev];
    const  DisjointBoxLayout& grids = levelCoords.grids();

    DataIterator dit = (*a_aCoef[lev]).dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& thisaCoef = (*a_aCoef[lev])[dit];
      Box thisBox = thisaCoef.box();

      FArrayBox& thisPHI = (*a_phi[lev])[dit];
      FArrayBox& thisHCAV = (*a_hcav[lev])[dit];

      const FArrayBox& thisZb = levelCoords.getTopography()[dit];
      const FArrayBox& thisZs = levelCoords.getSurfaceHeight()[dit];
      const FArrayBox& thisUb = (*a_cellVel[lev])[dit];
      FArrayBox& thisRunoff = (*a_runoff[lev])[dit];
      
      Real iceDensity = levelCoords.iceDensity();
      Real waterDensity = levelCoords.waterDensity();

      const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
      Real tol = 0.001;

      BoxIterator bit(thisBox);
      for (bit.begin(); bit.ok(); ++bit)
        {
          const IntVect& iv = bit();
          if (mask(iv) == GROUNDEDMASKVAL)
            {
              Real Nval    = iceDensity*m_g*thisZs(iv)+ m_g*(waterDensity - iceDensity)*thisZb(iv) - thisPHI(iv);
              Real sPrime  = m_sigma/waterDensity +m_h0*exp(-max(Nval, 0.0)/m_N0)/m_N0;
              Real g_2     = 2* m_Ahat* thisHCAV(iv)*pow(abs(max(Nval, 0.0)), m_n - 1)/(pow(m_n, m_n));
              thisaCoef(iv) = sPrime + a_dt*g_2;
            }
          else
            {
              thisaCoef(iv) = -m_BigNotGrounded;

            }
        }
      //impose effective boundary condition at ice margin if thickness<0 or thickness of neighbour<0
      for(int dir = 0; dir<SpaceDim; dir++)
        {
          Box facebox(grids[dit]);
          facebox.surroundingNodes(dir);
          for (BoxIterator bit(facebox); bit.ok(); ++bit)
            {
              const IntVect& iv = bit();
              IntVect ivm = iv - BASISV(dir);
              IntVect ivp = iv + BASISV(dir);
                 if (mask(iv) == GROUNDEDMASKVAL && (mask(ivp) !=GROUNDEDMASKVAL || mask(ivm) != GROUNDEDMASKVAL))
                {
                  if (m_g*waterDensity*thisZb(iv) -m_g* waterDensity*min(thisZb(iv), 0.0) - thisPHI(iv) < 0)
                    {
                      thisaCoef(iv) += m_BigaPhiCoef;
                    }
                }
            }

        }

    } // end loop over grids


  } // end loop over levels
}


void AMRSubglacial::compute_bCoefPhi(Vector<RefCountedPtr<LevelData<FluxBox>>> a_bCoef,
         Vector<LevelData<FArrayBox>* >&   a_phi,
         Vector<LevelData<FArrayBox>* >&   a_hcav,
 	 const Vector<RefCountedPtr<LevelSigmaCS> >& a_geometry,
         const Vector<LevelData<FArrayBox>* >&   a_cellVel,
         Vector<LevelData<FArrayBox>* >&   a_runoff,
         Real a_dt)
{
  Vector<RefCountedPtr<LevelData<FluxBox>>> edgePHI(m_finestLevel+1);
  Vector<RefCountedPtr<LevelData<FluxBox>>> edgeHCAV(m_finestLevel+1);
  Vector<RefCountedPtr<LevelData<FluxBox>>> edgeZb(m_finestLevel+1);
  Vector<RefCountedPtr<LevelData<FluxBox>>> edgeZs(m_finestLevel+1);
  for(int lev = 0 ; lev <= m_finestLevel; lev++)
  {
      DataIterator dit = (*a_bCoef[lev]).dataIterator();
      const LevelSigmaCS& levelCoords = *a_geometry[lev];
      const DisjointBoxLayout& grids = levelCoords.grids();
      
      Real iceDensity = levelCoords.iceDensity();
      Real waterDensity = levelCoords.waterDensity();
      
      for (dit.begin(); dit.ok(); ++dit)
        {
        FluxBox& thisBCoef = (*a_bCoef[lev])[dit];
        FArrayBox& thisPHI = (*a_phi[lev])[dit];
        FArrayBox& thisHCAV = (*a_hcav[lev])[dit];

        const FArrayBox& thisZb = levelCoords.getTopography()[dit];
        const FArrayBox& thisZs = levelCoords.getSurfaceHeight()[dit];


        //TODO: Check this doesn't cause a memory leak issue?
        edgePHI[lev] =  RefCountedPtr<LevelData<FluxBox> > (new LevelData<FluxBox>(m_grids[lev], 1, IntVect::Zero));
        edgeHCAV[lev]=  RefCountedPtr<LevelData<FluxBox> > (new LevelData<FluxBox>(m_grids[lev], 1, IntVect::Zero));
        edgeZb[lev] =  RefCountedPtr<LevelData<FluxBox> > (new LevelData<FluxBox>(m_grids[lev], 1, IntVect::Zero));
        edgeZs[lev]=  RefCountedPtr<LevelData<FluxBox> > (new LevelData<FluxBox>(m_grids[lev], 1, IntVect::Zero));


        CellToEdge(thisPHI, (*edgePHI[lev])[dit]);
        CellToEdge(thisHCAV, (*edgeHCAV[lev])[dit]);

        CellToEdge(thisZb, (*edgeZb[lev])[dit]);
        CellToEdge(thisZs, (*edgeZs[lev])[dit]);

        const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
        Real tol = 0.001;

        for (int dir=0; dir<SpaceDim; dir++)
          {
            FArrayBox& dirFlux = thisBCoef[dir];
            const Box& dirBox = dirFlux.box();

            FArrayBox& dirPHI = (*edgePHI[lev])[dit][dir];
            FArrayBox& dirHCAV = (*edgeHCAV[lev])[dit][dir];

            FArrayBox& dirZb = (*edgeZb[lev])[dit][dir];
            FArrayBox& dirZs = (*edgeZs[lev])[dit][dir];


          // this sets up a vector which is 0 in the dir
          // direct and 0.5 in the other (cell-centered) directions
          BoxIterator faceBit(dirBox);
          for(faceBit.begin(); faceBit.ok();  ++faceBit)
            {
              IntVect ivFace = faceBit();
              if (mask(ivFace) == GROUNDEDMASKVAL)
               {
                 // bcoef = dt*D(hcav, phi) = dt*K*h^3
                 Real Nval    =  (iceDensity*m_g*dirZs(ivFace, 0)+ m_g*(waterDensity - iceDensity)*dirZb(ivFace, 0) - dirPHI(ivFace, 0));
                 Real helVal  =  m_h0*exp(-max(Nval, 0.0)/m_N0);
                 dirFlux(ivFace, 0) =  a_dt*m_K*pow((helVal + dirHCAV(ivFace, 0.0)), m_alpha)/m_g/waterDensity;
                }
              else
                {
                  dirFlux(ivFace, 0) = 0.0;
                }
            } //end loop over faces in this direction

          Box facebox(grids[dit]);
          facebox.surroundingNodes(dir);
          for (BoxIterator bit(facebox); bit.ok(); ++bit)
            {
              const IntVect& iv = bit();
              IntVect ivm = iv - BASISV(dir);
              IntVect ivp = iv + BASISV(dir);
              //To cover all directions, check !grounded/grounded boudnaries above and below.
              //I think the `above' direction is covered, but include for now for completeness
              if (dirBox.contains(ivm))
               {
                  if (mask(iv) != GROUNDEDMASKVAL && mask(ivm) == GROUNDEDMASKVAL)
                   {
                      dirFlux(iv) = 0.0;
                   }
               }
              if (dirBox.contains(ivp))
               {
                  if (mask(iv) != GROUNDEDMASKVAL && mask(ivp) == GROUNDEDMASKVAL)
                   {
                      dirFlux(iv) = 0.0;
                   }
               }

            }

        } // end loop over directions
    } // end loop over grids
  } // end loop over levels
}



/********/
void AMRSubglacial::compute_RhsPhi(Vector<LevelData<FArrayBox>*>&      a_rhs,
          Vector<LevelData<FArrayBox>* >&   a_phi,
          Vector<LevelData<FArrayBox>* >&   a_hcav,
          Vector<LevelData<FArrayBox>* >&   a_phi_old,
          const Vector<RefCountedPtr<LevelSigmaCS> >& a_geometry,
          const Vector<LevelData<FArrayBox>* >&   a_cellVel,
          Vector<LevelData<FArrayBox>* >&   a_runoff,
          Real a_dt
        )
 {

  for (int lev = 0 ; lev<= m_finestLevel; lev++)
  {
    const LevelSigmaCS& levelCoords = *a_geometry[lev];
    const DisjointBoxLayout& grids = levelCoords.grids();
    for (DataIterator dit = (*a_rhs[lev]).dataIterator(); dit.ok(); ++dit)
      {
        FArrayBox& thisRHS = (*a_rhs[lev])[dit()];
        Box thisBox = thisRHS.box();

        FArrayBox& thisPHI = (*a_phi[lev])[dit()];
        FArrayBox& thisHCAV = (*a_hcav[lev])[dit()];
        FArrayBox& thisPHIOLD = (*a_phi_old[lev])[dit()];

        const FArrayBox& thisZb = levelCoords.getTopography()[dit];
        const FArrayBox& thisZs = levelCoords.getSurfaceHeight()[dit];
        const FArrayBox& thisUb = (*a_cellVel[lev])[dit];
        const FArrayBox& thisRunoff = (*a_runoff[lev])[dit];

        Real iceDensity = levelCoords.iceDensity();
        Real waterDensity = levelCoords.waterDensity();

        const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
        Real tol = 0.001;

        BoxIterator bit(thisBox);
        for (bit.begin(); bit.ok(); ++bit)
           {
             const IntVect& iv = bit();
             if (mask(iv) == GROUNDEDMASKVAL)
               {
                 //Now express as a function of iv
                 Real Nval    = iceDensity*m_g*thisZs(iv)+ m_g*(waterDensity - iceDensity)*thisZb(iv) - thisPHI(iv);
                 Real sPrime  = m_sigma/waterDensity +  m_h0*exp(-max(Nval, 0.0)/m_N0)/m_N0;
                 Real velMag = sqrt(pow(thisUb(iv, 0), 2) + pow(thisUb(iv, 1), 2));
                 Real g_1     = waterDensity*m_basalMelt/iceDensity + abs(velMag)*max((m_hr - thisHCAV(iv))/m_lr, 0.0)  - 2*m_Ahat*thisHCAV(iv)*pow(abs(max(Nval, 0.0)), m_n -1)*(iceDensity*m_g*thisZs(iv)+ (waterDensity*m_g - iceDensity*m_g)*thisZb(iv))/pow(m_n,(m_n));
                 thisRHS(iv) = sPrime*thisPHIOLD(iv) + a_dt*(m_basalMelt+ thisRunoff(iv)- g_1);
                 }
             else
               {
               //Set the value outside the domain of interest to conform to the boundary condition
                 if (thisZb(iv) <=0)
                   {thisRHS(iv) = 0.0;}
                 else
                   {thisRHS(iv) = -m_BigNotGrounded*m_g*iceDensity*thisZb(iv);}
               }
           }
         //Impose effective boundary condition
         for(int dir = 0; dir<SpaceDim; dir++)
           {
             Box facebox(grids[dit]);
             facebox.surroundingNodes(dir);
             for (BoxIterator bit(facebox); bit.ok(); ++bit)
               {
                 const IntVect& iv = bit();
                 IntVect ivm = iv - BASISV(dir);
                 IntVect ivp = iv + BASISV(dir);
                 if (mask(iv) == GROUNDEDMASKVAL && (mask(ivp) !=GROUNDEDMASKVAL || mask(ivm) != GROUNDEDMASKVAL))
                   {
                     //RHS = A = B (large) * phi_boundary_val
                     if (m_g*waterDensity*thisZb(iv) -m_g* waterDensity*min(thisZb(iv), 0.0) - thisPHI(iv) < 0)
                       {
                         thisRHS(iv) += m_BigaPhiCoef*(m_g*waterDensity*thisZb(iv) -m_g* waterDensity*min(thisZb(iv), 0.0));
                       }
                   }
//                    thisRHS(iv) += -m_BigaPhiCoef*(waterDensity*thisZb(iv) - waterDensity*min(thisZb(iv), 0.0));
               }

           }

    }
  }//Loop over levels
}

/*----------------------------------------------------------------*/
//            Set coefficients for h_cav
/*----------------------------------------------------------------*/


void AMRSubglacial::compute_aCoefHcav(Vector<RefCountedPtr<LevelData<FArrayBox>> > a_aCoef,
        Vector<LevelData<FArrayBox>* >&   a_phi,
        Vector<LevelData<FArrayBox>* >&   a_hcav,
        const Vector<RefCountedPtr<LevelSigmaCS> >& a_geometry,
        const Vector<LevelData<FArrayBox>* >&   a_cellVel,
        Vector<LevelData<FArrayBox>* >&   a_runoff,
        Real a_dt
        )
{
  for (int lev = 0 ; lev <= m_finestLevel; lev++){
    const LevelSigmaCS& levelCoords = *a_geometry[lev];
    const DisjointBoxLayout& grids = levelCoords.grids();

    Real iceDensity = levelCoords.iceDensity();
    Real waterDensity = levelCoords.waterDensity();
    
    DataIterator dit = (*a_aCoef[lev]).dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& thisaCoef = (*a_aCoef[lev])[dit];
      Box thisBox = thisaCoef.box();

      FArrayBox& thisPHI = (*a_phi[lev])[dit];
      FArrayBox& thisHCAV = (*a_hcav[lev])[dit];

      const FArrayBox& thisZb = levelCoords.getTopography()[dit];
      const FArrayBox& thisZs = levelCoords.getSurfaceHeight()[dit];
      const FArrayBox& thisUb = (*a_cellVel[lev])[dit];
      FArrayBox& thisRunoff = (*a_runoff[lev])[dit];

      const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];

      BoxIterator bit(thisBox);
      for (bit.begin(); bit.ok(); ++bit)
        {
          const IntVect& iv = bit();
          if (mask(iv) == GROUNDEDMASKVAL)
            {
              Real Nval    = iceDensity*m_g*thisZs(iv)+ m_g*(waterDensity - iceDensity)*thisZb(iv) - thisPHI(iv);
              Real g_4     =  - 2*m_Ahat*pow(abs(Nval), 2)*max(Nval, 0.0)/pow(m_n,m_n);
              thisaCoef(iv) = 1.0 - a_dt*g_4;
            }
          else
            {
              // Outside grounded ice set acoef large. This sets hcav to 0 (assuming rhs is set to 0)
              thisaCoef(iv) = -m_BigNotGrounded;
            }
         }


    } // end loop over grids
  } // end loop over levels
}



void AMRSubglacial::compute_bCoefHcav(Vector<RefCountedPtr<LevelData<FluxBox>> > a_bCoef,
         Vector<LevelData<FArrayBox>* >&   a_phi,
         Vector<LevelData<FArrayBox>* >&   a_hcav,
         const Vector<RefCountedPtr<LevelSigmaCS> >& a_geometry,
         const Vector<LevelData<FArrayBox>* >&   a_cellVel,
         Vector<LevelData<FArrayBox>* >&   a_runoff,
         Real a_dt)
{
  Vector<RefCountedPtr<LevelData<FluxBox>>> edgePHI(m_finestLevel+1);
  Vector<RefCountedPtr<LevelData<FluxBox>>> edgeHCAV(m_finestLevel+1);
  Vector<RefCountedPtr<LevelData<FluxBox>>> edgeZb(m_finestLevel+1);
  Vector<RefCountedPtr<LevelData<FluxBox>>> edgeZs(m_finestLevel+1);
  for(int lev = 0 ; lev <= m_finestLevel; lev++)
  {
      DataIterator dit = (*a_bCoef[lev]).dataIterator();
      const LevelSigmaCS& levelCoords = *a_geometry[lev];
      const DisjointBoxLayout& grids = levelCoords.grids();
      Real iceDensity = levelCoords.iceDensity();
      Real waterDensity = levelCoords.waterDensity();
      
      for (dit.begin(); dit.ok(); ++dit)
      {
        FluxBox& thisBCoef = (*a_bCoef[lev])[dit];
        FArrayBox& thisPHI = (*a_phi[lev])[dit];
        FArrayBox& thisHCAV = (*a_hcav[lev])[dit];

        const FArrayBox& thisZb = levelCoords.getTopography()[dit];
        const FArrayBox& thisZs = levelCoords.getSurfaceHeight()[dit];

        //TODO: Check, does this cause a memory leak issue?
        edgePHI[lev] =  RefCountedPtr<LevelData<FluxBox> > (new LevelData<FluxBox>(m_grids[lev], 1, IntVect::Zero));
        edgeHCAV[lev]=  RefCountedPtr<LevelData<FluxBox> > (new LevelData<FluxBox>(m_grids[lev], 1, IntVect::Zero));
        edgeZb[lev] =  RefCountedPtr<LevelData<FluxBox> > (new LevelData<FluxBox>(m_grids[lev], 1, IntVect::Zero));
        edgeZs[lev]=  RefCountedPtr<LevelData<FluxBox> > (new LevelData<FluxBox>(m_grids[lev], 1, IntVect::Zero));


        CellToEdge(thisPHI, (*edgePHI[lev])[dit]);
        CellToEdge(thisHCAV, (*edgeHCAV[lev])[dit]);

        CellToEdge(thisZb, (*edgeZb[lev])[dit]);
        CellToEdge(thisZs, (*edgeZs[lev])[dit]);

        const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
        Real tol = 0.001;

        for (int dir=0; dir<SpaceDim; dir++)
          {
            FArrayBox& dirFlux = thisBCoef[dir];
            const Box& dirBox = dirFlux.box();

            FArrayBox& dirPHI = (*edgePHI[lev])[dit][dir];
            FArrayBox& dirHCAV = (*edgeHCAV[lev])[dit][dir];

            FArrayBox& dirZb = (*edgeZb[lev])[dit][dir];
            FArrayBox& dirZs = (*edgeZs[lev])[dit][dir];


          // this sets up a vector which is 0 in the dir
          // direct and 0.5 in the other (cell-centered) directions
            BoxIterator faceBit(dirBox);
            for(faceBit.begin(); faceBit.ok();  ++faceBit)
            {
              IntVect ivFace = faceBit();
              if (mask(ivFace) == GROUNDEDMASKVAL)
                {
                  dirFlux(ivFace, 0) = a_dt*(m_epsilon_d);
                }
              else
                {
                  dirFlux(ivFace, 0) = 0.0;
                }
            } //end loop over faces in this direction

            Box facebox(grids[dit]);
            facebox.surroundingNodes(dir);
            for (BoxIterator bit(facebox); bit.ok(); ++bit)
              {
                const IntVect& iv = bit();
                IntVect ivm = iv - BASISV(dir);
                IntVect ivp = iv + BASISV(dir);
                if (mask(iv) != GROUNDEDMASKVAL)
                  {
                    if (dirBox.contains(ivm))
                     {
                        dirFlux(ivm) = 0.0;
                     }
                    if (dirBox.contains(ivp))
                     {
                        dirFlux(ivp) = 0.0;
                     }
                  }
              }

        } // end loop over directions
    } // end loop over grids
  } // end loop over levels
}

void AMRSubglacial::compute_RhsHcav(Vector<LevelData<FArrayBox>* >&      a_rhs,
         Vector<LevelData<FArrayBox>* >&   a_phi,
         Vector<LevelData<FArrayBox>* >&   a_hcav,
         Vector<LevelData<FArrayBox>* >&   a_hcav_old,
	 const Vector<RefCountedPtr<LevelSigmaCS> >& a_geometry,
         const Vector<LevelData<FArrayBox>* >&   a_cellVel,
         Vector<LevelData<FArrayBox>* >&   a_runoff,
         Real a_dt)
 {

  for (int lev = 0 ; lev<= m_finestLevel; lev++)
  {
    const LevelSigmaCS& levelCoords = *a_geometry[lev];
    const DisjointBoxLayout& grids = levelCoords.grids();
    Real iceDensity = levelCoords.iceDensity();
    Real waterDensity = levelCoords.waterDensity();

    for (DataIterator dit = (*a_rhs[lev]).dataIterator(); dit.ok(); ++dit)
      {
        FArrayBox& thisRHS = (*a_rhs[lev])[dit()];
        Box thisBox = thisRHS.box();

        FArrayBox& thisPHI = (*a_phi[lev])[dit()];
        FArrayBox& thisHCAV = (*a_hcav[lev])[dit()];
        FArrayBox& thisHCAVOLD = (*a_hcav_old[lev])[dit()];

        const FArrayBox& thisZb = levelCoords.getTopography()[dit];
        const FArrayBox& thisZs = levelCoords.getSurfaceHeight()[dit];
        const FArrayBox& thisUb = (*a_cellVel[lev])[dit];
        const FArrayBox& thisRunoff = (*a_runoff[lev])[dit];

        const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];

        BoxIterator bit(thisBox);
        for (bit.begin(); bit.ok(); ++bit)
           {
             const IntVect& iv = bit();
             if (mask(iv) == GROUNDEDMASKVAL )
               {
                 //Now express as a function of iv
                 Real velMag = sqrt(pow(thisUb(iv, 0), 2) + pow(thisUb(iv, 1), 2));
                 Real g_3     = waterDensity*m_basalMelt/iceDensity + abs(velMag)*max((m_hr - thisHCAV(iv))/m_lr, 0.0) ;
                 thisRHS(iv) = thisHCAVOLD(iv) + a_dt*g_3;
                
                 //thisRHS(iv) = thisHCAVOLD(iv) + a_dt*(m_basalMelt+ thisRunoff(iv));
               }
             else
               {
                 thisRHS(iv) = 0.0;
               }
           }

    } // end loop over grids
  } // end loop over levels
}


/*----------------------------------------------------------------*/
//            end: Set coefficients for h_cav
/*----------------------------------------------------------------*/


/*----------------------------------------------------------------*/
//            Utility functions
/*----------------------------------------------------------------*/

void AMRSubglacial::computeSinkTerm(
         Vector<LevelData<FArrayBox>* >&   a_source,
         Vector<LevelData<FArrayBox>* >&   a_phi,
         AmrIce& a_amrIce,
	 const Vector<RefCountedPtr<LevelSigmaCS> >& a_geometry,
         Real a_dt)
{
  for (int lev = 0 ; lev <= m_finestLevel; lev++)
    {
      
      RealVect dxLevel = a_amrIce.dx(lev);//dx of level
      const LevelSigmaCS& levelCoords = *a_geometry[lev];
      const DisjointBoxLayout& grids = levelCoords.grids();
      Real iceDensity = levelCoords.iceDensity();
      Real waterDensity = levelCoords.waterDensity();

      for (DataIterator dit(m_grids[lev]);dit.ok();++dit)
       {
          FArrayBox& gridSource = (*a_source[lev])[dit];
          FArrayBox& thisPHI = (*a_phi[lev])[dit()];
	  
          const FArrayBox& thick = levelCoords.getH()[dit];
       	  const FArrayBox& thisZs = levelCoords.getSurfaceHeight()[dit];
          const FArrayBox& thisZb = levelCoords.getTopography()[dit];
          const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
          
          for(int dir = 0; dir<SpaceDim; dir++)
           {
              Box facebox(grids[dit]);
              facebox.surroundingNodes(dir);
              for (BoxIterator bit(facebox); bit.ok(); ++bit)
               {
                  const IntVect& iv = bit();
                  IntVect ivm = iv - BASISV(dir);
                  IntVect ivp = iv + BASISV(dir);
                  Real tol = 0.001;
                  
                  Real source = 1.0;
                  
                  if (thisZs(iv) - thisZb(iv)<tol|| thisZs(ivm) - thisZb(ivm)<tol||thisZs(ivp) - thisZb(ivp) < tol)
                   {
                      //RHS = A = B (large) * phi_boundary_val
                      source =  m_BigaPhiCoef*( iceDensity*thisZb(iv) - waterDensity*min(thisZb(iv), 0.0) - thisPHI(iv));
                   }
                  if (thisZb(iv) <0.0 && ((thisZs(iv) + (waterDensity - iceDensity)/iceDensity*thisZb(iv))<tol || (thisZs(ivm) + (waterDensity - iceDensity)/iceDensity*thisZb(ivm))<tol || (thisZs(ivp) + (waterDensity - iceDensity)/iceDensity*thisZb(ivp))<tol) &&  mask(iv) == GROUNDEDMASKVAL )
                   {
                      source =  m_BigaPhiCoef *( waterDensity*thisZb(iv) - waterDensity*min(thisZb(iv), 0.0) - thisPHI(iv));
                   }
                  gridSource(iv) = source;
               }
           }
       }
   }
}
// Area integral of source terms near boundary--- equivalent to 
// boundary integral of hydraulic flux out of system-
Real AMRSubglacial::compute_fluxOutRate(
         Vector<LevelData<FArrayBox>* >&   a_phi,
         AmrIce& a_amrIce,
	 const Vector<RefCountedPtr<LevelSigmaCS> >& a_geometry,
         Vector<LevelData<FArrayBox>* >&   a_runoff,
         Real a_dt)
{

  Real fluxOut = 0.0;
  //for (int lev = 0 ; lev<= m_finestLevel; lev++)
  for (int lev = m_finestLevel ; lev<= (m_finestLevel); lev++)
  {
    RealVect dxLevel = a_amrIce.dx(lev);//dx of level
    
    const LevelSigmaCS& levelCoords = *a_geometry[lev];
    const DisjointBoxLayout& grids = levelCoords.grids();
    Real iceDensity = levelCoords.iceDensity();
    Real waterDensity = levelCoords.waterDensity();

    for (DataIterator dit = (*a_phi[lev]).dataIterator(); dit.ok(); ++dit)
      {
        FArrayBox& thisPHI = (*a_phi[lev])[dit()];
        Box thisBox = thisPHI.box();


        const FArrayBox& thisZb = levelCoords.getTopography()[dit];
        const FArrayBox& thisZs = levelCoords.getSurfaceHeight()[dit];

        const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];


         for(int dir = 0; dir<SpaceDim; dir++)
           {
             Box facebox(grids[dit]);
             facebox.surroundingNodes(dir);
             for (BoxIterator bit(facebox); bit.ok(); ++bit)
               {
                 const IntVect& iv = bit();
                 IntVect ivm = iv - BASISV(dir);
                 IntVect ivp = iv + BASISV(dir);
                 Real tol = 0.001;
                 if (mask(iv) == GROUNDEDMASKVAL && (mask(ivp) !=GROUNDEDMASKVAL || mask(ivm) != GROUNDEDMASKVAL))
                {
                  if (m_g*waterDensity*thisZb(iv) - m_g* waterDensity*min(thisZb(iv), 0.0) - thisPHI(iv) < 0)
                    {
                     Real source =  -m_BigaPhiCoef *(m_g*waterDensity*thisZb(iv) - m_g*waterDensity*min(thisZb(iv), 0.0) - thisPHI(iv));
                     fluxOut += dxLevel[0]*dxLevel[1]*source;
                    }
                }

//                 if (thisZs(iv) - thisZb(iv)<tol|| thisZs(ivm) - thisZb(ivm)<tol||thisZs(ivp) - thisZb(ivp) < tol)
//                   {
//                     //RHS = A = B (large) * phi_boundary_val
//                     Real source =  m_BigaPhiCoef *( iceDensity*thisZb(iv) - waterDensity*min(thisZb(iv), 0.0) - thisPHI(iv));
//                     fluxOut += dxLevel[0]*dxLevel[1]*source;
//                   }
//              if (thisZb(iv) <0.0 && ((thisZs(iv) + (waterDensity - iceDensity)/iceDensity*thisZb(iv))<tol || (thisZs(ivm) + (waterDensity - iceDensity)/iceDensity*thisZb(ivm))<tol || (thisZs(ivp) + (waterDensity - iceDensity)/iceDensity*thisZb(ivp))<tol) &&  mask(iv) == GROUNDEDMASKVAL )
//                   {
//                     Real source =  m_BigaPhiCoef *( waterDensity*thisZb(iv) - waterDensity*min(thisZb(iv), 0.0) - thisPHI(iv));
//                     fluxOut += dxLevel[0]*dxLevel[1]*source;
//                  //   pout()<<source<<endl;
//                  //   pout()<<dxLevel<<endl;
//                   }
               }

           }
    } // end loop over grids
  } // end loop over levels
  return fluxOut/a_dt;
}

// Area integral of total runoff

Real AMRSubglacial::compute_totalRunoffRate(
         Vector<LevelData<FArrayBox>* >&   a_phi,
         AmrIce& a_amrIce,
	 const Vector<RefCountedPtr<LevelSigmaCS> >& a_geometry,
         Vector<LevelData<FArrayBox>* >&   a_runoff,
         Real a_dt)
{
  Real totalRunoff = 0.0;
//  for (int lev = 0 ; lev<= m_finestLevel; lev++)
  for (int lev = 0 ; lev<= 0; lev++)
  {
    RealVect dxLevel = a_amrIce.dx(lev);//dx of level
    
    const LevelSigmaCS& levelCoords = *a_geometry[lev];
    const DisjointBoxLayout& grids = levelCoords.grids();
    Real iceDensity = levelCoords.iceDensity();
    Real waterDensity = levelCoords.waterDensity();

    for (DataIterator dit = (*a_runoff[lev]).dataIterator(); dit.ok(); ++dit)
      {
        FArrayBox& thisRunoff = (*a_runoff[lev])[dit()];
        Box thisBox = thisRunoff.box();

        const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];

        Real tol = 0.001;
        BoxIterator bit(thisBox);
        int counter = 0;
        for (bit.begin(); bit.ok(); ++bit)
           {
             const IntVect& iv = bit();
             if (mask(iv) == GROUNDEDMASKVAL )
               {
                 //Now express as a function of iv
                  Real temprunoff = thisRunoff(iv);
                  totalRunoff += dxLevel[0]*dxLevel[1]*thisRunoff(iv);
                  counter +=1;
               }
             else
               {
                 totalRunoff = totalRunoff;
               }
           }

    } // end loop over grids
  } // end loop over levels
  return totalRunoff;
}

/*----------------------------------------------------------------*/
//            end: Set coefficients for h_cav
/*----------------------------------------------------------------*/


int AMRSubglacial::problemSize()
{
    return m_hydraulicPotential.size();
}



#ifdef CH_USE_HDF5
void AMRSubglacial::addPlotVars(Vector<std::string>& a_vars)
{
  a_vars.push_back("HydraulicPotential");
  a_vars.push_back("CavityHeight");
  a_vars.push_back("EffectivePressure");
}

void AMRSubglacial::writePlotData(LevelData<FArrayBox>& a_data, int a_level)
{
  for (DataIterator dit(m_grids[a_level]);dit.ok();++dit)
    {
      int j = 0;
      a_data[dit].copy( (*m_hydraulicPotential[a_level])[dit],0,j++,1);
      a_data[dit].copy( (*m_cavityHeight[a_level])[dit],0,j++,1);
      a_data[dit].copy( (*m_effectivePressure[a_level])[dit],0,j++,1);
    }
}



/// fill a_var with the names of variables to add to the checkpoint file
void AMRSubglacial::addCheckVars(Vector<std::string>& a_vars)
{
  a_vars.push_back("Phi");
  a_vars.push_back("hcav");
  a_vars.push_back("N");
}

/// copy level a_level checkpoint data to  LevelData<FArrayBox>& a_data
void AMRSubglacial::writeCheckData(HDF5Handle& a_handle, int a_level)
{
  write(a_handle, *m_hydraulicPotential[a_level], "HydraulicPotential", m_hydraulicPotential[a_level]->ghostVect());
  write(a_handle, *m_cavityHeight[a_level], "CavityHeight", m_cavityHeight[a_level]->ghostVect());
  write(a_handle, *m_effectivePressure[a_level], "EffectivePressure", m_effectivePressure[a_level]->ghostVect());
}

///  read checkpoint data for a single  LevelData<FArrayBox> TODO MM: This is clearly not done yet!!
LevelData<FArrayBox>* AMRSubglacial::readCheckData(const std::string& a_varLabel,
                                               const std::string& a_dataLabel,
                                               HDF5Handle& a_handle,
                                               HDF5HeaderData&  a_header,
                                               const DisjointBoxLayout& a_grids)
{
  LevelData<FArrayBox>* field(NULL);
  bool contains(false);
  map<std::string, std::string>::const_iterator i;
  for (i = a_header.m_string.begin(); i!= a_header.m_string.end(); ++i)
    {
      contains |= (i->second == a_varLabel);
     }

   field = new LevelData<FArrayBox>(a_grids, DAMAGE_N_COMP,DAMAGE_N_GHOST*IntVect::Unit);
   if (contains)
     {
        int dataStatus =  read<FArrayBox>(a_handle, *field, a_dataLabel , a_grids);
        if (dataStatus != 0)
           {
              MayDay::Error("failed to read damage data from checkpoint, but the header indicated its presence. set to zero");
            }
      }
   else
     {
       MayDay::Warning("AMRDamage checkpoint data missing: setting to zero");
       for (DataIterator dit = field->dataIterator(); dit.ok(); ++dit)
         {
             (*field)[dit].setVal(0.0);
         }
      }
    return field;
}

//Read level a_level checkpoint data
void AMRSubglacial::readCheckData(HDF5Handle& a_handle, HDF5HeaderData&  a_header, int a_level, const DisjointBoxLayout& a_grids)
{
  LevelData<FArrayBox>* t = readCheckData("Phi", "HydraulicPotential", a_handle, a_header, a_grids);
  
  if (t)
    {
      if (m_hydraulicPotential.size() < a_level + 1)
        {
          m_hydraulicPotential.resize(a_level+1, NULL);
        }
      if (m_hydraulicPotential[a_level])
        {
          delete m_hydraulicPotential[a_level];
          m_hydraulicPotential[a_level] = NULL;
        }
      m_hydraulicPotential[a_level] = t;
    }

   t = readCheckData("hcav", "CavityHeight", a_handle, a_header, a_grids);
  
  if (t)
    {
      if (m_cavityHeight.size() < a_level + 1)
        {
          m_cavityHeight.resize(a_level+1, NULL);
        }
      if (m_cavityHeight[a_level])
        {
          delete m_cavityHeight[a_level];
          m_cavityHeight[a_level] = NULL;
        }
      m_cavityHeight[a_level] = t;
    }

   t = readCheckData("N", "EffectivePressure", a_handle, a_header, a_grids);
  
  if (t)
    {
      if (m_effectivePressure.size() < a_level + 1)
        {
          m_effectivePressure.resize(a_level+1, NULL);
        }
      if (m_effectivePressure[a_level])
        {
          delete m_effectivePressure[a_level];
          m_effectivePressure[a_level] = NULL;
        }
      m_effectivePressure[a_level] = t;
    }

}




///// read level a_level checkpoint data from  LevelData<FArrayBox>& a_data TODO: Michael Need to be able to do both phi and hcav
//void AMRSubglacial::readCheckData(HDF5Handle& a_handle, HDF5HeaderData&  a_header, int a_level, const DisjointBoxLayout& a_grids)
//{
//  bool containsSubglacialData(false);
//  map<std::string, std::string>::const_iterator i;
//  for (i = a_header.m_string.begin(); i!= a_header.m_string.end(); ++i)
//    {
//      containsSubglacialData |= (i->second == "Phi");
//    }
//
//  if (containsSubglacialData)
//    {
//
//      if (m_hydraulicPotential.size() <= a_level)
//	{
//	  m_hydraulicPotential.resize(a_level + 1, NULL);
//	  if (m_hydraulicPotential[a_level] != NULL)
//	    {
//	      delete m_hydraulicPotential[a_level];
//	    }
//
//	}
//      m_hydraulicPotential[a_level] = new LevelData<FArrayBox>(a_grids, SUBGLACIAL_N_COMP, SUBGLACIAL_N_GHOST * IntVect::Unit);
//      int dataStatus =  read<FArrayBox>(a_handle, *m_hydraulicPotential[a_level], "HydraulicPotential", a_grids);
//      if (dataStatus != 0)
//	{
//	  MayDay::Error("failed to read subglacial data from checkpoint, but the header indicated its presence");
//	}
//    }
//}
#endif
