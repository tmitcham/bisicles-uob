#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _AMRSUBGLACIAL_H_
#define _AMRSUBGLACIAL_H_

#include "AmrIce.H"
#include "LevelSigmaCS.H"
//Added by Michael
#include "VCAMRPoissonOp2.H"
#include "AMRMultiGrid.H"
//End added by Michael
#include "NamespaceHeader.H"

class SubglacialHydrologyModel;

#define SUBGLACIAL_N_GHOST 4
#define SUBGLACIAL_N_COMP 1



class SubglacialParameters
{
public:
  IntVect       nCells;
  int           maxGridSize;
  int           blockFactor;
  int           bufferSize;
  Real          fillRatio;
  Real          refineThresh;
  int           domBcType;
  int           maxLevel;
  int           numLevels;
  Vector<int>   refRatio;
  ProblemDomain coarsestDomain;
  Real          coarsestDx;
  RealVect      domainLength;
  RealVect      probLo;
  RealVect      probHi;
  Real          alpha=1.0;
  Real          beta=1.0;
  int           coefficient_average_type;
  int           verbosity;

};



  extern AMRLevelOpFactory<LevelData<FArrayBox> >*
  defineOperatorFactory_phi(
                        const Vector<DisjointBoxLayout>&             a_grids,
                        Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_aCoef,
                        Vector<RefCountedPtr<LevelData<FluxBox> > >& a_bCoef,
                        SubglacialParameters&                  a_params);

  extern AMRLevelOpFactory<LevelData<FArrayBox> >*
  defineOperatorFactory_hcav(
                        const Vector<DisjointBoxLayout>&               a_grids,
                        Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_aCoef,
                        Vector<RefCountedPtr<LevelData<FluxBox> > >&   a_bCoef,
                        SubglacialParameters&                     a_params);

  void ParseBC_phi(FArrayBox& a_state,
               const Box& a_valid,
               const ProblemDomain& a_domain,
               Real a_dx,
               bool a_homogeneous);

  void ParseValue(Real* pos,
                  int* dir,
                  Side::LoHiSide* side,
                  Real* a_values);

  void ParseBC_hcav(FArrayBox& a_state,
               const Box& a_valid,
               const ProblemDomain& a_domain,
               Real a_dx,
               bool a_homogeneous);


/// Subglacial hydrology model.
class AMRSubglacial
{

  friend class SubglacialHydrologyModel;


  /// AMR mesh hierarchy
  int m_finestLevel;
  Vector<DisjointBoxLayout> m_grids;
  Vector<int> m_ratio;
  Vector<RealVect> m_dx;

  //TODO: For now we will store the parameters here
//  Real m_Ahat = 3.1436e-18; //Pa^(-n)/yr
  Real m_Ahat; //Pa^(-n)/yr
  Real m_g;//9.7e15 m/yr^2
  Real m_n; //Glen's flow law exponent

  Real  m_sigma; // 10e-4
  Real  m_basalMelt;//basal melt m/yr
  Real  m_alpha; //hydraulic permeability exponent
  Real  m_epsilon_d;//m^2 /yr
  Real  m_h0;  // typical cavity thickness m
  Real  m_hr; // obstacale height m
  Real  m_K; //9.45e6 typical hydraulic permeability m/yr
  Real  m_lr; //obstacle length scale m
  Real  m_N0;//9e5typical hydraulic pressure scale Pa (time comp?)

  Real m_BigaPhiCoef; // Big number phi boundary condition
  Real m_BigNotGrounded; // Big number for outside domain conditions
  

public:
  ~AMRSubglacial();
  AMRSubglacial();
  
  //TODO Michael: replace with a getTime function
  Real m_time = 0.0;
  int m_time_step = 0;
  bool m_effectivePressureDefined = false;
  
  const LevelData<FArrayBox>* hydraulicPotential(int a_level) const;
  const LevelData<FArrayBox>* cavityHeight(int a_level) const;
  const LevelData<FArrayBox>* effectivePressure(int a_level) const;

  int finestLevel() const
  {
    return m_finestLevel;
  }

  
  const Vector<LevelData<FArrayBox>*> getEffectivePressure();

  void setInitialValues(
       AmrIce& a_amrIce,
       const Vector<RefCountedPtr<LevelSigmaCS> >& a_geometry);
  
  void define(const Vector<DisjointBoxLayout>& a_grids,
	      const Vector<int>& a_ratio,
	      int a_finestLevel,
	      const RealVect& a_crseDx,
              const Vector<RefCountedPtr<LevelSigmaCS> >& a_geometry);

  void runLargeTimestep(Real a_dt,
       AmrIce& a_amrIce,
			 const Vector<RefCountedPtr<LevelSigmaCS> >& a_geometry,
			 const Vector<LevelData<FArrayBox>* >& a_cellVel,
			 const Vector<LevelData<FluxBox>* >& a_faceVel);

  void subglacialAdvance(Real a_dt,
         Vector<RefCountedPtr<LevelData<FArrayBox> > > a_aCoef_phi,
         Vector<RefCountedPtr<LevelData<FluxBox> > > a_bCoef_phi,
         Vector<LevelData<FArrayBox> * >  a_rhs_phi,
         Vector<RefCountedPtr<LevelData<FArrayBox> > > a_aCoef_hcav,
         Vector<RefCountedPtr<LevelData<FluxBox> > > a_bCoef_hcav,
         Vector<LevelData<FArrayBox> * > a_rhs_hcav,
         Vector<LevelData<FArrayBox>* > a_phiOld,
         Vector<LevelData<FArrayBox>* > a_hcavOld,
         Vector<ProblemDomain> a_vectDomains,
         Vector<RealVect> a_vectDx,
         SubglacialParameters a_subParams,
         Vector<LevelData<FArrayBox>* > a_runoff,
  			 const Vector<RefCountedPtr<LevelSigmaCS> >& a_geometry,
  			 const Vector<LevelData<FArrayBox>* >& a_cellVel,
  			 const Vector<LevelData<FluxBox>* >& a_faceVel);

  int problemSize();
  
  void setObserverEffectivePressure(
         Vector<LevelData<FArrayBox>*>  a_observerEffectivePressure);

#ifdef CH_USE_HDF5
  //plot file methods

  /// add the names of variables to add to the plot file to a_var, in the order they will be written
  void addPlotVars(Vector<std::string>& a_vars);

  /// copy level a_level plot data to   LevelData<FArrayBox>& a_data
  void writePlotData(LevelData<FArrayBox>& a_data, int a_level);

  //check point methods

  /// fill a_var with the names of variables to add to the checkpoint file
  void addCheckVars(Vector<std::string>& a_vars);

  /// copy level a_level checkpoint data to  LevelData<FArrayBox>& a_data
  void writeCheckData(HDF5Handle& a_handle, int a_level);

  // read checkpoint data for a single LevelData<FArrayBox>
  LevelData<FArrayBox>* readCheckData(const std::string& a_varLabel,
                                               const std::string& a_dataLabel,
                                               HDF5Handle& a_handle,
                                               HDF5HeaderData&  a_header,
                                               const DisjointBoxLayout& a_grids);
  /// read level a_level checkpoint data from  LevelData<FArrayBox>& a_data
  void readCheckData(HDF5Handle& a_handle,  HDF5HeaderData& a_header, int a_level, const DisjointBoxLayout& a_grids);
#endif

protected:
  //cell-centered hydraulic potential and cavity height
  mutable Vector<LevelData<FArrayBox>* > m_hydraulicPotential;
  mutable Vector<LevelData<FArrayBox>* > m_cavityHeight;
  mutable Vector<LevelData<FArrayBox>*> m_effectivePressure;

private:

  //set values of surface runoff
  void setRunoff(Vector<LevelData<FArrayBox>* >& a_runoff,
			      const Vector<RefCountedPtr<LevelSigmaCS> >& a_geometry,
			      Real a_dt);
  // Set value of effective pressure
  void setEffectivePressure(
         const Vector<RefCountedPtr<LevelSigmaCS> >& a_geometry,
         const Real a_elapsedTime,
         const Real a_small_dt);

 // Take in the set the observers effective pressure to the calculated one

  //Suite of functions that set Helmholtz coefficients
  void compute_aCoefPhi(Vector<RefCountedPtr<LevelData<FArrayBox>>> a_aCoef,
         Vector<LevelData<FArrayBox>* >&   a_phi,
         Vector<LevelData<FArrayBox>* >&   a_hcav,
         const Vector<RefCountedPtr<LevelSigmaCS> >& a_geometry,
         const Vector<LevelData<FArrayBox>* >&   a_cellVel,
         Vector<LevelData<FArrayBox>* >&   a_runoff,
         Real a_dt
        );

  void compute_bCoefPhi(Vector<RefCountedPtr<LevelData<FluxBox>>> a_bCoef,
         Vector<LevelData<FArrayBox>* >&   a_phi,
         Vector<LevelData<FArrayBox>* >&   a_hcav,
 	 const Vector<RefCountedPtr<LevelSigmaCS> >& a_geometry,
         const Vector<LevelData<FArrayBox>* >&   a_cellVel,
         Vector<LevelData<FArrayBox>* >&   a_runoff,
         Real a_dt);

  void compute_RhsPhi(Vector<LevelData<FArrayBox>* >&      a_rhs,
          Vector<LevelData<FArrayBox>* >&   a_phi,
          Vector<LevelData<FArrayBox>* >&   a_hcav,
          Vector<LevelData<FArrayBox>* >&   a_phi_old,
          const Vector<RefCountedPtr<LevelSigmaCS> >& a_geometry,
          const Vector<LevelData<FArrayBox>* >&   a_cellVel,
          Vector<LevelData<FArrayBox>* >&   a_runoff,
          Real a_dt
        );

  void compute_aCoefHcav(Vector<RefCountedPtr<LevelData<FArrayBox>>> a_aCoef,
        Vector<LevelData<FArrayBox>* >&   a_phi,
        Vector<LevelData<FArrayBox>* >&   a_hcav,
	const Vector<RefCountedPtr<LevelSigmaCS> >& a_geometry,
        const Vector<LevelData<FArrayBox>* >&   a_cellVel,
        Vector<LevelData<FArrayBox>* >&   a_runoff,
        Real a_dt
        );

  void compute_bCoefHcav(Vector<RefCountedPtr<LevelData<FluxBox>> > a_bCoef,
         Vector<LevelData<FArrayBox>* >&   a_phi,
         Vector<LevelData<FArrayBox>* >&   a_hcav,
 	 const Vector<RefCountedPtr<LevelSigmaCS> >& a_geometry,
         const  Vector<LevelData<FArrayBox>* >&   a_cellVel,
         Vector<LevelData<FArrayBox>* >&   a_runoff,
         Real a_dt);

  void compute_RhsHcav(Vector<LevelData<FArrayBox>* >&      a_rhs,
         Vector<LevelData<FArrayBox>* >&   a_phi,
         Vector<LevelData<FArrayBox>* >&   a_hcav,
         Vector<LevelData<FArrayBox>* >&   a_hcav_old,
         const Vector<RefCountedPtr<LevelSigmaCS> >& a_geometry,
         const Vector<LevelData<FArrayBox>* >&   a_cellVel,
         Vector<LevelData<FArrayBox>* >&   a_runoff,
         Real a_dt);

  void computeSinkTerm(
         Vector<LevelData<FArrayBox>* >&   a_source,
         Vector<LevelData<FArrayBox>* >&   a_phi,
         AmrIce& a_amrIce,
	 const Vector<RefCountedPtr<LevelSigmaCS> >& a_geometry,
         Real a_dt);

  Real compute_fluxOutRate(
         Vector<LevelData<FArrayBox>* >&   a_phi,
         AmrIce& a_amrIce,
	 const Vector<RefCountedPtr<LevelSigmaCS> >& a_geometry,
         Vector<LevelData<FArrayBox>* >&   a_runoff,
         Real a_dt);


  Real compute_totalRunoffRate(
         Vector<LevelData<FArrayBox>* >&   a_phi,
         AmrIce& a_amrIce,
	 const Vector<RefCountedPtr<LevelSigmaCS> >& a_geometry,
         Vector<LevelData<FArrayBox>* >&   a_runoff,
         Real a_dt);

};





/// Contains an AMRSubglacial, responds to changes in the ice sheet state
class SubglacialIceObserver : public AmrIce::Observer
{
 //Michael TODO: might need multiple pointers, one to each field
  AMRSubglacial* m_subglacialPtr;

public:

  SubglacialIceObserver();

  ~SubglacialIceObserver();

  //Each of these is a function that returns the field
  AMRSubglacial& subglacial() const;

  /// General purpose notify method: do something in response to a signal from an AmrIce
  void notify(AmrIce::Observer::Notification a_n, AmrIce& a_amrIce);

  // output effective pressure

  void setEffectivePressure(Vector<LevelData<FArrayBox>*> a_effectivePressure);
 //TODO: delete 
  bool m_effectivePressureDefined = false;

  // plot file methods

  /// add the names of variables to add to the plot file to a_var, in the order they will be written
  void addPlotVars(Vector<std::string>& a_vars);

  /// copy level a_level plot data to   LevelData<FArrayBox>& a_data
  void writePlotData(LevelData<FArrayBox>& a_data, int a_level);

  //check point methods

  /// fill a_var with the names of variables to add to the checkpoint file
  void addCheckVars(Vector<std::string>& a_vars);

  /// copy level a_level checkpoint data to  LevelData<FArrayBox>& a_data
  void writeCheckData(HDF5Handle& a_handle, int a_level);

  /// read level a_level checkpoint data from  LevelData<FArrayBox>& a_data
  void readCheckData(HDF5Handle& a_handle,  HDF5HeaderData& a_header, int a_level, const DisjointBoxLayout& a_grids);

};





#include "NamespaceFooter.H"

#endif
