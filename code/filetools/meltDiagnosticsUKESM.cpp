#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//===========================================================================
// diagnostics.cpp
// read in a bisicles plotfile 
// and an optional mask, and write out a bunch of diagnostics about the ice sheet. 
// These are
// 0. time
// 1. Volume of ice
// 2. Volume of ice above flotation
// dh/dt
// SMB
// BMB
// Volume of calved ice
// Calving flux  
//===========================================================================

#include <iostream>
#include "ParmParse.H"
#include "AMRIO.H"
#include "LevelSigmaCS.H"
#include "IceConstants.H"
#include "computeSum.H"
#include "CellToEdge.H"
#include "amrIceF_F.H"
#include "AdvectPhysics.H"
#include "PatchGodunov.H"
#include "SigmaCSF_F.H"
#include "PiecewiseLinearFillPatch.H"
#include "FillFromReference.H"
#include "IceThicknessIBC.H"
#include "LevelDataIBC.H"
#include "DomainDiagnosticData.H"
#include <functional>
 

/**
   convert single level mask data, to which stores integer mask labels (as reals...) 
   between a_mask_label_start and a_mask_no, to multi-level real fractional coverage areas for each 
   distinct mask label. Resulting a_mask_fraction FABS will have a_mask_no components
*/
void maskToAMR(Vector<LevelData<FArrayBox>* >& a_mask_fraction,
	       const Vector<Real>& a_mask_fraction_dx,
	       const LevelData<FArrayBox>& a_mask,
	       int a_mask_id_start,
	       int a_n_mask_ids,
	       const Real& a_mask_dx)
{

  
  Real tol = 1.0e-10; // tolerance in Abs(x - y) < tol test
  
  for (int comp = 0; comp < a_n_mask_ids ; comp++)
    {
      Real mask_id = Real(comp + a_mask_id_start);
      
      // create field of mask areas on the mask level (all cells 1 or 0)
      const DisjointBoxLayout mask_grids = a_mask.disjointBoxLayout();
      LevelData<FArrayBox> mask_fraction(mask_grids,1,IntVect::Zero);
      for (DataIterator dit(mask_grids); dit.ok(); ++dit)
	{
	  for (BoxIterator bit(mask_grids[dit]); bit.ok(); ++bit)
	    {
	      const IntVect& iv = bit();
	      mask_fraction[dit](iv) = (Abs(a_mask[dit](iv) - mask_id) < tol)?1.0:0.0;  
	    }
	}
      // coarsen/refine the mask areas, store in component mask_no - a_mask_no_start)
      for (int lev = 0; lev < a_mask_fraction.size(); lev++)
	{
	  LevelData<FArrayBox>& dest = *a_mask_fraction[lev];
	  LevelData<FArrayBox> tmp(dest.disjointBoxLayout(),1,IntVect::Zero);
	  RealVect src_dx = RealVect::Unit * a_mask_dx;
	  RealVect dest_dx = RealVect::Unit * a_mask_fraction_dx[lev];
	  FillFromReference(tmp, mask_fraction, dest_dx, src_dx, false);
	  tmp.localCopyTo(Interval(0,0), dest, Interval(comp, comp));
	}
    }
}


struct NameUnitValue
{
  std::string name;
  std::string unit;
  Real value;
  NameUnitValue(std::string n, std::string u, Real v)
    :name(n), unit(u), value(v){}

};

std::ostream& json(std::ostream& o, NameUnitValue& a)
  {
    o << "\"" << a.name << "\":{\"unit\":\"" << a.unit << "\", \"value\":" << a.value << "}";
    return o;
  }

std::ostream& operator<<(std::ostream& o, NameUnitValue& a)
  {
    o << a.name << "," << a.unit << "," << a.value;
    return o;
  }


void reportConservationInside(Vector<NameUnitValue>& report,
			      function<bool(Real h, Real f, int mask)> inside,
			      const Vector<RefCountedPtr<LevelSigmaCS > >& coords,
			      const Vector<LevelData<FArrayBox>* >& NEMOBasalThicknessSource, 
			      const Vector<LevelData<FArrayBox>* >& basalThicknessSource,
			      const Vector<LevelData<FArrayBox>* >& activeBasalThicknessSource,
			      const Vector<LevelData<FArrayBox>* >& topography,
			      const Vector<LevelData<FArrayBox>* >& thickness,
			      const Vector<LevelData<FArrayBox>* >& iceFrac, 
			      const Vector<LevelData<FArrayBox>* >& sectorMaskFraction,
			      const Vector<Real>& dx, const Vector<int>& ratio, 
			      int finestLevel, int maskNo, int maskComp)
{

  CH_TIME("reportConservationInside");
  // just to make the args less reptitive...

  auto sumScalar = [inside,coords,topography, thickness, iceFrac, sectorMaskFraction, dx, ratio, finestLevel, maskNo, maskComp ]
    (const Vector<LevelData<FArrayBox>* >& scalar)
		   {
		     Real sumScalar;
		     MaskedIntegration::integrateScalarInside
		       (sumScalar, inside ,coords, scalar, topography, thickness,
			iceFrac,sectorMaskFraction, dx, ratio, finestLevel,  maskNo, maskComp);
		     return sumScalar;
		   };

  std::string dhunit("m3/a");

  Real NEMO = sumScalar(NEMOBasalThicknessSource);
  report.push_back(NameUnitValue("NEMO",dhunit,NEMO));

  Real BMB = sumScalar(basalThicknessSource);
  report.push_back(NameUnitValue("BMB",dhunit,BMB));
  
  Real activeBMB = sumScalar(activeBasalThicknessSource);
  report.push_back(NameUnitValue("activeBMB",dhunit,activeBMB));
  
}


void extractThicknessSource(Vector<LevelData<FArrayBox>* >& NEMOBasalThicknessSource, 
			    Vector<LevelData<FArrayBox>* >& basalThicknessSource,
			    Vector<LevelData<FArrayBox>* >& activeBasalThicknessSource,
			    Vector<LevelData<FArrayBox>* >& topography,
			    const Vector<Real>& dx, const Vector<int>& ratio, 
			    const Vector<std::string>& name, 
			    const Vector<LevelData<FArrayBox>* >& data)

{
  int numLevels = data.size();

  for (int lev = 0; lev < numLevels; lev++)
    {

      const DisjointBoxLayout& grids = topography[lev]->disjointBoxLayout();
      for (DataIterator dit(grids); dit.ok(); ++dit)
	{
          (*NEMOBasalThicknessSource[lev])[dit].setVal(0.0);
          (*basalThicknessSource[lev])[dit].setVal(0.0);
          (*activeBasalThicknessSource[lev])[dit].setVal(0.0);
	}

      for (int j = 0; j < name.size(); j++)
	{
	  if (name[j] == "melt_water")
	    {
	      data[lev]->copyTo(Interval(j,j),*NEMOBasalThicknessSource[lev],Interval(0,0));
	    }
	  else if (name[j] == "basalThicknessSource")
	    {
	      data[lev]->copyTo(Interval(j,j),*basalThicknessSource[lev],Interval(0,0));
	    }	  
	  else if (name[j] == "activeBasalThicknessSource")
	    {
	      data[lev]->copyTo(Interval(j,j),*activeBasalThicknessSource[lev],Interval(0,0));
	    }
	}
     

      if (lev > 0)
	{
	  
	  const DisjointBoxLayout& crseGrids = topography[lev-1]->disjointBoxLayout();
	  PiecewiseLinearFillPatch filler(grids , crseGrids, 1, 
					  crseGrids.physDomain(), ratio[lev-1], 1);
	  Real time_interp_coeff = 0.0;

	  filler.fillInterp(*NEMOBasalThicknessSource[lev],*NEMOBasalThicknessSource[lev-1],
			    *NEMOBasalThicknessSource[lev-1],time_interp_coeff,0, 0, 1);
	  filler.fillInterp(*basalThicknessSource[lev],*basalThicknessSource[lev-1],
			    *basalThicknessSource[lev-1],time_interp_coeff,0, 0, 1);
	  filler.fillInterp(*activeBasalThicknessSource[lev],*activeBasalThicknessSource[lev-1],
			    *activeBasalThicknessSource[lev-1],time_interp_coeff,0, 0, 1);
	}
      NEMOBasalThicknessSource[lev]->exchange();
      basalThicknessSource[lev]->exchange();
      activeBasalThicknessSource[lev]->exchange();
    }
} 

void stateDiagnostics(std::ostream& sout, bool append, std::string plot_file,
		      const Vector<LevelData<FArrayBox>* >& sectorMaskFraction,
		      int maskNoStart, int maskNoEnd,
		      const Vector<LevelData<FArrayBox>* >& data,
		      const Vector<DisjointBoxLayout >& grids,
		      const Vector<std::string>& name,
		      int numLevels,
		      Vector<Real>& dx, // cant be const due to createSigmaCS
		      Real dt, Real time,
		      Vector<int>& ratio,// cant be const due to createSigmaCS
		      Real iceDensity,
		      Real waterDensity,
		      Real gravity,
		      Real h_min, Real f_min)
{

  CH_TIME("stateDiagnostics");

  Vector<LevelData<FArrayBox>* > thickness(numLevels,NULL);
  Vector<LevelData<FArrayBox>* > topography(numLevels,NULL);
  Vector<LevelData<FArrayBox>* > NEMOBasalThicknessSource(numLevels,NULL);
  Vector<LevelData<FArrayBox>* > basalThicknessSource(numLevels,NULL);
  Vector<LevelData<FArrayBox>* > activeBasalThicknessSource(numLevels,NULL);
  Vector<LevelData<FArrayBox>* > iceFrac(numLevels,NULL);


  for (int lev=0;lev<numLevels;++lev)
    {
      // Extra ghost cells are needed to calculate advection using PPM. 
      thickness[lev] = new LevelData<FArrayBox>(grids[lev],1,4*IntVect::Unit);
      topography[lev] = new LevelData<FArrayBox>(grids[lev],1,4*IntVect::Unit);
      NEMOBasalThicknessSource[lev] = new LevelData<FArrayBox>(grids[lev],1,IntVect::Unit);
      basalThicknessSource[lev] = new LevelData<FArrayBox>(grids[lev],1,IntVect::Unit);
      activeBasalThicknessSource[lev] = new LevelData<FArrayBox>(grids[lev],1,IntVect::Unit);
      iceFrac[lev] = new LevelData<FArrayBox>(grids[lev],1,IntVect::Unit);
    }

    
  sout.setf(ios_base::scientific,ios_base::floatfield); 
  sout.precision(12);

  extractThicknessSource(NEMOBasalThicknessSource, basalThicknessSource, activeBasalThicknessSource,
			 topography, dx, ratio, name, data);
  
  
  // CSV style output of diagnostics
  if (!append)
    {
      sout << "csvheader,filename,time,maskNo,region,quantity,unit,value" << endl;
    }
  
  for (int maskNo = maskNoStart; maskNo <= maskNoEnd; ++maskNo)
    {
      int maskComp = maskNo - maskNoStart; // component of sectorMaskFraction FABs corresponding to maskNo
       
      typedef std::map<std::string, function<bool(Real, Real, int)> > MapSF;
      MapSF regions;
      
      auto entire = [h_min](Real h, Real f, int mask){ return true;} ;
      regions["entire"] = entire;
      
      for (MapSF::const_iterator mit = regions.begin(); mit != regions.end(); ++mit)
	{
	  Vector<NameUnitValue> report;
	 
	  reportConservationInside(report, mit->second,coords, NEMOBasalThicknessSource,
				   basalThicknessSource, activeBasalThicknessSource,topography, thickness, iceFrac,
				   sectorMaskFraction, dx, ratio, numLevels-1, maskNo, maskComp);
	  
	  for (int i = 0; i < report.size(); i++)
	    {
	      sout << "csvdata," << plot_file << "," << time << "," << maskNo << "," << mit->first << "," << report[i] << endl;
	    }
	}
    }
    
  for (int lev=0;lev<numLevels;++lev)
    {
      //if (sectorMask[lev] != NULL) delete sectorMask[lev];
      if (thickness[lev] != NULL) delete thickness[lev];
      if (topography[lev] != NULL) delete topography[lev];
      if (NEMOBasalThicknessSource[lev] != NULL) delete NEMOBasalThicknessSource[lev];
      if (basalThicknessSource[lev] != NULL) delete basalThicknessSource[lev];
      if (activeBasalThicknessSource[lev] != NULL) delete activeBasalThicknessSource[lev];
      if (iceFrac[lev] != NULL) delete iceFrac[lev];
    }
}


int main(int argc, char* argv[]) {

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif 

  { // Begin nested scope

#ifdef CH_MPI
    MPI_Barrier(Chombo_MPI::comm);
#endif
    int rank, number_procs;
#ifdef CH_MPI
    MPI_Comm_rank(Chombo_MPI::comm, &rank);
    MPI_Comm_size(Chombo_MPI::comm, &number_procs);
#else
    rank=0;
    number_procs=1;
#endif

    if(argc < 2) 
      {
	std::cerr << " usage: " << argv[0] << "plot_file=plot_file [-options] [key=value args]" << endl;
	std::cerr << " example: " << argv[0] << "plot_file=test.2d.hdf5 out_file=test.csv -append ice_density=918.0 water_density=1028.0"  << endl;
	std::cerr << "        computes diagnostic integrals from test.2d.hdf5, appends results to test.csv " << endl;
	std::cerr << " example: " << argv[0] << "plot_file=test.2d.hdf5 out_file=test.csv -append ice_density=918.0 water_density=1028.0 mask_file=mask2.d.hdf5 mask_no_start=0 mask_no_end=4"  << endl;
	std::cerr << "        computes diagnostic integrals from test.2d.hdf5, for the whole domain and each subdomain 1-4, appends results to test.csv, mask.2d.hdf5 indicates subdomains with a grid of values 1-4 " << endl;
	exit(1);
      }
  
    ParmParse pp(argc-1,argv+1,NULL,NULL);

    // the one mandatory arg : a plot file
    std::string plot_file; pp.get("plot_file",plot_file);

    // out_file : if not specified use stdio
    std::string out_file(""); pp.query("out_file", out_file);
    std::ofstream fout;
    bool append = pp.contains("append");
    if (out_file != "")
	{
	  fout.open(out_file, append?std::ofstream::app:std::ofstream::trunc);
	}
    std::ostream& sout = (out_file == "")?pout():fout;
      
    
    Real ice_density(917.0); pp.query("ice_density", ice_density);
    Real water_density(1028.0); pp.query("water_density", water_density);
    Real gravity(9.81); pp.query("gravity", gravity);

    std::string mask_file("");
    pp.query("mask_file",mask_file);
    bool maskFile = mask_file.size() > 0;

    int mask_no_start(0); pp.query("mask_no_start", mask_no_start);
    int mask_no_end(mask_no_start); pp.query("mask_no_end", mask_no_end);

    Real h_min(100.0); pp.query("h_min", h_min);
    Real f_min(1.0e-1); pp.query("f_min", f_min); 
    

    Box domainBox;
    Vector<std::string> name;
    Vector<LevelData<FArrayBox>* > data;
    Vector<DisjointBoxLayout > grids;
    Vector<int > ratio;
    int numLevels;

    Real dt ,crseDx, time;
   
    ReadAMRHierarchyHDF5(std::string(plot_file), grids, data, name , 
			 domainBox, crseDx, dt, time, ratio, numLevels);

    Vector<ProblemDomain> domain(numLevels,domainBox);
    Vector<RealVect> vdx(numLevels,RealVect::Unit*crseDx);
    Vector<Real> dx(numLevels,crseDx);
    for (int lev=1;lev<numLevels;++lev)
      {
	dx[lev] = dx[lev-1] / Real(ratio[lev-1]);
	vdx[lev] = vdx[lev-1] / Real(ratio[lev-1]);
	domain[lev] = domain[lev-1];
	domain[lev].refine(ratio[lev-1]);
      }

    Vector<LevelData<FArrayBox>* > sector_mask_fraction(numLevels);
    if (maskFile)
      {
	//load the sector mask, if it exists    
	Box mdomainBox;
	Vector<std::string> mname;
	Vector<LevelData<FArrayBox>* > mdata;
	Vector<DisjointBoxLayout > mgrids;
	Vector<int > mratio;
	int mnumLevels;
	Real mdt ,mcrseDx, mtime;
	ReadAMRHierarchyHDF5(std::string(mask_file), mgrids, mdata, mname , 
			     mdomainBox, mcrseDx, mdt, mtime, mratio, mnumLevels);

	if (mnumLevels > 1) MayDay::Error("mask files assumed to have one AMR level");

	// create storage for sector mask fractions - one FAB component for each mask label
	int n_mask_nos = 1 + mask_no_end - mask_no_start;
	
	for (int lev = 0; lev < data.size(); lev++)
	  {
	    sector_mask_fraction[lev] = new LevelData<FArrayBox>(grids[lev],n_mask_nos,IntVect::Unit);
	  }
	
	maskToAMR(sector_mask_fraction, dx, *mdata[0], mask_no_start, n_mask_nos, mcrseDx);

	//free heap
	for (int lev = 0; lev < mdata.size(); lev++)
	  {
	    if (mdata[lev]) delete mdata[lev];
	  }
	
      }
    

    
    
    stateDiagnostics(sout, append, plot_file, sector_mask_fraction, mask_no_start, mask_no_end,
		     data, grids, name, numLevels,  dx, dt, time,  ratio,
		     ice_density, water_density, gravity, h_min, f_min);


    // free heap
    for (int lev = 0; lev < data.size(); lev++)
      {
	if (sector_mask_fraction[lev]) delete sector_mask_fraction[lev];
	if (data[lev]) delete data[lev];
      }
   
    		  
  }  // end nested scope
  CH_TIMER_REPORT();

#ifdef CH_MPI
  MPI_Finalize();
#endif
  
  return 0;
}
