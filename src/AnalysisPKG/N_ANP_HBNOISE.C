//-------------------------------------------------------------------------
//   Copyright (c) 2025 Meysam Bahmanian
//   Heinz Nixdorf Institute, University of Paderborn, Germany
//
//   This file is part of the Xyce(TM) Parallel Electrical Simulator.
//
//   Xyce(TM) is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   Xyce(TM) is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with Xyce(TM).
//   If not, see <http://www.gnu.org/licenses/>.
//-------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Purpose       : This is the Harmonic Balance Noise Analysis implementation
// Special Notes :
// Creator       : Meysam Bahmanian
// Creation Date : 5/13/2025
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iomanip>

#include <fstream>

#include <N_ANP_HBNOISE.h>
#include "N_ANP_NoiseData.h"

#include <N_ANP_AnalysisManager.h>
#include <N_ANP_OutputMgrAdapter.h>
#include <N_ANP_SweepParam.h>
#include <N_ANP_SweepParamFreeFunctions.h>
#include <N_IO_CircuitBlock.h>
#include <N_IO_CmdParse.h>
#include <N_IO_OptionBlock.h>
#include <N_IO_PkgOptionsMgr.h>
#include <N_IO_SpiceSeparatedFieldTool.h>
#include <N_LAS_BlockMatrix.h>
#include <N_LAS_BlockSystemHelpers.h>
#include <N_LAS_BlockVector.h>
#include <N_LAS_Builder.h>
#include <N_LAS_Graph.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_MultiVector.h>
#include <N_LAS_System.h>
#include <N_LAS_SystemHelpers.h>
#include <N_LAS_Solver.h>
#include <N_LAS_Problem.h>
#include <N_LAS_TranSolverFactory.h>
#include <N_LOA_Loader.h>
#include <N_NLS_Manager.h>

#include <N_TIA_DataStore.h>
#include <N_TIA_fwd.h>
#include <N_TIA_StepErrorControl.h>
#include <N_TIA_WorkingIntegrationMethod.h>

#include <N_UTL_Diagnostic.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_LogStream.h>
#include <N_UTL_Math.h>
#include <N_UTL_NodeSymbols.h>
#include <N_UTL_SaveIOSState.h>
#include <N_UTL_Timer.h>

#include <N_PDS_Comm.h>
#include <N_PDS_MPI.h>
#include <N_PDS_Manager.h>
#include <N_PDS_Serial.h>

#include <N_UTL_ExtendedString.h>
#include <N_UTL_FFTInterface.hpp>
#include <N_UTL_MachDepParams.h>

#include <N_LAS_HBBuilder.h>
#include <N_LAS_HBPrecondFactory.h>
#include <N_LAS_HBSolverFactory.h>
#include <N_LAS_PrecondFactory.h>

#include <N_DEV_DeviceMgr.h>
#include <N_LOA_HBLoader.h>
#include <N_LOA_NonlinearEquationLoader.h>
#include <N_UTL_APFT.h>

#include <Teuchos_BLAS.hpp>
#include <Teuchos_Utils.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialDenseVector.hpp>
#include <Teuchos_SerialDenseHelpers.hpp>
#include <Teuchos_SerialDenseSolver.hpp>


#include <N_PDS_ParMap.h>

#include <N_TOP_Topology.h>

#include "N_ANP_HB.h"

#include <Teuchos_RCP.hpp>
using Teuchos::rcp;
using Teuchos::RCP;
using Teuchos::rcp_dynamic_cast;

// putting these here for now (from spice3's noisedefs.h file)
#define N_MINLOG          1E-38       /* the smallest number we can take the log of */
#define N_MINGAIN         1E-20  // the smallest input-output gain we can tolerate
                                 // (to calculate input-referred noise we divide
                                 // the output noise by the gain)

#define N_INTFTHRESH   1E-10       // the largest slope (of a log-log noise spectral
                                   //    density vs. freq plot) at which the noise
                                   //    spectum is still considered flat. (no need for
                                   //    log curve fitting)
#define N_INTUSELOG      1E-10       // decides which expression to use for the integral of
                                     //    x**k.  If k is -1, then we must use a 'ln' form.
                                     //    Otherwise, we use a 'power' form.  This
                                     //    parameter is the region around (k=) -1 for which
                                     //    use the 'ln' form.

namespace Xyce {
namespace Analysis {

//-----------------------------------------------------------------------------
// Function      : HBNOISE::setLinSol
// Purpose       : Save linear solver options
// Special Notes :
// Scope         : public
// Creator       : Meysam Bahmanian
// Creation Date : 5/13/2025
//-----------------------------------------------------------------------------
bool HBNOISE::setLinSol(const Util::OptionBlock &option_block)
{
  linSolOptionBlock_ = option_block;
  return true;
}


//-----------------------------------------------------------------------------
// Function      : HBNOISE::convertDataToSweepParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Meysam Bahmanian
// Creation Date : 5/14/2025
//-----------------------------------------------------------------------------
bool HBNOISE::convertDataToSweepParams()
{
  return convertData(hbnoiseSweepVector_, dataNamesMap_, dataTablesMap_);
}

//-----------------------------------------------------------------------------
// Function      : NOISE::setDataStatements
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 9/5/18
//-----------------------------------------------------------------------------
bool HBNOISE::setDataStatements(const Util::OptionBlock & paramsBlock)
{
  return processDataStatements(paramsBlock, dataNamesMap_, dataTablesMap_);
}

//-----------------------------------------------------------------------------
// Function      : HBNOISE::HBNOISE
// Purpose       : Constructor
// Special Notes :
// Scope         : public
// Creator       : Meysam Bahmanian
// Creation Date : 5/13/2025
//-----------------------------------------------------------------------------
HBNOISE::HBNOISE(
  AnalysisManager &                     analysis_manager,
  Linear::System &                      linear_system,
  Nonlinear::Manager &                  nonlinear_manager,
  Loader::Loader &                      loader,
  Device::DeviceMgr &                   device_manager,
  Topo::Topology &                      topology)
  : AnalysisBase(analysis_manager, "HBNOISE"),
    StepEventListener(&analysis_manager),
    analysisManager_(analysis_manager),
    loader_(loader),
    linearSystem_(linear_system),
    nonlinearManager_(nonlinear_manager),
    deviceManager_(device_manager),
    topology_(topology),
    outputManagerAdapter_(analysis_manager.getOutputManagerAdapter()),
    pdsMgrPtr_(0),
    currentAnalysisObject_(0),
    hbLoaderPtr_(0),
    hbBuilderPtr_(0),
    builderPtr_(0),
    hbLinearSystem_(0),
    outputNodeSingle_(true),
    outputNode1_(""),
    outputNode2_(""),
    harmonicNumber_(1.0),
    type_("DEC"),
    np_(10.0),
    fOffsetStart_(1.0),
    fOffsetStop_(1.0),
    stepMult_(0.0),
    fstep_(0.0),
    pts_per_summary_(0),
    dataSpecification_(false),
    hbnoiseLoopSize_(0),
    hbAnalysis_(0),
    bNoiseVecPtr(linearSystem_.builder().createVector()),
    // bNoiseVecRealPtr(linearSystem_.builder().createVector()),
    // bNoiseVecImagPtr(linearSystem_.builder().createVector()),
    calcNoiseIntegrals_(true),
    totalAMNoiseDens_(0.0),
    totalPMNoiseDens_(0.0),
    freq_(0.0),
    size_(0),
    period_(0.0)
{
  bNoiseVecPtr->putScalar(0.0);
  // bNoiseVecRealPtr->putScalar(0.0);
  // bNoiseVecImagPtr->putScalar(0.0);

  outputManagerAdapter_.setDotHBNOISESpecified(true);
  
  pdsMgrPtr_ = analysisManager_.getPDSManager();

  // noiseDataVec holds the total noise results and the noise
  // results for each device
  int numNoiseDevices = loader_.getNumNoiseDevices();
  noiseDataVecI_.resize(numNoiseDevices);
  for (int i=0;i<numNoiseDevices;++i)
  {
    noiseDataVecI_[i] = new NoiseData();
  }

  // Note: setting up the noise sources and putting the relevant entries
  // into the symbol table must be done after the devices are created,
  // but before the creation of the DNI() and DNO() operators.

  // set up the noise sources in each device
  loader_.setupNoiseSources(noiseDataVecI_);

  // Put an entry for each device (with noise source(s)) into the
  // symbol table owned by Topology.  This will be used during
  // operator creation for DNI() and DNO().  DNO or DNI operators come
  // in two forms, DNO(deviceName) or DNO(deviceName,noiseSource)
  Util::SymbolTable& symbol_table = topology_.getNodeSymbols();
  std::multimap<std::string,int> noiseNamesMap;
  for (int i=0;i<numNoiseDevices;++i)
  {
    // noiseDataVec_[i]->deviceName is the individual device name (e.g., Q1)
    addSymbol(symbol_table, Util::NOISE_DEVICE_SYMBOL, i, noiseDataVecI_[i]->deviceName +"_ND");

    // Account for duplicate names in noiseDataVec_[i]->noiseNames, which happens with
    // some of the ADMS device, by placing them into a multimap first before adding them
    // to the symbol table.
    for (int j=0;j<noiseDataVecI_[i]->noiseNames.size();++j)
    {
      std::string prefix = "noise_" + noiseDataVecI_[i]->deviceName;
      if (prefix == noiseDataVecI_[i]->noiseNames[j])
      {
        // noiseDataVec_[i]->noiseNames[j] are the names of the noise types (e.g., rc, rb
        // re, ic, ib and fn for a Q device).  Don't add entries if noiseName[j] is equal
        // to the string "noise_ + deviceName" (e.g., for R devices).  Those entries are
        // superfluous since (for example) DNO(R1,R1) doesn't work by design.  Only DNO(R1)
        // works.  This block should be changed from a "no op" if that design decision
        // changes.
      }
      else if (prefix.length() < noiseDataVecI_[i]->noiseNames[j].length())
      {
        // For more complex devices, noiseNames[j] will be (for example) noise_Q1_RC .
        // For some ADMS device, the noise type (RC in this example) may have inconvenient
        // characters like ( or ) in it.  The ExtendedString method removeBadChars()
        // will remove them before insertion into noiseNamesMap.
        ExtendedString noiseType(noiseDataVecI_[i]->noiseNames[j].substr(prefix.length()+1));
	std::string noiseName = prefix + "_" + noiseType.removeBadChars();
        noiseNamesMap.insert(std::pair<std::string,int>(noiseName, j));
      }
    }

    std::string prevName="";
    int suffix = 0;
    for (std::multimap<std::string,int>::iterator it=noiseNamesMap.begin(); it!=noiseNamesMap.end(); ++it)
    {
      // For ADMS devices, that may have duplicate entries for a given noise type, the entries
      // are "suffixed" with _0, _1, _2, etc.  If there are no duplicate entries (e.g., for
      // the Q device) then just the _0 suffix is used.
      (*it).first != prevName ? suffix=0 : ++suffix;
      std::ostringstream s;
      s << suffix;
      addSymbol(symbol_table, Util::NOISE_TYPE_SYMBOL, (*it).second, (*it).first + "_" + s.str());
      prevName = (*it).first;
    }

    noiseNamesMap.clear();
  }
}

//-----------------------------------------------------------------------------
// Function      : HBNOISE::~HBNOISE
// Purpose       : Destructor
// Special Notes :
// Scope         : public
// Creator       : Meysam Bahmanian
// Creation Date : 5/13/2025
//-----------------------------------------------------------------------------
HBNOISE::~HBNOISE()
{
  for (size_t i = 0; i < noiseDataVecI_.size(); ++i) {
    delete noiseDataVecI_[i];
  }
  noiseDataVecI_.clear();
}

//-----------------------------------------------------------------------------
// Function      : HBNOISE::getHBAnalysis
// Purpose       : Get the HB analysis object from the analysis manager
// Special Notes :
// Scope         : private
// Creator       : Meysam Bahmanian
// Creation Date : 5/16/2025
//-----------------------------------------------------------------------------
Analysis::HB* HBNOISE::getHBAnalysis()
{
  std::vector<ProcessorBase *>& analyses = analysisManager_.getAnalysisVector();
  for (std::vector<ProcessorBase *>::const_iterator it = analyses.begin(); it != analyses.end(); ++it) {
    if (Analysis::HB* hb = dynamic_cast<Analysis::HB*>(*it)) {
      return hb;
    }
  }
  return nullptr;
}

//-----------------------------------------------------------------------------
// Function      : HBNOISE::notify
// Purpose       : Handle step events
// Special Notes :
// Scope         : public
// Creator       : Meysam Bahmanian
// Creation Date : 5/13/2025
//-----------------------------------------------------------------------------
void HBNOISE::notify(const StepEvent &event) 
{
  if (event.state_ == StepEvent::STEP_STARTED)
  {
    AnalysisBase::resetForStepAnalysis();
  }
}

//-----------------------------------------------------------------------------
// Function      : HBNOISE::doInit
// Purpose       : Initialize the analysis
// Special Notes :
// Scope         : protected
// Creator       : Meysam Bahmanian
// Creation Date : 5/13/2025
//-----------------------------------------------------------------------------
bool HBNOISE::doInit()
{
  bool bsuccess = true;

  // get the HB analysis object first
  hbAnalysis_ = getHBAnalysis();
  if (!hbAnalysis_)
  {
    Report::UserError0() << "HBNOISE analysis requires an HB analysis to be defined first";
    return false;
  }
  hbAnalysis_->hbNoise_ = true;
  // check if the HB analysis has a single frequency
  std::vector<double> freqs = hbAnalysis_->freqs_;
  if (freqs.size() > 1)
  {
    Report::UserError0() << "HBNOISE analysis requires a single frequency to be specified";
    return false;
  }
  // get the frequency from the HB analysis
  size_ = hbAnalysis_->size_;
  numHarms_ = (size_-1)/2;
  freq_ = freqs[0];
  omega_ = 2.0*M_PI*freq_;
  period_ = 1.0/freq_;
  times_.resize(size_);
  for( int i = 0; i < size_; ++i )
    times_[i] = i*period_/(size_-1);

  // TODO: check maximum offset frequency to be below the Nyquist frequency

  // Noise data vector for each harmonic
  // int numHarms = (size_-1)/2; reminder
  // for baseband we have 1 vector, for harmonics we have 2 vectors (in-phase and quadrature) for each harmonic
  int numNoiseDevices = noiseDataVecI_.size();

  // noiseDataVecI_ is already set up in the constructor
  noiseDataVecQ_.resize(numNoiseDevices);
  for (int i=0;i<numNoiseDevices;++i)
  {
    noiseDataVecQ_[i] = new NoiseData();
  }
  loader_.setupNoiseSources(noiseDataVecQ_);

  // set up the noise sources for each harmonic
  noiseDataVecVecI_.resize(size_);
  noiseDataVecVecQ_.resize(size_);
  for (int i = 0; i < size_; ++i) {
    noiseDataVecVecI_[i].resize(numNoiseDevices);
    noiseDataVecVecQ_[i].resize(numNoiseDevices);
    for (int j = 0; j < numNoiseDevices; ++j) {
      noiseDataVecVecI_[i][j] = new Analysis::NoiseData();
      noiseDataVecVecQ_[i][j] = new Analysis::NoiseData();
    }
    loader_.setupNoiseSources(noiseDataVecVecI_[i]);
    loader_.setupNoiseSources(noiseDataVecVecQ_[i]);
  }

  // check if the "DATA" specification was used.  If so, create a new vector of
  // SweepParams, in the "TABLE" style.
  if (dataSpecification_)
  {
    if (!convertDataToSweepParams())
    {
      Report::UserFatal() << "Invalid data=<name> parameter on .HBNOISE line.";
      return false;
    }

    std::vector<SweepParam>::iterator begin = hbnoiseSweepVector_.begin();
    std::vector<SweepParam>::iterator end = hbnoiseSweepVector_.end();
    std::vector<SweepParam>::iterator it = begin;
    for ( ; it != end; ++it)
    {
      SweepParam &sweep_param = (*it);
      std::string name = (*it).name; Util::toUpper(name);
      if (name == "FREQ" || name == "HERTZ")
      {
        // used to check that the specified frequencies are monotonically
        // increasing, to determine whether the noise integrals can be
        // calculated when DATA=<name> is used on the .HBNOISE line
        double prevFreq=-1.0;

        // frequency values for .HBNOISE must be > 0
        for (int i=0; i<(*it).valList.size(); ++i)
	{
          if ( (*it).valList[i] <= 0 )
	  {
            Report::UserFatal() << "Frequency values in .DATA for .HBNOISE analysis must be > 0";
            return false;
          }
          if ( calcNoiseIntegrals_ && ((*it).valList[i] <= prevFreq) )
	  {
            calcNoiseIntegrals_ = false;
	    Report::UserWarning0() << "Total Noise Integrals will not be calculated, "
		<< "since frequencies in .DATA table are not monotonically increasing";
          }
          prevFreq = (*it).valList[i];
        }
      }
      else
      {
        loader_.getParamAndReduce(analysisManager_.getComm(), sweep_param.name);
      }
    }

    // now set up the looping, etc
    hbnoiseLoopSize_ = setSweepLoopVals(begin, end);
  }
  else
  {
    hbnoiseLoopSize_ = setupSweepParam_();
  }

  processOutputNodes ();

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : HBNOISE::doRun
// Purpose       : Run the HBNOISE analysis
// Special Notes :
// Scope         : public
// Creator       : Meysam Bahmanian
// Creation Date : 5/13/2025
//-----------------------------------------------------------------------------
bool HBNOISE::doRun()
{
  return doInit() && doLoopProcess() && doFinish();
}

//-----------------------------------------------------------------------------
// Function      : HBNOISE::doLoopProcess
// Purpose       : Process the main analysis loop
// Special Notes :
// Scope         : protected
// Creator       : Meysam Bahmanian
// Creation Date : 5/13/2025
//-----------------------------------------------------------------------------
bool HBNOISE::doLoopProcess()
{
  bool bsuccess = true;

  // analysisManager_.pushActiveAnalysis(hbAnalysis_);
  bsuccess = hbAnalysis_->run();
  analysisManager_.pushActiveAnalysis(this);
  // add if fails ...

  // save the mag & phase of solution for AM/PM noise calculations
  

  // some of these pointers are set after HB::run(), so I have to copy them at this point.
  hbBuilderPtr_ = hbAnalysis_->hbBuilderPtr_;
  builderPtr_ = &(hbAnalysis_->builder_);
  hbLoaderPtr_ = hbAnalysis_->hbLoaderPtr_;

  updateLinearTimeVariantSystem_C_and_G_();
  createHarmonicSpaceLinearSystem_();

  // static_cast<Xyce::Util::Notifier<AnalysisEvent> &>(analysisManager_).publish
  //   (AnalysisEvent(AnalysisEvent::INITIALIZE, AnalysisEvent::NOISE));

  clearNoiseIntegrals_();

  setupAdjointRHS_();

  Parallel::Manager &pds_manager = *analysisManager_.getPDSManager();
  Parallel::Communicator &comm = *(pds_manager.getPDSComm());
  int myPID = comm.procID();

  // ///////////////////////////////////////////////////////////////////////////
  // // frequency loop
  // // loop over all the specified frequencies, do an AC solve and a NOISE solve at each.
  for (int currentStep = 0; currentStep < hbnoiseLoopSize_; ++currentStep)
  {
  //   // solve the AC system, to get up-to-date currents and voltages
    if (dataSpecification_)
    {
      updateDataParams_(currentStep);
    }
    else
    {
      updateCurrentFreq_(currentStep);
    }

  //   static_cast<Xyce::Util::Notifier<AnalysisEvent> &>(analysisManager_).publish
  //     (AnalysisEvent(AnalysisEvent::STEP_STARTED, AnalysisEvent::NOISE, currentFreq_, currentStep));

  //   updateACLinearSystem_C_and_G_();
  //   updateACLinearSystemFreq_();
    updateHarmonicSpaceFreq_();
  //   updateACLinearSystemMagAndPhase_();

  //   bool stepAttemptStatus;
  //   {
  //     Stats::StatTop _ACsolveStat("AC Linear Solve");
  //     Stats::TimeBlock _AC_Timer(_ACsolveStat);

  //     stepAttemptStatus = solveACLinearSystem_();
  //   }

  //   // save a copy of X_ (the AC solution, already computed), for output purposes, etc.
  //   *saved_AC_X_ = *X_;

  //   // Compute AC gain.
  //   Linear::Vector & Xreal = X_->block( 0 );
  //   Linear::Vector & Ximag = X_->block( 1 );

  //   double v1r = 0.0;
  //   double v1i = 0.0;
  //   double v2r = 0.0;
  //   double v2i = 0.0;

  //   //comm.barrier();
  //   int root=-1;
  //   if (outputVarGIDs_.size()>0)
  //   {
  //     if (outputVarGIDs_[0] > -1)
  //     {
  //       v1r = Xreal.getElementByGlobalIndex(outputVarGIDs_[0]);
  //       v1i = Ximag.getElementByGlobalIndex(outputVarGIDs_[0]);
  //       root = myPID;
  //     }
  //     Xyce::Parallel::AllReduce(comm.comm(), MPI_MAX, &root, 1);
  //     comm.bcast( &v1r, 1, root );
  //     comm.bcast( &v1i, 1, root );
  //   }

  //   root=-1;
  //   if (outputVarGIDs_.size()>1)
  //   {
  //     if (outputVarGIDs_[1] > -1)
  //     {
  //       v2r = Xreal.getElementByGlobalIndex(outputVarGIDs_[1]);
  //       v2i = Ximag.getElementByGlobalIndex(outputVarGIDs_[1]);
  //       root = myPID;
  //     }
  //     Xyce::Parallel::AllReduce(comm.comm(), MPI_MAX, &root, 1);
  //     comm.bcast( &v2r, 1, root );
  //     comm.bcast( &v2i, 1, root );
  //   }

  //   double realVal = v1r-v2r;
  //   double imagVal = v1i-v2i;
  //   GainSqInv_ = 1.0 / std::max(((realVal*realVal) + (imagVal*imagVal)),N_MINGAIN);
  //   lnGainInv_ = std::log(GainSqInv_);

  //   // save previous (last) lnNoise densities
  //   for (int i=0;i<noiseDataVec_.size();++i)
  //   {
  //     int numNoiseThisDevice = noiseDataVec_[i]->numSources;
  //     for (int j=0;j<numNoiseThisDevice;++j)
  //     {
  //       noiseDataVec_[i]->lastLnNoiseDens[j] = noiseDataVec_[i]->lnNoiseDens[j];
  //     }
  //   }

    // do NOISE analysis for this frequency.
    resetAdjointHBNOISELinearSystem_();
    solveAdjointHBNOISE_();

  //   // Perform total noise integrals, if the specified frequency values are
  //   // monotonically increasing.  This is always true if DATA=<name> is NOT
  //   // used on the .NOISE line.
  //   if (currentStep != 0 && calcNoiseIntegrals_)
  //   {
  //     for (int i=0;i<noiseDataVec_.size();++i)
  //     {
  //       int numNoiseThisDevice = noiseDataVec_[i]->numSources;
  //       for (int j=0;j<numNoiseThisDevice;++j)
  //       {
  //         double noizDens = noiseDataVec_[i]->outputNoiseDens[j];
  //         double lnDens = noiseDataVec_[i]->lnNoiseDens[j];
  //         double lnlastDens = noiseDataVec_[i]->lastLnNoiseDens[j];

  //         double tempOutNoise = noiseIntegral( noizDens, lnDens, lnlastDens,
  //                delLnFreq_, delFreq_, lnFreq_, lnLastFreq_);

  //         double tempInNoise = noiseIntegral(
  //                noizDens * GainSqInv_,
  //                lnDens + lnGainInv_,
  //                lnlastDens + lnGainInv_,
  //                delLnFreq_, delFreq_, lnFreq_, lnLastFreq_);


  //         noiseDataVec_[i]->outputNoiseTotal[j] += tempOutNoise;
  //         noiseDataVec_[i]->inputNoiseTotal[j] += tempInNoise;

  //         totalOutputNoise_+=tempOutNoise;
  //         totalInputNoise_+=tempInNoise;
  //       }
  //     }
  //   }

  //   // process success/failure
  //   if (stepAttemptStatus)
  //   {
  //     static_cast<Xyce::Util::Notifier<AnalysisEvent> &>(analysisManager_).publish
  //       (AnalysisEvent(AnalysisEvent::STEP_SUCCESSFUL, AnalysisEvent::NOISE, currentFreq_, currentStep));
  //     doProcessSuccessfulStep();
  //   }
  //   else // stepAttemptStatus  (ie do this if the step FAILED)
  //   {
  //     static_cast<Xyce::Util::Notifier<AnalysisEvent> &>(analysisManager_).publish
  //       (AnalysisEvent(AnalysisEvent::STEP_FAILED, AnalysisEvent::NOISE, currentFreq_, currentStep));
  //     doProcessFailedStep();
  //   }
      doProcessSuccessfulStep();
  }

  // Xyce::Parallel::AllReduce(comm.comm(), MPI_SUM, &totalOutputNoise_, 1);
  // Xyce::Parallel::AllReduce(comm.comm(), MPI_SUM, &totalInputNoise_, 1);

  // if (calcNoiseIntegrals_)
  // {
  //   // Outputs to the screen
  //   noiseOutputToScreen_( Xyce::lout() );
  // }

  // static_cast<Xyce::Util::Notifier<AnalysisEvent> &>(analysisManager_).publish
  //   (AnalysisEvent(AnalysisEvent::FINISH, AnalysisEvent::NOISE));

  return true;
}

//-----------------------------------------------------------------------------
// Function      : NOISE::clearNoiseIntegrals
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Meysam Bahmanian
// Creation Date : 6/5/2025
//-----------------------------------------------------------------------------
void HBNOISE::clearNoiseIntegrals_() 
{
  int numNoiseDevices = loader_.getNumNoiseDevices();
  // clear out the integral arrays

  for (int i=0;i<numNoiseDevices;++i)
  {
    int numNoiseThisDevice = noiseDataVecI_[i]->numSources;
    for (int j=0;j<numNoiseThisDevice;++j)
    {
      noiseDataVecI_[i]->inputNoiseTotal[j] = 0.0;
      noiseDataVecI_[i]->outputNoiseTotal[j] = 0.0;
    }
  }

  for (int i = 0; i < size_; ++i) {
    for (int j = 0; j < numNoiseDevices; ++j) {
      int numNoiseThisDevice = noiseDataVecVecI_[i][j]->numSources;
      for (int k = 0; k < numNoiseThisDevice; ++k) {
        noiseDataVecVecI_[i][j]->inputNoiseTotal[k] = 0.0;
        noiseDataVecVecI_[i][j]->outputNoiseTotal[k] = 0.0;
      }
    }
  }

  if (harmonicNumber_ != 0)
  {
    for (int i=0;i<numNoiseDevices;++i)
    {
      int numNoiseThisDevice = noiseDataVecQ_[i]->numSources;
      for (int j=0;j<numNoiseThisDevice;++j)
      {
        noiseDataVecQ_[i]->inputNoiseTotal[j] = 0.0;
        noiseDataVecQ_[i]->outputNoiseTotal[j] = 0.0;
      }
    }
  }

  for (int i = 0; i < size_; ++i) {
    for (int j = 0; j < numNoiseDevices; ++j) {
      int numNoiseThisDevice = noiseDataVecVecQ_[i][j]->numSources;
      for (int k = 0; k < numNoiseThisDevice; ++k) {
        noiseDataVecVecQ_[i][j]->inputNoiseTotal[k] = 0.0;
        noiseDataVecVecQ_[i][j]->outputNoiseTotal[k] = 0.0;
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : NOISE::setupAdjointRHS_
// Purpose       : Sets up stuff that only needs to be set up once,
// Special Notes :
// Scope         : public
// Creator       : Meysam Bahmanian
// Creation Date : 5/31/2025
// Originally written by Eric Keiter for the NOISE analysis class
//-----------------------------------------------------------------------------
// This function (copied from the NOISE analysis class) is still usable for HBNOISE. 
// The GIDs are set properly in the bNoiseVecRealPtr vector.
// We need to place this vector properly in right block in the RHS, 
// taking into account the harmonic number and in-phase/quadrature components
void HBNOISE::setupAdjointRHS_()
{
  Parallel::Manager &pds_manager = *analysisManager_.getPDSManager();
  Parallel::Communicator &comm = *(pds_manager.getPDSComm());

  bNoiseVecPtr->putScalar(0.0);
  // bNoiseVecRealPtr->putScalar(0.0);
  // bNoiseVecImagPtr->putScalar(0.0);

  int numOutVars = outputVarNames_.size();
  for (int iout=0;iout<numOutVars;++iout)
  {
    int tmpGID=outputVarGIDs_[iout];
    if (tmpGID > -1)
    {
      double val=1.0;
      if (iout>0) val=-1.0;
      bNoiseVecPtr->setElementByGlobalIndex( tmpGID, val, 0);
    }
  }
  bNoiseVecPtr->fillComplete();
  if (DEBUG_ANALYSIS)
  {
    Xyce::dout() << "bNoiseVecPtr:" << std::endl;
    bNoiseVecPtr->print(Xyce::dout());
    Xyce::dout() << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : NOISE::processOutputNodes
// Purpose       : determines the GIDs for the nodes specified in the first argument
//                 of the .NOISE line.
// Special Notes :
// Scope         : public
// Creator       : Meysam Bahmanian
// Creation Date : 5/31/2025
// Originally written by Eric Keiter for the NOISE analysis class
//-----------------------------------------------------------------------------
void HBNOISE::processOutputNodes ()
{
  // setup the names:
  outputVarNames_.clear();
  outputVarNames_.push_back(outputNode1_);
  if (!outputNodeSingle_)
  {
    outputVarNames_.push_back(outputNode2_);
  }
  int numOutVars = outputVarNames_.size();

  // set up the gid's:
  int found(0);
  int found2(0);
  bool foundLocal(false);
  bool foundLocal2(false);

  Parallel::Manager &pds_manager = *analysisManager_.getPDSManager();
  Parallel::Communicator &comm = *(pds_manager.getPDSComm());

  outputVarGIDs_.resize( numOutVars, -1 );
  for (int iout = 0; iout < numOutVars; ++iout)
  {
    std::vector<int> svGIDList1, dummyList;
    char type1;
    foundLocal = topology_.getNodeSVarGIDs(NodeID(outputVarNames_[iout], Xyce::_VNODE),
        svGIDList1, dummyList, type1);

    found = static_cast<int>(foundLocal);
    Xyce::Parallel::AllReduce(comm.comm(), MPI_LOR, &found, 1);

    foundLocal2 = false;
    if (!found)// if looking for this as a voltage node failed, try a "device" (i.e. current) node.
    {
      foundLocal2 = topology_.getNodeSVarGIDs(NodeID(outputVarNames_[iout], Xyce::_DNODE),
          svGIDList1, dummyList, type1);
    }
    found2 = static_cast<int>(foundLocal2);
    Xyce::Parallel::AllReduce(comm.comm(), MPI_LOR, &found2, 1);

    if (!found && !found2)
    {
      Report::UserError() << "Output function variable " << outputVarNames_[iout] << " not found";
    }

    if (found || found2)
    {
      int tmpGID=-1;
      if(svGIDList1.size()==1)
      {
        tmpGID = svGIDList1.front();
      }
      outputVarGIDs_[iout] = tmpGID;
    }
  }
}


//-----------------------------------------------------------------------------
// Function      : HBNOISE::updateHarmonicSpaceFreq_()
// Purpose       :
// Special Notes :
//
// Scope         : public
// Creator       : Meysam Bahmanian
// Creation Date : 6/31/2025
//-----------------------------------------------------------------------------
bool HBNOISE::updateHarmonicSpaceFreq_()
{
  // First diagonal block
  harmonicSpaceMatrix_->put( 0.0 ); // Zero out whole matrix
  harmonicSpaceMatrix_->add(*harmonicSpaceMatrixConstant_);

  double omega =  2.0 * M_PI * currentFreq_;

  harmonicSpaceMatrix_omegamC_2_->put( 0.0 );
  harmonicSpaceMatrix_omegamC_2_->add(*harmonicSpaceMatrix_C_2_);
  harmonicSpaceMatrix_omegamC_2_->scale(omega);

  harmonicSpaceMatrix_->add(*harmonicSpaceMatrix_omegamC_2_);

  harmonicSpaceMatrix_->assembleGlobalMatrix();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : HBNOISE::resetAdjointHBNOISELinearSystem_
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Meysam Bahmanian
// Creation Date : 6/3/2025
//-----------------------------------------------------------------------------
// Now we need to be careful here! We have to set up the adjoint RHS properly.
// our harmonic space matrix is cosine absolute-referenced.
// What does it mean? The harmonic space matrix has no information about the phase of the desired carrier and its harmonic number.
// However, we want to find the sensitivities (in this context the transfer functions) with respect to an output at a given harmonic with the phase provided by the solution.
// This means, we need to properly project the adjoint RHS to the solution vector.
// Unfortunately, I have to again write some math so we can follow what is going on.
// The HB solution at the desired output node has the form of:
// a = A_I * cos(k*wc*t) - A_Q * sin(k*wc*t)
// The solution ITSELF should be considered as in-phase. The coefficients a_I and a_Q appear beacuse our mathematical approach is cosine absolute-referenced.
// Now we reconstruct the solution as
// a = A * cos(k*wc*t + phi)
// where :
// A = sqrt(A_I^2 + A_Q^2)
// cos(phi) = A_I/A
// sin(phi) = A_Q/A
// What does this mean for the Adjoint RHS? Let's look at the form of our desired output-referenced (not absolute-referenced) "objective functions" for the adjoint solve:
// J = x_I * cos(k*wc*t+phi) - x_Q * sin(k*wc*t+phi)
// We want to find the sensitivities of x_I and x_Q with respect to the noise sources.
// So we need to solve two adjoint problems, one for the in-phase component x_I, and one for the quadrature component x_Q.
// Now we need to derive the RHS for each case and translate J to the cosine absolute-referenced form:
// J = (x_I * cos(phi) - x_Q * sin(phi)) * cos(k*wc*t) - (x_I * sin(phi) + x_Q * cos(phi)) * sin(k*wc*t)
// Therefore, for the in-phase component, we need to solve the following adjoint problem:
// J_I = + x_I * (A_I/A) * cos(k*wc*t) - x_I * (A_Q/A) * sin(k*wc*t)
// and for the quadrature component, we need to solve the following adjoint problem:
// J_Q = - x_Q * (A_Q/A) * cos(k*wc*t) - x_Q * (A_I/A) * sin(k*wc*t)
// We clearly see that J_I and J_Q both have in-phase and quadrature components, since we are now referencing our objective functions to the desired output.

void HBNOISE::resetAdjointHBNOISELinearSystem_()
{
  // clear out the harmonicSpaceX_ vector used by the linear solver.
  harmonicSpaceXI_->putScalar(0.0);
  harmonicSpaceXQ_->putScalar(0.0);

  // setup the harmonicSpaceB_ vector RHS for the adjoint solve:
  harmonicSpaceBI_->putScalar( 0.0 );
  harmonicSpaceBQ_->putScalar( 0.0 );
  if (harmonicNumber_ == 0)
  { // baseband noise
    // baseband has just in-phase component
    harmonicSpaceBI_->block( 0 ).update( 1.0, *bNoiseVecPtr);
  }
  else
  {
    // in-phase
    // RHS is positive
    harmonicSpaceBI_->block( 4*harmonicNumber_ - 2 ).update( +outputValCosPhi_, *bNoiseVecPtr);
    // RHS is negative
    harmonicSpaceBI_->block( 4*harmonicNumber_     ).update( +outputValSinPhi_, *bNoiseVecPtr);

    // quadrature
    // RHS is positive
    harmonicSpaceBQ_->block( 4*harmonicNumber_ - 2 ).update( -outputValSinPhi_, *bNoiseVecPtr);
    // RHS is negative
    harmonicSpaceBQ_->block( 4*harmonicNumber_     ).update( +outputValCosPhi_, *bNoiseVecPtr);
  }
  if (DEBUG_ANALYSIS)
  {
    Xyce::dout()<<"adjoint noise B vector (in-phase):"<<std::endl;
    harmonicSpaceBI_->print(Xyce::dout());
    Xyce::dout()<<"adjoint noise B vector (quadrature):"<<std::endl;
    harmonicSpaceBQ_->print(Xyce::dout());
  }

  if (DEBUG_ANALYSIS)
  {
    Xyce::dout()<<"harmonic space matrix:"<<std::endl;
    harmonicSpaceMatrix_->print(Xyce::dout());
  }
}

//-----------------------------------------------------------------------------
// Function      : HBNOISE::solveAdjointHBNOISE_
// Purpose       : Solves for HBNOISE using the adjoint method.
// Special Notes : This is the production solve.
// Scope         : public
// Creator       : Meysam Bahmanian
// Creation Date : 6/04/2025
//-----------------------------------------------------------------------------
bool HBNOISE::solveAdjointHBNOISE_()
{
  bool bsuccess = true;

  // always solve for in-phase
  int linearStatus = blockSolverI_->solveTranspose();
  if (linearStatus != 0)
  {
    Xyce::dout() << "Linear solve for in-phase exited with error: " << linearStatus << std::endl;
    bsuccess = false;
  }

  if (harmonicNumber_ != 0) // if not baseband, solve for quadrature
  {
  linearStatus = blockSolverQ_->solveTranspose();
  if (linearStatus != 0)
  {
      Xyce::dout() << "Linear solve for quadrature exited with error: " << linearStatus << std::endl;
      bsuccess = false;
    }
  }

  Xyce::dout() << "Preparing output vectors for in-phase adjoint solve: " << std::endl;
  prepareHBNOISEOutputVectors_(harmonicSpaceXI_, noiseDataVecVecI_, noiseDataVecI_, totalAMNoiseDens_);

  if (harmonicNumber_ != 0) {
    Xyce::dout() << "Preparing output vectors for quadrature adjoint solve: " << std::endl;
    prepareHBNOISEOutputVectors_(harmonicSpaceXQ_, noiseDataVecVecQ_, noiseDataVecQ_, totalPMNoiseDens_);
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : HBNOISE::prepareHBNOISEOutputVectors_
// Purpose       : Prepare the output vectors for the adjoint solve.
//                 This is put to a separate function to be used for both in-phase and quadrature adjoint solves.
// Special Notes :
// Scope         : private
// Creator       : Meysam Bahmanian
// Creation Date : 6/5/2025
//-----------------------------------------------------------------------------
// In this function we face a new mathematical challenge.
// First, we need to convert upper and lower sideband noise PSDs to in-phase and quadrature components.
// Howevenr, if the noise PSD for LSB/USB is not flat, the IQ noise PSDs will be correlated. These correlation
// factors have to be taken into account for the output noise vectors.
// Let's say we have a band-limited noise n(t) in the neighborhood of the carrier frequency w, with PSD for LSB and USB as follows:
// n(t) = n_L(t) + n_U(t)
// where n_L(t) has no frequency components above w, and n_U(t) has no frequency components below w.
// Let's say we have a signal 
// s(t) = cos(wc*t + wm*t)
// which enters the system with gains of
// for AC we have the source b and the output v as
// b = b_C * cos(wc*t) - b_S * sin(wc*t)
// v = v_C * cos(wc*t) - v_S * sin(wc*t)
// the relation between b and v is:
// v = b_C * H_R * cos(wc*t) - b_S * H_R * sin(wc*t) - b_C * H_I * sin(wc*t) - b_S * H_I * cos(wc*t)
// so we have:
// v_C = + b_C * H_R - b_S * H_I
// v_S = + b_C * H_I + b_S * H_R
// in adjoint formulation only v_C is used:
// H_R = +dv_C/db_C
// H_I = -dv_C/db_S
// Now for the noise sources, we have:

// form of b_ in harmonic space:
// b = b_C0 * cos(wm*t) - b_S0 * sin(wm*t) + b_CLk * cos(k*wc*t-wm*t) - b_SLk * sin(k*wc*t-wm*t) + b_CUk * cos(k*wc*t+wm*t) - b_SUk * sin(k*wc*t+wm*t)
// C/S is cosine/since and L/U is lower/upper sideband

void HBNOISE::prepareHBNOISEOutputVectors_(
  Linear::BlockVector *           harmonicSpaceX,
  std::vector<std::vector<Xyce::Analysis::NoiseData*> > &noiseDataVecVec,
  std::vector<Xyce::Analysis::NoiseData*> &noiseDataVec,
  double &totalRelativeNoiseDens
  )
{
  double omega =  2.0 * M_PI * currentFreq_;

  int numHarms = (size_-1)/2;
  int numNoiseDevices = loader_.getNumNoiseDevices();

  // baseband
  for (int i=0;i<numNoiseDevices;++i)
  {
    (noiseDataVecVec[0][i])->omega = omega;
    (noiseDataVecVec[0][i])->freq = currentFreq_;
  }
  loader_.getNoiseSources(noiseDataVecVec[0]);

  // harmonics
  // the noise data vector for the harmonics: The device noise PSDs are loaded for the upper and lower sidebands around the harmonics,
  for (int i = 1; i <= numHarms; ++i) {
    for (int j = 0; j < numNoiseDevices; ++j) {
      // lower sideband
      (noiseDataVecVec[2*i-1][j])->omega = i*omega_ - omega;
      (noiseDataVecVec[2*i-1][j])->freq = i*freq_ - currentFreq_;

      // upper sideband
      (noiseDataVecVec[2*i  ][j])->omega = i*omega_ + omega;
      (noiseDataVecVec[2*i  ][j])->freq = i*freq_ + currentFreq_;
    }
    loader_.getNoiseSources(noiseDataVecVec[2*i-1]);
    loader_.getNoiseSources(noiseDataVecVec[2*i  ]);
  }

  // for baseband we had just two vectors real/imag. 
  // for harmonic space we have 2 vectors for the baseband and 4 vectors for each harmonic
  // storage for the solve vectors
  std::vector< Teuchos::RCP<Linear::Vector> > outputVectors;
  for (int i = 0; i < 4*numHarms+2; ++i) 
  {
    outputVectors.push_back( rcp(linearSystem_.builder().createVector()) );
  }

  copyFromBlockVector( *harmonicSpaceX, outputVectors);

  if (DEBUG_ANALYSIS)
  {
    // baseband components
    Xyce::dout() << "d(o)/d(b_ real) at baseband adjoint solve:" << std::endl;
    outputVectors[0]->print( Xyce::dout() );
    Xyce::dout() << "d(o)/d(b_ imag) at baseband adjoint solve:" << std::endl;
    outputVectors[1]->print( Xyce::dout() );

    // harmonics components
    for (int i = 1; i <= numHarms; ++i) 
    {
      Xyce::dout() << "d(o)/d(b_ LSB real) at harmonic " << i << " adjoint solve:" << std::endl;
      outputVectors[4*i-2]->print( Xyce::dout() );
      Xyce::dout() << "d(o)/d(b_ LSB imag) at harmonic " << i << " adjoint solve:" << std::endl;
      outputVectors[4*i-1]->print( Xyce::dout() );
      Xyce::dout() << "d(o)/d(b_ USB real) at harmonic " << i << " adjoint solve:" << std::endl;
      outputVectors[4*i  ]->print( Xyce::dout() );
      Xyce::dout() << "d(o)/d(b_ USB imag) at harmonic " << i << " adjoint solve:" << std::endl;
      outputVectors[4*i+1]->print( Xyce::dout() );

      Xyce::dout() << std::endl;
    }
  }

  Parallel::Manager &pds_manager = *analysisManager_.getPDSManager();
  Parallel::Communicator &comm = *(pds_manager.getPDSComm());

  for (int i = 0; i < numNoiseDevices; ++i) 
  { // select a noise device
    for (int k = 0; k < 2*numHarms+1; ++k) 
    { // sweeping over noise harmonics
      evalDeviceNoiseDensities(*(noiseDataVecVec[k][i]), *(outputVectors[2*k]), *(outputVectors[2*k+1]));
    }
  }

  // now we sum over the noise sidebands of all harmonics and store the results in the noiseDataVec
  totalRelativeNoiseDens = 0.0;
  for (int i = 0; i < numNoiseDevices; ++i) 
  {
    noiseDataVec[i]->totalOutputNoise = 0.0;
    noiseDataVec[i]->totalRelativeNoiseDens = 0.0;
    for (int j = 0; j < noiseDataVecVec[0][i]->numSources; ++j) 
    {
      for (int k = 0; k < 2*numHarms+1; ++k) 
      { // sum of all noise sidebands for each noise source
        noiseDataVec[i]->outputNoiseDens[j] += noiseDataVecVec[k][i]->outputNoiseDens[j];
      }
      noiseDataVec[i]->relativeNoiseDens[j] = noiseDataVec[i]->outputNoiseDens[j] / outputValSqr_;
    }

    for (int k = 0; k < 2*numHarms+1; ++k) 
    { // sum of all noise sidebands for total device noise
      noiseDataVec[i]->totalOutputNoise += noiseDataVecVec[k][i]->totalOutputNoise;
    }
    noiseDataVec[i]->totalRelativeNoiseDens = noiseDataVec[i]->totalOutputNoise / outputValSqr_;
    totalRelativeNoiseDens += noiseDataVec[i]->totalRelativeNoiseDens;
  }

  // reduce the total noise density over all processors
  Xyce::Parallel::AllReduce(comm.comm(), MPI_SUM, &totalRelativeNoiseDens, 1);

  // if (comm.isSerial() )
  // {
  //   // FIX:  replace this output call!
  //   hackTecplotOutput();
  // }
}

//-----------------------------------------------------------------------------
// Function      : HBNOISE::evalDeviceNoiseDensities
// Purpose       : Evaluates the noise densities for a given device
// Special Notes : This function is used to evaluate the noise densities for a given device
// Scope         : This function is used to evaluate the noise densities for a given device
// Creator       : Meysam Bahmanian
// Creation Date : 6/7/2025
//-----------------------------------------------------------------------------
inline void HBNOISE::evalDeviceNoiseDensities(
  Xyce::Analysis::NoiseData&    noiseData, 
  Linear::Vector&               Xreal, 
  Linear::Vector&               Ximag) 
  {
  int numNoiseThisDevice = noiseData.numSources;

  noiseData.totalNoise = 0.0;                  // not used
  noiseData.totalOutputNoise = 0.0;            // total output noise density of the device
  noiseData.totalRelativeNoiseDens = 0.0;      // total AM or PM noise density of the device
  for (int j=0;j<numNoiseThisDevice;++j)
  {
    int li_Pos = noiseData.li_Pos[j];
    int li_Neg = noiseData.li_Neg[j];
    int li_PosCorl = noiseData.li_PosCorl[j];
    int li_NegCorl = noiseData.li_NegCorl[j];

    // if two sets of nodes available, calculate correlated gain.
    // otherwise, calculate uncorrelated gain.
    if ((li_PosCorl != -1) && (li_NegCorl != -1))
    {
      double realVal1 = ((li_Pos!=-1)?Xreal[li_Pos]:0) - ((li_Neg!=-1)?Xreal[li_Neg]:0);
      double imagVal1 = ((li_Pos!=-1)?Ximag[li_Pos]:0) - ((li_Neg!=-1)?Ximag[li_Neg]:0);
      double realVal2 = ((li_PosCorl!=-1)?Xreal[li_PosCorl]:0) - ((li_NegCorl!=-1)?Xreal[li_NegCorl]:0);
      double imagVal2 = ((li_PosCorl!=-1)?Ximag[li_PosCorl]:0) - ((li_NegCorl!=-1)?Ximag[li_NegCorl]:0);
      double realOut = noiseData.T0 * realVal1 + noiseData.T2 * realVal2 - noiseData.T3 * imagVal2;
      double imagOut = noiseData.T0 * imagVal1 + noiseData.T2 * imagVal2 + noiseData.T3 * realVal2;
      noiseData.gainSqr[j] = (realOut*realOut) + (imagOut*imagOut);
    }
    else
    {
      double realVal = ((li_Pos!=-1)?Xreal[li_Pos]:0) - ((li_Neg!=-1)?Xreal[li_Neg]:0);
      double imagVal = ((li_Pos!=-1)?Ximag[li_Pos]:0) - ((li_Neg!=-1)?Ximag[li_Neg]:0);
      noiseData.gainSqr[j] = (realVal*realVal) + (imagVal*imagVal);
    }

    noiseData.totalNoise += fabs(noiseData.noiseDens[j]); // sum of all noise sources of the device (not sure what this is for!)
    noiseData.outputNoiseDens[j] = noiseData.gainSqr[j] * fabs(noiseData.noiseDens[j]); // output noise density of the j'th noise source
    noiseData.lnNoiseDens[j] = std::log(std::max( noiseData.outputNoiseDens[j],N_MINLOG) ); // log of output noise density of the j'th noise source
    noiseData.relativeNoiseDens[j] = noiseData.outputNoiseDens[j] / outputValSqr_; // AM or PM noise density
    noiseData.totalOutputNoise += noiseData.outputNoiseDens[j]; // total output noise density of the device
  }
  noiseData.totalRelativeNoiseDens = noiseData.totalOutputNoise / outputValSqr_; // total AM or PM noise density of the device
}




//-----------------------------------------------------------------------------
// Function      : HBNOISE::createHarmonicSpaceLinearSystem_
// Purpose       : Creates the harmonic coupled matrix based on LTV system
//                 Ct_ and Gt_. Similar to NOISE::createACLinearSystem_() 
//                 but handles all the harmonics processed in HB
// Special Notes : The set of functions updateLinearTimeVariantSystem_C_and_G_() and 
//                 createHarmonicSpaceLinearSystem_() are computationally suboptimal
//                 and they apply Fourier tranform to a sparse matrix. They most likely
//                 will also fail to handle Touchstone files.
//                 But they are easier to understand for defining our Gold Standard
//                 They will be replaced later with a more efficient approach.
//                 similar to how HBLoader::loadDAEVectors() works
// Scope         : private
// Creator       : Meysam Bahmanian
// Creation Date : 5/22/2025
//-----------------------------------------------------------------------------
//                                Formulation in time domain
//                     for AC system (with Einstein convention):
// G_ij * v_j + C_ij * dv_j/dt = b_i
// Now assume linear system with single-tone excitation:
// v_j = V_jC*cos(w*t) - V_jS*sin(w*t)
// b_j = B_jC*cos(w*t) - B_jS*sin(w*t)
// where C means cos and S means sin.
// Then we have:
// G_ij * V_jC*cos(w*t) - G_ij * V_jS*sin(w*t) + C_ij * (-w*V_jC*sin(w*t) - w*V_jS*cos(w*t)) = B_iC*cos(w*t) - B_iS*sin(w*t)
// This gives us the famous matrix equation:
// G_ij * V_jC - w * C_ij * V_jS = B_iC
// w * C_ij * V_jC + G_ij * V_jS  = B_iS
//-----------------------------------------------------------------------------
//                     for Harmonic space system (with Einstein convention):
// G_ij(t) = G_ij0 + G_ijfI*cos(f*wc*t) - G_ijfQ*sin(f*wc*t) // sum over f
// where f is the harmonic number, wc is the carrier frequency, I means in-phase and Q means quadrature.
// The form of b_i is
// b_i = B_iC*cos(wm*t) - B_iS*sin(wm*t) + B_ieLC*cos(e*wc*t-wm*t) - B_ieLS*sin(e*wc*t-wm*t) + B_ieUC*cos(e*wc*t+wm*t) - B_ieUS*sin(e*wc*t+wm*t) // sum over e
// where e is the harmonic number, wc is the carrier frequency, L means lower sideband and U means upper sideband.
// The form of v_j is
// v_j = V_jC*cos(wm*t) - V_jS*sin(wm*t) + (V_jeIC*cos(wm*t) - V_jeIS*sin(wm*t))*cos(e*wc*t) - (V_jeQC*cos(wm*t) - V_jeQS*sin(wm*t))*sin(e*wc*t) // sum over e
// where e is the harmonic number, wc is the carrier frequency, I means in-phase and Q means quadrature.
// and wm is the modulation frequency.
// Now we need to multiply G_ij(t) with v_j and equate it to b_i(t)
// This gives us the Harmonic space matrix
// we first restructure v_j as lower and upper sidebands:
// baseband:
// V_jC*cos(wm*t) - V_jS*sin(wm*t) 
// LSB terms:
// + 0.5*(+V_jeIC+V_jeQS)*cos(e*wc*t-wm*t) - 0.5*(-V_jeIS+V_jeQC)*sin(e*wc*t-wm*t)
// USB terms:
// + 0.5*(+V_jeIC-V_jeQS)*cos(e*wc*t+wm*t) - 0.5*(+V_jeIS+V_jeQC)*sin(e*wc*t+wm*t)

// instead of writing the matrix equation, we will write the contributions:
// for all harmonics k, G_ij0 represents the AC linear system without any frequency transformation
// All harmonics are linearly transformed. But we need to formulate them as upper and lower sidebands.
// So we have:
// G_ij0 * (V_jC*cos(wm*t) - V_jS*sin(wm*t)) = B_iC*cos(wm*t) - B_iS*sin(wm*t)
// for LSBs and USBs:
// G_ij0 * 0.5 * [ (+V_jeIC+V_jeQS)*cos(e*wc*t-wm*t) - (-V_jeIS+V_jeQC)*sin(e*wc*t-wm*t) ] = B_ieLC*cos(e*wc*t-wm*t) - B_ieLS*sin(e*wc*t-wm*t)
// G_ij0 * 0.5 * [ (+V_jeIC-V_jeQS)*cos(e*wc*t+wm*t) - (+V_jeIS+V_jeQC)*sin(e*wc*t+wm*t) ] = B_ieUC*cos(e*wc*t+wm*t) - B_ieUS*sin(e*wc*t+wm*t)

// Now we have to find the contrinbutions of terms
// [ G_ijfI*cos(f*wc*t) - G_ijfQ*sin(f*wc*t) ] * v_j
//
// a) The baseband terms: V_jC*cos(wm*t) - V_jS*sin(wm*t)
// are directly transformed to f'th harmonic upper and lower sidebands, resulting the terms:
// + G_ijfI*V_jC*cos(wm*t)*cos(f*wc*t) // f'th harmonic upper sideband cosine term
// - G_ijfI*V_jS*sin(wm*t)*cos(f*wc*t) // f'th harmonic upper sideband sine term
// - G_ijfQ*V_jC*cos(wm*t)*sin(f*wc*t) // f'th harmonic quadrature cosine term
// + G_ijfQ*V_jS*sin(wm*t)*sin(f*wc*t) // f'th harmonic quadrature sine term
//
// b) equal harmonic mixing terms (e=f): 
// c) non-equal harmonic mixing terms (e!=f):
// Calculations of these terms are done in the function updateHarmonicSpaceMatrix_omegaC_0_():

// You may wonder what is the benefit of sine and cosine terms in the above equations.
// They show their significance when we have caps (surprisingly, even linear caps), these terms model the transfer functions in harmonic space.

// Now let's find the caps contribution
// we need to first correct our interpretation of the caps
// i = d/dt( C(t) * v(t) ) // this can also be shown that is true for small signal v(t). But let's accept it for now, I don't want to prove it here.
// i = dC(t)/dt * v(t) + C(t) * dv/dt
// So dC(t)/dt acts similar to conductance matrix and be added to G
// I call dC(t)/dt PART 0
// We simply multiply each Cf_ with its frequency. (see details in PART 0 below)
// Let's focus on other parts genereated by C(t) * dv/dt
// for C we have:
// C_ij(t) = C_ij0 + C_ijfI*cos(f*wc*t) - C_ijfQ*sin(f*wc*t) // sum over f
// v_j(t) has the form of
// v_j = V_jC*cos(wm*t) - V_jS*sin(wm*t) + 0.5*(+V_jeIC+V_jeQS)*cos(e*wc*t-wm*t) - 0.5*(-V_jeIS+V_jeQC)*sin(e*wc*t-wm*t) + 0.5*(+V_jeIC-V_jeQS)*cos(e*wc*t+wm*t) - 0.5*(+V_jeIS+V_jeQC)*sin(e*wc*t+wm*t)
// we divide the derivative in 2 parts, one that scales with wm and one that scales with wc
// the wm part is offset freq and has to be swept in a loop and has to be evaluated for every loop iteration
// PART 1 scales with wc and performs cross-harmonic quadrature transformation:
// dv_j/dt |1 = 0.5*(+V_jeIS-V_jeQC)*cos(e*wc*t-wm*t) - 0.5*(+V_jeIC+V_jeQS)*sin(e*wc*t-wm*t) + 0.5*(-V_jeIS-V_jeQC)*cos(e*wc*t+wm*t) - 0.5*(+V_jeIC-V_jeQS)*sin(e*wc*t-wm*t) // all * e*wc
// PART 2 scales with wm and performs coss-harmonic shaping:
// dv_j/dt |2 = - V_jS*cos(wm*t) - V_jC*sin(wm*t) + 0.5*(-V_jeIS+V_jeQC)*cos(e*wc*t-wm*t) - 0.5*(-V_jeIC-V_jeQS)*sin(e*wc*t-wm*t) + 0.5*(-V_jeIS-V_jeQC)*cos(e*wc*t+wm*t) - 0.5*(+V_jeIC-V_jeQS)*sin(e*wc*t-wm*t) // all * wm

// This approach gives us the desired format:
// ( [G] + [omega*C, PART 0] + [omega*C, PART 1] + wm*[C, PART 2] ) * v = b
// The matrices will be calculated once, for each offset frequency PART 2 will be scaled with wm and added to the matrix.


bool HBNOISE::createHarmonicSpaceLinearSystem_()
{
  // first take the Fourier Transform of Ct_ and Gt_

  int BlockCount = Ct_[0]->blockCount(); // number of time points
  int BlockSize = Ct_[0]->blockSize(); // number of GIDs

  for (int i=0; i<BlockSize; i++){
    Cf_.push_back(hbBuilderPtr_->createExpandedRealFormTransposeBlockVector());
    Gf_.push_back(hbBuilderPtr_->createExpandedRealFormTransposeBlockVector());
    Cf_[i]->putScalar(0.0);
    Gf_[i]->putScalar(0.0);
    hbLoaderPtr_->permutedFFT2(*(Ct_[i]), &*(Cf_[i]));
    hbLoaderPtr_->permutedFFT2(*(Gt_[i]), &*(Gf_[i]));
  }
  if (DEBUG_HBNOISE)
  {
    Xyce::dout() << "Reporting Gf_ Matrices, each array element is a row of the matrix, each block is a node" << std::endl;
    for (int i=0; i<BlockSize; i++){
      Xyce::dout() << "Gf_[" << i << "]: " << std::endl;
      Gf_[i]->print(Xyce::dout());
      Xyce::dout() << std::endl;
    }
    // Xyce::dout() << "Reporting Cf_ Matrices, each array element is a row of the matrix, each block is a node" << std::endl;
    // for (int i=0; i<BlockSize; i++){
    //   Xyce::dout() << "Cf_[" << i << "]: " << std::endl;
    //   Cf_[i]->print(Xyce::dout());
    //   Xyce::dout() << std::endl;
    // }
  }

  Parallel::Manager &pds_manager = *analysisManager_.getPDSManager();

  RCP<Parallel::ParMap> baseMap = rcp(pds_manager.getParallelMap( Parallel::SOLUTION ), false);
  const Linear::Graph* baseFullGraph = pds_manager.getMatrixGraph(Parallel::JACOBIAN);

  int numHarms = (size_-1)/2;
  int numBlocks = 2 + 4 * numHarms; // 2 is for our old AC part (real/imag) and every harmonic has 4 terms (real/imag in-phase and real/imag quadrature)
  int offset = baseMap->maxGlobalEntity() + 1;  // Use this offset to create a contiguous gid map for direct solvers.

  RCP<Parallel::ParMap> blockMap = Linear::createBlockParMap(numBlocks, *baseMap, 0, 0, offset);
  harmonicSpaceBI_ = Xyce::Linear::createBlockVector(numBlocks, blockMap, baseMap);
  harmonicSpaceBI_->putScalar(0.0);
  harmonicSpaceBQ_ = Xyce::Linear::createBlockVector(numBlocks, blockMap, baseMap);
  harmonicSpaceBQ_->putScalar(0.0);

  std::vector<std::vector<int> > blockPattern(numBlocks);
  for (int i=0; i<numBlocks; i++){
    blockPattern[i].resize(numBlocks);
    for (int j=0; j<numBlocks; j++){
      blockPattern[i][j] = j;
    }
  }

  RCP<Linear::Graph> blockGraph = Linear::createBlockGraph( offset, blockPattern, *blockMap, *baseFullGraph);

  // There will be 4 fundamental Matrices to be created:
  // harmonicSpaceMatrix_G_        : G matrix in harmonic space
  // harmonicSpaceMatrix_omegaC_0_ : wc*dC/dt matrix in harmonic space
  // harmonicSpaceMatrix_omegaC_1_ : wc*C*dv/dt matrix in harmonic space
  // harmonicSpaceMatrix_C_2_      : wm*C*dv/dt matrix in harmonic space
  // The master equation is: 
  // ( [G] + [omegaC, PART 0] + [omegaC, PART 1] + wm*[C, PART 2] ) * v = b
  // The first 3 matrices on LHS are not functions of offset frequency, they are constant.
  // The last matrix is a function of offset frequency, it is the only matrix that changes with offset frequency.
  // So I define a new Matrix to sum the first 3 matrices: 
  // harmonicSpaceMatrixConstant_ = [G] + [omegaC, PART 0] + [omegaC, PART 1]
  // the master equation becomes:
  // ( harmonicSpaceMatrixConstant_ + wm*[C, PART 2] ) * v = b
  // the sum of these two matrices is the overall matrix to be solved for:
  // harmonicSpaceMatrix_ = harmonicSpaceMatrixConstant_ + wm*[C, PART 2]

  harmonicSpaceMatrix_G_ = Xyce::Linear::createBlockMatrix( numBlocks, offset, blockPattern, blockGraph.get(), baseFullGraph);
  updateHarmonicSpaceMatrix_G_();

  harmonicSpaceMatrix_omegaC_0_ = Xyce::Linear::createBlockMatrix( numBlocks, offset, blockPattern, blockGraph.get(), baseFullGraph);
  updateHarmonicSpaceMatrix_omegaC_0_();

  harmonicSpaceMatrix_omegaC_1_ = Xyce::Linear::createBlockMatrix( numBlocks, offset, blockPattern, blockGraph.get(), baseFullGraph);
  updateHarmonicSpaceMatrix_omegaC_1_();

  harmonicSpaceMatrix_C_2_ = Xyce::Linear::createBlockMatrix( numBlocks, offset, blockPattern, blockGraph.get(), baseFullGraph);
  updateHarmonicSpaceMatrix_C_2_();

  // this will be used in the frequency loop
  harmonicSpaceMatrix_omegamC_2_ = Xyce::Linear::createBlockMatrix( numBlocks, offset, blockPattern, blockGraph.get(), baseFullGraph);
  harmonicSpaceMatrix_omegamC_2_->put( 0.0 );
  
  if (DEBUG_HBNOISE)
  {
    Xyce::dout() << "Reporting harmonicSpaceMatrix_G_:" << std::endl;
    harmonicSpaceMatrix_G_->print(Xyce::dout());
    // Xyce::dout() << "Reporting harmonicSpaceMatrix_omegaC_0_:" << std::endl;
    // harmonicSpaceMatrix_omegaC_0_->print(Xyce::dout());
    // Xyce::dout() << "Reporting harmonicSpaceMatrix_omegaC_1_:" << std::endl;
    // harmonicSpaceMatrix_omegaC_1_->print(Xyce::dout());
    // Xyce::dout() << "Reporting harmonicSpaceMatrix_C_2_:" << std::endl;
    // harmonicSpaceMatrix_C_2_->print(Xyce::dout());
    // Xyce::dout() << std::endl;
  }

  // This will be the overal matrix to be solved for.
  harmonicSpaceMatrix_ = Xyce::Linear::createBlockMatrix( numBlocks, offset, blockPattern, blockGraph.get(), baseFullGraph);
  harmonicSpaceMatrix_->put( 0.0 ); 

  // This is the constant part of the matrix, it is the sum of the first 3 matrices:
  harmonicSpaceMatrixConstant_ = Xyce::Linear::createBlockMatrix( numBlocks, offset, blockPattern, blockGraph.get(), baseFullGraph);
  harmonicSpaceMatrixConstant_->put( 0.0 ); 
  harmonicSpaceMatrixConstant_->add( *harmonicSpaceMatrix_G_ );
  harmonicSpaceMatrixConstant_->add( *harmonicSpaceMatrix_omegaC_0_ );
  harmonicSpaceMatrixConstant_->add( *harmonicSpaceMatrix_omegaC_1_ );


  harmonicSpaceXI_ = Xyce::Linear::createBlockVector (numBlocks, blockMap, baseMap);
  harmonicSpaceXI_->putScalar( 0.0 );

  harmonicSpaceXQ_ = Xyce::Linear::createBlockVector (numBlocks, blockMap, baseMap);
  harmonicSpaceXQ_->putScalar( 0.0 );

  harmonicSpace_SavedXI_ = Xyce::Linear::createBlockVector (numBlocks, blockMap, baseMap);
  harmonicSpace_SavedXI_->putScalar( 0.0 );

  harmonicSpace_SavedXQ_ = Xyce::Linear::createBlockVector (numBlocks, blockMap, baseMap);
  harmonicSpace_SavedXQ_->putScalar( 0.0 );

  blockProblemI_ = Xyce::Linear::createProblem( harmonicSpaceMatrix_, harmonicSpaceXI_, harmonicSpaceBI_ );
  blockProblemQ_ = Xyce::Linear::createProblem( harmonicSpaceMatrix_, harmonicSpaceXQ_, harmonicSpaceBQ_ );

  Linear::TranSolverFactory factory;
  blockSolverI_ = factory.create( linSolOptionBlock_, *blockProblemI_, analysisManager_.getCommandLine() );
  blockSolverQ_ = factory.create( linSolOptionBlock_, *blockProblemQ_, analysisManager_.getCommandLine() );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : HBNOISE::updateHarmonicSpaceMatrix_G_()
// Purpose       : Updates G matrix in harmonic space
// Scope         : private
// Creator       : Meysam Bahmanian
// Creation Date : 5/30/2025
//-----------------------------------------------------------------------------
// I currently find the harmonic space matrices in a loop iteration. But I think it will be easier using tensors. Fourier transform of G_
// is a tensor of rank 3. We can use a tranformation tensor to get the harmonic space G matrix.
// Ct_ is also a tensor:
//
//     Gt_ (or Ct_) tensor                     Fourier matrix                       Gf_ (or Cf_) tensor    
//
//           ---------|                          \    /                                  \-------\           
//  time->  /       / |         freq expansion->  \  / <-time compression                |\       \  <-freq
//         /-------/  |    *                       \/                       =            \ \-------\  
//  node-> | Ct/Gt | /                                                            node->  \| Cf/Gf |   
//         |-------|/                                                                      \-------|   
//             ^                                                                               ^        
//           node                                                                             node       
//
// I think an interesting approach would be to find a way to transform the Gf_ tensor to a flat matrix for the harmonic space.
// If we properly write the upper and lower rank indices, we have:
// Gf_{jkf} : j,k are members of {0-N}, f is member of {0-numHarms}, where N = GIDmax
// G_{j'k'} : j',k' are members of {0-N'} where N'=GIDmax*(2+4*numHarms)
// From a purely mathematical point of view, it would be interesting to find the form of harmonic space tranformation tensor, HST:
// HST_{j',k'}^{j,k,f}
// which has lower rank of 2 and upper rank of 3. Multiplying Gf tensor with HST tensor directly gives harmonic space matrix
// The elements of HST are most likey just 0, +-1, +-0.5 but this requires further formulation.

// Let's start with G matrix, it is easier!
// Gf_ strucure: std::vector with GID elements (each element corresponds to a row of Gf matrix), each element has GID blocks, each block has real/imag freq points
// the frequency points are in total 2*size_ = 2 * (2*numberOfHarmonics + 1)
// 0: dc real
// 1: dc imag
// 2: 1st harmonic real 
// 3: 1st harmonic imag etc.
// summary:
// G_ij(t) = G_ij0 + G_ijfI*cos(f*wc*t) - G_ijfQ*sin(f*wc*t) // sum over f
// b_i = B_iC*cos(wm*t) - B_iS*sin(wm*t) + B_ieLC*cos(e*wc*t-wm*t) - B_ieLS*sin(e*wc*t-wm*t) + B_ieUC*cos(e*wc*t+wm*t) - B_ieUS*sin(e*wc*t+wm*t) // sum over e
// v_j = V_jC*cos(wm*t) - V_jS*sin(wm*t) + (V_jeIC*cos(wm*t) - V_jeIS*sin(wm*t))*cos(e*wc*t) - (V_jeQC*cos(wm*t) - V_jeQS*sin(wm*t))*sin(e*wc*t) // sum over e
// v_j = V_jC*cos(wm*t) - V_jS*sin(wm*t) + 0.5*(+V_jeIC+V_jeQS)*cos(e*wc*t-wm*t) - 0.5*(-V_jeIS+V_jeQC)*sin(e*wc*t-wm*t) + 0.5*(+V_jeIC-V_jeQS)*cos(e*wc*t+wm*t) - 0.5*(+V_jeIS+V_jeQC)*sin(e*wc*t+wm*t)


// *************** Syntax decryption ***************
// G_ijfA : i/j are row/column indices; f is harmonic number; A is I or Q and means cosine or sine term
// for dc we have just G_ij0
// V_jeAB: j is the row index, e is the harmonic number, A is I (In-phase) or Q (Quadrature), B is C (Cosine) or S (Sine)
// for dc we have just V_jC and V_jS
// ***************      example      ***************
// - 0.5 * G_ijfQ * V_jeIS * sin(f*wc*t+wm*t)
// This equation defines setting a harmonic space matrix element.
// sin(f*wc*t+wm*t) : this is USB sine term and has the row block-address of 4*f+1 in the RHS and consequently the column block-address of 4*f+1 in the LHS Matrix
// V_jeIS : this is the j-th row of the LHS block-vector and consequently determinies the column index of the block-matrix
// eIS determines the block address: eIC (4*e-2); eIS (4*e-1); eQC (4*e); eQS (4*e+1)

bool HBNOISE::updateHarmonicSpaceMatrix_G_()
{
  harmonicSpaceMatrix_G_->put( 0.0 );
  int numHarms = (size_-1)/2;
  int numBlocks = 2 + 4 * numHarms;
  int numRows = Gf_.size();

  for (int i=0; i<numRows; i++)
  { // selectring a row of G Matrix
    // we should select only non-zero columns, But vectors are not treated like matrices and I can't get a view of augmented indices.
    // So I currently do a zero checking. This approach should be optimized.
    for (int j=0; j<numRows; j++)
    { // selecting a column of G Matrix
      for (int f=0; f<=numHarms; f++)
      { // selecting a frequency point of Gf_
        // zero checking
        if ( Gf_[i]->block(j)[2*f] == 0.0 && Gf_[i]->block(j)[2*f+1] == 0.0 )
        { continue; }
        // now we fill the matrix
        // first the diagonal of G matrix
        if (f==0)
        { // linear transformation
          // G_ij0 * [ V_jC*cos(wm*t) - V_jS*sin(wm*t) + 0.5*(+V_jeIC+V_jeQS)*cos(e*wc*t-wm*t) - 0.5*(-V_jeIS+V_jeQC)*sin(e*wc*t-wm*t) + 0.5*(+V_jeIC-V_jeQS)*cos(e*wc*t+wm*t) - 0.5*(+V_jeIS+V_jeQC)*sin(e*wc*t+wm*t) ]
          for (int e=0; e<=numHarms; e++){
            if (e==0) 
            { // G_ij0 * [ V_jC*cos(wm*t) - V_jS*sin(wm*t) ]

              // + G_ij0 * V_jC * cos(wm*t)
              // RHS has positive sign
              setMatrixElement(harmonicSpaceMatrix_G_->block(0,0), i, j, Gf_[i]->block(j)[0]);

              // - G_ij0 * V_jS*sin(wm*t)
              // RHS has negative sign
              setMatrixElement(harmonicSpaceMatrix_G_->block(1,1), i, j, Gf_[i]->block(j)[0]);
            } else 
            { // translation of I/Q format to LSB/USB

              // LSB cosine
              // + G_ij0 * 0.5 * (+V_jeIC+V_jeQS)*cos(e*wc*t-wm*t)
              // RHS has positive sign
              // + G_ij0 * 0.5 * V_jeIC * cos(e*wc*t-wm*t)
              setMatrixElement(harmonicSpaceMatrix_G_->block(4*e-2,4*e-2), i, j, +0.5*Gf_[i]->block(j)[0]);
              // + G_ij0 * 0.5 * V_jeQS * cos(e*wc*t-wm*t)
              setMatrixElement(harmonicSpaceMatrix_G_->block(4*e-2,4*e+1), i, j, +0.5*Gf_[i]->block(j)[0]);
              
              // LSB sine
              // + G_ij0 * [ - 0.5*(-V_jeIS+V_jeQC)*sin(e*wc*t-wm*t) ]
              // RHS has negative sign
              // + G_ij0 * 0.5 * V_jeIS * sin(e*wc*t-wm*t)
              setMatrixElement(harmonicSpaceMatrix_G_->block(4*e-1,4*e-1), i, j, -0.5*Gf_[i]->block(j)[0]);
              // - G_ij0 * 0.5 * V_jeQC * sin(e*wc*t-wm*t)
              setMatrixElement(harmonicSpaceMatrix_G_->block(4*e-1,4*e  ), i, j, +0.5*Gf_[i]->block(j)[0]);

              // USB cosine
              // G_ij0 * [ + 0.5*(+V_jeIC-V_jeQS)*cos(e*wc*t+wm*t) ]
              // RHS has positive sign
              // + G_ij0 * 0.5 * V_jeIC * cos(e*wc*t+wm*t)
              setMatrixElement(harmonicSpaceMatrix_G_->block(4*e  ,4*e-2), i, j, +0.5*Gf_[i]->block(j)[0]);
              // - G_ij0 * 0.5 * V_jeQS * cos(e*wc*t+wm*t)
              setMatrixElement(harmonicSpaceMatrix_G_->block(4*e  ,4*e+1), i, j, -0.5*Gf_[i]->block(j)[0]);

              // USB sine
              // + G_ij0 * [- 0.5*(+V_jeIS+V_jeQC)*sin(e*wc*t+wm*t) ]
              // RHS has negative sign
              // - G_ij0 * 0.5 * V_jeIS * sin(e*wc*t+wm*t)
              setMatrixElement(harmonicSpaceMatrix_G_->block(4*e+1,4*e-1), i, j, +0.5*Gf_[i]->block(j)[0]);
              // - G_ij0 * 0.5 * V_jeQC * sin(e*wc*t+wm*t)
              setMatrixElement(harmonicSpaceMatrix_G_->block(4*e+1,4*e  ), i, j, +0.5*Gf_[i]->block(j)[0]);
            }
          }
        } else
        { // now the mixing parts (harmonic coupling)
          for (int e=0; e<=numHarms; e++)
          {
            if (e==0) 
            { // baseband modulation
              // [ G_ijfI*cos(f*wc*t) - G_ijfQ*sin(f*wc*t) ] * [ V_jC*cos(wm*t) - V_jS*sin(wm*t) ]

              // LSB cosine
              // + G_ijfI*cos(f*wc*t) * V_jC*cos(wm*t) + G_ijfQ*sin(f*wc*t) * V_jS*sin(wm*t)
              // RHS has positive sign
              // + 0.5 * G_ijfI * V_jC * cos(f*wc*t-wm*t)
              setMatrixElement(harmonicSpaceMatrix_G_->block(4*f-2,0), i, j, +0.5*2*Gf_[i]->block(j)[2*f  ]);
              // + 0.5 * G_ijfQ * V_jS * cos(f*wc*t-wm*t)
              setMatrixElement(harmonicSpaceMatrix_G_->block(4*f-2,1), i, j, +0.5*2*Gf_[i]->block(j)[2*f+1]);

              // LSB sine
              // - G_ijfI*cos(f*wc*t) * V_jS*sin(wm*t) - G_ijfQ*sin(f*wc*t) * V_jC*cos(wm*t)
              // RHS has negative sign
              // + 0.5 * G_ijfI * V_jS * sin(f*wc*t-wm*t)
              setMatrixElement(harmonicSpaceMatrix_G_->block(4*f-1,1), i, j, -0.5*2*Gf_[i]->block(j)[2*f  ]);
              // - 0.5 * G_ijfQ * V_jC * sin(f*wc*t-wm*t)
              setMatrixElement(harmonicSpaceMatrix_G_->block(4*f-1,0), i, j, +0.5*2*Gf_[i]->block(j)[2*f+1]);

              // USB cosine
              // + G_ijfI*cos(f*wc*t) * V_jC*cos(wm*t) + G_ijfQ*sin(f*wc*t) * V_jS*sin(wm*t)
              // RHS has positive sign
              // + 0.5 * G_ijfI * V_jC * cos(f*wc*t+wm*t)
              setMatrixElement(harmonicSpaceMatrix_G_->block(4*f  ,0), i, j, +0.5*2*Gf_[i]->block(j)[2*f  ]);
              // - 0.5 * G_ijfQ * V_jS * cos(f*wc*t+wm*t)
              setMatrixElement(harmonicSpaceMatrix_G_->block(4*f  ,1), i, j, -0.5*2*Gf_[i]->block(j)[2*f+1]);

              // USB sine
              // - G_ijfI*cos(f*wc*t) * V_jS*sin(wm*t) - G_ijfQ*sin(f*wc*t) * V_jC*cos(wm*t)
              // RHS has negative sign
              // - 0.5 * G_ijfI * V_jS * sin(f*wc*t+wm*t)
              setMatrixElement(harmonicSpaceMatrix_G_->block(4*f+1,1), i, j, +0.5*2*Gf_[i]->block(j)[2*f  ]);
              // - 0.5 * G_ijfQ * V_jC * sin(f*wc*t+wm*t)
              setMatrixElement(harmonicSpaceMatrix_G_->block(4*f+1,0), i, j, +0.5*2*Gf_[i]->block(j)[2*f+1]);
            } else 
            { // now f>0 and e>0
              // 0.5 * [ G_ijfI*cos(f*wc*t) - G_ijfQ*sin(f*wc*t) ] * [ (+V_jeIC+V_jeQS)*cos(e*wc*t-wm*t) - (-V_jeIS+V_jeQC)*sin(e*wc*t-wm*t) + (+V_jeIC-V_jeQS)*cos(e*wc*t+wm*t) - (+V_jeIS+V_jeQC)*sin(e*wc*t+wm*t) ]

              int sigma = e+f;

              // case sigma
              if (sigma<=numHarms) 
              { // f+e should not be larget than numHarms, otherwise ignore it
                // in finite-harmonics space the system still shows nonlineary and some mixing products have to be ignored

                // LSB terms:
                // + 0.5 * [ G_ijfI*cos(f*wc*t) - G_ijfQ*sin(f*wc*t) ] * [ (+V_jeIC+V_jeQS)*cos(e*wc*t-wm*t) - (-V_jeIS+V_jeQC)*sin(e*wc*t-wm*t) ]

                // LSB cosine
                // + 0.5 *G_ijfI*cos(f*wc*t) * (+V_jeIC+V_jeQS)*cos(e*wc*t-wm*t) + 0.5 * G_ijfQ*sin(f*wc*t) * (-V_jeIS+V_jeQC)*sin(e*wc*t-wm*t)
                // at sigma:
                // + 0.25 * G_ijfI * (+V_jeIC+V_jeQS) * cos(sigma*wc*t-wm*t) + 0.25 * G_ijfQ * (+V_jeIS-V_jeQC) * cos(sigma*wc*t-wm*t)
                // RHS has positive sign
                // + 0.25 * G_ijfI * V_jeIC * cos(sigma*wc*t-wm*t)
                setMatrixElement(harmonicSpaceMatrix_G_->block(4*sigma-2, 4*e-2), i, j, +0.25*2*Gf_[i]->block(j)[2*f  ]);
                // + 0.25 * G_ijfI * V_jeQS * cos(sigma*wc*t-wm*t)
                setMatrixElement(harmonicSpaceMatrix_G_->block(4*sigma-2, 4*e+1), i, j, +0.25*2*Gf_[i]->block(j)[2*f  ]);
                // + 0.25 * G_ijfQ * V_jeIS * cos(sigma*wc*t-wm*t)
                setMatrixElement(harmonicSpaceMatrix_G_->block(4*sigma-2, 4*e-1), i, j, +0.25*2*Gf_[i]->block(j)[2*f+1]);
                // - 0.25 * G_ijfQ * V_jeQC * cos(sigma*wc*t-wm*t)
                setMatrixElement(harmonicSpaceMatrix_G_->block(4*sigma-2, 4*e  ), i, j, -0.25*2*Gf_[i]->block(j)[2*f+1]);

                // LSB sine
                // + 0.5 * G_ijfI * cos(f*wc*t) * (+V_jeIS-V_jeQC)*sin(e*wc*t-wm*t) + 0.5 * G_ijfQ * sin(f*wc*t) * (-V_jeIC-V_jeQS)*cos(e*wc*t-wm*t)
                // at sigma:
                // + 0.25 * G_ijfI * (+V_jeIS-V_jeQC)*sin(sigma*wc*t-wm*t) + 0.25 * G_ijfQ * (-V_jeIC-V_jeQS)*sin(sigma*wc*t-wm*t)
                // RHS has negative sign
                // + 0.25 * G_ijfI * V_jeIS * sin(sigma*wc*t-wm*t)
                setMatrixElement(harmonicSpaceMatrix_G_->block(4*sigma-1, 4*e-1), i, j, -0.25*2*Gf_[i]->block(j)[2*f  ]);
                // - 0.25 * G_ijfI * V_jeQC * sin(sigma*wc*t-wm*t)
                setMatrixElement(harmonicSpaceMatrix_G_->block(4*sigma-1, 4*e  ), i, j, +0.25*2*Gf_[i]->block(j)[2*f  ]);
                // - 0.25 * G_ijfQ * V_jeIC * sin(sigma*wc*t-wm*t)
                setMatrixElement(harmonicSpaceMatrix_G_->block(4*sigma-1, 4*e-2), i, j, +0.25*2*Gf_[i]->block(j)[2*f+1]);
                // - 0.25 * G_ijfQ * V_jeQS * sin(sigma*wc*t-wm*t)
                setMatrixElement(harmonicSpaceMatrix_G_->block(4*sigma-1, 4*e+1), i, j, +0.25*2*Gf_[i]->block(j)[2*f+1]);


                // USB terms:
                // + 0.5 * [ G_ijfI*cos(f*wc*t) - G_ijfQ*sin(f*wc*t) ] * [ (+V_jeIC-V_jeQS)*cos(e*wc*t+wm*t) - (+V_jeIS+V_jeQC)*sin(e*wc*t+wm*t) ]

                // USB cosine
                // + 0.5 * G_ijfI*cos(f*wc*t) * (+V_jeIC-V_jeQS)*cos(e*wc*t+wm*t) + 0.5 * G_ijfQ*sin(f*wc*t) * (+V_jeIS+V_jeQC)*sin(e*wc*t+wm*t)
                // at sigma:
                // + 0.25 * G_ijfI * (+V_jeIC-V_jeQS) * cos(sigma*wc*t+wm*t) + 0.25 * G_ijfQ * (-V_jeIS-V_jeQC) * cos(sigma*wc*t+wm*t)
                // RHS has positive sign
                // + 0.25 * G_ijfI * V_jeIC * cos(sigma*wc*t+wm*t)
                setMatrixElement(harmonicSpaceMatrix_G_->block(4*sigma  , 4*e-2), i, j, +0.25*2*Gf_[i]->block(j)[2*f  ]);
                // - 0.25 * G_ijfI * V_jeQS * cos(sigma*wc*t+wm*t)
                setMatrixElement(harmonicSpaceMatrix_G_->block(4*sigma  , 4*e+1), i, j, -0.25*2*Gf_[i]->block(j)[2*f  ]);
                // - 0.25 * G_ijfQ * V_jeIS * cos(sigma*wc*t+wm*t)
                setMatrixElement(harmonicSpaceMatrix_G_->block(4*sigma  , 4*e-1), i, j, -0.25*2*Gf_[i]->block(j)[2*f+1]);
                // - 0.25 * G_ijfQ * V_jeQC * cos(sigma*wc*t+wm*t)
                setMatrixElement(harmonicSpaceMatrix_G_->block(4*sigma  , 4*e  ), i, j, -0.25*2*Gf_[i]->block(j)[2*f+1]);

                // USB sine
                // + 0.5 * G_ijfI * cos(f*wc*t) * (-V_jeIS-V_jeQC)*sin(e*wc*t+wm*t) + 0.5 * G_ijfQ * sin(f*wc*t) * (-V_jeIC+V_jeQS)*cos(e*wc*t+wm*t)
                // at sigma:
                // + 0.25 * G_ijfI * (-V_jeIS-V_jeQC)*sin(sigma*wc*t+wm*t) + 0.25 * G_ijfQ * (-V_jeIC+V_jeQS)*sin(sigma*wc*t+wm*t)
                // RHS has negative sign
                // - 0.25 * G_ijfI * V_jeIS * sin(sigma*wc*t+wm*t)
                setMatrixElement(harmonicSpaceMatrix_G_->block(4*sigma+1, 4*e-1), i, j, +0.25*2*Gf_[i]->block(j)[2*f  ]);
                // - 0.25 * G_ijfI * V_jeQC * sin(sigma*wc*t+wm*t)
                setMatrixElement(harmonicSpaceMatrix_G_->block(4*sigma+1, 4*e  ), i, j, +0.25*2*Gf_[i]->block(j)[2*f  ]);
                // - 0.25 * G_ijfQ * V_jeIC * sin(sigma*wc*t+wm*t)
                setMatrixElement(harmonicSpaceMatrix_G_->block(4*sigma+1, 4*e-2), i, j, +0.25*2*Gf_[i]->block(j)[2*f+1]);
                // + 0.25 * G_ijfQ * V_jeQS * sin(sigma*wc*t+wm*t)
                setMatrixElement(harmonicSpaceMatrix_G_->block(4*sigma+1, 4*e+1), i, j, -0.25*2*Gf_[i]->block(j)[2*f+1]);
              }

              // int delta = e-f; this line is just for the sake of understanding the code. Equations are with reference to delta, not deltaAbs.
              int deltaAbs = std::abs(e-f);

              // case delta
              // now I have to deal with the index issue! if delta==0 then the indices become negative!
              if (deltaAbs==0)
              {
                // transformation from harmonics to baseband
                // 0.5 * [ G_ijfI*cos(f*wc*t) - G_ijfQ*sin(f*wc*t) ] * [ (+V_jeIC+V_jeQS)*cos(e*wc*t-wm*t) - (-V_jeIS+V_jeQC)*sin(e*wc*t-wm*t) + (+V_jeIC-V_jeQS)*cos(e*wc*t+wm*t) - (+V_jeIS+V_jeQC)*sin(e*wc*t+wm*t) ]

                // cosine terms:
                // + 0.5 * G_ijfI*cos(f*wc*t) * (+V_jeIC+V_jeQS)*cos(e*wc*t-wm*t) 
                // + 0.5 * G_ijfI*cos(f*wc*t) * (+V_jeIC-V_jeQS)*cos(e*wc*t+wm*t)
                // + 0.5 * G_ijfQ*sin(f*wc*t) * (-V_jeIS+V_jeQC)*sin(e*wc*t-wm*t)
                // + 0.5 * G_ijfQ*sin(f*wc*t) * (+V_jeIS+V_jeQC)*sin(e*wc*t+wm*t)

                // at delta=0
                // + 0.25 * G_ijfI * (+V_jeIC+V_jeQS) * cos(wm*t)
                // + 0.25 * G_ijfI * (+V_jeIC-V_jeQS) * cos(wm*t)
                // + 0.25 * G_ijfQ * (-V_jeIS+V_jeQC) * cos(wm*t)
                // + 0.25 * G_ijfQ * (+V_jeIS+V_jeQC) * cos(wm*t)

                // final terms:
                // + 0.5 * G_ijfI * V_jeIC * cos(wm*t) + 0.5 * G_ijfQ * V_jeQC * cos(wm*t)

                // we could have directly used I/Q forms too!
                // + G_ijfI*cos(f*wc*t) * V_jeIC*cos(wm*t) * cos(e*wc*t) + G_ijfQ*sin(f*wc*t) * V_jeQC*cos(wm*t) * sin(e*wc*t)
                // + 0.5*G_ijfI * V_jeIC*cos(wm*t) + 0.5*G_ijfQ * V_jeQC*cos(wm*t)
                // RHS has positive sign
                // + 0.5*G_ijfI * V_jeIC*cos(wm*t)
                setMatrixElement(harmonicSpaceMatrix_G_->block(0, 4*e-2), i, j, +0.5*2*Gf_[i]->block(j)[2*f  ]); // baseband cosine translated from in-phase cosine by in-phase G
                // + 0.5*G_ijfQ * V_jeQC*cos(wm*t)
                setMatrixElement(harmonicSpaceMatrix_G_->block(0, 4*e  ), i, j, +0.5*2*Gf_[i]->block(j)[2*f+1]); // baseband cosine translated from quadrature cosine by quadrature G

                // - G_ijfI*cos(f*wc*t) * V_jeIS*sin(wm*t) * cos(e*wc*t) - G_ijfQ*sin(f*wc*t) * V_jeQS*sin(wm*t) * sin(e*wc*t)
                // - 0.5*G_ijfI * V_jeIS*sin(wm*t) - 0.5*G_ijfQ * V_jeQS*sin(wm*t)
                // RHS has negative sign
                // - 0.5*G_ijfI * V_jeIS*sin(wm*t)
                setMatrixElement(harmonicSpaceMatrix_G_->block(1, 4*e-1), i, j, +0.5*2*Gf_[i]->block(j)[2*f  ]); // baseband sine translated from in-phase sine by in-phase G
                // - 0.5*G_ijfQ * V_jeQS*sin(wm*t)
                setMatrixElement(harmonicSpaceMatrix_G_->block(1, 4*e+1), i, j, +0.5*2*Gf_[i]->block(j)[2*f+1]); // baseband sine translated from quadrature sine by quadrature G
              }
              else
              {
                // here LSB and USB depend on the sign of delta, so I use the same structure as sigma case, but I call it conditional LSB and USB
                // a conditional LSB can be both LSB and USB!
                // a conditional USB can be both LSB and USB!
                int sign = e>f ? 1 : -1;
                int rowCosine = sign==1 ? 4*deltaAbs-2 : 4*deltaAbs  ;
                int rowSine   = sign==1 ? 4*deltaAbs-1 : 4*deltaAbs+1;

                // conditional LSB terms:
                // + 0.5 * [ G_ijfI*cos(f*wc*t) - G_ijfQ*sin(f*wc*t) ] * [ (+V_jeIC+V_jeQS)*cos(e*wc*t-wm*t) - (-V_jeIS+V_jeQC)*sin(e*wc*t-wm*t) ]

                // conditional LSB cosine
                // + 0.5 *G_ijfI*cos(f*wc*t) * (+V_jeIC+V_jeQS)*cos(e*wc*t-wm*t) + 0.5 * G_ijfQ*sin(f*wc*t) * (-V_jeIS+V_jeQC)*sin(e*wc*t-wm*t)
                // at delta:
                // + 0.25 * G_ijfI * (+V_jeIC+V_jeQS) * cos(delta*wc*t-wm*t) + 0.25 * G_ijfQ * (-V_jeIS+V_jeQC) * cos(delta*wc*t-wm*t)
                // RHS has positive sign
                // + 0.25 * G_ijfI * V_jeIC * cos(delta*wc*t-wm*t)
                setMatrixElement(harmonicSpaceMatrix_G_->block(rowCosine, 4*e-2), i, j, +0.25*2*Gf_[i]->block(j)[2*f  ]);
                // + 0.25 * G_ijfI * V_jeQS * cos(delta*wc*t-wm*t)
                setMatrixElement(harmonicSpaceMatrix_G_->block(rowCosine, 4*e+1), i, j, +0.25*2*Gf_[i]->block(j)[2*f  ]);
                // - 0.25 * G_ijfQ * V_jeIS * cos(delta*wc*t-wm*t)
                setMatrixElement(harmonicSpaceMatrix_G_->block(rowCosine, 4*e-1), i, j, -0.25*2*Gf_[i]->block(j)[2*f+1]);
                // + 0.25 * G_ijfQ * V_jeQC * cos(delta*wc*t-wm*t)
                setMatrixElement(harmonicSpaceMatrix_G_->block(rowCosine, 4*e  ), i, j, +0.25*2*Gf_[i]->block(j)[2*f+1]);

                // conditional LSB sine
                // + 0.5 * G_ijfI * cos(f*wc*t) * (+V_jeIS-V_jeQC)*sin(e*wc*t-wm*t) + 0.5 * G_ijfQ * sin(f*wc*t) * (-V_jeIC-V_jeQS)*cos(e*wc*t-wm*t)
                // at delta:
                // + 0.25 * G_ijfI * (+V_jeIS-V_jeQC)*sin(delta*wc*t-wm*t) + 0.25 * G_ijfQ * (+V_jeIC+V_jeQS)*sin(delta*wc*t-wm*t)
                // RHS has negative sign
                // + 0.25 * G_ijfI * V_jeIS * sin(delta*wc*t-wm*t)
                setMatrixElement(harmonicSpaceMatrix_G_->block(rowSine, 4*e-1), i, j, -sign*0.25*2*Gf_[i]->block(j)[2*f  ]);
                // - 0.25 * G_ijfI * V_jeQC * sin(delta*wc*t-wm*t)
                setMatrixElement(harmonicSpaceMatrix_G_->block(rowSine, 4*e  ), i, j, +sign*0.25*2*Gf_[i]->block(j)[2*f  ]);
                // + 0.25 * G_ijfQ * V_jeIC * sin(delta*wc*t-wm*t)
                setMatrixElement(harmonicSpaceMatrix_G_->block(rowSine, 4*e-2), i, j, -sign*0.25*2*Gf_[i]->block(j)[2*f+1]);
                // + 0.25 * G_ijfQ * V_jeQS * sin(delta*wc*t-wm*t)
                setMatrixElement(harmonicSpaceMatrix_G_->block(rowSine, 4*e+1), i, j, -sign*0.25*2*Gf_[i]->block(j)[2*f+1]);

                // conditional USB terms:
                // + 0.5 * [ G_ijfI*cos(f*wc*t) - G_ijfQ*sin(f*wc*t) ] * [ (+V_jeIC-V_jeQS)*cos(e*wc*t+wm*t) - (+V_jeIS+V_jeQC)*sin(e*wc*t+wm*t) ]

                // conditional USB cosine
                // + 0.5 * G_ijfI*cos(f*wc*t) * (+V_jeIC-V_jeQS)*cos(e*wc*t+wm*t) + 0.5 * G_ijfQ*sin(f*wc*t) * (+V_jeIS+V_jeQC)*sin(e*wc*t+wm*t)
                // at delta:
                // + 0.25 * G_ijfI * (+V_jeIC-V_jeQS) * cos(delta*wc*t+wm*t) + 0.25 * G_ijfQ * (+V_jeIS+V_jeQC) * cos(delta*wc*t+wm*t)
                // RHS has positive sign
                // + 0.25 * G_ijfI * V_jeIC * cos(delta*wc*t+wm*t)
                setMatrixElement(harmonicSpaceMatrix_G_->block(rowCosine, 4*e-2), i, j, +0.25*2*Gf_[i]->block(j)[2*f  ]);
                // - 0.25 * G_ijfI * V_jeQS * cos(delta*wc*t+wm*t)
                setMatrixElement(harmonicSpaceMatrix_G_->block(rowCosine, 4*e+1), i, j, -0.25*2*Gf_[i]->block(j)[2*f  ]);
                // + 0.25 * G_ijfQ * V_jeIS * cos(delta*wc*t+wm*t)
                setMatrixElement(harmonicSpaceMatrix_G_->block(rowCosine, 4*e-1), i, j, +0.25*2*Gf_[i]->block(j)[2*f+1]);
                // + 0.25 * G_ijfQ * V_jeQC * cos(delta*wc*t+wm*t)
                setMatrixElement(harmonicSpaceMatrix_G_->block(rowCosine, 4*e  ), i, j, +0.25*2*Gf_[i]->block(j)[2*f+1]);

                // conditional USB sine
                // + 0.5 * G_ijfI * cos(f*wc*t) * (-V_jeIS-V_jeQC)*sin(e*wc*t+wm*t) + 0.5 * G_ijfQ * sin(f*wc*t) * (-V_jeIC+V_jeQS)*cos(e*wc*t+wm*t)
                // at delta:
                // + 0.25 * G_ijfI * (-V_jeIS-V_jeQC)*sin(delta*wc*t+wm*t) + 0.25 * G_ijfQ * (+V_jeIC-V_jeQS)*sin(delta*wc*t+wm*t)
                // RHS has negative sign
                // - 0.25 * G_ijfI * V_jeIS * sin(delta*wc*t+wm*t)
                setMatrixElement(harmonicSpaceMatrix_G_->block(rowSine, 4*e-1), i, j, +0.25*sign*2*Gf_[i]->block(j)[2*f  ]);
                // - 0.25 * G_ijfI * V_jeQC * sin(delta*wc*t+wm*t)
                setMatrixElement(harmonicSpaceMatrix_G_->block(rowSine, 4*e  ), i, j, +0.25*sign*2*Gf_[i]->block(j)[2*f  ]);
                // + 0.25 * G_ijfQ * V_jeIC * sin(delta*wc*t+wm*t)
                setMatrixElement(harmonicSpaceMatrix_G_->block(rowSine, 4*e-2), i, j, -0.25*sign*2*Gf_[i]->block(j)[2*f+1]);
                // - 0.25 * G_ijfQ * V_jeQS * sin(delta*wc*t+wm*t)
                setMatrixElement(harmonicSpaceMatrix_G_->block(rowSine, 4*e+1), i, j, +0.25*sign*2*Gf_[i]->block(j)[2*f+1]);
              }
            } // end of f>0 and e>0
          } // end of e loop
        } // end of f!=0
      } // end of f loop of Gf_
    } // end of column loop
  } // end of row loop
  // And this was just for G matrix. Now we have to do the same for three parts of C matrix and we have to deal with derivatives too!
  return true;
}

//-----------------------------------------------------------------------------
// Function      : HBNOISE::updateHarmonicSpaceMatrix_omegaC_0_()
// Purpose       : Updates omegaC_0 matrix in harmonic space
// Scope         : private
// Creator       : Meysam Bahmanian
// Creation Date : 5/30/2025
//-----------------------------------------------------------------------------
// **************************************************
//           C Matrix PART 0: dC(t)/dt
// **************************************************
// This part behaves like G matrix. We just need to use time-derivative of C matrix instead of G matrix.
// We do the time-derivative in frequency domain.
// we have
// C_ij(t) = C_ij0 + C_ijfI*cos(f*wc*t) - C_ijfQ*sin(f*wc*t)            // sum over f
// dC_ij(t)/dt = - f*wc*C_ijfQ*cos(f*wc*t) - f*wc*C_ijfI*sin(f*wc*t)     // sum over f
// So I have to the same thing as G, but use 
// Gf_[i]->block(j)[2*f  ] -> -f*wc*Cf_[i]->block(j)[2*f+1]
// Gf_[i]->block(j)[2*f+1] -> +f*wc*Cf_[i]->block(j)[2*f  ]
// The only good news is that the linear transformation is now zero (well, I mean that tiny loop!)
// summary:
// C_ij(t) = C_ij0 + C_ijfI*cos(f*wc*t) - C_ijfQ*sin(f*wc*t)            // sum over f
// dC_ij(t)/dt = - f*wc*C_ijfQ*cos(f*wc*t) - f*wc*C_ijfI*sin(f*wc*t)     // sum over f
// b_i = B_iC*cos(wm*t) - B_iS*sin(wm*t) + B_ieLC*cos(e*wc*t-wm*t) - B_ieLS*sin(e*wc*t-wm*t) + B_ieUC*cos(e*wc*t+wm*t) - B_ieUS*sin(e*wc*t+wm*t) // sum over e
// v_j = V_jC*cos(wm*t) - V_jS*sin(wm*t) + (V_jeIC*cos(wm*t) - V_jeIS*sin(wm*t))*cos(e*wc*t) - (V_jeQC*cos(wm*t) - V_jeQS*sin(wm*t))*sin(e*wc*t) // sum over e
// v_j = V_jC*cos(wm*t) - V_jS*sin(wm*t) + 0.5*(+V_jeIC+V_jeQS)*cos(e*wc*t-wm*t) - 0.5*(-V_jeIS+V_jeQC)*sin(e*wc*t-wm*t) + 0.5*(+V_jeIC-V_jeQS)*cos(e*wc*t+wm*t) - 0.5*(+V_jeIS+V_jeQC)*sin(e*wc*t+wm*t)
bool HBNOISE::updateHarmonicSpaceMatrix_omegaC_0_()
{
  harmonicSpaceMatrix_omegaC_0_->put( 0.0 );
  int numHarms = (size_-1)/2;
  int numBlocks = 2 + 4 * numHarms;
  int numRows = Cf_.size();
  for (int i=0; i<numRows; i++)
  { // selectring a row of G Matrix
    for (int j=0; j<numRows; j++)
    { // selecting a column of G Matrix
      for (int f=0; f<=numHarms; f++)
      { // selecting a frequency point of Gf_
        if ( Cf_[i]->block(j)[2*f] == 0.0 && Cf_[i]->block(j)[2*f+1] == 0.0 )
        { continue; }
        // now we fill the matrix
        // first the diagon of omega*C matrix
        if (f==0)
        { // linear transformation
          // these terms are zero for dC(t)/dt
        } else
        { // now the mixing parts (harmonic coupling)
          for (int e=0; e<=numHarms; e++)
          {
            // ************************************************************************
            // dC_ij(t)/dt = - f*wc*C_ijfQ*cos(f*wc*t) - f*wc*C_ijfI*sin(f*wc*t)
            // I'm going to make these replacements:
            // G_ijfI -> -f*wc*C_ijfQ
            // G_ijfQ -> +f*wc*C_ijfI
            // and do not change the original comments written for G matrix.
            // ************************************************************************
            if (e==0) 
            { // baseband modulation
              // [ G_ijfI*cos(f*wc*t) - G_ijfQ*sin(f*wc*t) ] * [ V_jC*cos(wm*t) - V_jS*sin(wm*t) ]

              // LSB cosine
              // + G_ijfI*cos(f*wc*t) * V_jC*cos(wm*t) + G_ijfQ*sin(f*wc*t) * V_jS*sin(wm*t)
              // RHS has positive sign
              // + 0.5 * G_ijfI * V_jC * cos(f*wc*t-wm*t)
              // setMatrixElement(harmonicSpaceMatrix_G_->block(4*f-2,0), i, j, +0.5*Gf_[i]->block(j)[2*f  ]);
              setMatrixElement(harmonicSpaceMatrix_omegaC_0_->block(4*f-2,0), i, j, +0.5*(-f*omega_)*2*Cf_[i]->block(j)[2*f+1]);
              // + 0.5 * G_ijfQ * V_jS * cos(f*wc*t-wm*t)
              // setMatrixElement(harmonicSpaceMatrix_G_->block(4*f-2,1), i, j, +0.5*Gf_[i]->block(j)[2*f+1]);
              setMatrixElement(harmonicSpaceMatrix_omegaC_0_->block(4*f-2,1), i, j, +0.5*(+f*omega_)*2*Cf_[i]->block(j)[2*f]);

              // LSB sine
              // - G_ijfI*cos(f*wc*t) * V_jS*sin(wm*t) - G_ijfQ*sin(f*wc*t) * V_jC*cos(wm*t)
              // RHS has negative sign
              // + 0.5 * G_ijfI * V_jS * sin(f*wc*t-wm*t)
              // setMatrixElement(harmonicSpaceMatrix_G_->block(4*f-1,1), i, j, -0.5*Gf_[i]->block(j)[2*f  ]);
              setMatrixElement(harmonicSpaceMatrix_omegaC_0_->block(4*f-1,1), i, j, -0.5*(-f*omega_)*2*Cf_[i]->block(j)[2*f+1]);
              // - 0.5 * G_ijfQ * V_jC * sin(f*wc*t-wm*t)
              // setMatrixElement(harmonicSpaceMatrix_G_->block(4*f-1,0), i, j, +0.5*Gf_[i]->block(j)[2*f+1]);
              setMatrixElement(harmonicSpaceMatrix_omegaC_0_->block(4*f-1,0), i, j, +0.5*(+f*omega_)*2*Cf_[i]->block(j)[2*f  ]);

              // USB cosine
              // + G_ijfI*cos(f*wc*t) * V_jC*cos(wm*t) + G_ijfQ*sin(f*wc*t) * V_jS*sin(wm*t)
              // RHS has positive sign
              // + 0.5 * G_ijfI * V_jC * cos(f*wc*t+wm*t)
              // setMatrixElement(harmonicSpaceMatrix_G_->block(4*f  ,0), i, j, +0.5*Gf_[i]->block(j)[2*f  ]);
              setMatrixElement(harmonicSpaceMatrix_omegaC_0_->block(4*f  ,0), i, j, +0.5*(-f*omega_)*2*Cf_[i]->block(j)[2*f+1]);
              // - 0.5 * G_ijfQ * V_jS * cos(f*wc*t+wm*t)
              // setMatrixElement(harmonicSpaceMatrix_G_->block(4*f  ,1), i, j, -0.5*Gf_[i]->block(j)[2*f+1]);
              setMatrixElement(harmonicSpaceMatrix_omegaC_0_->block(4*f  ,1), i, j, -0.5*(+f*omega_)*2*Cf_[i]->block(j)[2*f  ]);

              // USB sine
              // - G_ijfI*cos(f*wc*t) * V_jS*sin(wm*t) - G_ijfQ*sin(f*wc*t) * V_jC*cos(wm*t)
              // RHS has negative sign
              // - 0.5 * G_ijfI * V_jS * sin(f*wc*t+wm*t)
              // setMatrixElement(harmonicSpaceMatrix_G_->block(4*f+1,1), i, j, +0.5*Gf_[i]->block(j)[2*f  ]);
              setMatrixElement(harmonicSpaceMatrix_omegaC_0_->block(4*f+1,1), i, j, +0.5*(-f*omega_)*2*Cf_[i]->block(j)[2*f+1]);
              // - 0.5 * G_ijfQ * V_jC * sin(f*wc*t+wm*t)
              // setMatrixElement(harmonicSpaceMatrix_G_->block(4*f+1,0), i, j, +0.5*Gf_[i]->block(j)[2*f+1]);
              setMatrixElement(harmonicSpaceMatrix_omegaC_0_->block(4*f+1,0), i, j, +0.5*(+f*omega_)*2*Cf_[i]->block(j)[2*f  ]);
            } else 
            { // now f>0 and e>0
              // 0.5 * [ G_ijfI*cos(f*wc*t) - G_ijfQ*sin(f*wc*t) ] * [ (+V_jeIC+V_jeQS)*cos(e*wc*t-wm*t) - (-V_jeIS+V_jeQC)*sin(e*wc*t-wm*t) + (+V_jeIC-V_jeQS)*cos(e*wc*t+wm*t) - (+V_jeIS+V_jeQC)*sin(e*wc*t+wm*t) ]

              int sigma = e+f;

              // case sigma
              if (sigma<=numHarms) 
              { // f+e should not be larget than numHarms, otherwise ignore it
                // in finite-harmonics space the system still shows nonlineary and some mixing products have to be ignored

                // LSB terms:
                // + 0.5 * [ G_ijfI*cos(f*wc*t) - G_ijfQ*sin(f*wc*t) ] * [ (+V_jeIC+V_jeQS)*cos(e*wc*t-wm*t) - (-V_jeIS+V_jeQC)*sin(e*wc*t-wm*t) ]

                // LSB cosine
                // + 0.5 *G_ijfI*cos(f*wc*t) * (+V_jeIC+V_jeQS)*cos(e*wc*t-wm*t) + 0.5 * G_ijfQ*sin(f*wc*t) * (-V_jeIS+V_jeQC)*sin(e*wc*t-wm*t)
                // at sigma:
                // + 0.25 * G_ijfI * (+V_jeIC+V_jeQS) * cos(sigma*wc*t-wm*t) + 0.25 * G_ijfQ * (+V_jeIS-V_jeQC) * cos(sigma*wc*t-wm*t)
                // RHS has positive sign
                // + 0.25 * G_ijfI * V_jeIC * cos(sigma*wc*t-wm*t)
                // setMatrixElement(harmonicSpaceMatrix_G_->block(4*sigma-2, 4*e-2), i, j, +0.25*Gf_[i]->block(j)[2*f  ]);
                setMatrixElement(harmonicSpaceMatrix_omegaC_0_->block(4*sigma-2, 4*e-2), i, j, +0.25*(-f*omega_)*2*Cf_[i]->block(j)[2*f+1]);
                // + 0.25 * G_ijfI * V_jeQS * cos(sigma*wc*t-wm*t)
                // setMatrixElement(harmonicSpaceMatrix_G_->block(4*sigma-2, 4*e+1), i, j, +0.25*Gf_[i]->block(j)[2*f  ]);
                setMatrixElement(harmonicSpaceMatrix_omegaC_0_->block(4*sigma-2, 4*e+1), i, j, +0.25*(-f*omega_)*2*Cf_[i]->block(j)[2*f+1]);
                // + 0.25 * G_ijfQ * V_jeIS * cos(sigma*wc*t-wm*t)
                // setMatrixElement(harmonicSpaceMatrix_G_->block(4*sigma-2, 4*e-1), i, j, +0.25*Gf_[i]->block(j)[2*f+1]);
                setMatrixElement(harmonicSpaceMatrix_omegaC_0_->block(4*sigma-2, 4*e-1), i, j, +0.25*(+f*omega_)*2*Cf_[i]->block(j)[2*f  ]);
                // - 0.25 * G_ijfQ * V_jeQC * cos(sigma*wc*t-wm*t)
                // setMatrixElement(harmonicSpaceMatrix_G_->block(4*sigma-2, 4*e  ), i, j, -0.25*Gf_[i]->block(j)[2*f+1]);
                setMatrixElement(harmonicSpaceMatrix_omegaC_0_->block(4*sigma-2, 4*e  ), i, j, -0.25*(+f*omega_)*2*Cf_[i]->block(j)[2*f  ]);

                // LSB sine
                // + 0.5 * G_ijfI * cos(f*wc*t) * (+V_jeIS-V_jeQC)*sin(e*wc*t-wm*t) + 0.5 * G_ijfQ * sin(f*wc*t) * (-V_jeIC-V_jeQS)*cos(e*wc*t-wm*t)
                // at sigma:
                // + 0.25 * G_ijfI * (+V_jeIS-V_jeQC)*sin(sigma*wc*t-wm*t) + 0.25 * G_ijfQ * (-V_jeIC-V_jeQS)*sin(sigma*wc*t-wm*t)
                // RHS has negative sign
                // + 0.25 * G_ijfI * V_jeIS * sin(sigma*wc*t-wm*t)
                // setMatrixElement(harmonicSpaceMatrix_G_->block(4*sigma-1, 4*e-1), i, j, -0.25*Gf_[i]->block(j)[2*f  ]);
                setMatrixElement(harmonicSpaceMatrix_omegaC_0_->block(4*sigma-1, 4*e-1), i, j, -0.25*(-f*omega_)*2*Cf_[i]->block(j)[2*f+1]);
                // - 0.25 * G_ijfI * V_jeQC * sin(sigma*wc*t-wm*t)
                // setMatrixElement(harmonicSpaceMatrix_G_->block(4*sigma-1, 4*e  ), i, j, +0.25*Gf_[i]->block(j)[2*f  ]);
                setMatrixElement(harmonicSpaceMatrix_omegaC_0_->block(4*sigma-1, 4*e  ), i, j, +0.25*(-f*omega_)*2*Cf_[i]->block(j)[2*f+1]);
                // - 0.25 * G_ijfQ * V_jeIC * sin(sigma*wc*t-wm*t)
                // setMatrixElement(harmonicSpaceMatrix_G_->block(4*sigma-1, 4*e-2), i, j, +0.25*Gf_[i]->block(j)[2*f+1]);
                setMatrixElement(harmonicSpaceMatrix_omegaC_0_->block(4*sigma-1, 4*e-2), i, j, +0.25*(+f*omega_)*2*Cf_[i]->block(j)[2*f  ]);
                // - 0.25 * G_ijfQ * V_jeQS * sin(sigma*wc*t-wm*t)
                // setMatrixElement(harmonicSpaceMatrix_G_->block(4*sigma-1, 4*e+1), i, j, +0.25*Gf_[i]->block(j)[2*f+1]);
                setMatrixElement(harmonicSpaceMatrix_omegaC_0_->block(4*sigma-1, 4*e+1), i, j, +0.25*(+f*omega_)*2*Cf_[i]->block(j)[2*f  ]);


                // USB terms:
                // + 0.5 * [ G_ijfI*cos(f*wc*t) - G_ijfQ*sin(f*wc*t) ] * [ (+V_jeIC-V_jeQS)*cos(e*wc*t+wm*t) - (+V_jeIS+V_jeQC)*sin(e*wc*t+wm*t) ]

                // USB cosine
                // + 0.5 * G_ijfI*cos(f*wc*t) * (+V_jeIC-V_jeQS)*cos(e*wc*t+wm*t) + 0.5 * G_ijfQ*sin(f*wc*t) * (+V_jeIS+V_jeQC)*sin(e*wc*t+wm*t)
                // at sigma:
                // + 0.25 * G_ijfI * (+V_jeIC-V_jeQS) * cos(sigma*wc*t+wm*t) + 0.25 * G_ijfQ * (-V_jeIS-V_jeQC) * cos(sigma*wc*t+wm*t)
                // RHS has positive sign
                // + 0.25 * G_ijfI * V_jeIC * cos(sigma*wc*t+wm*t)
                // setMatrixElement(harmonicSpaceMatrix_G_->block(4*sigma  , 4*e-2), i, j, +0.25*Gf_[i]->block(j)[2*f  ]);
                setMatrixElement(harmonicSpaceMatrix_omegaC_0_->block(4*sigma  , 4*e-2), i, j, +0.25*(-f*omega_)*2*Cf_[i]->block(j)[2*f+1]);
                // - 0.25 * G_ijfI * V_jeQS * cos(sigma*wc*t+wm*t)
                // setMatrixElement(harmonicSpaceMatrix_G_->block(4*sigma  , 4*e+1), i, j, -0.25*Gf_[i]->block(j)[2*f  ]);
                setMatrixElement(harmonicSpaceMatrix_omegaC_0_->block(4*sigma  , 4*e+1), i, j, -0.25*(-f*omega_)*2*Cf_[i]->block(j)[2*f+1]);
                // - 0.25 * G_ijfQ * V_jeIS * cos(sigma*wc*t+wm*t)
                // setMatrixElement(harmonicSpaceMatrix_G_->block(4*sigma  , 4*e-1), i, j, -0.25*Gf_[i]->block(j)[2*f+1]);
                setMatrixElement(harmonicSpaceMatrix_omegaC_0_->block(4*sigma  , 4*e-1), i, j, -0.25*(+f*omega_)*2*Cf_[i]->block(j)[2*f  ]);
                // - 0.25 * G_ijfQ * V_jeQC * cos(sigma*wc*t+wm*t)
                // setMatrixElement(harmonicSpaceMatrix_G_->block(4*sigma  , 4*e  ), i, j, -0.25*Gf_[i]->block(j)[2*f+1]);
                setMatrixElement(harmonicSpaceMatrix_omegaC_0_->block(4*sigma  , 4*e  ), i, j, -0.25*(+f*omega_)*2*Cf_[i]->block(j)[2*f  ]);

                // USB sine
                // + 0.5 * G_ijfI * cos(f*wc*t) * (-V_jeIS-V_jeQC)*sin(e*wc*t+wm*t) + 0.5 * G_ijfQ * sin(f*wc*t) * (-V_jeIC+V_jeQS)*cos(e*wc*t+wm*t)
                // at sigma:
                // + 0.25 * G_ijfI * (-V_jeIS-V_jeQC)*sin(sigma*wc*t+wm*t) + 0.25 * G_ijfQ * (-V_jeIC+V_jeQS)*sin(sigma*wc*t+wm*t)
                // RHS has negative sign
                // - 0.25 * G_ijfI * V_jeIS * sin(sigma*wc*t+wm*t)
                // setMatrixElement(harmonicSpaceMatrix_G_->block(4*sigma+1, 4*e-1), i, j, +0.25*Gf_[i]->block(j)[2*f  ]);
                setMatrixElement(harmonicSpaceMatrix_omegaC_0_->block(4*sigma+1, 4*e-1), i, j, +0.25*(-f*omega_)*2*Cf_[i]->block(j)[2*f+1]);
                // - 0.25 * G_ijfI * V_jeQC * sin(sigma*wc*t+wm*t)
                // setMatrixElement(harmonicSpaceMatrix_G_->block(4*sigma+1, 4*e  ), i, j, +0.25*Gf_[i]->block(j)[2*f  ]);
                setMatrixElement(harmonicSpaceMatrix_omegaC_0_->block(4*sigma+1, 4*e  ), i, j, +0.25*(-f*omega_)*2*Cf_[i]->block(j)[2*f+1]);
                // - 0.25 * G_ijfQ * V_jeIC * sin(sigma*wc*t+wm*t)
                // setMatrixElement(harmonicSpaceMatrix_G_->block(4*sigma+1, 4*e-2), i, j, +0.25*Gf_[i]->block(j)[2*f+1]);
                setMatrixElement(harmonicSpaceMatrix_omegaC_0_->block(4*sigma+1, 4*e-2), i, j, +0.25*(+f*omega_)*2*Cf_[i]->block(j)[2*f  ]);
                // + 0.25 * G_ijfQ * V_jeQS * sin(sigma*wc*t+wm*t)
                // setMatrixElement(harmonicSpaceMatrix_G_->block(4*sigma+1, 4*e+1), i, j, -0.25*Gf_[i]->block(j)[2*f+1]);
                setMatrixElement(harmonicSpaceMatrix_omegaC_0_->block(4*sigma+1, 4*e+1), i, j, -0.25*(+f*omega_)*2*Cf_[i]->block(j)[2*f  ]);
              }

              // int delta = e-f; this line is just for the sake of understanding the code. Equations are with reference to delta, not deltaAbs.
              int deltaAbs = std::abs(e-f);

              // case delta
              // now I have to deal with the index issue! if delta==0 then the indices become negative!
              if (deltaAbs==0)
              {
                // transformation from harmonics to baseband
                // 0.5 * [ G_ijfI*cos(f*wc*t) - G_ijfQ*sin(f*wc*t) ] * [ (+V_jeIC+V_jeQS)*cos(e*wc*t-wm*t) - (-V_jeIS+V_jeQC)*sin(e*wc*t-wm*t) + (+V_jeIC-V_jeQS)*cos(e*wc*t+wm*t) - (+V_jeIS+V_jeQC)*sin(e*wc*t+wm*t) ]

                // cosine terms:
                // + 0.5 * G_ijfI*cos(f*wc*t) * (+V_jeIC+V_jeQS)*cos(e*wc*t-wm*t) 
                // + 0.5 * G_ijfI*cos(f*wc*t) * (+V_jeIC-V_jeQS)*cos(e*wc*t+wm*t)
                // + 0.5 * G_ijfQ*sin(f*wc*t) * (-V_jeIS+V_jeQC)*sin(e*wc*t-wm*t)
                // + 0.5 * G_ijfQ*sin(f*wc*t) * (+V_jeIS+V_jeQC)*sin(e*wc*t+wm*t)

                // at delta=0
                // + 0.25 * G_ijfI * (+V_jeIC+V_jeQS) * cos(wm*t)
                // + 0.25 * G_ijfI * (+V_jeIC-V_jeQS) * cos(wm*t)
                // + 0.25 * G_ijfQ * (-V_jeIS+V_jeQC) * cos(wm*t)
                // + 0.25 * G_ijfQ * (+V_jeIS+V_jeQC) * cos(wm*t)

                // final terms:
                // + 0.5 * G_ijfI * V_jeIC * cos(wm*t) + 0.5 * G_ijfQ * V_jeQC * cos(wm*t)

                // we could have directly used I/Q forms too!
                // + G_ijfI*cos(f*wc*t) * V_jeIC*cos(wm*t) * cos(e*wc*t) + G_ijfQ*sin(f*wc*t) * V_jeQC*cos(wm*t) * sin(e*wc*t)
                // + 0.5*G_ijfI * V_jeIC*cos(wm*t) + 0.5*G_ijfQ * V_jeQC*cos(wm*t)
                // RHS has positive sign
                // setMatrixElement(harmonicSpaceMatrix_G_->block(0, 4*e-2), i, j, +0.5*Gf_[i]->block(j)[2*f  ]); // baseband cosine translated from in-phase cosine by in-phase G
                setMatrixElement(harmonicSpaceMatrix_omegaC_0_->block(0, 4*e-2), i, j, +0.5*(-f*omega_)*2*Cf_[i]->block(j)[2*f+1]); 
                // setMatrixElement(harmonicSpaceMatrix_G_->block(0, 4*e  ), i, j, +0.5*Gf_[i]->block(j)[2*f+1]); // baseband cosine translated from quadrature cosine by quadrature G
                setMatrixElement(harmonicSpaceMatrix_omegaC_0_->block(0, 4*e  ), i, j, +0.5*(+f*omega_)*2*Cf_[i]->block(j)[2*f  ]); 

                // - G_ijfI*cos(f*wc*t) * V_jeIS*sin(wm*t) * cos(e*wc*t) - G_ijfQ*sin(f*wc*t) * V_jeQS*sin(wm*t) * sin(e*wc*t)
                // - 0.5*G_ijfI * V_jeIS*sin(wm*t) - 0.5*G_ijfQ * V_jeQS*sin(wm*t)
                // RHS has negative sign
                // setMatrixElement(harmonicSpaceMatrix_G_->block(1, 4*e-1), i, j, +0.5*Gf_[i]->block(j)[2*f  ]); // baseband sine translated from in-phase sine by in-phase G
                setMatrixElement(harmonicSpaceMatrix_omegaC_0_->block(1, 4*e-1), i, j, +0.5*(-f*omega_)*2*Cf_[i]->block(j)[2*f+1]);
                // setMatrixElement(harmonicSpaceMatrix_G_->block(1, 4*e+1), i, j, +0.5*Gf_[i]->block(j)[2*f+1]); // baseband sine translated from quadrature sine by quadrature G
                setMatrixElement(harmonicSpaceMatrix_omegaC_0_->block(1, 4*e+1), i, j, +0.5*(+f*omega_)*2*Cf_[i]->block(j)[2*f  ]);
              }
              else
              {
                // here LSB and USB depend on the sign of delta, so I use the same structure as sigma case, but I call it conditional LSB and USB
                // a conditional LSB can be both LSB and USB!
                // a conditional USB can be both LSB and USB!
                int sign = e>f ? 1 : -1;
                int rowCosine = sign==1 ? 4*deltaAbs-2 : 4*deltaAbs  ;
                int rowSine   = sign==1 ? 4*deltaAbs-1 : 4*deltaAbs+1;

                // conditional LSB terms:
                // + 0.5 * [ G_ijfI*cos(f*wc*t) - G_ijfQ*sin(f*wc*t) ] * [ (+V_jeIC+V_jeQS)*cos(e*wc*t-wm*t) - (-V_jeIS+V_jeQC)*sin(e*wc*t-wm*t) ]

                // conditional LSB cosine
                // + 0.5 *G_ijfI*cos(f*wc*t) * (+V_jeIC+V_jeQS)*cos(e*wc*t-wm*t) + 0.5 * G_ijfQ*sin(f*wc*t) * (-V_jeIS+V_jeQC)*sin(e*wc*t-wm*t)
                // at delta:
                // + 0.25 * G_ijfI * (+V_jeIC+V_jeQS) * cos(delta*wc*t-wm*t) + 0.25 * G_ijfQ * (-V_jeIS+V_jeQC) * cos(delta*wc*t-wm*t)
                // RHS has positive sign
                // + 0.25 * G_ijfI * V_jeIC * cos(delta*wc*t-wm*t)
                // setMatrixElement(harmonicSpaceMatrix_G_->block(rowCosine, 4*e-2), i, j, +0.25*Gf_[i]->block(j)[2*f  ]);
                setMatrixElement(harmonicSpaceMatrix_omegaC_0_->block(rowCosine, 4*e-2), i, j, +0.25*(-f*omega_)*2*Cf_[i]->block(j)[2*f+1]);
                // + 0.25 * G_ijfI * V_jeQS * cos(delta*wc*t-wm*t)
                //setMatrixElement(harmonicSpaceMatrix_G_->block(rowCosine, 4*e+1), i, j, +0.25*Gf_[i]->block(j)[2*f  ]);
                setMatrixElement(harmonicSpaceMatrix_omegaC_0_->block(rowCosine, 4*e+1), i, j, +0.25*(-f*omega_)*2*Cf_[i]->block(j)[2*f+1]);
                // - 0.25 * G_ijfQ * V_jeIS * cos(delta*wc*t-wm*t)
                // setMatrixElement(harmonicSpaceMatrix_G_->block(rowCosine, 4*e-1), i, j, -0.25*Gf_[i]->block(j)[2*f+1]);
                setMatrixElement(harmonicSpaceMatrix_omegaC_0_->block(rowCosine, 4*e-1), i, j, -0.25*(+f*omega_)*2*Cf_[i]->block(j)[2*f  ]);
                // + 0.25 * G_ijfQ * V_jeQC * cos(delta*wc*t-wm*t)
                // setMatrixElement(harmonicSpaceMatrix_G_->block(rowCosine, 4*e  ), i, j, +0.25*Gf_[i]->block(j)[2*f+1]);
                setMatrixElement(harmonicSpaceMatrix_omegaC_0_->block(rowCosine, 4*e  ), i, j, +0.25*(+f*omega_)*2*Cf_[i]->block(j)[2*f  ]);

                // conditional LSB sine
                // + 0.5 * G_ijfI * cos(f*wc*t) * (+V_jeIS-V_jeQC)*sin(e*wc*t-wm*t) + 0.5 * G_ijfQ * sin(f*wc*t) * (-V_jeIC-V_jeQS)*cos(e*wc*t-wm*t)
                // at delta:
                // + 0.25 * G_ijfI * (+V_jeIS-V_jeQC)*sin(delta*wc*t-wm*t) + 0.25 * G_ijfQ * (+V_jeIC+V_jeQS)*sin(delta*wc*t-wm*t)
                // RHS has negative sign
                // + 0.25 * G_ijfI * V_jeIS * sin(delta*wc*t-wm*t)
                // setMatrixElement(harmonicSpaceMatrix_G_->block(rowSine, 4*e-1), i, j, -sign*0.25*Gf_[i]->block(j)[2*f  ]);
                setMatrixElement(harmonicSpaceMatrix_omegaC_0_->block(rowSine, 4*e-1), i, j, -sign*0.25*(-f*omega_)*2*Cf_[i]->block(j)[2*f+1]);
                // - 0.25 * G_ijfI * V_jeQC * sin(delta*wc*t-wm*t)
                //setMatrixElement(harmonicSpaceMatrix_G_->block(rowSine, 4*e  ), i, j, +sign*0.25*Gf_[i]->block(j)[2*f  ]);
                setMatrixElement(harmonicSpaceMatrix_omegaC_0_->block(rowSine, 4*e  ), i, j, +sign*0.25*(-f*omega_)*2*Cf_[i]->block(j)[2*f+1]);
                // + 0.25 * G_ijfQ * V_jeIC * sin(delta*wc*t-wm*t)
                // setMatrixElement(harmonicSpaceMatrix_G_->block(rowSine, 4*e-2), i, j, -sign*0.25*Gf_[i]->block(j)[2*f+1]);
                setMatrixElement(harmonicSpaceMatrix_omegaC_0_->block(rowSine, 4*e-2), i, j, -sign*0.25*(+f*omega_)*2*Cf_[i]->block(j)[2*f  ]);
                // + 0.25 * G_ijfQ * V_jeQS * sin(delta*wc*t-wm*t)
                // setMatrixElement(harmonicSpaceMatrix_G_->block(rowSine, 4*e+1), i, j, -sign*0.25*Gf_[i]->block(j)[2*f+1]);
                setMatrixElement(harmonicSpaceMatrix_omegaC_0_->block(rowSine, 4*e+1), i, j, -sign*0.25*(+f*omega_)*2*Cf_[i]->block(j)[2*f  ]);

                // conditional USB terms:
                // + 0.5 * [ G_ijfI*cos(f*wc*t) - G_ijfQ*sin(f*wc*t) ] * [ (+V_jeIC-V_jeQS)*cos(e*wc*t+wm*t) - (+V_jeIS+V_jeQC)*sin(e*wc*t+wm*t) ]

                // conditional USB cosine
                // + 0.5 * G_ijfI*cos(f*wc*t) * (+V_jeIC-V_jeQS)*cos(e*wc*t+wm*t) + 0.5 * G_ijfQ*sin(f*wc*t) * (+V_jeIS+V_jeQC)*sin(e*wc*t+wm*t)
                // at delta:
                // + 0.25 * G_ijfI * (+V_jeIC-V_jeQS) * cos(delta*wc*t+wm*t) + 0.25 * G_ijfQ * (+V_jeIS+V_jeQC) * cos(delta*wc*t+wm*t)
                // RHS has positive sign
                // + 0.25 * G_ijfI * V_jeIC * cos(delta*wc*t+wm*t)
                // setMatrixElement(harmonicSpaceMatrix_G_->block(rowCosine, 4*e-2), i, j, +0.25*Gf_[i]->block(j)[2*f  ]);
                setMatrixElement(harmonicSpaceMatrix_omegaC_0_->block(rowCosine, 4*e-2), i, j, +0.25*(-f*omega_)*2*Cf_[i]->block(j)[2*f+1]);
                // - 0.25 * G_ijfI * V_jeQS * cos(delta*wc*t+wm*t)
                // setMatrixElement(harmonicSpaceMatrix_G_->block(rowCosine, 4*e+1), i, j, -0.25*Gf_[i]->block(j)[2*f  ]);
                setMatrixElement(harmonicSpaceMatrix_omegaC_0_->block(rowCosine, 4*e+1), i, j, -0.25*(-f*omega_)*2*Cf_[i]->block(j)[2*f+1]);
                // + 0.25 * G_ijfQ * V_jeIS * cos(delta*wc*t+wm*t)
                // setMatrixElement(harmonicSpaceMatrix_G_->block(rowCosine, 4*e-1), i, j, +0.25*Gf_[i]->block(j)[2*f+1]);
                setMatrixElement(harmonicSpaceMatrix_omegaC_0_->block(rowCosine, 4*e-1), i, j, +0.25*(+f*omega_)*2*Cf_[i]->block(j)[2*f  ]);
                // + 0.25 * G_ijfQ * V_jeQC * cos(delta*wc*t+wm*t)
                // setMatrixElement(harmonicSpaceMatrix_G_->block(rowCosine, 4*e  ), i, j, +0.25*Gf_[i]->block(j)[2*f+1]);
                setMatrixElement(harmonicSpaceMatrix_omegaC_0_->block(rowCosine, 4*e  ), i, j, +0.25*(+f*omega_)*2*Cf_[i]->block(j)[2*f  ]);

                // conditional USB sine
                // + 0.5 * G_ijfI * cos(f*wc*t) * (-V_jeIS-V_jeQC)*sin(e*wc*t+wm*t) + 0.5 * G_ijfQ * sin(f*wc*t) * (-V_jeIC+V_jeQS)*cos(e*wc*t+wm*t)
                // at delta:
                // + 0.25 * G_ijfI * (-V_jeIS-V_jeQC)*sin(delta*wc*t+wm*t) + 0.25 * G_ijfQ * (+V_jeIC-V_jeQS)*sin(delta*wc*t+wm*t)
                // RHS has negative sign
                // - 0.25 * G_ijfI * V_jeIS * sin(delta*wc*t+wm*t)
                // setMatrixElement(harmonicSpaceMatrix_G_->block(rowSine, 4*e-1), i, j, +0.25*sign*Gf_[i]->block(j)[2*f  ]);
                setMatrixElement(harmonicSpaceMatrix_omegaC_0_->block(rowSine, 4*e-1), i, j, +0.25*sign*(-f*omega_)*2*Cf_[i]->block(j)[2*f+1]);
                // - 0.25 * G_ijfI * V_jeQC * sin(delta*wc*t+wm*t)
                // setMatrixElement(harmonicSpaceMatrix_G_->block(rowSine, 4*e  ), i, j, +0.25*sign*Gf_[i]->block(j)[2*f  ]);
                setMatrixElement(harmonicSpaceMatrix_omegaC_0_->block(rowSine, 4*e  ), i, j, +0.25*sign*(-f*omega_)*2*Cf_[i]->block(j)[2*f+1]);
                // + 0.25 * G_ijfQ * V_jeIC * sin(delta*wc*t+wm*t)
                // setMatrixElement(harmonicSpaceMatrix_G_->block(rowSine, 4*e-2), i, j, -0.25*sign*Gf_[i]->block(j)[2*f+1]);
                setMatrixElement(harmonicSpaceMatrix_omegaC_0_->block(rowSine, 4*e-2), i, j, -0.25*sign*(+f*omega_)*2*Cf_[i]->block(j)[2*f  ]);
                // - 0.25 * G_ijfQ * V_jeQS * sin(delta*wc*t+wm*t)
                // setMatrixElement(harmonicSpaceMatrix_G_->block(rowSine, 4*e+1), i, j, +0.25*sign*Gf_[i]->block(j)[2*f+1]);
                setMatrixElement(harmonicSpaceMatrix_omegaC_0_->block(rowSine, 4*e+1), i, j, +0.25*sign*(+f*omega_)*2*Cf_[i]->block(j)[2*f  ]);
              }
            } // end of f>0 and e>0
          } // end of e loop
        } // end of f!=0
      } // end of f loop of Gf_
    } // end of column loop
  } // end of row loop
  return true;
}

//-----------------------------------------------------------------------------
// Function      : HBNOISE::updateHarmonicSpaceMatrix_omegaC_1_()
// Purpose       : Updates omegaC_1 matrix in harmonic space
// Scope         : private
// Creator       : Meysam Bahmanian
// Creation Date : 5/30/2025
//-----------------------------------------------------------------------------
// These patterns are really ugly! I wish I could put them into a function. Probably it needs someone smarter than me to do it!
// Now the real ugly equations start! C Matrix PART 1 and C Matrix PART 2
// The time derivative makes very confusing equations, an additinal quadrature transformation, both for carrier and modulated components.
// **************************************************
//           C Matrix PART 1: 
// **************************************************
// Reminder:
// C_ij(t) = C_ij0 + C_ijfI*cos(f*wc*t) - C_ijfQ*sin(f*wc*t) // sum over f
// v_j(t) has the form of
// v_j = V_jC*cos(wm*t) - V_jS*sin(wm*t) + (V_jeIC*cos(wm*t) - V_jeIS*sin(wm*t))*cos(e*wc*t) - (V_jeQC*cos(wm*t) - V_jeQS*sin(wm*t))*sin(e*wc*t) // sum over e
// PART 1 scales with wc and performs cross-harmonic quadrature transformation:
// dv_j/dt |1 = - ( V_jeIC*cos(wm*t) - V_jeIS*sin(wm*t) )*e*wc*sin(e*wc*t) - ( V_jeQC*cos(wm*t) - V_jeQS*sin(wm*t) )*e*wc*cos(e*wc*t)
// dv_j/dt |1 = [ - ( V_jeQC*cos(wm*t) - V_jeQS*sin(wm*t) )*cos(e*wc*t) - ( V_jeIC*cos(wm*t) - V_jeIS*sin(wm*t) )*sin(e*wc*t) ]*e*wc
// First observation: dv_j/dt |1 has no baseband component => the first two column blocks of PART 1 are 0.
// This equation has also some similarities with the equation of G Matrix, in the sense that the wm parts are intact. 
// However, we have " quadrature <=> in-phase " transformation because of derivative with respect to carrier component.
// So we can use the same pattern of filling G matrix, but this time we need to swap in-phase and quadrature components.
// for instance, row 0 goes to row 1 and row 1 goes to row 0.
// summary:
// C_ij(t) = C_ij0 + C_ijfI*cos(f*wc*t) - C_ijfQ*sin(f*wc*t)            // sum over f
// b_i = B_iC*cos(wm*t) - B_iS*sin(wm*t) + B_ieLC*cos(e*wc*t-wm*t) - B_ieLS*sin(e*wc*t-wm*t) + B_ieUC*cos(e*wc*t+wm*t) - B_ieUS*sin(e*wc*t+wm*t) // sum over e
// v_j = V_jC*cos(wm*t) - V_jS*sin(wm*t) + ( V_jeIC*cos(wm*t) - V_jeIS*sin(wm*t) ) * cos(e*wc*t) - ( V_jeQC*cos(wm*t) - V_jeQS*sin(wm*t) ) * sin(e*wc*t) // sum over e
// v_j = V_jC*cos(wm*t) - V_jS*sin(wm*t) + 0.5*(+V_jeIC+V_jeQS)*cos(e*wc*t-wm*t) - 0.5*(-V_jeIS+V_jeQC)*sin(e*wc*t-wm*t) + 0.5*(+V_jeIC-V_jeQS)*cos(e*wc*t+wm*t) - 0.5*(+V_jeIS+V_jeQC)*sin(e*wc*t+wm*t)
// IQ form:
// dv_j/dt |1 = [ (- V_jeQC*cos(wm*t) + V_jeQS*sin(wm*t) )*cos(e*wc*t) - ( V_jeIC*cos(wm*t) - V_jeIS*sin(wm*t) )*sin(e*wc*t) ] *e*wc
// LSB/USB form:
// dv_j/dt |1 = [ + 0.5*(+V_jeIS-V_jeQC)*cos(e*wc*t-wm*t) - 0.5*(+V_jeIC+V_jeQS)*sin(e*wc*t-wm*t) + 0.5*(-V_jeIS-V_jeQC)*cos(e*wc*t+wm*t) - 0.5*(+V_jeIC-V_jeQS)*sin(e*wc*t+wm*t) ]*e*wc
bool HBNOISE::updateHarmonicSpaceMatrix_omegaC_1_()
{
  harmonicSpaceMatrix_omegaC_1_->put( 0.0 );
  int numHarms = (size_-1)/2;
  int numBlocks = 2 + 4 * numHarms;
  int numRows = Cf_.size();
  for (int i=0; i<numRows; i++)
  { // selectring a row of C Matrix
    for (int j=0; j<numRows; j++)
    { // selecting a column of C Matrix
      for (int f=0; f<=numHarms; f++)
      { // selecting a frequency point of Cf_
        if ( Cf_[i]->block(j)[2*f] == 0.0 && Cf_[i]->block(j)[2*f+1] == 0.0 )
        { continue; }
        // now we fill the matrix
        if (f==0)
        { // linear transformation
          // C_ij0 * [ + 0.5*(+V_jeIS-V_jeQC)*cos(e*wc*t-wm*t) - 0.5*(+V_jeIC+V_jeQS)*sin(e*wc*t-wm*t) + 0.5*(-V_jeIS-V_jeQC)*cos(e*wc*t+wm*t) - 0.5*(+V_jeIC-V_jeQS)*sin(e*wc*t+wm*t) ] * e*wc
          for (int e=1; e<numHarms; e++)
          { // the first two blocks are baseband and zero (this matrix is not proportional to wm), hence starting from e=1

            // LSB cosine
            // C_ij0 * [ + 0.5*(+V_jeIS-V_jeQC)*cos(e*wc*t-wm*t) ] *e*wc
            // RHS has positive sign
            // + 0.5 * C_ij0 * V_jeIS * cos(e*wc*t-wm*t) *e*wc
            setMatrixElement(harmonicSpaceMatrix_omegaC_1_->block(4*e-2, 4*e-1), i, j, +0.5*(e*omega_)*Cf_[i]->block(j)[0]);
            // - 0.5 * C_ij0 * V_jeQC * cos(e*wc*t-wm*t) *e*wc
            setMatrixElement(harmonicSpaceMatrix_omegaC_1_->block(4*e-2, 4*e  ), i, j, -0.5*(e*omega_)*Cf_[i]->block(j)[0]);

            // LSB sine
            // C_ij0 * [ - 0.5*(+V_jeIC+V_jeQS)*sin(e*wc*t-wm*t) ] *e*wc
            // RHS has negative sign
            // - 0.5 * C_ij0 * V_jeIC * sin(e*wc*t-wm*t) *e*wc
            setMatrixElement(harmonicSpaceMatrix_omegaC_1_->block(4*e-1, 4*e-2), i, j, +0.5*(e*omega_)*Cf_[i]->block(j)[0]);
            // - 0.5 * C_ij0 * V_jeQS * sin(e*wc*t-wm*t) *e*wc
            setMatrixElement(harmonicSpaceMatrix_omegaC_1_->block(4*e-1, 4*e+1), i, j, +0.5*(e*omega_)*Cf_[i]->block(j)[0]);

            // USB cosine
            // C_ij0 * [ + 0.5*(-V_jeIS-V_jeQC)*cos(e*wc*t+wm*t) ] *e*wc
            // RHS has positive sign
            // - 0.5 * C_ij0 * V_jeIS * cos(e*wc*t+wm*t) *e*wc
            setMatrixElement(harmonicSpaceMatrix_omegaC_1_->block(4*e  , 4*e-1), i, j, -0.5*(e*omega_)*Cf_[i]->block(j)[0]);
            // - 0.5 * C_ij0 * V_jeQC * cos(e*wc*t+wm*t) *e*wc
            setMatrixElement(harmonicSpaceMatrix_omegaC_1_->block(4*e  , 4*e  ), i, j, -0.5*(e*omega_)*Cf_[i]->block(j)[0]);

            // USB sine
            // C_ij0 * [- 0.5*(+V_jeIC-V_jeQS)*sin(e*wc*t+wm*t) ] *e*wc
            // RHS has negative sign
            // - 0.5 * C_ij0 * V_jeIC * sin(e*wc*t+wm*t) *e*wc
            setMatrixElement(harmonicSpaceMatrix_omegaC_1_->block(4*e+1, 4*e-2), i, j, +0.5*(e*omega_)*Cf_[i]->block(j)[0]);
            // + 0.5 * C_ij0 * V_jeQS * sin(e*wc*t+wm*t) *e*wc
            setMatrixElement(harmonicSpaceMatrix_omegaC_1_->block(4*e+1, 4*e+1), i, j, -0.5*(e*omega_)*Cf_[i]->block(j)[0]);

          }
        } else
        { // now the mixing parts (harmonic coupling)
          for (int e=1; e<=numHarms; e++)
          { // now f>0 and e>0
            // [ C_ijfI*cos(f*wc*t) - C_ijfQ*sin(f*wc*t) ] * [ + 0.5*(+V_jeIS-V_jeQC)*cos(e*wc*t-wm*t) - 0.5*(+V_jeIC+V_jeQS)*sin(e*wc*t-wm*t) + 0.5*(-V_jeIS-V_jeQC)*cos(e*wc*t+wm*t) - 0.5*(+V_jeIC-V_jeQS)*sin(e*wc*t+wm*t) ] * e*wc

            // case sigma
            int sigma = f+e;
            if (sigma<=numHarms) 
            { // f+e should not be larget than numHarms, otherwise ignore it
              // in finite-harmonics space the system still shows nonlineary and some mixing products have to be ignored

              // LSB terms
              // [ C_ijfI*cos(f*wc*t) - C_ijfQ*sin(f*wc*t) ] * [ + 0.5*(+V_jeIS-V_jeQC)*cos(e*wc*t-wm*t) - 0.5*(+V_jeIC+V_jeQS)*sin(e*wc*t-wm*t) ] * e*wc

              // LSB cosine
              // + 0.5 *C_ijfI*cos(f*wc*t) * (+V_jeIS-V_jeQC)*cos(e*wc*t-wm*t) + 0.5 * C_ijfQ*sin(f*wc*t) * (+V_jeIC+V_jeQS)*sin(e*wc*t-wm*t) // all *e*wc
              // at sigma:
              // + 0.25 * C_ijfI * (+V_jeIS-V_jeQC) * cos(sigma*wc*t-wm*t) + 0.25 * C_ijfQ * (-V_jeIC-V_jeQS) * cos(sigma*wc*t-wm*t) // all *e*wc
              // RHS has positive sign
              // + 0.25 * C_ijfI * V_jeIS * cos(sigma*wc*t-wm*t) *e*wc
              setMatrixElement(harmonicSpaceMatrix_omegaC_1_->block(4*sigma-2, 4*e-1), i, j, +0.25*(e*omega_)*2*Cf_[i]->block(j)[2*f  ]);
              // - 0.25 * C_ijfI * V_jeQC * cos(sigma*wc*t-wm*t) *e*wc
              setMatrixElement(harmonicSpaceMatrix_omegaC_1_->block(4*sigma-2, 4*e  ), i, j, -0.25*(e*omega_)*2*Cf_[i]->block(j)[2*f  ]);
              // - 0.25 * C_ijfQ * V_jeIC * cos(sigma*wc*t-wm*t) *e*wc
              setMatrixElement(harmonicSpaceMatrix_omegaC_1_->block(4*sigma-2, 4*e-2), i, j, -0.25*(e*omega_)*2*Cf_[i]->block(j)[2*f+1]);
              // - 0.25 * C_ijfQ * V_jeQS * cos(sigma*wc*t-wm*t) *e*wc
              setMatrixElement(harmonicSpaceMatrix_omegaC_1_->block(4*sigma-2, 4*e+1), i, j, -0.25*(e*omega_)*2*Cf_[i]->block(j)[2*f+1]);

              // LSB sine
              // + 0.5 *C_ijfI*cos(f*wc*t) * (-V_jeIC-V_jeQS)*sin(e*wc*t-wm*t) + 0.5 * C_ijfQ*sin(f*wc*t) * (-V_jeIS+V_jeQC)*cos(e*wc*t-wm*t) // all *e*wc
              // at sigma:
              // + 0.25 * C_ijfI * (-V_jeIC-V_jeQS) * sin(sigma*wc*t-wm*t) + 0.25 * C_ijfQ * (-V_jeIS+V_jeQC) * sin(sigma*wc*t-wm*t) // all *e*wc
              // RHS has negative sign
              // - 0.25 * C_ijfI * V_jeIC * sin(sigma*wc*t-wm*t) *e*wc
              setMatrixElement(harmonicSpaceMatrix_omegaC_1_->block(4*sigma-1, 4*e-2), i, j, +0.25*(e*omega_)*2*Cf_[i]->block(j)[2*f  ]);
              // - 0.25 * C_ijfI * V_jeQS * sin(sigma*wc*t-wm*t) *e*wc
              setMatrixElement(harmonicSpaceMatrix_omegaC_1_->block(4*sigma-1, 4*e+1), i, j, +0.25*(e*omega_)*2*Cf_[i]->block(j)[2*f  ]);
              // - 0.25 * C_ijfQ * V_jeIS * sin(sigma*wc*t-wm*t) *e*wc
              setMatrixElement(harmonicSpaceMatrix_omegaC_1_->block(4*sigma-1, 4*e-1), i, j, +0.25*(e*omega_)*2*Cf_[i]->block(j)[2*f+1]);
              // + 0.25 * C_ijfQ * V_jeQC * sin(sigma*wc*t-wm*t) *e*wc
              setMatrixElement(harmonicSpaceMatrix_omegaC_1_->block(4*sigma-1, 4*e  ), i, j, -0.25*(e*omega_)*2*Cf_[i]->block(j)[2*f+1]);

              // USB terms
              // [ C_ijfI*cos(f*wc*t) - C_ijfQ*sin(f*wc*t) ] * [ + 0.5*(-V_jeIS-V_jeQC)*cos(e*wc*t+wm*t) - 0.5*(+V_jeIC-V_jeQS)*sin(e*wc*t+wm*t) ] * e*wc

              // USB cosine
              // + 0.5 *C_ijfI*cos(f*wc*t) * (-V_jeIS-V_jeQC)*cos(e*wc*t+wm*t) + 0.5 * C_ijfQ*sin(f*wc*t) * (+V_jeIC-V_jeQS)*sin(e*wc*t+wm*t) // all *e*wc
              // at sigma:
              // + 0.25 * C_ijfI * (-V_jeIS-V_jeQC) * cos(sigma*wc*t+wm*t) + 0.25 * C_ijfQ * (-V_jeIC+V_jeQS) * cos(sigma*wc*t+wm*t) // all *e*wc
              // RHS has positive sign
              // - 0.25 * C_ijfI * V_jeIS * cos(sigma*wc*t+wm*t) *e*wc
              setMatrixElement(harmonicSpaceMatrix_omegaC_1_->block(4*sigma  , 4*e-1), i, j, -0.25*(e*omega_)*2*Cf_[i]->block(j)[2*f  ]);
              // - 0.25 * C_ijfI * V_jeQC * cos(sigma*wc*t+wm*t) *e*wc
              setMatrixElement(harmonicSpaceMatrix_omegaC_1_->block(4*sigma  , 4*e  ), i, j, -0.25*(e*omega_)*2*Cf_[i]->block(j)[2*f  ]);
              // - 0.25 * C_ijfQ * V_jeIC * cos(sigma*wc*t+wm*t) *e*wc
              setMatrixElement(harmonicSpaceMatrix_omegaC_1_->block(4*sigma  , 4*e-2), i, j, -0.25*(e*omega_)*2*Cf_[i]->block(j)[2*f+1]);
              // + 0.25 * C_ijfQ * V_jeQS * cos(sigma*wc*t+wm*t) *e*wc
              setMatrixElement(harmonicSpaceMatrix_omegaC_1_->block(4*sigma  , 4*e+1), i, j, +0.25*(e*omega_)*2*Cf_[i]->block(j)[2*f+1]);

              // USB sine
              // + 0.5 *C_ijfI*cos(f*wc*t) * (-V_jeIC+V_jeQS)*sin(e*wc*t+wm*t) + 0.5 * C_ijfQ*sin(f*wc*t) * (+V_jeIS+V_jeQC)*cos(e*wc*t+wm*t) // all *e*wc
              // at sigma:
              // + 0.25 * C_ijfI * (-V_jeIC+V_jeQS) * sin(sigma*wc*t+wm*t) + 0.25 * C_ijfQ * (+V_jeIS+V_jeQC) * sin(sigma*wc*t+wm*t) // all *e*wc
              // RHS has negative sign
              // - 0.25 * C_ijfI * V_jeIC * sin(sigma*wc*t+wm*t) *e*wc
              setMatrixElement(harmonicSpaceMatrix_omegaC_1_->block(4*sigma+1, 4*e-2), i, j, +0.25*(e*omega_)*2*Cf_[i]->block(j)[2*f  ]);
              // + 0.25 * C_ijfI * V_jeQS * sin(sigma*wc*t+wm*t) *e*wc
              setMatrixElement(harmonicSpaceMatrix_omegaC_1_->block(4*sigma+1, 4*e+1), i, j, -0.25*(e*omega_)*2*Cf_[i]->block(j)[2*f  ]);
              // + 0.25 * C_ijfQ * V_jeIS * sin(sigma*wc*t+wm*t) *e*wc
              setMatrixElement(harmonicSpaceMatrix_omegaC_1_->block(4*sigma+1, 4*e-1), i, j, -0.25*(e*omega_)*2*Cf_[i]->block(j)[2*f+1]);
              // + 0.25 * C_ijfQ * V_jeQC * sin(sigma*wc*t+wm*t) *e*wc
              setMatrixElement(harmonicSpaceMatrix_omegaC_1_->block(4*sigma+1, 4*e  ), i, j, -0.25*(e*omega_)*2*Cf_[i]->block(j)[2*f+1]);
            }

            // int delta = e-f; this line is just for the sake of understanding the code. Equations are with reference to delta, not deltaAbs.
            int deltaAbs = std::abs(e-f);

            // case delta
            // now I have to deal with the index issue! if delta==0 then the indices become negative!
            if (deltaAbs==0)
            {
              // transformation from harmonics to baseband
              // I use I/Q form of dv_j/dt here, instead of LSB/USB forms. This will be easier for the baseband transformation.
              // [ C_ijfI*cos(f*wc*t) - C_ijfQ*sin(f*wc*t) ] * [(- V_jeQC*cos(wm*t) + V_jeQS*sin(wm*t) )*cos(e*wc*t) - ( V_jeIC*cos(wm*t) - V_jeIS*sin(wm*t) )*sin(e*wc*t)] * e*wc

              // baseband cosine term:
              // [ - C_ijfI*cos(f*wc*t) * V_jeQC * cos(wm*t) * cos(e*wc*t) + C_ijfQ*sin(f*wc*t) * V_jeIC * cos(wm*t) * sin(e*wc*t) ] *e*wc
              // at delta=0
              // - 0.5 * C_ijfI * V_jeQC * cos(wm*t) + 0.5 * C_ijfQ * V_jeIC * cos(wm*t) // all *e*wc
              // RHS has positive sign
              // - 0.5 * C_ijfI * V_jeQC * cos(wm*t) *e*wc
              setMatrixElement(harmonicSpaceMatrix_omegaC_0_->block(0, 4*e  ), i, j, -0.5*(e*omega_)*2*Cf_[i]->block(j)[2*f  ]); 
              // + 0.5 * C_ijfQ * V_jeIC * cos(wm*t) *e*wc
              setMatrixElement(harmonicSpaceMatrix_omegaC_0_->block(0, 4*e-2), i, j, +0.5*(e*omega_)*2*Cf_[i]->block(j)[2*f+1]); 
              
              // baseband sine term:
              // [ + C_ijfI*cos(f*wc*t) * V_jeQS * sin(wm*t) * cos(e*wc*t) - C_ijfQ*sin(f*wc*t) * V_jeIS * sin(wm*t) * sin(e*wc*t) ] *e*wc
              // at delta=0
              // + 0.5 * C_ijfI * V_jeQS * sin(wm*t) - 0.5 * C_ijfQ * V_jeIS * sin(wm*t) // all *e*wc
              // RHS has negative sign
              // + 0.5 * C_ijfI * V_jeQS * sin(wm*t) *e*wc
              setMatrixElement(harmonicSpaceMatrix_omegaC_0_->block(1, 4*e+1), i, j, -0.5*(e*omega_)*2*Cf_[i]->block(j)[2*f  ]);
              // - 0.5 * C_ijfQ * V_jeIS * sin(wm*t) *e*wc
              setMatrixElement(harmonicSpaceMatrix_omegaC_0_->block(1, 4*e-1), i, j, +0.5*(e*omega_)*2*Cf_[i]->block(j)[2*f+1]); 

            }
            else
            {
              // here LSB and USB depend on the sign of delta, so I use the same structure as sigma case, but I call it conditional LSB and USB
              // a conditional LSB can be both LSB and USB!
              // a conditional USB can be both LSB and USB!
              int sign = e>=f ? 1 : -1;
              int rowCosine = sign==1 ? 4*deltaAbs-2 : 4*deltaAbs  ;
              int rowSine   = sign==1 ? 4*deltaAbs-1 : 4*deltaAbs+1;

              // conditional LSB terms:
              // [ C_ijfI*cos(f*wc*t) - C_ijfQ*sin(f*wc*t) ] * [ + 0.5*(+V_jeIS-V_jeQC)*cos(e*wc*t-wm*t) - 0.5*(+V_jeIC+V_jeQS)*sin(e*wc*t-wm*t) ] * e*wc

              // conditional LSB cosine
              // + 0.5 *C_ijfI*cos(f*wc*t) * (+V_jeIS-V_jeQC)*cos(e*wc*t-wm*t) + 0.5 * C_ijfQ*sin(f*wc*t) * (+V_jeIC+V_jeQS)*sin(e*wc*t-wm*t) // all *e*wc
              // at delta:
              // + 0.25 * C_ijfI * (+V_jeIS-V_jeQC) * cos(delta*wc*t-wm*t) + 0.25 * C_ijfQ * (+V_jeIC+V_jeQS) * cos(delta*wc*t-wm*t) // all *e*wc
              // RHS has positive sign
              // + 0.25 * C_ijfI * V_jeIS * cos(delta*wc*t-wm*t) *e*wc
              setMatrixElement(harmonicSpaceMatrix_omegaC_1_->block(rowCosine, 4*e-1), i, j, +0.25*(e*omega_)*2*Cf_[i]->block(j)[2*f  ]);
              // - 0.25 * C_ijfI * V_jeQC * cos(delta*wc*t-wm*t) *e*wc
              setMatrixElement(harmonicSpaceMatrix_omegaC_1_->block(rowCosine, 4*e  ), i, j, -0.25*(e*omega_)*2*Cf_[i]->block(j)[2*f  ]);
              // + 0.25 * C_ijfQ * V_jeIC * cos(delta*wc*t-wm*t) *e*wc
              setMatrixElement(harmonicSpaceMatrix_omegaC_1_->block(rowCosine, 4*e-2), i, j, +0.25*(e*omega_)*2*Cf_[i]->block(j)[2*f+1]);
              // + 0.25 * C_ijfQ * V_jeQS * cos(delta*wc*t-wm*t) *e*wc
              setMatrixElement(harmonicSpaceMatrix_omegaC_1_->block(rowCosine, 4*e+1), i, j, +0.25*(e*omega_)*2*Cf_[i]->block(j)[2*f+1]);

              // conditional LSB sine
              // + 0.5 *C_ijfI*cos(f*wc*t) * (-V_jeIC-V_jeQS)*sin(e*wc*t-wm*t) + 0.5 * C_ijfQ*sin(f*wc*t) * (-V_jeIS+V_jeQC)*cos(e*wc*t-wm*t) // all *e*wc
              // at delta:
              // + 0.25 * C_ijfI * (-V_jeIC-V_jeQS) * sin(delta*wc*t-wm*t) + 0.25 * C_ijfQ * (+V_jeIS-V_jeQC) * sin(delta*wc*t-wm*t) // all *e*wc
              // RHS has negative sign
              // - 0.25 * C_ijfI * V_jeIC * sin(delta*wc*t-wm*t) *e*wc
              setMatrixElement(harmonicSpaceMatrix_omegaC_1_->block(rowSine, 4*e-2), i, j, +0.25*sign*(e*omega_)*2*Cf_[i]->block(j)[2*f  ]);
              // - 0.25 * C_ijfI * V_jeQS * sin(delta*wc*t-wm*t) *e*wc
              setMatrixElement(harmonicSpaceMatrix_omegaC_1_->block(rowSine, 4*e+1), i, j, +0.25*sign*(e*omega_)*2*Cf_[i]->block(j)[2*f  ]);
              // + 0.25 * C_ijfQ * V_jeIS * sin(delta*wc*t-wm*t) *e*wc
              setMatrixElement(harmonicSpaceMatrix_omegaC_1_->block(rowSine, 4*e-1), i, j, -0.25*sign*(e*omega_)*2*Cf_[i]->block(j)[2*f+1]);
              // - 0.25 * C_ijfQ * V_jeQC * sin(delta*wc*t-wm*t) *e*wc
              setMatrixElement(harmonicSpaceMatrix_omegaC_1_->block(rowSine, 4*e  ), i, j, +0.25*sign*(e*omega_)*2*Cf_[i]->block(j)[2*f+1]);

              // conditional USB terms
              // [ C_ijfI*cos(f*wc*t) - C_ijfQ*sin(f*wc*t) ] * [ + 0.5*(-V_jeIS-V_jeQC)*cos(e*wc*t+wm*t) - 0.5*(+V_jeIC-V_jeQS)*sin(e*wc*t+wm*t) ] * e*wc

              // conditional USB cosine
              // + 0.5 *C_ijfI*cos(f*wc*t) * (-V_jeIS-V_jeQC)*cos(e*wc*t+wm*t) + 0.5 * C_ijfQ*sin(f*wc*t) * (+V_jeIC-V_jeQS)*sin(e*wc*t+wm*t) // all *e*wc
              // at delta:
              // + 0.25 * C_ijfI * (-V_jeIS-V_jeQC) * cos(delta*wc*t+wm*t) + 0.25 * C_ijfQ * (+V_jeIC-V_jeQS) * cos(delta*wc*t+wm*t) // all *e*wc
              // RHS has positive sign
              // - 0.25 * C_ijfI * V_jeIS * cos(delta*wc*t+wm*t) *e*wc
              setMatrixElement(harmonicSpaceMatrix_omegaC_1_->block(rowCosine, 4*e-1), i, j, -0.25*(e*omega_)*2*Cf_[i]->block(j)[2*f  ]);
              // - 0.25 * C_ijfI * V_jeQC * cos(delta*wc*t+wm*t) *e*wc
              setMatrixElement(harmonicSpaceMatrix_omegaC_1_->block(rowCosine, 4*e  ), i, j, -0.25*(e*omega_)*2*Cf_[i]->block(j)[2*f  ]);
              // + 0.25 * C_ijfQ * V_jeIC * cos(delta*wc*t+wm*t) *e*wc
              setMatrixElement(harmonicSpaceMatrix_omegaC_1_->block(rowCosine, 4*e-2), i, j, +0.25*(e*omega_)*2*Cf_[i]->block(j)[2*f+1]);
              // - 0.25 * C_ijfQ * V_jeQS * cos(delta*wc*t+wm*t) *e*wc
              setMatrixElement(harmonicSpaceMatrix_omegaC_1_->block(rowCosine, 4*e+1), i, j, -0.25*(e*omega_)*2*Cf_[i]->block(j)[2*f+1]);

              // conditional USB sine
              // + 0.5 *C_ijfI*cos(f*wc*t) * (-V_jeIC+V_jeQS)*sin(e*wc*t+wm*t) + 0.5 * C_ijfQ*sin(f*wc*t) * (+V_jeIS+V_jeQC)*cos(e*wc*t+wm*t) // all *e*wc
              // at delta:
              // + 0.25 * C_ijfI * (-V_jeIC+V_jeQS) * sin(delta*wc*t+wm*t) + 0.25 * C_ijfQ * (-V_jeIS-V_jeQC) * sin(delta*wc*t+wm*t) // all *e*wc
              // RHS has negative sign
              // - 0.25 * C_ijfI * V_jeIC * sin(delta*wc*t+wm*t) *e*wc
              setMatrixElement(harmonicSpaceMatrix_omegaC_1_->block(rowSine, 4*e-2), i, j, +0.25*sign*(e*omega_)*2*Cf_[i]->block(j)[2*f  ]);
              // + 0.25 * C_ijfI * V_jeQS * sin(delta*wc*t+wm*t) *e*wc
              setMatrixElement(harmonicSpaceMatrix_omegaC_1_->block(rowSine, 4*e+1), i, j, -0.25*sign*(e*omega_)*2*Cf_[i]->block(j)[2*f  ]);
              // - 0.25 * C_ijfQ * V_jeIS * sin(delta*wc*t+wm*t) *e*wc
              setMatrixElement(harmonicSpaceMatrix_omegaC_1_->block(rowSine, 4*e-1), i, j, +0.25*sign*(e*omega_)*2*Cf_[i]->block(j)[2*f+1]);
              // - 0.25 * C_ijfQ * V_jeQC * sin(delta*wc*t+wm*t) *e*wc
              setMatrixElement(harmonicSpaceMatrix_omegaC_1_->block(rowSine, 4*e  ), i, j, +0.25*sign*(e*omega_)*2*Cf_[i]->block(j)[2*f+1]);
            }
          } // end of e loop
        } // end of f!=0
      } // end of f loop of Gf_
    } // end of column loop
  } // end of row loop
  return true;
}

//-----------------------------------------------------------------------------
// Function      : HBNOISE::updateHarmonicSpaceMatrix_C_2_()
// Purpose       : Updates C_2 matrix in harmonic space
// Scope         : private
// Creator       : Meysam Bahmanian
// Creation Date : 5/30/2025
//-----------------------------------------------------------------------------
// **************************************************
//           C Matrix PART 2: 
// **************************************************
// Reminder:
// C_ij(t) = C_ij0 + C_ijfI*cos(f*wc*t) - C_ijfQ*sin(f*wc*t) // sum over f
// v_j(t) has the form of
// v_j = V_jC*cos(wm*t) - V_jS*sin(wm*t) + (V_jeIC*cos(wm*t) - V_jeIS*sin(wm*t))*cos(e*wc*t) - (V_jeQC*cos(wm*t) - V_jeQS*sin(wm*t))*sin(e*wc*t) // sum over e
// PART 2 scales with wm and performs same-harmonic transformation of sine<=>cosine components. You can think of it as real<=>imaginary transformation.
// dv_j/dt |2 = - wm*V_jC*sin(wm*t) - wm*V_jS*cos(wm*t) + wm * ( - V_jeIC*sin(wm*t) - V_jeIS*cos(wm*t) )*cos(e*wc*t) - wm * ( - V_jeQC*sin(wm*t) - V_jeQS*cos(wm*t) )*sin(e*wc*t) 
// (1/wm)*dv_j/dt |2 = - V_jS*cos(wm*t) - V_jC*sin(wm*t) + ( - V_jeIS*cos(wm*t) - V_jeIC*sin(wm*t) )*cos(e*wc*t) - ( - V_jeQS*cos(wm*t) - V_jeQC*sin(wm*t) )*sin(e*wc*t) 
// summary:
// C_ij(t) = C_ij0 + C_ijfI*cos(f*wc*t) - C_ijfQ*sin(f*wc*t)            // sum over f
// b_i = B_iC*cos(wm*t) - B_iS*sin(wm*t) + B_ieLC*cos(e*wc*t-wm*t) - B_ieLS*sin(e*wc*t-wm*t) + B_ieUC*cos(e*wc*t+wm*t) - B_ieUS*sin(e*wc*t+wm*t) // sum over e
// v_j = V_jC*cos(wm*t) - V_jS*sin(wm*t) + (V_jeIC*cos(wm*t) - V_jeIS*sin(wm*t))*cos(e*wc*t) - (V_jeQC*cos(wm*t) - V_jeQS*sin(wm*t))*sin(e*wc*t) // sum over e
// v_j = V_jC*cos(wm*t) - V_jS*sin(wm*t) + 0.5*(+V_jeIC+V_jeQS)*cos(e*wc*t-wm*t) - 0.5*(-V_jeIS+V_jeQC)*sin(e*wc*t-wm*t) + 0.5*(+V_jeIC-V_jeQS)*cos(e*wc*t+wm*t) - 0.5*(+V_jeIS+V_jeQC)*sin(e*wc*t+wm*t)
// IQ form:
// (1/wm)*dv_j/dt |2 = - V_jS*cos(wm*t) - V_jC*sin(wm*t) + ( - V_jeIS*cos(wm*t) - V_jeIC*sin(wm*t) )*cos(e*wc*t) - ( - V_jeQS*cos(wm*t) - V_jeQC*sin(wm*t) )*sin(e*wc*t) 
// LSB/USB form:
// (1/wm)*dv_j/dt |2 = - V_jS*cos(wm*t) - V_jC*sin(wm*t) + 0.5*(-V_jeIS+V_jeQC)*cos(e*wc*t-wm*t) - 0.5*(-V_jeIC-V_jeQS)*sin(e*wc*t-wm*t) + 0.5*(-V_jeIS-V_jeQC)*cos(e*wc*t+wm*t) - 0.5*(+V_jeIC-V_jeQS)*sin(e*wc*t+wm*t)

bool HBNOISE::updateHarmonicSpaceMatrix_C_2_() 
{
  harmonicSpaceMatrix_C_2_->put( 0.0 );
  int numHarms = (size_-1)/2;
  int numBlocks = 2 + 4 * numHarms;
  int numRows = Cf_.size();
  for (int i=0; i<numRows; i++)
  { // selectring a row of G Matrix
    for (int j=0; j<numRows; j++)
    { // selecting a column of G Matrix
      for (int f=0; f<=numHarms; f++)
      { // selecting a frequency point of Gf_
        if ( Cf_[i]->block(j)[2*f] == 0.0 && Cf_[i]->block(j)[2*f+1] == 0.0 )
        { continue; }
        // now we fill the matrix
        // first the diagon of omega*C matrix
        if (f==0)
        { // linear transformation
          // C_ij0 * [ - V_jS*cos(wm*t) - V_jC*sin(wm*t) + 0.5*(-V_jeIS+V_jeQC)*cos(e*wc*t-wm*t) - 0.5*(-V_jeIC-V_jeQS)*sin(e*wc*t-wm*t) + 0.5*(-V_jeIS-V_jeQC)*cos(e*wc*t+wm*t) - 0.5*(+V_jeIC-V_jeQS)*sin(e*wc*t+wm*t) ]
          for (int e=0; e<=numHarms; e++){
            if (e==0) 
            { // our old friend: AC analysis!
              // C_ij0 * [ - V_jS*cos(wm*t) - V_jC*sin(wm*t) ]

              // - C_ij0 * V_jS * cos(wm*t)
              // RHS has positive sign
              setMatrixElement(harmonicSpaceMatrix_C_2_->block(0,1), i, j, -Cf_[i]->block(j)[0]);

              // - C_ij0 * V_jC * sin(wm*t)
              // RHS has negative sign
              setMatrixElement(harmonicSpaceMatrix_C_2_->block(1,0), i, j, +Cf_[i]->block(j)[0]);

            } else 
            { // translation of I/Q format to LSB/USB
              // C_ij0 * [ + 0.5*(-V_jeIS+V_jeQC)*cos(e*wc*t-wm*t) - 0.5*(-V_jeIC-V_jeQS)*sin(e*wc*t-wm*t) + 0.5*(-V_jeIS-V_jeQC)*cos(e*wc*t+wm*t) - 0.5*(+V_jeIC-V_jeQS)*sin(e*wc*t+wm*t) ]

              // LSB cosine
              // C_ij0 * [ + 0.5*(-V_jeIS+V_jeQC)*cos(e*wc*t-wm*t) ]
              // RHS has positive sign
              // - 0.5 * C_ij0 * V_jeIS * cos(e*wc*t-wm*t)
              setMatrixElement(harmonicSpaceMatrix_C_2_->block(4*e-2,4*e-1), i, j, -0.5*Cf_[i]->block(j)[0]);
              // + 0.5 * C_ij0 * V_jeQC * cos(e*wc*t-wm*t)
              setMatrixElement(harmonicSpaceMatrix_C_2_->block(4*e-2,4*e  ), i, j, +0.5*Cf_[i]->block(j)[0]);
              
              // LSB sine
              // C_ij0 * [ - 0.5*(-V_jeIC-V_jeQS)*sin(e*wc*t-wm*t) ]
              // RHS has negative sign
              // + 0.5 * C_ij0 * V_jeIC * sin(e*wc*t-wm*t)
              setMatrixElement(harmonicSpaceMatrix_C_2_->block(4*e-1,4*e-2), i, j, -0.5*Cf_[i]->block(j)[0]);
              // + 0.5 * C_ij0 * V_jeQS * sin(e*wc*t-wm*t)
              setMatrixElement(harmonicSpaceMatrix_C_2_->block(4*e-1,4*e+1), i, j, -0.5*Cf_[i]->block(j)[0]);

              // USB cosine
              // C_ij0 * [ + 0.5*(-V_jeIS-V_jeQC)*cos(e*wc*t+wm*t) ]
              // RHS has positive sign
              // - 0.5 * C_ij0 * V_jeIS * cos(e*wc*t+wm*t)
              setMatrixElement(harmonicSpaceMatrix_C_2_->block(4*e  ,4*e-1), i, j, -0.5*Cf_[i]->block(j)[0]);
              // - 0.5 * C_ij0 * V_jeQC * cos(e*wc*t+wm*t)
              setMatrixElement(harmonicSpaceMatrix_C_2_->block(4*e  ,4*e  ), i, j, -0.5*Cf_[i]->block(j)[0]);

              // USB sine
              // C_ij0 * [ - 0.5*(+V_jeIC-V_jeQS)*sin(e*wc*t+wm*t) ]
              // RHS has negative sign
              // - 0.5 * C_ij0 * V_jeIC * sin(e*wc*t+wm*t)
              setMatrixElement(harmonicSpaceMatrix_C_2_->block(4*e+1,4*e-2), i, j, +0.5*Cf_[i]->block(j)[0]);
              // + 0.5 * C_ij0 * V_jeQS * sin(e*wc*t+wm*t)
              setMatrixElement(harmonicSpaceMatrix_C_2_->block(4*e+1,4*e+1), i, j, -0.5*Cf_[i]->block(j)[0]);

            }
          }
        } else
        { // now the mixing parts (harmonic coupling)
          // [ C_ijfI*cos(f*wc*t) - C_ijfQ*sin(f*wc*t) ] * [ - V_jS*cos(wm*t) - V_jC*sin(wm*t) + 0.5*(-V_jeIS+V_jeQC)*cos(e*wc*t-wm*t) - 0.5*(-V_jeIC-V_jeQS)*sin(e*wc*t-wm*t) + 0.5*(-V_jeIS-V_jeQC)*cos(e*wc*t+wm*t) - 0.5*(+V_jeIC-V_jeQS)*sin(e*wc*t+wm*t) ]
          for (int e=0; e<=numHarms; e++)
          {
            if (e==0) 
            { // baseband modulation
              // [ C_ijfI*cos(f*wc*t) - C_ijfQ*sin(f*wc*t) ] * [ - V_jS*cos(wm*t) - V_jC*sin(wm*t) ]

              // LSB cosine
              // - C_ijfI*cos(f*wc*t) * V_jS*cos(wm*t) + C_ijfQ*sin(f*wc*t) * V_jC*sin(wm*t)
              // RHS has positive sign
              // - 0.5 * C_ijfI * V_jS * cos(f*wc*t-wm*t)
              setMatrixElement(harmonicSpaceMatrix_C_2_->block(4*f-2,1), i, j, -0.5*2*Cf_[i]->block(j)[2*f  ]);
              // + 0.5 * C_ijfQ * V_jC * cos(f*wc*t-wm*t)
              setMatrixElement(harmonicSpaceMatrix_C_2_->block(4*f-2,0), i, j, +0.5*2*Cf_[i]->block(j)[2*f+1]);

              // LSB sine
              // - C_ijfI*cos(f*wc*t) * V_jC*sin(wm*t) + C_ijfQ*sin(f*wc*t) * V_jS*cos(wm*t)
              // RHS has negative sign
              // + 0.5 * C_ijfI * V_jC * sin(f*wc*t-wm*t)
              setMatrixElement(harmonicSpaceMatrix_C_2_->block(4*f-1,0), i, j, -0.5*2*Cf_[i]->block(j)[2*f  ]);
              // + 0.5 * C_ijfQ * V_jS * sin(f*wc*t-wm*t)
              setMatrixElement(harmonicSpaceMatrix_C_2_->block(4*f-1,1), i, j, -0.5*2*Cf_[i]->block(j)[2*f+1]);

              // USB cosine
              // - C_ijfI*cos(f*wc*t) * V_jS*cos(wm*t) + C_ijfQ*sin(f*wc*t) * V_jC*sin(wm*t)
              // RHS has positive sign
              // - 0.5 * C_ijfI * V_jS * cos(f*wc*t+wm*t)
              setMatrixElement(harmonicSpaceMatrix_C_2_->block(4*f  ,1), i, j, -0.5*2*Cf_[i]->block(j)[2*f  ]);
              // - 0.5 * C_ijfQ * V_jC * cos(f*wc*t+wm*t)
              setMatrixElement(harmonicSpaceMatrix_C_2_->block(4*f  ,0), i, j, -0.5*2*Cf_[i]->block(j)[2*f+1]);

              // USB sine
              // - C_ijfI*cos(f*wc*t) * V_jC*sin(wm*t) + C_ijfQ*sin(f*wc*t) * V_jS*cos(wm*t)
              // RHS has negative sign
              // - 0.5 * C_ijfI * V_jC * sin(f*wc*t+wm*t)
              setMatrixElement(harmonicSpaceMatrix_C_2_->block(4*f+1,0), i, j, +0.5*2*Cf_[i]->block(j)[2*f  ]);
              // + 0.5 * C_ijfQ * V_jS * sin(f*wc*t+wm*t)
              setMatrixElement(harmonicSpaceMatrix_C_2_->block(4*f+1,1), i, j, -0.5*2*Cf_[i]->block(j)[2*f+1]);
            } else 
            { // now f>0 and e>0
              // [ C_ijfI*cos(f*wc*t) - C_ijfQ*sin(f*wc*t) ] * [ + 0.5*(-V_jeIS+V_jeQC)*cos(e*wc*t-wm*t) - 0.5*(-V_jeIC-V_jeQS)*sin(e*wc*t-wm*t) + 0.5*(-V_jeIS-V_jeQC)*cos(e*wc*t+wm*t) - 0.5*(+V_jeIC-V_jeQS)*sin(e*wc*t+wm*t) ]

              int sigma = e+f;

              // case sigma
              if (sigma<=numHarms) 
              { // f+e should not be larger than numHarms, otherwise ignore it
                // in finite-harmonics space the system still shows nonlineary and some mixing products have to be ignored

                // LSB terms:
                // [ C_ijfI*cos(f*wc*t) - C_ijfQ*sin(f*wc*t) ] * [ + 0.5*(-V_jeIS+V_jeQC)*cos(e*wc*t-wm*t) - 0.5*(-V_jeIC-V_jeQS)*sin(e*wc*t-wm*t) ]

                // LSB cosine
                // + 0.5 * C_ijfI*cos(f*wc*t) * (-V_jeIS+V_jeQC)*cos(e*wc*t-wm*t) + 0.5 * C_ijfQ*sin(f*wc*t) * (-V_jeIC-V_jeQS)*sin(e*wc*t-wm*t)
                // at sigma:
                // + 0.25 * C_ijfI * (-V_jeIS+V_jeQC) * cos(sigma*wc*t-wm*t) + 0.25 * C_ijfQ * (+V_jeIC+V_jeQS) * cos(sigma*wc*t-wm*t)
                // RHS has positive sign
                // - 0.25 * C_ijfI * V_jeIS * cos(sigma*wc*t-wm*t)
                setMatrixElement(harmonicSpaceMatrix_C_2_->block(4*sigma-2, 4*e-1), i, j, -0.25*2*Cf_[i]->block(j)[2*f  ]);
                // + 0.25 * C_ijfI * V_jeQC * cos(sigma*wc*t-wm*t)
                setMatrixElement(harmonicSpaceMatrix_C_2_->block(4*sigma-2, 4*e  ), i, j, +0.25*2*Cf_[i]->block(j)[2*f  ]);
                // + 0.25 * C_ijfQ * V_jeIC * cos(sigma*wc*t-wm*t)
                setMatrixElement(harmonicSpaceMatrix_C_2_->block(4*sigma-2, 4*e-2), i, j, +0.25*2*Cf_[i]->block(j)[2*f+1]);
                // + 0.25 * C_ijfQ * V_jeQS * cos(sigma*wc*t-wm*t)
                setMatrixElement(harmonicSpaceMatrix_C_2_->block(4*sigma-2, 4*e+1), i, j, +0.25*2*Cf_[i]->block(j)[2*f+1]);

                // LSB sine
                // + 0.5 * C_ijfI*cos(f*wc*t) * (+V_jeIC+V_jeQS)*sin(e*wc*t-wm*t) + 0.5 * C_ijfQ*sin(f*wc*t) * (+V_jeIS-V_jeQC)*cos(e*wc*t-wm*t)
                // at sigma:
                // + 0.25 * C_ijfI * (+V_jeIC+V_jeQS)*sin(sigma*wc*t-wm*t) + 0.25 * C_ijfQ * (+V_jeIS-V_jeQC)*sin(sigma*wc*t-wm*t)
                // RHS has negative sign
                // + 0.25 * C_ijfI * V_jeIC * sin(sigma*wc*t-wm*t)
                setMatrixElement(harmonicSpaceMatrix_C_2_->block(4*sigma-1, 4*e-2), i, j, -0.25*2*Cf_[i]->block(j)[2*f  ]);
                // + 0.25 * C_ijfI * V_jeQS * sin(sigma*wc*t-wm*t)
                setMatrixElement(harmonicSpaceMatrix_C_2_->block(4*sigma-1, 4*e+1), i, j, -0.25*2*Cf_[i]->block(j)[2*f  ]);
                // + 0.25 * C_ijfQ * V_jeIS * sin(sigma*wc*t-wm*t)
                setMatrixElement(harmonicSpaceMatrix_C_2_->block(4*sigma-1, 4*e-1), i, j, -0.25*2*Cf_[i]->block(j)[2*f+1]);
                // - 0.25 * C_ijfQ * V_jeQC * sin(sigma*wc*t-wm*t)
                setMatrixElement(harmonicSpaceMatrix_C_2_->block(4*sigma-1, 4*e  ), i, j, +0.25*2*Cf_[i]->block(j)[2*f+1]);


                // USB terms:
                // [ C_ijfI*cos(f*wc*t) - C_ijfQ*sin(f*wc*t) ] * [ + 0.5*(-V_jeIS-V_jeQC)*cos(e*wc*t+wm*t) - 0.5*(+V_jeIC-V_jeQS)*sin(e*wc*t+wm*t) ]

                // USB cosine
                // + 0.5 * C_ijfI*cos(f*wc*t) * (-V_jeIS-V_jeQC)*cos(e*wc*t+wm*t) + 0.5 * C_ijfQ*sin(f*wc*t) * (+V_jeIC-V_jeQS)*sin(e*wc*t+wm*t)
                // at sigma:
                // + 0.25 * C_ijfI * (-V_jeIS-V_jeQC) * cos(sigma*wc*t+wm*t) + 0.25 * C_ijfQ * (-V_jeIC+V_jeQS) * cos(sigma*wc*t+wm*t)
                // RHS has positive sign
                // - 0.25 * C_ijfI * V_jeIS * cos(sigma*wc*t+wm*t)
                setMatrixElement(harmonicSpaceMatrix_C_2_->block(4*sigma  , 4*e-1), i, j, -0.25*2*Cf_[i]->block(j)[2*f  ]);
                // - 0.25 * C_ijfI * V_jeQC * cos(sigma*wc*t+wm*t)
                setMatrixElement(harmonicSpaceMatrix_C_2_->block(4*sigma  , 4*e  ), i, j, -0.25*2*Cf_[i]->block(j)[2*f  ]);
                // - 0.25 * C_ijfQ * V_jeIC * sin(sigma*wc*t+wm*t)
                setMatrixElement(harmonicSpaceMatrix_C_2_->block(4*sigma  , 4*e-2), i, j, -0.25*2*Cf_[i]->block(j)[2*f+1]);
                // + 0.25 * C_ijfQ * V_jeQS * sin(sigma*wc*t+wm*t)
                setMatrixElement(harmonicSpaceMatrix_C_2_->block(4*sigma  , 4*e+1), i, j, +0.25*2*Cf_[i]->block(j)[2*f+1]);

                // USB sine
                // + 0.5 * C_ijfI*cos(f*wc*t) * (-V_jeIC+V_jeQS)*sin(e*wc*t+wm*t) + 0.5 * C_ijfQ*sin(f*wc*t) * (+V_jeIS+V_jeQC)*cos(e*wc*t+wm*t)
                // at sigma:
                // + 0.25 * C_ijfI * (-V_jeIC+V_jeQS)*sin(sigma*wc*t+wm*t) + 0.25 * C_ijfQ * (+V_jeIS+V_jeQC)*sin(sigma*wc*t+wm*t)
                // RHS has negative sign
                // - 0.25 * C_ijfI * V_jeIC * sin(sigma*wc*t+wm*t)
                setMatrixElement(harmonicSpaceMatrix_C_2_->block(4*sigma+1, 4*e-2), i, j, +0.25*2*Cf_[i]->block(j)[2*f  ]);
                // + 0.25 * C_ijfI * V_jeQS * sin(sigma*wc*t+wm*t)
                setMatrixElement(harmonicSpaceMatrix_C_2_->block(4*sigma+1, 4*e+1), i, j, -0.25*2*Cf_[i]->block(j)[2*f  ]);
                // + 0.25 * C_ijfQ * V_jeIS * sin(sigma*wc*t+wm*t)
                setMatrixElement(harmonicSpaceMatrix_C_2_->block(4*sigma+1, 4*e-1), i, j, -0.25*2*Cf_[i]->block(j)[2*f+1]);
                // + 0.25 * C_ijfQ * V_jeQC * sin(sigma*wc*t+wm*t)
                setMatrixElement(harmonicSpaceMatrix_C_2_->block(4*sigma+1, 4*e  ), i, j, -0.25*2*Cf_[i]->block(j)[2*f+1]);
              }

              // int delta = e-f; this line is just for the sake of understanding the code. Equations are with reference to delta, not deltaAbs.
              int deltaAbs = std::abs(e-f);

              // case delta
              // now I have to deal with the index issue! if delta==0 then the indices become negative!
              if (deltaAbs==0)
              { // transformation from harmonics to baseband
                // I use I/Q form of (1/wm)*dv_j/dt here, instead of LSB/USB forms. This will be easier for the baseband transformation.
                // [ C_ijfI*cos(f*wc*t) - C_ijfQ*sin(f*wc*t) ] * [ ( - V_jeIS*cos(wm*t) - V_jeIC*sin(wm*t) )*cos(e*wc*t) - ( - V_jeQS*cos(wm*t) - V_jeQC*sin(wm*t) )*sin(e*wc*t) ]

                // baseband cosine term:
                // - C_ijfI * cos(f*wc*t) * V_jeIS * cos(wm*t) * cos(e*wc*t) - C_ijfQ * sin(f*wc*t) * V_jeQS * cos(wm*t) * sin(e*wc*t)
                // at delta=0
                // - 0.5 * C_ijfI * V_jeIS * cos(wm*t) - 0.5 * C_ijfQ * V_jeQS * cos(wm*t)
                // RHS has positive sign
                // - 0.5 * C_ijfI * V_jeIS * cos(wm*t)
                setMatrixElement(harmonicSpaceMatrix_C_2_->block(0, 4*e-1), i, j, -0.5*2*Cf_[i]->block(j)[2*f  ]);
                // - 0.5 * C_ijfQ * V_jeQS * cos(wm*t)
                setMatrixElement(harmonicSpaceMatrix_C_2_->block(0, 4*e+1), i, j, -0.5*2*Cf_[i]->block(j)[2*f+1]); 
                
                // baseband sine term:
                // - C_ijfI * cos(f*wc*t) * V_jeIC * sin(wm*t) * cos(e*wc*t) - C_ijfQ * sin(f*wc*t) * V_jeQC * sin(wm*t) * sin(e*wc*t)
                // at delta=0
                // - 0.5 * C_ijfI * V_jeIC * sin(wm*t) - 0.5 * C_ijfQ * V_jeQC * sin(wm*t)
                // RHS has negative sign
                // - 0.5 * C_ijfI * V_jeIC * sin(wm*t)
                setMatrixElement(harmonicSpaceMatrix_C_2_->block(1, 4*e-2), i, j, +0.5*2*Cf_[i]->block(j)[2*f  ]);
                // - 0.5 * C_ijfQ * V_jeQC * sin(wm*t)
                setMatrixElement(harmonicSpaceMatrix_C_2_->block(1, 4*e  ), i, j, +0.5*2*Cf_[i]->block(j)[2*f+1]); 

              }
              else
              {
                // here LSB and USB depend on the sign of delta, so I use the same structure as sigma case, but I call it conditional LSB and USB
                // a conditional LSB can be both LSB and USB!
                // a conditional USB can be both LSB and USB!
                int sign = e>=f ? 1 : -1;
                int rowCosine = sign==1 ? 4*deltaAbs-2 : 4*deltaAbs  ;
                int rowSine   = sign==1 ? 4*deltaAbs-1 : 4*deltaAbs+1;

                // conditional LSB terms:
                // [ C_ijfI*cos(f*wc*t) - C_ijfQ*sin(f*wc*t) ] * [ + 0.5*(-V_jeIS+V_jeQC)*cos(e*wc*t-wm*t) - 0.5*(-V_jeIC-V_jeQS)*sin(e*wc*t-wm*t) ]

                // conditional LSB cosine
                // + 0.5 * C_ijfI*cos(f*wc*t) * (-V_jeIS+V_jeQC)*cos(e*wc*t-wm*t) + 0.5 * C_ijfQ*sin(f*wc*t) * (-V_jeIC-V_jeQS)*sin(e*wc*t-wm*t)
                // at delta:
                // + 0.25 * C_ijfI * (-V_jeIS+V_jeQC) * cos(delta*wc*t-wm*t) + 0.25 * C_ijfQ * (-V_jeIC-V_jeQS) * cos(delta*wc*t-wm*t)
                // RHS has positive sign
                // - 0.25 * C_ijfI * V_jeIS * cos(delta*wc*t-wm*t)
                setMatrixElement(harmonicSpaceMatrix_C_2_->block(rowCosine, 4*e-1), i, j, -0.25*2*Cf_[i]->block(j)[2*f  ]);
                // + 0.25 * C_ijfI * V_jeQC * cos(delta*wc*t-wm*t)
                setMatrixElement(harmonicSpaceMatrix_C_2_->block(rowCosine, 4*e  ), i, j, +0.25*2*Cf_[i]->block(j)[2*f  ]);
                // - 0.25 * C_ijfQ * V_jeIC * cos(delta*wc*t-wm*t)
                setMatrixElement(harmonicSpaceMatrix_C_2_->block(rowCosine, 4*e-2), i, j, -0.25*2*Cf_[i]->block(j)[2*f+1]);
                // - 0.25 * C_ijfQ * V_jeQS * cos(delta*wc*t-wm*t)
                setMatrixElement(harmonicSpaceMatrix_C_2_->block(rowCosine, 4*e+1), i, j, -0.25*2*Cf_[i]->block(j)[2*f+1]);

                // conditional LSB sine
                // + 0.5 * C_ijfI*cos(f*wc*t) * (+V_jeIC+V_jeQS)*sin(e*wc*t-wm*t) + 0.5 * C_ijfQ*sin(f*wc*t) * (+V_jeIS-V_jeQC)*cos(e*wc*t-wm*t)
                // at delta:
                // + 0.25 * C_ijfI * (+V_jeIC+V_jeQS)*sin(delta*wc*t-wm*t) + 0.25 * C_ijfQ * (-V_jeIS+V_jeQC)*sin(delta*wc*t-wm*t)
                // RHS has negative sign
                // + 0.25 * C_ijfI * V_jeIC * sin(delta*wc*t-wm*t)
                setMatrixElement(harmonicSpaceMatrix_C_2_->block(rowSine, 4*e-2), i, j, -0.25*sign*2*Cf_[i]->block(j)[2*f  ]);
                // + 0.25 * C_ijfI * V_jeQS * sin(delta*wc*t-wm*t)
                setMatrixElement(harmonicSpaceMatrix_C_2_->block(rowSine, 4*e+1), i, j, -0.25*sign*2*Cf_[i]->block(j)[2*f  ]);
                // - 0.25 * C_ijfQ * V_jeIS * sin(delta*wc*t-wm*t)
                setMatrixElement(harmonicSpaceMatrix_C_2_->block(rowSine, 4*e-1), i, j, +0.25*sign*2*Cf_[i]->block(j)[2*f+1]);
                // + 0.25 * C_ijfQ * V_jeQC * sin(delta*wc*t-wm*t)
                setMatrixElement(harmonicSpaceMatrix_C_2_->block(rowSine, 4*e  ), i, j, -0.25*sign*2*Cf_[i]->block(j)[2*f+1]);

                // conditional USB terms:
                // [ C_ijfI*cos(f*wc*t) - C_ijfQ*sin(f*wc*t) ] * [ + 0.5*(-V_jeIS-V_jeQC)*cos(e*wc*t+wm*t) - 0.5*(+V_jeIC-V_jeQS)*sin(e*wc*t+wm*t) ]

                // conditional USB cosine
                // + 0.5 * C_ijfI*cos(f*wc*t) * (-V_jeIS-V_jeQC)*cos(e*wc*t+wm*t) + 0.5 * C_ijfQ*sin(f*wc*t) * (+V_jeIC-V_jeQS)*sin(e*wc*t+wm*t)
                // at delta:
                // + 0.25 * C_ijfI * (-V_jeIS-V_jeQC) * cos(delta*wc*t+wm*t) + 0.25 * C_ijfQ * (+V_jeIC-V_jeQS) * cos(delta*wc*t+wm*t)
                // RHS has positive sign
                // - 0.25 * C_ijfI * V_jeIS * cos(delta*wc*t+wm*t)
                setMatrixElement(harmonicSpaceMatrix_C_2_->block(rowCosine, 4*e-1), i, j, -0.25*2*Cf_[i]->block(j)[2*f  ]);
                // - 0.25 * C_ijfI * V_jeQC * cos(delta*wc*t+wm*t)
                setMatrixElement(harmonicSpaceMatrix_C_2_->block(rowCosine, 4*e  ), i, j, -0.25*2*Cf_[i]->block(j)[2*f  ]);
                // + 0.25 * C_ijfQ * V_jeIC * cos(delta*wc*t+wm*t)
                setMatrixElement(harmonicSpaceMatrix_C_2_->block(rowCosine, 4*e-2), i, j, +0.25*2*Cf_[i]->block(j)[2*f+1]);
                // - 0.25 * C_ijfQ * V_jeQS * cos(delta*wc*t+wm*t)
                setMatrixElement(harmonicSpaceMatrix_C_2_->block(rowCosine, 4*e+1), i, j, -0.25*2*Cf_[i]->block(j)[2*f+1]);

                // conditional USB sine
                // + 0.5 * C_ijfI*cos(f*wc*t) * (-V_jeIC+V_jeQS)*sin(e*wc*t+wm*t) + 0.5 * C_ijfQ*sin(f*wc*t) * (+V_jeIS+V_jeQC)*cos(e*wc*t+wm*t)
                // at delta:
                // + 0.25 * C_ijfI * (-V_jeIC+V_jeQS)*sin(delta*wc*t+wm*t) + 0.25 * C_ijfQ * (-V_jeIS-V_jeQC)*sin(delta*wc*t+wm*t)
                // RHS has negative sign
                // - 0.25 * C_ijfI * V_jeIC * sin(delta*wc*t+wm*t)
                setMatrixElement(harmonicSpaceMatrix_C_2_->block(rowSine, 4*e-2), i, j, +0.25*sign*2*Cf_[i]->block(j)[2*f  ]);
                // + 0.25 * C_ijfI * V_jeQS * sin(delta*wc*t+wm*t)
                setMatrixElement(harmonicSpaceMatrix_C_2_->block(rowSine, 4*e+1), i, j, -0.25*sign*2*Cf_[i]->block(j)[2*f  ]);
                // - 0.25 * C_ijfQ * V_jeIS * sin(delta*wc*t+wm*t)
                setMatrixElement(harmonicSpaceMatrix_C_2_->block(rowSine, 4*e-1), i, j, +0.25*sign*2*Cf_[i]->block(j)[2*f+1]);
                // - 0.25 * C_ijfQ * V_jeQC * sin(delta*wc*t+wm*t)
                setMatrixElement(harmonicSpaceMatrix_C_2_->block(rowSine, 4*e  ), i, j, +0.25*sign*2*Cf_[i]->block(j)[2*f+1]);
              }
            } // end of f>0 and e>0
          } // end of e loop
        } // end of f!=0
      } // end of f loop of Cf_
    } // end of column loop
  } // end of row loop
  return true;
}

//-----------------------------------------------------------------------------
// Function      : HBNOISE::updateLinearTimeVariantSystem_C_and_G_()
// Purpose       : Updates C and G matrices for all time points in HB
//                 Similar to NOISE::updateACLinearSystem_C_and_G_()
//                 but handles multiple time points for HBNOISE
// Special Notes : The set of functions updateLinearTimeVariantSystem_C_and_G_() and 
//                 createHarmonicSpaceLinearSystem_() are computationally suboptimal
//                 and they apply Fourier tranform to a sparse matrix. They most likely
//                 will also fail to handle Touchstone files.
//                 But they are easier to understand for defining our Gold Standard.
//                 They will be replaced later with a more efficient approach,
//                 similar to how HBLoader::loadDAEVectors() works
// Scope         : private
// Creator       : Meysam Bahmanian
// Creation Date : 5/20/2025
//-----------------------------------------------------------------------------
bool HBNOISE::updateLinearTimeVariantSystem_C_and_G_()
{
  // Get the frequency domain HB solution from DataStore
  // b means to "block"
  Linear::Vector *Xf = analysisManager_.getDataStore()->currSolutionPtr;
  Linear::BlockVector & bXf = *dynamic_cast<Linear::BlockVector*>(Xf);

  // first store the output node results
  {
    double v1r = 0.0;
    double v1i = 0.0;
    double v2r = 0.0;
    double v2i = 0.0;
    int lid;
    if (outputVarGIDs_.size()>0)
    {
      if (outputVarGIDs_[0] > -1)
      {
        // TODO: correct for parallel
        // v1r = bXf.getElementByGlobalIndex(outputVarGIDs_[0], 2*harmonicNumber_  );
        // v1i = bXf.getElementByGlobalIndex(outputVarGIDs_[0], 2*harmonicNumber_+1);
        v1r = bXf.block(outputVarGIDs_[0])[2*harmonicNumber_  ];
        v1i = bXf.block(outputVarGIDs_[0])[2*harmonicNumber_+1];
      }
    }

    if (outputVarGIDs_.size()>1)
    {
      if (outputVarGIDs_[0] > -1)
      {
        // v2r = bXf.getElementByGlobalIndex(outputVarGIDs_[1], 2*harmonicNumber_  );
        // v2i = bXf.getElementByGlobalIndex(outputVarGIDs_[1], 2*harmonicNumber_+1);
        v2r = bXf.block(outputVarGIDs_[1])[2*harmonicNumber_  ];
        v2i = bXf.block(outputVarGIDs_[1])[2*harmonicNumber_+1];
      }
    }
    outputValReal_ = v1r - v2r;
    outputValImag_ = v1i - v2i;
    outputValSqr_ = outputValReal_*outputValReal_ + outputValImag_*outputValImag_;
    outputValCosPhi_ = outputValReal_ / sqrt(outputValSqr_);
    outputValSinPhi_ = outputValImag_ / sqrt(outputValSqr_);
  }

  Teuchos::RCP<Linear::BlockVector> bXtPtr_ = hbBuilderPtr_->createTimeDomainBlockVector();
  bXtPtr_->putScalar(0.0);

  hbLoaderPtr_->permutedIFT(bXf, &*bXtPtr_);

  Linear::BlockVector & bX = *bXtPtr_;

  int BlockCount = bX.blockCount(); // number of time points
  int BlockSize = bX.blockSize(); // number of GIDs

  if (DEBUG_HBNOISE)
  {
    for (int i = 0; i < BlockCount; ++i)
    {
      Xyce::dout() << "Solution time domain, block (" << i << "): each block is a time point" << std::endl;
      bX.block(i).print( Xyce::dout() );
      Xyce::dout() << std::endl;
    }
    for (int i = 0; i < BlockSize; ++i)
    {
      Xyce::dout() << "Solution frequency domain, block (" << i << "): each block is a node" << std::endl;
      bXf.block(i).print( Xyce::dout() );
      Xyce::dout() << std::endl;
    }
  }

  // Solutions:
  Linear::Vector * currSolutionPtr = builderPtr_->createVector();

  for (int i=0; i<BlockSize; i++){
    Ct_.push_back(hbBuilderPtr_->createTimeDomainBlockVector());
    Gt_.push_back(hbBuilderPtr_->createTimeDomainBlockVector());
    Ct_[i]->putScalar(0.0);
    Gt_[i]->putScalar(0.0);
  }

  Linear::Vector * tmpQ = builderPtr_->createVector();
  Linear::Vector * tmpF = builderPtr_->createVector();
  Linear::Vector * tmpB = builderPtr_->createVector();
  Linear::Matrix * tmpC;
  Linear::Matrix * tmpG;
  // Linear::Matrix * dQdxMatrixPtr = builderPtr_->createMatrix();
  // Linear::Matrix * dFdxMatrixPtr = builderPtr_->createMatrix();
  Teuchos::RCP<Linear::Matrix> dQdxMatrixPtr = rcp(builderPtr_->createMatrix());
  Teuchos::RCP<Linear::Matrix> dFdxMatrixPtr = rcp(builderPtr_->createMatrix());

  // now we store dFdx and dQdx matrices
  for (int i = 0; i < BlockCount; ++i)
  {
    deviceManager_.setFastTime(hbAnalysis_->fastTimes_[i]);
    loader_.updateSources();  // this is here to handle "fast" sources.
    // *appVecPtr_ = bX.block(i);
    *currSolutionPtr = bXtPtr_->block(i);
    analysisManager_.getDataStore()->daeQVectorPtr->putScalar(0.0);
    analysisManager_.getDataStore()->daeFVectorPtr->putScalar(0.0);

    analysisManager_.getDataStore()->dFdxdVpVectorPtr->putScalar(0.0);
    analysisManager_.getDataStore()->dQdxdVpVectorPtr->putScalar(0.0);

    analysisManager_.getDataStore()->daeQVectorPtr->putScalar(0.0);
    analysisManager_.getDataStore()->daeFVectorPtr->putScalar(0.0);
    analysisManager_.getDataStore()->daeBVectorPtr->putScalar(0.0);

    dQdxMatrixPtr->put(0.0);
    dFdxMatrixPtr->put(0.0);

    loader_.updateState(
                (currSolutionPtr),
                (currSolutionPtr),
                (currSolutionPtr),
                (analysisManager_.getDataStore()->nextStatePtr),
                (analysisManager_.getDataStore()->currStatePtr),
                (analysisManager_.getDataStore()->lastStatePtr),
                (analysisManager_.getDataStore()->nextStorePtr),
                (analysisManager_.getDataStore()->currStorePtr),
                (analysisManager_.getDataStore()->lastStorePtr),
                Xyce::Device::NONLINEAR_FREQ
                );

    loader_.loadDAEVectors(
                (currSolutionPtr),
                (currSolutionPtr),
                (currSolutionPtr),
                (analysisManager_.getDataStore()->nextStatePtr),
                (analysisManager_.getDataStore()->currStatePtr),
                (analysisManager_.getDataStore()->lastStatePtr),
                (analysisManager_.getDataStore()->nextStateDerivPtr),
                (analysisManager_.getDataStore()->nextStorePtr),
                (analysisManager_.getDataStore()->currStorePtr),
                (analysisManager_.getDataStore()->lastStorePtr),
                (analysisManager_.getDataStore()->nextLeadCurrentPtr),
                (analysisManager_.getDataStore()->nextLeadCurrentQPtr),
                (analysisManager_.getDataStore()->nextLeadDeltaVPtr),
                (analysisManager_.getDataStore()->daeQVectorPtr),
                (analysisManager_.getDataStore()->daeFVectorPtr),
                (analysisManager_.getDataStore()->daeBVectorPtr),
                (analysisManager_.getDataStore()->dFdxdVpVectorPtr),
                (analysisManager_.getDataStore()->dQdxdVpVectorPtr),
                Xyce::Device::NONLINEAR_FREQ
                );

    loader_.loadBVectorsforSources();
    analysisManager_.getDataStore()->daeBVectorPtr->fillComplete();

    loader_.loadDAEMatrices(
                (currSolutionPtr),
                analysisManager_.getDataStore()->nextStatePtr, 
                analysisManager_.getDataStore()->nextStateDerivPtr, 
                analysisManager_.getDataStore()->nextStorePtr, 
                &*dQdxMatrixPtr,  
                &*dFdxMatrixPtr,
                Xyce::Device::NONLINEAR_FREQ
                );

    loader_.loadDAEMatrices(
                (currSolutionPtr),
                analysisManager_.getDataStore()->nextStatePtr, 
                analysisManager_.getDataStore()->nextStateDerivPtr, 
                analysisManager_.getDataStore()->nextStorePtr, 
                &*dQdxMatrixPtr,  
                &*dFdxMatrixPtr,
                Xyce::Device::LINEAR_FREQ
                );
    
    tmpG = &*dFdxMatrixPtr;
    tmpC = &*dQdxMatrixPtr;

    int numEntries;
    std::vector<double> coeffs(BlockSize); 
    std::vector<int> colIndices(BlockSize);

    for (int j=0; j<BlockSize; j++) {
      tmpC->getLocalRowCopy(j, BlockSize, numEntries, coeffs.data(), colIndices.data());
      for (int k = 0; k < numEntries; k++) {
        Ct_[j]->block(i)[colIndices[k]] = coeffs[k];
      }
    }

    for (int j=0; j<BlockSize; j++) {
      tmpG->getLocalRowCopy(j, BlockSize, numEntries, coeffs.data(), colIndices.data());
      for (int k = 0; k < numEntries; k++) {
        Gt_[j]->block(i)[colIndices[k]] = coeffs[k];
      }
    }

    if (DEBUG_HBNOISE)
    {
      // print conductance matrix
      Xyce::dout() << "dFdxMatrixPtr time point of " << i << ":" << std::endl;
      dFdxMatrixPtr->print( Xyce::dout() );
      Xyce::dout() << std::endl;

      // print capacitance matrix
      // Xyce::dout() << "dQdxMatrixPtr time point of " << i << ":" << std::endl;
      // dQdxMatrixPtr->print( Xyce::dout() );
      // Xyce::dout() << std::endl;
    }
  }

  if (DEBUG_HBNOISE)
  {
    Xyce::dout() << "Reporting Gt_ Matrices, each array element is a row of the matrix, each block is a time point" << std::endl;
    for (int i=0; i<BlockSize; i++){
      Xyce::dout() << "Gt_[" << i << "]: " << std::endl;
      Gt_[i]->print(Xyce::dout());
      Xyce::dout() << std::endl;
    }
    // Xyce::dout() << "Reporting Ct_ Matrices, each array element is a row of the matrix, each block is a time point" << std::endl;
    // for (int i=0; i<BlockSize; i++){
    //   Xyce::dout() << "Ct_[" << i << "]: " << std::endl;
    //   Ct_[i]->print(Xyce::dout());
    //   Xyce::dout() << std::endl;
    // }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : HBNOISE::setMatrixElement
// Purpose       : Set a matrix element
// Special Notes : This helper function should be moved to a more appropriate place
// Scope         : private
// Creator       : Meysam Bahmanian
// Creation Date : 6/4/2025
//-----------------------------------------------------------------------------
void HBNOISE::setMatrixElement(Linear::Matrix& mat, int row, int col, double value) {
    static std::vector<double> values{0.0};
    static std::vector<int> indices{0};
    values[0] = value;
    indices[0] = col;
    mat.putLocalRow(row, 1, values.data(), indices.data());
}

//-----------------------------------------------------------------------------
// Function      : HBNOISE::doProcessSuccessfulStep
// Purpose       : Process a successful step
// Special Notes :
// Scope         : protected
// Creator       : Meysam Bahmanian
// Creation Date : 5/13/2025
//-----------------------------------------------------------------------------
bool HBNOISE::doProcessSuccessfulStep()
{
  outputManagerAdapter_.outputHBNoise (
    currentFreq_, 
    fOffsetStart_, fOffsetStop_, 
    harmonicSpaceXI_->block(0),  // I'm not sure what these paramenters in NOISE class do. They are needed for get values function, So I'm gonna pass them for now
    harmonicSpaceXI_-> block(1), // I'm not sure what these paramenters in NOISE class do. They are needed for get values function, So I'm gonna pass them for now
    totalAMNoiseDens_, 
    totalPMNoiseDens_, 
    noiseDataVecI_, 
    noiseDataVecQ_,
    noiseDataVecVecI_,
    noiseDataVecVecQ_);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : HBNOISE::doProcessFailedStep
// Purpose       : Process a failed step
// Special Notes :
// Scope         : protected
// Creator       : Meysam Bahmanian
// Creation Date : 5/13/2025
//-----------------------------------------------------------------------------
bool HBNOISE::doProcessFailedStep()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : HBNOISE::doFinish
// Purpose       : Finish the analysis
// Special Notes :
// Scope         : protected
// Creator       : Meysam Bahmanian
// Creation Date : 5/13/2025
//-----------------------------------------------------------------------------
bool HBNOISE::doFinish()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : HBNOISE::doHandlePredictor
// Purpose       : Handle predictor step
// Special Notes :
// Scope         : protected
// Creator       : Meysam Bahmanian
// Creation Date : 5/13/2025
//-----------------------------------------------------------------------------
bool HBNOISE::doHandlePredictor()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : HBNOISE::processSuccessfulDCOP
// Purpose       : Process successful DC operating point
// Special Notes :
// Scope         : public
// Creator       : Meysam Bahmanian
// Creation Date : 5/13/2025
//-----------------------------------------------------------------------------
bool HBNOISE::processSuccessfulDCOP()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : HBNOISE::processFailedDCOP
// Purpose       : Process failed DC operating point
// Special Notes :
// Scope         : public
// Creator       : Meysam Bahmanian
// Creation Date : 5/13/2025
//-----------------------------------------------------------------------------
bool HBNOISE::processFailedDCOP()
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : HBNOISE::getTIAParams
// Purpose       : Get Time Integration Analysis Parameters (non-const version)
// Special Notes :
// Scope         : public
// Creator       : Meysam Bahmanian
// Creation Date : 5/13/2025
//-----------------------------------------------------------------------------
TimeIntg::TIAParams &HBNOISE::getTIAParams()
{
  static TimeIntg::TIAParams dummyTIAParams;
  return dummyTIAParams;
}

//-----------------------------------------------------------------------------
// Function      : HBNOISE::getTIAParams
// Purpose       : Get Time Integration Analysis Parameters (const version)
// Special Notes :
// Scope         : public
// Creator       : Meysam Bahmanian
// Creation Date : 5/13/2025
//-----------------------------------------------------------------------------
const TimeIntg::TIAParams &HBNOISE::getTIAParams() const
{
  static TimeIntg::TIAParams dummyTIAParams;
  return dummyTIAParams;
}

//-----------------------------------------------------------------------------
// Function      : HBNOISE::finalExpressionBasedSetup
// Purpose       : Setup final expressions
// Special Notes :
// Scope         : public
// Creator       : Meysam Bahmanian
// Creation Date : 5/13/2025
//-----------------------------------------------------------------------------
void HBNOISE::finalExpressionBasedSetup()
{
  return;
}

//-----------------------------------------------------------------------------
// Function      : HBNOISE::getDCOPFlag
// Purpose       : Get DC Operating Point Flag
// Special Notes :
// Scope         : public
// Creator       : Meysam Bahmanian
// Creation Date : 5/13/2025
//-----------------------------------------------------------------------------
bool HBNOISE::getDCOPFlag() const
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : NOISE::updateDataParams_
// Purpose       :
// Special Notes : Used for AC analysis classes, when .DATA is used
// Scope         : private
// Creator       : Meysam Bahmanian
// Creation Date : 5/31/2025
// Originally written by Eric Keiter for the NOISE analysis class
//-----------------------------------------------------------------------------
bool HBNOISE::updateDataParams_ (int stepNumber)
{
  bool reset = updateSweepParams(stepNumber, hbnoiseSweepVector_.begin(), hbnoiseSweepVector_.end());

  bool nonFreqPresent=false;
  lastFreq_ = currentFreq_;
  for (int iac=0;iac<hbnoiseSweepVector_.size();++iac)
  {
    std::string name = hbnoiseSweepVector_[iac].name; Util::toUpper(name);
    double val = hbnoiseSweepVector_[iac].currentVal;
    if (name == "FREQ" || name == "HERTZ")
    {
      currentFreq_ = val;

      delFreq_ = currentFreq_ - lastFreq_;
      lnFreq_     = std::log(std::max(currentFreq_,N_MINLOG));
      lnLastFreq_ = std::log(std::max(lastFreq_,N_MINLOG));
      delLnFreq_  = lnFreq_ - lnLastFreq_;
    }
    else
    {
      nonFreqPresent=true;
      loader_.setParam(name, val, true);
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : HBNOISE::updateCurrentFreq_
// Purpose       :
// Special Notes : Used for AC analysis classes, when .DATA is used
// Scope         : private
// Creator       : Meysam Bahmanian
// Creation Date : 5/31/2025
// Originally written by Eric Keiter for the NOISE analysis class
//-----------------------------------------------------------------------------
bool HBNOISE::updateCurrentFreq_(int stepNumber)
{
  lastFreq_ = currentFreq_;
  if (type_ == "LIN")
  {

    currentFreq_  = fOffsetStart_ + static_cast<double>(stepNumber)*fstep_;
  }
  else if(type_ == "DEC" || type_ == "OCT")
  {

    currentFreq_ = fOffsetStart_*pow(stepMult_, static_cast<double>(stepNumber) );
  }
  else
  {
    Report::DevelFatal().in("HBNOISE::updateCurrentFreq_")
      << "HBNOISE::updateCurrentFreq_: unsupported STEP type";
  }

  delFreq_ = currentFreq_ - lastFreq_;

  lnFreq_     = std::log(std::max(currentFreq_,N_MINLOG));
  lnLastFreq_ = std::log(std::max(lastFreq_,N_MINLOG));
  delLnFreq_  = lnFreq_ - lnLastFreq_;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : HBNOISE::setAnalysisParams
// Purpose       : Sets the HBNOISE analysis parameters
// Special Notes :
// Scope         : public
// Creator       : Meysam Bahmanian
// Creation Date : 5/13/2025
//-----------------------------------------------------------------------------
bool HBNOISE::setAnalysisParams(const Util::OptionBlock & paramsBlock)
{
  bool retval=true;

  // Check for DATA first.  If DATA is present, then use the sweep functions,
  // rather than the HBNOISE specific built-in ones.  This also supports the case
  // of having multiple .HBNOISE lines in the netlist, wherein only the last .HBNOISE
  // line is used.
  if (isDataSpecified(paramsBlock))
  {
    dataSpecification_ = true;
    type_="TYPE";
    hbnoiseSweepVector_.push_back(parseSweepParams(paramsBlock.begin(), paramsBlock.end()));
  }

  for (Util::ParamList::const_iterator it = paramsBlock.begin(),
      end = paramsBlock.end(); it != end; ++it)
  {
    if ((*it).uTag() == "V")
    {
      if ((*it).getImmutableValue<double>()==1.0)
      {
        outputNodeSingle_ = true;
        Util::ParamList::const_iterator itNode = it;
        itNode++;
        outputNode1_ = (*itNode).uTag();
      }
      else if ((*it).getImmutableValue<double>()==2.0)
      {
        outputNodeSingle_ = false;
        Util::ParamList::const_iterator itNode = it;
        itNode++;
        outputNode1_ = (*itNode).uTag();
        itNode++;
        outputNode2_ = (*itNode).uTag();
      }
    }
    else if ((*it).uTag() == "HARMONIC")
    {
      harmonicNumber_ = (*it).getImmutableValue<double>();
      ExtendedString npStr((*it).stringValue());
      if ( !npStr.isInt() )
      {
        Report::UserError0() << "Harmonic parameter on .HBNOISE line must be an integer";
        retval = false;
      }
    }
    else if ((*it).uTag() == "TYPE" && !dataSpecification_)
    {
      type_ = (*it).stringValue();
    }
    else if ((*it).uTag() == "NP")
    {
      np_ = (*it).getImmutableValue<double>();
      ExtendedString npStr((*it).stringValue());
      if ( !npStr.isInt() )
      {
        Report::UserError0() << "Points Value parameter on .HBNOISE line must be an integer";
        retval = false;
      }
    }
    else if ((*it).uTag() == "FSTART")
    {
      fOffsetStart_ = (*it).getImmutableValue<double>();
    }
    else if ((*it).uTag() == "FSTOP")
    {
      fOffsetStop_ = (*it).getImmutableValue<double>();
    }
    else if ((*it).uTag() == "PTS_PER_SUMMARY")
    {
      pts_per_summary_ = (*it).getImmutableValue<int>();
    }
  }

  // exit from here if DATA=<name> is used on the .NOISE line
  if (dataSpecification_) return retval;

  // debug output, when DATA=<name> is not used
  if (DEBUG_ANALYSIS && isActive(Diag::TIME_PARAMETERS))
  {
    dout() << section_divider << std::endl
           << "HBNOISE simulation parameters"
           << std::endl;

    if (outputNodeSingle_)
    {
      dout() << "Output Node: V(" << outputNode1_ << ")" <<std::endl;
    }
    else
    {
      dout() << "Output Node: V(" << outputNode1_ << ","<<outputNode2_<<")" <<std::endl;
    }

    dout() << "harmonic number = " << harmonicNumber_ << std::endl
           << "number of points  = " << np_ << std::endl
           << "start offset frequency = " << fOffsetStart_ << std::endl
           << "stop offset frequency = " << fOffsetStop_ << std::endl
           << "pts_per_summary = " << pts_per_summary_
             << std::endl;
  }

  // error checking of parameters, when DATA=<name> is not used
  if ( np_ < 1 )
  {
    Report::UserError0() << "Points Value parameter on .HBNOISE line must be >= 1";
    retval = false;
  }
  if ( (fOffsetStart_ <=0) || (fOffsetStop_ <= 0) )
  {
    Report::UserError0() << "Illegal values for start or end offset frequencies on .HBNOISE line. " <<
       "Both values must be > 0";
    retval = false;
  }
  if ( fOffsetStop_ < fOffsetStart_ )
  {
    Report::UserError0() << "End offset frequency must not be less than start offset frequency on .HBNOISE line";
    retval = false;
  }

  return retval;
}

//-----------------------------------------------------------------------------
// Function      : HBNOISE::setupSweepParam_
// Purpose       : Processes sweep parameters.
// Special Notes : Used for HBNOISE analysis classes.
// Scope         : private
// Creator       : Meysam Bahmanian
// Creation Date : 5/13/2025
//-----------------------------------------------------------------------------
int HBNOISE::setupSweepParam_()
{
  double fstart, fstop;
  double fcount = 0.0;

  fstart = fOffsetStart_;
  fstop = fOffsetStop_;

  if (DEBUG_ANALYSIS && isActive(Diag::TIME_PARAMETERS))
  {
    Xyce::dout() << std::endl << std::endl;
    Xyce::dout() << section_divider << std::endl;
    Xyce::dout() << "HBNOISE::setupSweepParam_" << std::endl;
  }

    if (type_ == "LIN")
    {
      int np = static_cast<int>(np_);

      if ( np == 1)
        fstep_ = 0;
      else
        fstep_  = (fstop - fstart)/(np_ - 1.0);

      fcount = np_;
      if (DEBUG_ANALYSIS && isActive(Diag::TIME_PARAMETERS))
      {
        Xyce::dout() << "fstep   = " << fstep_  << std::endl;
      }
    }
    else if (type_ == "DEC")
    {
      stepMult_ = std::pow(10.0, 1.0/np_);
      fcount   = floor(fabs(std::log10(fstart) - std::log10(fstop)) * np_ + 1.0);
      if (DEBUG_ANALYSIS && isActive(Diag::TIME_PARAMETERS))
      {
        Xyce::dout() << "stepMult_ = " << stepMult_  << std::endl;
      }
    }
    else if (type_ == "OCT")
    {
      stepMult_ = std::pow(2.0, 1.0/np_);

      // changed to remove dependence on "log2" function, which apparently
      // doesn't exist in the math libraries of FreeBSD or the mingw
      // cross-compilation suite.   Log_2(x)=log_e(x)/log_e(2.0)
      double ln2 = std::log(2.0);
      fcount   = floor(fabs(std::log(fstart) - std::log(fstop))/ln2 * np_ + 1.0);
      if (DEBUG_ANALYSIS && isActive(Diag::TIME_PARAMETERS))
      {
        Xyce::dout() << "stepMult_ = " << stepMult_  << std::endl;
      }
    }
    else
    {
      Report::UserFatal0() << "Unsupported NOISE sweep type: " << type_;
    }

  // At this point, pinterval equals the total number of steps
  // for the step loop.
  return static_cast<int> (fcount);
}

namespace {

typedef Util::Factory<AnalysisBase, HBNOISE> HBNOISEFactoryBase;

//-----------------------------------------------------------------------------
// Class         : HBNOISEFactory
// Purpose       : Factory for creating HBNOISE analysis
// Special Notes :
// Scope         : private
// Creator       : Meysam Bahmanian
// Creation Date : 5/13/2025
//-----------------------------------------------------------------------------
class HBNOISEFactory : public HBNOISEFactoryBase
{
//-----------------------------------------------------------------------------
// Function      : HBNOISEFactory
// Purpose       : Constructor
// Special Notes :
// Scope         : public
// Creator       : Meysam Bahmanian
// Creation Date : 5/13/2025
//-----------------------------------------------------------------------------
public:
  HBNOISEFactory(
    Analysis::AnalysisManager &          analysis_manager,
    Linear::System &                     linear_system,
    Nonlinear::Manager &                 nonlinear_manager,
    Loader::Loader &                     loader,
    Device::DeviceMgr &                  device_manager,
    Topo::Topology &                     topology,
    IO::InitialConditionsManager &       initial_conditions_manager)
    : HBNOISEFactoryBase(),
      analysisManager_(analysis_manager),
      linearSystem_(linear_system),
      nonlinearManager_(nonlinear_manager),
      loader_(loader),
      deviceManager_(device_manager),
      topology_(topology)
  {}

  virtual ~HBNOISEFactory()
  {}

  //-----------------------------------------------------------------------------
  // Function      : create
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : Meysam Bahmanian
  // Creation Date : 5/14/2025
  //-----------------------------------------------------------------------------
  ///
  /// Create a new HBNOISE analysis and applies the analysis and time integrator option blocks.
  ///
  /// @return new HBNOISE analysis object
  ///

  HBNOISE *create() const
  {
    analysisManager_.setAnalysisMode(ANP_MODE_HBNOISE);

    HBNOISE *hbnoise = new HBNOISE(analysisManager_, linearSystem_,
                                  nonlinearManager_, loader_, deviceManager_,
                                  topology_);

    hbnoise->setAnalysisParams(hbnoiseAnalysisOptionBlock_);
    hbnoise->setLinSol(linSolOptionBlock_);

    // Process data statements
    for (std::vector<Util::OptionBlock>::const_iterator it = dataOptionBlockVec_.begin(), end = dataOptionBlockVec_.end(); it != end; ++it)
    {
      hbnoise->setDataStatements(*it);
    }

    return hbnoise;
  }


  //-----------------------------------------------------------------------------
  // Function      : setHBNOISEAnalysisOptionBlock
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : Meysam Bahmanian
  // Creation Date : Thu Jan 29 13:00:14 2015
  //-----------------------------------------------------------------------------
  ///
  /// Saves the analysis parsed options block in the factory.
  ///
  /// @invariant Overwrites any previously specified analysis option block.
  ///
  /// @param option_block parsed option block
  ///
  void setHBNOISEAnalysisOptionBlock(const Util::OptionBlock &option_block)
  {
    hbnoiseAnalysisOptionBlock_ = option_block;
  }

  bool setDotDataBlock(const Util::OptionBlock &option_block)
  {
    dataOptionBlockVec_.push_back(option_block);
    return true;
  }

  // Move these to public section
  AnalysisManager &                     analysisManager_;
  Linear::System &                      linearSystem_;
  Nonlinear::Manager &                  nonlinearManager_;
  Loader::Loader &                      loader_;
  Device::DeviceMgr &                   deviceManager_;
  Topo::Topology &                      topology_;

private:
  Util::OptionBlock     hbnoiseAnalysisOptionBlock_;
  Util::OptionBlock     linSolOptionBlock_;
  std::vector<Util::OptionBlock>        dataOptionBlockVec_;
};

// .HBNOISE
struct HBNOISEAnalysisReg : public IO::PkgOptionsReg
{
  HBNOISEAnalysisReg(
    HBNOISEFactory &   factory )
    : factory_(factory)
  {}

  bool operator()(const Util::OptionBlock &option_block)
  {
    factory_.setHBNOISEAnalysisOptionBlock(option_block);
    factory_.deviceManager_.setBlockAnalysisFlag(true);

    factory_.analysisManager_.addAnalysis(&factory_);

    return true;
  }

  HBNOISEFactory &         factory_;
};

//-----------------------------------------------------------------------------
// Function      : extractHBNOISEData
// Purpose       : Extract the parameters from a netlist .HBNOISE line
// Special Notes :
// Scope         : public
// Creator       : Meysam Bahmanian
// Creation Date : 5/13/2025
//-----------------------------------------------------------------------------
bool extractHBNOISEData(
  IO::PkgOptionsMgr &           options_manager,
  IO::CircuitBlock &            circuit_block,
  const std::string &           netlist_filename,
  const IO::TokenVector &       parsed_line)
{
  // Create option block for HBNOISE analysis
  Util::OptionBlock option_block("HBNOISE", Util::OptionBlock::NO_EXPRESSIONS, netlist_filename, parsed_line[0].lineNumber_);

  int numFields = parsed_line.size();

  // check for "DATA" first.
  bool dataFound=false;
  int pos1=1;
  while ( pos1 < numFields )
  {
    ExtendedString stringVal ( parsed_line[pos1].string_ );
    stringVal.toUpper ();

    if (stringVal == "DATA")
    {
      dataFound=true;
      break;
    }
    ++pos1;
  }

  // Check minimum required fields: .HBNOISE V(OUT) harmonic_number
  if ((!dataFound && numFields < 10) || (dataFound && numFields < 9))
  {
    Report::UserError0().at(netlist_filename, parsed_line[0].lineNumber_)
      << ".HBNOISE line has an unexpected number of fields.  NumFields = " << numFields;
    return false;
  }

  int linePosition = 1;   // Start of parameters

  Util::Param parameter("", "");

  // output node(s):
  ExtendedString stringVal("");
  std::ostringstream msg;
  int p_err=0;
  if( parsed_line[linePosition].string_ == "V" || parsed_line[linePosition].string_ == "v")
  {
    if( parsed_line[linePosition+3].string_ == ")" )
    {
      stringVal = parsed_line[linePosition].string_;
      stringVal.toUpper();
      parameter.setTag(stringVal);
      parameter.setVal( 1.0 );
      option_block.addParam( parameter );

      stringVal = parsed_line[linePosition+2].string_;
      stringVal.toUpper();
      parameter.setTag( stringVal );
      parameter.setVal( 0.0 );
      option_block.addParam( parameter );

      linePosition += 4;
    }
    else if( parsed_line[linePosition+5].string_ == ")" )
    {
      stringVal = parsed_line[linePosition].string_;
      stringVal.toUpper();
      parameter.setTag(stringVal);
      parameter.setVal( 2.0 );
      option_block.addParam( parameter );

      stringVal = parsed_line[linePosition+2].string_;
      stringVal.toUpper();
      parameter.setTag( stringVal );
      parameter.setVal( 0.0 );
      option_block.addParam( parameter );

      stringVal = parsed_line[linePosition+4].string_;
      stringVal.toUpper();
      parameter.setTag( stringVal );
      parameter.setVal( 0.0 );
      option_block.addParam( parameter );

      linePosition += 6;
    }
    else
    {
      msg << "Unrecognized parenthetical specification for HBNOISE output ";
      p_err = linePosition;
    }
  }
  else
  {
    msg << "Incorrect format for HBNOISE output.";
    p_err = linePosition;
  }

  // harmonic numberis required
  stringVal = parsed_line[linePosition].string_;
  stringVal.toUpper();
  parameter.setTag( "HARMONIC" );
  parameter.setVal(std::string(stringVal));
  option_block.addParam( parameter );
  ++linePosition;     // Advance to next parameter.

  // sweep type is required
  stringVal = parsed_line[linePosition].string_;
  stringVal.toUpper();
  parameter.setTag( "TYPE" );
  parameter.setVal(std::string(stringVal));
  option_block.addParam( parameter );
  ++linePosition;     // Advance to next parameter.

  if (dataFound)
  {
    // handle DATA=<name> format
    ++linePosition;  // skip over the = sign
    parameter.setTag( "DATASET" );
    parameter.setVal( parsed_line[ linePosition ].string_ );
    option_block.addParam( parameter );
  }
  else
  {
    // handle format of <points value> <start frequency value> <end frequency value>
    // np is required
    parameter.setTag( "NP" );
    parameter.setVal( parsed_line[linePosition].string_ );
    option_block.addParam( parameter );
    ++linePosition;     // Advance to next parameter.

    // fstart is required
    parameter.setTag( "FSTART" );
    parameter.setVal( parsed_line[linePosition].string_ );
    option_block.addParam( parameter );
    ++linePosition;     // Advance to next parameter.

    // fstop is required
    parameter.setTag( "FSTOP" );
    parameter.setVal( parsed_line[linePosition].string_ );
    option_block.addParam( parameter );
  }

  // pts_per_summary is optional.  If value is negative, assume it wasn't set.
  parameter.setTag( "PTS_PER_SUMMARY" );
  if ( (!dataFound && numFields >= 11) || (dataFound && numFields >= 10) )
  {
    ++linePosition;     // Advance to next parameter.
    parameter.setVal( parsed_line[linePosition].string_ );
  }
  else
  {
    parameter.setVal( std::string("-1") );
  }

  option_block.addParam( parameter );
  circuit_block.addOptions(option_block);

  return true;
}

} // namespace <unnamed>

//-----------------------------------------------------------------------------
// Function      : registerHBNOISEFactory
// Purpose       : Registers the HBNOISE factory with the factory block
// Special Notes :
// Scope         : public
// Creator       : Meysam Bahmanian
// Creation Date : 5/13/2025
//-----------------------------------------------------------------------------
bool registerHBNOISEFactory(FactoryBlock &factory_block)
{
  HBNOISEFactory *factory = new HBNOISEFactory(factory_block.analysisManager_,
    factory_block.linearSystem_, factory_block.nonlinearManager_,
    factory_block.loader_, factory_block.deviceManager_,
    factory_block.topology_, factory_block.initialConditionsManager_);

  addAnalysisFactory(factory_block, factory);

  // Register the command parser and processor
  factory_block.optionsManager_.addCommandParser(".HBNOISE", extractHBNOISEData);
  factory_block.optionsManager_.addCommandProcessor("HBNOISE", new HBNOISEAnalysisReg(*factory));

  factory_block.optionsManager_.addCommandProcessor("DATA",
    IO::createRegistrationOptions(*factory, &HBNOISEFactory::setDotDataBlock) );


  return true;
}

} // namespace Analysis
} // namespace Xyce 