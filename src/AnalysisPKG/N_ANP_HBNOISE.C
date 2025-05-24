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
#include <N_ANP_AnalysisManager.h>
#include <N_ANP_DCSweep.h>
#include <N_ANP_OutputMgrAdapter.h>
#include <N_ANP_Report.h>
#include <N_ANP_Transient.h>
#include <N_ANP_SweepParam.h>
#include <N_ANP_SweepParamFreeFunctions.h>
#include <N_DEV_DeviceMgr.h>
#include <N_IO_ActiveOutput.h>
#include <N_IO_CircuitBlock.h>
#include <N_IO_CmdParse.h>
#include <N_IO_OptionBlock.h>
#include <N_IO_PkgOptionsMgr.h>
#include <N_IO_SpiceSeparatedFieldTool.h>
#include <N_IO_PrintTypes.h>
#include <N_LAS_BlockMatrix.h>
#include <N_LAS_BlockSystemHelpers.h>
#include <N_LAS_BlockVector.h>
#include <N_LAS_HBBuilder.h>
#include <N_LAS_HBPrecondFactory.h>
#include <N_LAS_HBSolverFactory.h>
#include <N_LAS_PrecondFactory.h>
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
#include <N_LOA_HBLoader.h>
#include <N_LOA_NonlinearEquationLoader.h>
#include <N_NLS_Manager.h>
#include <N_PDS_ParMap.h>
#include <N_TIA_DataStore.h>
#include <N_TIA_fwd.h>
#include <N_TIA_StepErrorControl.h>
#include <N_TIA_WorkingIntegrationMethod.h>
#include <N_UTL_APFT.h>
#include <N_UTL_Diagnostic.h>
#include <N_UTL_ExtendedString.h>
#include <N_UTL_FFTInterface.hpp>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_MachDepParams.h>
#include <N_UTL_Math.h>
#include <N_UTL_Timer.h>

#include <Teuchos_BLAS.hpp>
#include <Teuchos_Utils.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialDenseVector.hpp>
#include <Teuchos_SerialDenseHelpers.hpp>
#include <Teuchos_SerialDenseSolver.hpp>

#include <N_TOP_Topology.h>
#include <N_PDS_Comm.h>

#include "N_ANP_HB.h"
#include "N_ANP_NoiseData.h"

using Teuchos::rcp;
using Teuchos::RCP;
using Teuchos::rcp_dynamic_cast;

namespace Xyce {
namespace Analysis {

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
  Topo::Topology &                      topology,
  IO::InitialConditionsManager &        initial_conditions_manager)
  : AnalysisBase(analysis_manager, "HBNOISE"),
    StepEventListener(&analysis_manager),
    analysisManager_(analysis_manager),
    loader_(loader),
    linearSystem_(linear_system),
    nonlinearManager_(nonlinear_manager),
    deviceManager_(device_manager),
    topology_(topology),
    initialConditionsManager_(initial_conditions_manager),
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
    bVecRealPtr(linearSystem_.builder().createVector()),
    bVecImagPtr(linearSystem_.builder().createVector()),
    bNoiseVecRealPtr(linearSystem_.builder().createVector()),
    bNoiseVecImagPtr(linearSystem_.builder().createVector()),
    calcNoiseIntegrals_(true),
    freq_(0.0),
    size_(0),
    period_(0.0)
{
  bVecRealPtr->putScalar(0.0);
  bVecImagPtr->putScalar(0.0);
  bNoiseVecRealPtr->putScalar(0.0);
  bNoiseVecImagPtr->putScalar(0.0);
  
  pdsMgrPtr_ = analysisManager_.getPDSManager();
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
  for (size_t i = 0; i < noiseDataVec_.size(); ++i) {
    delete noiseDataVec_[i];
  }
  noiseDataVec_.clear();
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
  freq_ = freqs[0];
  period_ = 1.0/freq_;
  times_.resize(size_);
  for( int i = 0; i < size_; ++i )
    times_[i] = i*period_/(size_-1);

  // Setup noiseDataVec_ (similar to NOISE constructor/init)
  // Ensure loader_ is ready and devices are instantiated.
  // This might need to happen after hbAnalysis_->doInit() if that's what sets up the relevant loader state.
  // For now, assume loader_ is the correct one and ready.
  int numNoiseDevices = loader_.getNumNoiseDevices();
  noiseDataVec_.resize(numNoiseDevices);
  for (int i = 0; i < numNoiseDevices; ++i) {
    noiseDataVec_[i] = new Analysis::NoiseData();
  }
  loader_.setupNoiseSources(noiseDataVec_); // This populates deviceName, noiseNames, li_Pos, li_Neg, etc.


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

  // after dataSpecification_ validation, run the HB analysis
  analysisManager_.pushActiveAnalysis(hbAnalysis_);
  bsuccess = hbAnalysis_->run();
  hbBuilderPtr_ = hbAnalysis_->hbBuilderPtr_;
  builderPtr_ = &(hbAnalysis_->builder_);
  hbLoaderPtr_ = hbAnalysis_->hbLoaderPtr_;
  updateLinearTimeVariantSystem_C_and_G_();
  createHarmonicSpaceLinearSystem_();

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
  return true;
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
bool HBNOISE::createHarmonicSpaceLinearSystem_(){
  // first take the Fourier Transform of Ct_ and Gt_
  for (int i=0; i<BlockSize_; i++){
    Cf_.push_back(hbBuilderPtr_->createExpandedRealFormTransposeBlockVector());
    Gf_.push_back(hbBuilderPtr_->createExpandedRealFormTransposeBlockVector());
    Cf_[i]->putScalar(0.0);
    Gf_[i]->putScalar(0.0);
    hbLoaderPtr_->permutedFFT2(*(Ct_[i]), &*(Cf_[i]));
    hbLoaderPtr_->permutedFFT2(*(Gt_[i]), &*(Gf_[i]));
  }
  if (DEBUG_HBNOISE)
  {
    Xyce::dout() << "Reporting Gf_ Matrices, each block is a node" << std::endl;
    for (int i=0; i<BlockSize_; i++){
      Xyce::dout() << "Gf_[" << i << "]: " << std::endl;
      Gf_[i]->print(Xyce::dout());
      Xyce::dout() << std::endl;
    }
    //for (int i=0; i<BlockSize_; i++){
    //  Xyce::dout() << "Cf_[" << i << "]: " << std::endl;
    //  Cf_[i]->print(Xyce::dout());
    //  Xyce::dout() << std::endl;
    //}
  }



  // now take the Fourier Transform of Ct_ and Gt_
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

  Teuchos::RCP<Linear::BlockVector> bXtPtr_ = hbBuilderPtr_->createTimeDomainBlockVector();
  bXtPtr_->putScalar(0.0);

  hbLoaderPtr_->permutedIFT(bXf, &*bXtPtr_);

  Linear::BlockVector & bX = *bXtPtr_;
  BlockCount_ = bX.blockCount(); // number of time points
  BlockSize_ = bX.blockSize(); // number of GIDs

  if (DEBUG_HBNOISE)
  {
    for (int i = 0; i < BlockCount_; ++i)
    {
      Xyce::dout() << "Solution time domain, block (" << i << "): each block is a time point" << std::endl;
      bX.block(i).print( Xyce::dout() );
      Xyce::dout() << std::endl;
    }
    for (int i = 0; i < BlockSize_; ++i)
    {
      Xyce::dout() << "Solution frequency domain, block (" << i << "): each block is a node" << std::endl;
      bXf.block(i).print( Xyce::dout() );
      Xyce::dout() << std::endl;
    }
  }

  // Solutions:
  Linear::Vector * currSolutionPtr = builderPtr_->createVector();

  for (int i=0; i<BlockSize_; i++){
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
  for (int i = 0; i < BlockCount_; ++i)
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
    std::vector<double> coeffs(BlockSize_); 
    std::vector<int> colIndices(BlockSize_);

    for (int j=0; j<BlockSize_; j++) {
      tmpC->getLocalRowCopy(j, BlockSize_, numEntries, coeffs.data(), colIndices.data());
      for (int k = 0; k < numEntries; k++) {
        Ct_[j]->block(i)[colIndices[k]] = coeffs[k];
      }
    }

    for (int j=0; j<BlockSize_; j++) {
      tmpG->getLocalRowCopy(j, BlockSize_, numEntries, coeffs.data(), colIndices.data());
      for (int k = 0; k < numEntries; k++) {
        Gt_[j]->block(i)[colIndices[k]] = coeffs[k];
      }
    }

    if (DEBUG_HBNOISE)
    {
      // print conductance matrix
      Xyce::dout() << "dFdxMatrixPtr block(" << i << "):" << std::endl;
      dFdxMatrixPtr->print( Xyce::dout() );
      Xyce::dout() << std::endl;

      // print capacitance matrix
      Xyce::dout() << "dQdxMatrixPtr block(" << i << "):" << std::endl;
      dQdxMatrixPtr->print( Xyce::dout() );
      Xyce::dout() << std::endl;
    }
  }

  if (DEBUG_HBNOISE)
  {
    Xyce::dout() << "Reporting Gt_ Matrices, each block is a time point" << std::endl;
    for (int i=0; i<BlockSize_; i++){
      Xyce::dout() << "Gt_[" << i << "]: " << std::endl;
      Gt_[i]->print(Xyce::dout());
      Xyce::dout() << std::endl;
    }
  }

  return true;
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
// Function      : HBNOISE::setLinSol
// Purpose       : Save linear solver options
// Special Notes :
// Scope         : public
// Creator       : Meysam Bahmanian
// Creation Date : 5/13/2025
//-----------------------------------------------------------------------------
bool HBNOISE::setLinSol(const Util::OptionBlock & OB)
{
  saved_lsOB_ = OB;
  return true;
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
      topology_(topology),
      initialConditionsManager_(initial_conditions_manager)
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
                                  topology_, initialConditionsManager_);

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
  IO::InitialConditionsManager &        initialConditionsManager_;

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