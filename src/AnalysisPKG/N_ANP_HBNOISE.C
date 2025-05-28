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
  omega_ = 2.0*M_PI*freq_;
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
// G_ij * V_jS + w * C_ij * V_jC = B_iS
//-----------------------------------------------------------------------------
//                     for Harmonic space system (with Einstein convention):
// G_ij(t) = G_ij0 + G_ijfI*cos(f*wc*t) - G_ijfQ*sin(f*wc*t) // sum over f
// where f is the harmonic number, wc is the carrier frequency, I means in-phase and Q means quadrature.
// The form of b_j now is
// b_j = B_jC*cos(wm*t) - B_jS*sin(wm*t) + (B_jeIC*cos(wm*t) - B_jeIS*sin(wm*t))*cos(e*wc*t) - (B_jeQC*cos(wm*t) - B_jeQS*sin(wm*t))*sin(e*wc*t) // sum over e
// The form of v_j now is
// v_j = V_jC*cos(wm*t) - V_jS*sin(wm*t) + (V_jeIC*cos(wm*t) - V_jeIS*sin(wm*t))*cos(e*wc*t) - (V_jeQC*cos(wm*t) - V_jeQS*sin(wm*t))*sin(e*wc*t) // sum over e
// where e is the harmonic number, wc is the carrier frequency, I means in-phase and Q means quadrature.
// where wm is the modulation frequency
// Now we need to multiply G_ij(t) with v_j and equate it to b_j(t)
// This gives us the Harmonic space matrix
// instead of writing the matrix equation, we will write the contributions:
// for all harmonics k, G_ij0 represennt the AC linear system without any frequency translation
// All harmonics are linearly transformed
// Now we have to find the contrinbutions of terms
// [ G_ijfI*cos(f*wc*t) - G_ijfQ*sin(f*wc*t) ] * v_j
//
// a) The baseband terms: B_jC*cos(wm*t) - B_jS*sin(wm*t)
// are directly transformed to f'th harmonic in-phase and quadrature components resulting the terms:
// + G_ijfI*V_jC*cos(wm*t)*cos(f*wc*t) // f'th harmonic in-phase cosine term
// - G_ijfI*V_jS*sin(wm*t)*cos(f*wc*t) // f'th harmonic in-phase sine term
// - G_ijfQ*V_jC*cos(wm*t)*sin(f*wc*t) // f'th harmonic quadrature cosine term
// + G_ijfQ*V_jS*sin(wm*t)*sin(f*wc*t) // f'th harmonic quadrature sine term
//
// b) equal harmonic mixing terms (e=f): 
// [G_ijfI*cos(f*wc*t) - G_ijfQ*sin(f*wc*t)] * [(V_jeIC*cos(wm*t) - V_jeIS*sin(wm*t))*cos(e*wc*t) - (V_jeQC*cos(wm*t) - V_jeQS*sin(wm*t))*sin(e*wc*t)]
// This generate 2 baseband terms and 4 harmonic terms at f'th (or e'th) harmonic.
// The overal baseband terms are:
// + 0.5 * G_ijfI * (V_jeIC*cos(wm*t) - V_jeIS*sin(wm*t)) + 0.5 * G_ijfQ * (V_jeQC*cos(wm*t) - V_jeQS*sin(wm*t))
// So the baseband terms are:
// + 0.5 * G_ijfI * V_jeIC*cos(wm*t) + 0.5 * G_ijfQ * V_jeQC*cos(wm*t) // cosine term
// - 0.5 * G_ijfI * V_jeIS*sin(wm*t) - 0.5 * G_ijfQ * V_jeQS*sin(wm*t) // sine term
// The overal in-phase terms are:
// [ + 0.5 * G_ijfI * (V_jeIC*cos(wm*t) - V_jeIS*sin(wm*t)) - 0.5 * G_ijfQ * (V_jeQC*cos(wm*t) - V_jeQS*sin(wm*t)) ] * cos(2*e*wc*t)
// The overal quadrature terms are:
// [ - 0.5 * G_ijfI * (V_jeQC*cos(wm*t) - V_jeQS*sin(wm*t)) - 0.5 * G_ijfQ * (V_jeIC*cos(wm*t) - V_jeIS*sin(wm*t)) ] * sin(2*e*wc*t)
// So the 2*e'th (2*f'th) harmonic terms are:
// [ + 0.5 * G_ijfI * V_jeIC - 0.5 * G_ijfQ * V_jeQC ] * cos(wm*t) * cos(2*e*wc*t) // 2*e'th harmonic in-phase cosine term
// [ - 0.5 * G_ijfI * V_jeIS + 0.5 * G_ijfQ * V_jeQS ] * sin(wm*t) * cos(2*e*wc*t) // 2*e'th harmonic in-phase sine term
// [ - 0.5 * G_ijfI * V_jeQC - 0.5 * G_ijfQ * V_jeIC ] * cos(wm*t) * sin(2*e*wc*t) // 2*e'th harmonic quadrature cosine term
// [ + 0.5 * G_ijfI * V_jeQS + 0.5 * G_ijfQ * V_jeIS ] * sin(wm*t) * sin(2*e*wc*t) // 2*e'th harmonic quadrature sine term

// b) non-equal harmonic mixing terms (e!=f):
// [G_ijfI*cos(f*wc*t) - G_ijfQ*sin(f*wc*t)] * [(V_jeIC*cos(wm*t) - V_jeIS*sin(wm*t))*cos(e*wc*t) - (V_jeQC*cos(wm*t) - V_jeQS*sin(wm*t))*sin(e*wc*t)]
// This gives us 4 terms at delta=f-e and another 4 terms at sigma=f+e
// delta'th harmonic terms:
// [ + 0.5 * G_ijfI * V_jeIC + 0.5 * G_ijfQ * V_jeQC ] * cos(wm*t) * cos(delta*wc*t) // delta'th harmonic in-phase cosine term
// [ - 0.5 * G_ijfI * V_jeIS - 0.5 * G_ijfQ * V_jeQS ] * sin(wm*t) * cos(delta*wc*t) // delta'th harmonic in-phase sine term
// [ + 0.5 * G_ijfI * V_jeQC - 0.5 * G_ijfQ * V_jeIC ] * cos(wm*t) * sin(delta*wc*t) // delta'th harmonic quadrature cosine term
// [ - 0.5 * G_ijfI * V_jeQS + 0.5 * G_ijfQ * V_jeIS ] * sin(wm*t) * sin(delta*wc*t) // delta'th harmonic quadrature sine term
// sigma'th harmonic terms:
// [ + 0.5 * G_ijfI * V_jeIC - 0.5 * G_ijfQ * V_jeQC ] * cos(wm*t) * cos(sigma*wc*t) // sigma'th harmonic in-phase cosine term
// [ - 0.5 * G_ijfI * V_jeIS + 0.5 * G_ijfQ * V_jeQS ] * sin(wm*t) * cos(sigma*wc*t) // sigma'th harmonic in-phase sine term
// [ - 0.5 * G_ijfI * V_jeQC - 0.5 * G_ijfQ * V_jeIC ] * cos(wm*t) * sin(sigma*wc*t) // sigma'th harmonic quadrature cosine term
// [ + 0.5 * G_ijfI * V_jeQS + 0.5 * G_ijfQ * V_jeIS ] * sin(wm*t) * sin(sigma*wc*t) // sigma'th harmonic quadrature sine term
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
// v_j = V_jC*cos(wm*t) - V_jS*sin(wm*t) + (V_jeIC*cos(wm*t) - V_jeIS*sin(wm*t))*cos(e*wc*t) - (V_jeQC*cos(wm*t) - V_jeQS*sin(wm*t))*sin(e*wc*t) // sum over e
// we divide the derivative in 2 parts, one that scales with wm and one that scales with wc
// the wm part is offset freq and has to be swept in a loop and has to be evaluated for every loop iteration
// PART 1 scales with wc and performs cross-harmonic quadrature transformation:
// dv_j/dt |1 = - ( V_jeIC*cos(wm*t) - V_jeIS*sin(wm*t) )*e*wc*sin(e*wc*t) - ( V_jeQC*cos(wm*t) - V_jeQS*sin(wm*t) )*e*wc*cos(e*wc*t)
// PART 2 scales with wm and performs coss-harmonic shaping:
// dv_j/dt |2 = - wm*V_jC*sin(wm*t) - wm*V_jS*cos(wm*t) + wm * ( - V_jeIC*sin(wm*t) - V_jeIS*cos(wm*t) )*cos(e*wc*t) - wm * ( - V_jeQC*sin(wm*t) - V_jeQS*cos(wm*t))*sin(e*wc*t) 

// This approach gives us the desired format:
// ( [G] + [omega*C, PART 0] + [omega*C, PART 1] + wm*[C, PART 2] ) * v = b
// The matrices will be calculated once, of for each offset frequency PART 2 will be scaled with wm and added to the matrix.


bool HBNOISE::createHarmonicSpaceLinearSystem_(){
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
    Xyce::dout() << "Reporting Gf_ Matrices, each block is a node" << std::endl;
    for (int i=0; i<BlockSize; i++){
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

  Parallel::Manager &pds_manager = *analysisManager_.getPDSManager();

  RCP<Parallel::ParMap> baseMap = rcp(pds_manager.getParallelMap( Parallel::SOLUTION ), false);
  const Linear::Graph* baseFullGraph = pds_manager.getMatrixGraph(Parallel::JACOBIAN);

  int numHarms = (size_-1)/2;
  int numBlocks = 2 + 4 * numHarms; // 2 is for our old AC part (real/imag) and every harmonic has 4 terms (real/imag in-phase and real/imag quadrature)
  int offset = baseMap->maxGlobalEntity() + 1;  // Use this offset to create a contiguous gid map for direct solvers.

  RCP<Parallel::ParMap> blockMap = Linear::createBlockParMap(numBlocks, *baseMap, 0, 0, offset);
  harmonicSpaceB_ = Xyce::Linear::createBlockVector(numBlocks, blockMap, baseMap);

  std::vector<std::vector<int> > blockPattern(numBlocks);
  for (int i=0; i<numBlocks; i++){
    blockPattern[i].resize(numBlocks);
    for (int j=0; j<numBlocks; j++){
      blockPattern[i][j] = j;
    }
  }

  RCP<Linear::Graph> blockGraph = Linear::createBlockGraph( offset, blockPattern, *blockMap, *baseFullGraph);
  harmonicSpaceMatrix_ = Xyce::Linear::createBlockMatrix( numBlocks, offset, blockPattern, blockGraph.get(), baseFullGraph);
  harmonicSpaceMatrix_->put( 0.0 ); 

  // I currently do it in a loop iteration. But I think it will be easier using tensors. Fourier transform of G_
  // is a tensor of rank 3. We can use a tranformation tensor to get the harmonic space G matrix.
  // Ct_ is also a tensor:
  //
  //     Gt_ (or Ct_) tensor            Fourier matrix             Gf_ (or Cf_) tensor    
  //
  //           ---------|                 \    /                        \-------\           
  //  time->  /       / |          freq->  \  / <-time                  |\       \  <-freq
  //         /-------/  |    *              \/             =            \ \-------\  
  //  node-> | Ct/Gt | /                                         node->  \| Cf/Gf |   
  //         |-------|/                                                   \-------|   
  //             ^                                                            ^        
  //           node                                                          node       
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
  // Gf_ strucure: std::vector with GID elements, each element has GID blocks, each block has real/imag freq points
  // the frequency points are in total 2*size_ = 2 * (2*numberOfHarmonics + 1)
  // 0: dc real
  // 1: dc imag
  // 2: 1st harmonic real 
  // 3: 1st harmonic imag etc.

  // Harmonic space structure:
  // TBD
  // 
  harmonicSpaceMatrix_G_ = Xyce::Linear::createBlockMatrix( numBlocks, offset, blockPattern, blockGraph.get(), baseFullGraph);
  harmonicSpaceMatrix_G_->put( 0.0 );
  int numRows = Gf_.size();

  for (int i=0; i<numRows; i++)
  { // selectring a row of G Matrix
    for (int j=0; j<numRows; j++)
    { // selecting a column of G Matrix
      for (int f=0; f<=numHarms; f++)
      { // selecting a frequency point of Gf_
        // now we fill the matrix
        // first the diagon of G matrix
        if (f==0)
        { // linear transformation
          for (int l=0; l<numBlocks; l++){
            harmonicSpaceMatrix_G_->block(l,l)[i][j] = Gf_[i]->block(j)[0];
          }
        } else
        { // now the mixing parts (harmonic coupling)
          for (int e=0; e<=numHarms; e++)
          { // harmonic index of v
            if (e==0) 
            {
              // Maybe I should have made two dummy blocks for baseband so the k-indices were not so confusing!
              // baseband modulation
              // [ G_ijfI*cos(f*wc*t) - G_ijfQ*sin(f*wc*t) ] * [ V_jC*cos(wm*t) - V_jS*sin(wm*t) ]

              // + G_ijfI*cos(f*wc*t) * V_jC*cos(wm*t)
              // RHS has positive sign
              harmonicSpaceMatrix_G_->block(4*f-2,0)[i][j] = +2*Gf_[i]->block(j)[2*f]; // in-phase cosine

              // - G_ijfI*cos(f*wc*t) * V_jS*sin(wm*t)
              // RHS has negative sign
              harmonicSpaceMatrix_G_->block(4*f-1,0)[i][j] = +2*Gf_[i]->block(j)[2*f]; // in-phase sine

              // - G_ijfQ*sin(f*wc*t) * V_jC*cos(wm*t)
              // RHS has positive sign
              harmonicSpaceMatrix_G_->block(4*f  ,0)[i][j] = -2*Gf_[i]->block(j)[2*f+1]; // quadrature cosine

              // + G_ijfQ*sin(f*wc*t) * V_jS*sin(wm*t)
              // RHS has negative sign
              harmonicSpaceMatrix_G_->block(4*f+1,0)[i][j] = -2*Gf_[i]->block(j)[2*f+1]; // quadrature sine
            } else 
            { // now f>0 and e>0
              // [ G_ijfI*cos(f*wc*t) - G_ijfQ*sin(f*wc*t) ] * [ (V_jeIC*cos(wm*t) - V_jeIS*sin(wm*t))*cos(e*wc*t) - (V_jeQC*cos(wm*t) - V_jeQS*sin(wm*t))*sin(e*wc*t) ]
              int sigma = f+e;
              // int delta = f-e; this line is just for the sake of understanding the code. Equations are with reference to delta, not deltaAbs.
              int deltaAbs = std::abs(f-e);
              int sign = f>=e ? 1 : -1;

              // case sigma
              if (sigma<=numHarms) 
              { // f+e should not be larget than numHarms, otherwise ignore it
                // in finite-harmonics space the system stil shows nonlineary and some mixing products have to be ignored

                // + G_ijfI*cos(f*wc*t) * V_jeIC*cos(wm*t) * cos(e*wc*t) + G_ijfQ*sin(f*wc*t) * V_jeQC*cos(wm*t) * sin(e*wc*t)
                // + 0.5*G_ijfI*cos(sigma*wc*t) * V_jeIC*cos(wm*t) - 0.5*G_ijfQ*cos(sigma*wc*t) * V_jeQC*cos(wm*t)
                // RHS has positive sign
                harmonicSpaceMatrix_G_->block(4*sigma-2, 4*e-2)[i][j] = +Gf_[i]->block(j)[2*f  ]; // in-phase cosine translated from in-phase cosine by in-phase G
                harmonicSpaceMatrix_G_->block(4*sigma-2, 4*e  )[i][j] = -Gf_[i]->block(j)[2*f+1]; // in-phase cosine translated from quadrature cosine by quadrature G

                // - G_ijfI*cos(f*wc*t) * V_jeIS*sin(wm*t) * cos(e*wc*t) - G_ijfQ*sin(f*wc*t) * V_jeQS*sin(wm*t) * sin(e*wc*t)
                // - 0.5*G_ijfI*cos(sigma*wc*t) * V_jeIS*sin(wm*t) + 0.5*G_ijfQ*cos(sigma*wc*t) * V_jeQS*sin(wm*t)
                // RHS has negative sign
                harmonicSpaceMatrix_G_->block(4*sigma-1,4*e-1)[i][j] = +Gf_[i]->block(j)[2*f  ]; // in-phase sine translated from in-phase sine by in-phase G
                harmonicSpaceMatrix_G_->block(4*sigma-1,4*e+1)[i][j] = -Gf_[i]->block(j)[2*f+1]; // in-phase sine translated from quadrature sine by quadrature G

                // - G_ijfI*cos(f*wc*t) * V_jeQC*cos(wm*t) * sin(e*wc*t) - G_ijfQ*sin(f*wc*t) * V_jeIC*cos(wm*t) * cos(e*wc*t)
                // - 0.5*G_ijfI*sin(sigma*wc*t) * V_jeQC*cos(wm*t) - 0.5*G_ijfQ*sin(sigma*wc*t) * V_jeIC*cos(wm*t)
                // RHS has positive sign
                harmonicSpaceMatrix_G_->block(4*sigma  ,4*e  )[i][j] = -Gf_[i]->block(j)[2*f  ]; // quadrature cosine translated from quadrature cosine by in-phase G
                harmonicSpaceMatrix_G_->block(4*sigma  ,4*e-2)[i][j] = -Gf_[i]->block(j)[2*f+1]; // quadrature cosine translated from in-phase cosine by quadrature G

                // + G_ijfI*cos(f*wc*t) * V_jeQS*sin(wm*t) * sin(e*wc*t) + G_ijfQ*sin(f*wc*t) * V_jeIS*sin(wm*t) * cos(e*wc*t)
                // + 0.5*G_ijfI*sin(sigma*wc*t) * V_jeQS*sin(wm*t) + 0.5*G_ijfQ*sin(sigma*wc*t) * V_jeIS*sin(wm*t)
                // RHS has negative sign
                harmonicSpaceMatrix_G_->block(4*sigma+1,4*e+1)[i][j] = -Gf_[i]->block(j)[2*f  ]; // quadrature sine translated from quadrature sine by in-phase G
                harmonicSpaceMatrix_G_->block(4*sigma+1,4*e-1)[i][j] = -Gf_[i]->block(j)[2*f+1]; // quadrature sine translated from in-phase sine by quadrature G
              }

              // case delta
              // now I have to deal with the index issue! if delta==0 then the indices become negative!
              if (deltaAbs==0)
              {
                // translation from harmonics to baseband
                // + G_ijfI*cos(f*wc*t) * V_jeIC*cos(wm*t) * cos(e*wc*t) + G_ijfQ*sin(f*wc*t) * V_jeQC*cos(wm*t) * sin(e*wc*t)
                // + 0.5*G_ijfI * V_jeIC*cos(wm*t) + 0.5*G_ijfQ * V_jeQC*cos(wm*t)
                // RHS has positive sign
                harmonicSpaceMatrix_G_->block(0, 4*e-2)[i][j] = +Gf_[i]->block(j)[2*f  ]; // baseband cosine translated from in-phase cosine by in-phase G
                harmonicSpaceMatrix_G_->block(1, 4*e  )[i][j] = +Gf_[i]->block(j)[2*f+1]; // baseband cosine translated from quadrature cosine by quadrature G

                // - G_ijfI*cos(f*wc*t) * V_jeIS*sin(wm*t) * cos(e*wc*t) - G_ijfQ*sin(f*wc*t) * V_jeQS*sin(wm*t) * sin(e*wc*t)
                // - 0.5*G_ijfI * V_jeIS*sin(wm*t) - 0.5*G_ijfQ * V_jeQS*sin(wm*t)
                // RHS has negative sign
                harmonicSpaceMatrix_G_->block(0, 4*e-1)[i][j] = +Gf_[i]->block(j)[2*f  ]; // baseband sine translated from in-phase sine by in-phase G
                harmonicSpaceMatrix_G_->block(1, 4*e+1)[i][j] = +Gf_[i]->block(j)[2*f+1]; // baseband sine translated from quadrature sine by quadrature G
              }
              else
              {
                // + G_ijfI*cos(f*wc*t) * V_jeIC*cos(wm*t) * cos(e*wc*t) + G_ijfQ*sin(f*wc*t) * V_jeQC*cos(wm*t) * sin(e*wc*t)
                // + 0.5*G_ijfI*cos(delta*wc*t) * V_jeIC*cos(wm*t) + 0.5*G_ijfQ*cos(delta*wc*t) * V_jeQC*cos(wm*t)
                // RHS has positive sign
                harmonicSpaceMatrix_G_->block(4*deltaAbs-2, 4*e-2)[i][j] = +Gf_[i]->block(j)[2*f  ]; // in-phase cosine translated from in-phase cosine by in-phase G
                harmonicSpaceMatrix_G_->block(4*deltaAbs-2, 4*e  )[i][j] = +Gf_[i]->block(j)[2*f+1]; // in-phase cosine translated from quadrature cosine by quadrature G

                // - G_ijfI*cos(f*wc*t) * V_jeIS*sin(wm*t) * cos(e*wc*t) - G_ijfQ*sin(f*wc*t) * V_jeQS*sin(wm*t) * sin(e*wc*t)
                // - 0.5*G_ijfI*cos(delta*wc*t) * V_jeIS*sin(wm*t) - 0.5*G_ijfQ*cos(delta*wc*t) * V_jeQS*sin(wm*t)
                // RHS has negative sign
                harmonicSpaceMatrix_G_->block(4*deltaAbs-1,4*e-1)[i][j] = +Gf_[i]->block(j)[2*f  ]; // in-phase sine translated from in-phase sine by in-phase G
                harmonicSpaceMatrix_G_->block(4*deltaAbs-1,4*e+1)[i][j] = +Gf_[i]->block(j)[2*f+1]; // in-phase sine translated from quadrature sine by quadrature G

                // - G_ijfI*cos(f*wc*t) * V_jeQC*cos(wm*t) * sin(e*wc*t) - G_ijfQ*sin(f*wc*t) * V_jeIC*cos(wm*t) * cos(e*wc*t)
                // + 0.5*G_ijfI*sin(delta*wc*t) * V_jeQC*cos(wm*t) - 0.5*G_ijfQ*sin(delta*wc*t) * V_jeIC*cos(wm*t)
                // RHS has positive sign
                harmonicSpaceMatrix_G_->block(4*deltaAbs  ,4*e  )[i][j] = +sign*Gf_[i]->block(j)[2*f  ]; // quadrature cosine translated from quadrature cosine by in-phase G
                harmonicSpaceMatrix_G_->block(4*deltaAbs  ,4*e-2)[i][j] = -sign*Gf_[i]->block(j)[2*f+1]; // quadrature cosine translated from in-phase cosine by quadrature G

                // + G_ijfI*cos(f*wc*t) * V_jeQS*sin(wm*t) * sin(e*wc*t) + G_ijfQ*sin(f*wc*t) * V_jeIS*sin(wm*t) * cos(e*wc*t)
                // - 0.5*G_ijfI*sin(delta*wc*t) * V_jeQS*sin(wm*t) + 0.5*G_ijfQ*sin(delta*wc*t) * V_jeIS*sin(wm*t)
                // RHS has negative sign
                harmonicSpaceMatrix_G_->block(4*deltaAbs+1,4*e+1)[i][j] = +sign*Gf_[i]->block(j)[2*f  ]; // quadrature sine translated from quadrature sine by in-phase G
                harmonicSpaceMatrix_G_->block(4*deltaAbs+1,4*e-1)[i][j] = -sign*Gf_[i]->block(j)[2*f+1]; // quadrature sine translated from in-phase sine by quadrature G

              }
            } // end of f>0 and e>0
          } // end of e loop
        } // end of f!=0
      } // end of f loop of Gf_
    } // end of column loop
  } // end of row loop
  // And this was just for G matrix. Now we have to do the same for three parts of C matrix and we have to deal with derivatives too!
  
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
  harmonicSpaceMatrix_omegaC_0_ = Xyce::Linear::createBlockMatrix( numBlocks, offset, blockPattern, blockGraph.get(), baseFullGraph);
  harmonicSpaceMatrix_omegaC_0_->put( 0.0 );

  for (int i=0; i<numRows; i++)
  { // selectring a row of G Matrix
    for (int j=0; j<numRows; j++)
    { // selecting a column of G Matrix
      for (int f=0; f<=numHarms; f++)
      { // selecting a frequency point of Gf_
        // now we fill the matrix
        // first the diagon of omega*C matrix
        if (f==0)
        { // linear transformation
          // these terms are zero for dC(t)/dt
          // for (int l=0; l<numBlocks; l++){
          //   harmonicSpaceMatrix_omegaC0_->block(l,l)[i][j] = 0*Cf_[i]->block(j)[0];
          // }
        } else
        { // now the mixing parts (harmonic coupling)
          for (int e=0; e<=numHarms; e++)
          { // harmonic index of v
            if (e==0) 
            {
              // Maybe I should have made two dummy blocks for baseband so the k-indices were not so confusing!
              // baseband modulation
              // [ G_ijfI*cos(f*wc*t) - G_ijfQ*sin(f*wc*t) ] * [ V_jC*cos(wm*t) - V_jS*sin(wm*t) ]

              // + G_ijfI*cos(f*wc*t) * V_jC*cos(wm*t)
              // RHS has positive sign
              harmonicSpaceMatrix_omegaC_0_->block(4*f-2,0)[i][j] = +2*(-omega_*f)*Cf_[i]->block(j)[2*f+1]; // in-phase cosine

              // - G_ijfI*cos(f*wc*t) * V_jS*sin(wm*t)
              // RHS has negative sign
              harmonicSpaceMatrix_omegaC_0_->block(4*f-1,0)[i][j] = +2*(-omega_*f)*Cf_[i]->block(j)[2*f+1]; // in-phase sine

              // - G_ijfQ*sin(f*wc*t) * V_jC*cos(wm*t)
              // RHS has positive sign
              harmonicSpaceMatrix_omegaC_0_->block(4*f  ,0)[i][j] = -2*(+omega_*f)*Cf_[i]->block(j)[2*f]; // quadrature cosine

              // + G_ijfQ*sin(f*wc*t) * V_jS*sin(wm*t)
              // RHS has negative sign
              harmonicSpaceMatrix_omegaC_0_->block(4*f+1,0)[i][j] = -2*(+omega_*f)*Cf_[i]->block(j)[2*f]; // quadrature sine
            } else 
            { // now f>0 and e>0
              // [ G_ijfI*cos(f*wc*t) - G_ijfQ*sin(f*wc*t) ] * [ (V_jeIC*cos(wm*t) - V_jeIS*sin(wm*t))*cos(e*wc*t) - (V_jeQC*cos(wm*t) - V_jeQS*sin(wm*t))*sin(e*wc*t) ]
              int sigma = f+e;
              // int delta = f-e; this line is just for the sake of understanding the code. Equations are with reference to delta, not deltaAbs.
              int deltaAbs = std::abs(f-e);
              int sign = f>=e ? 1 : -1;

              // case sigma
              if (sigma<=numHarms) 
              { // f+e should not be larget than numHarms, otherwise ignore it
                // in finite-harmonics space the system stil shows nonlineary and some mixing products have to be ignored

                // + G_ijfI*cos(f*wc*t) * V_jeIC*cos(wm*t) * cos(e*wc*t) + G_ijfQ*sin(f*wc*t) * V_jeQC*cos(wm*t) * sin(e*wc*t)
                // + 0.5*G_ijfI*cos(sigma*wc*t) * V_jeIC*cos(wm*t) - 0.5*G_ijfQ*cos(sigma*wc*t) * V_jeQC*cos(wm*t)
                // RHS has positive sign
                harmonicSpaceMatrix_omegaC_0_->block(4*sigma-2, 4*e-2)[i][j] = +(-omega_*f)*Cf_[i]->block(j)[2*f+1]; // in-phase cosine translated from in-phase cosine by in-phase G
                harmonicSpaceMatrix_omegaC_0_->block(4*sigma-2, 4*e  )[i][j] = -(+omega_*f)*Cf_[i]->block(j)[2*f]; // in-phase cosine translated from quadrature cosine by quadrature G

                // - G_ijfI*cos(f*wc*t) * V_jeIS*sin(wm*t) * cos(e*wc*t) - G_ijfQ*sin(f*wc*t) * V_jeQS*sin(wm*t) * sin(e*wc*t)
                // - 0.5*G_ijfI*cos(sigma*wc*t) * V_jeIS*sin(wm*t) + 0.5*G_ijfQ*cos(sigma*wc*t) * V_jeQS*sin(wm*t)
                // RHS has negative sign
                harmonicSpaceMatrix_omegaC_0_->block(4*sigma-1,4*e-1)[i][j] = +(-omega_*f)*Cf_[i]->block(j)[2*f+1]; // in-phase sine translated from in-phase sine by in-phase G
                harmonicSpaceMatrix_omegaC_0_->block(4*sigma-1,4*e+1)[i][j] = -(+omega_*f)*Cf_[i]->block(j)[2*f]; // in-phase sine translated from quadrature sine by quadrature G

                // - G_ijfI*cos(f*wc*t) * V_jeQC*cos(wm*t) * sin(e*wc*t) - G_ijfQ*sin(f*wc*t) * V_jeIC*cos(wm*t) * cos(e*wc*t)
                // - 0.5*G_ijfI*sin(sigma*wc*t) * V_jeQC*cos(wm*t) - 0.5*G_ijfQ*sin(sigma*wc*t) * V_jeIC*cos(wm*t)
                // RHS has positive sign
                harmonicSpaceMatrix_omegaC_0_->block(4*sigma  ,4*e  )[i][j] = -(-omega_*f)*Cf_[i]->block(j)[2*f+1]; // quadrature cosine translated from quadrature cosine by in-phase G
                harmonicSpaceMatrix_omegaC_0_->block(4*sigma  ,4*e-2)[i][j] = -(+omega_*f)*Cf_[i]->block(j)[2*f]; // quadrature cosine translated from in-phase cosine by quadrature G

                // + G_ijfI*cos(f*wc*t) * V_jeQS*sin(wm*t) * sin(e*wc*t) + G_ijfQ*sin(f*wc*t) * V_jeIS*sin(wm*t) * cos(e*wc*t)
                // + 0.5*G_ijfI*sin(sigma*wc*t) * V_jeQS*sin(wm*t) + 0.5*G_ijfQ*sin(sigma*wc*t) * V_jeIS*sin(wm*t)
                // RHS has negative sign
                harmonicSpaceMatrix_omegaC_0_->block(4*sigma+1,4*e+1)[i][j] = -(-omega_*f)*Cf_[i]->block(j)[2*f+1]; // quadrature sine translated from quadrature sine by in-phase G
                harmonicSpaceMatrix_omegaC_0_->block(4*sigma+1,4*e-1)[i][j] = -(+omega_*f)*Cf_[i]->block(j)[2*f]; // quadrature sine translated from in-phase sine by quadrature G
              }

              // case delta
              // now I have to deal with the index issue! if delta==0 then the indices become negative!
              if (deltaAbs==0)
              {
                // translation from harmonics to baseband
                // + G_ijfI*cos(f*wc*t) * V_jeIC*cos(wm*t) * cos(e*wc*t) + G_ijfQ*sin(f*wc*t) * V_jeQC*cos(wm*t) * sin(e*wc*t)
                // + 0.5*G_ijfI * V_jeIC*cos(wm*t) + 0.5*G_ijfQ * V_jeQC*cos(wm*t)
                // RHS has positive sign
                harmonicSpaceMatrix_omegaC_0_->block(0, 4*e-2)[i][j] = +(-omega_*f)*Cf_[i]->block(j)[2*f+1]; // baseband cosine translated from in-phase cosine by in-phase G
                harmonicSpaceMatrix_omegaC_0_->block(1, 4*e  )[i][j] = +(+omega_*f)*Cf_[i]->block(j)[2*f]; // baseband cosine translated from quadrature cosine by quadrature G

                // - G_ijfI*cos(f*wc*t) * V_jeIS*sin(wm*t) * cos(e*wc*t) - G_ijfQ*sin(f*wc*t) * V_jeQS*sin(wm*t) * sin(e*wc*t)
                // - 0.5*G_ijfI * V_jeIS*sin(wm*t) - 0.5*G_ijfQ * V_jeQS*sin(wm*t)
                // RHS has negative sign
                harmonicSpaceMatrix_omegaC_0_->block(0, 4*e-1)[i][j] = +(-omega_*f)*Cf_[i]->block(j)[2*f+1]; // baseband sine translated from in-phase sine by in-phase G
                harmonicSpaceMatrix_omegaC_0_->block(1, 4*e+1)[i][j] = +(+omega_*f)*Cf_[i]->block(j)[2*f]; // baseband sine translated from quadrature sine by quadrature G
              }
              else
              {
                // + G_ijfI*cos(f*wc*t) * V_jeIC*cos(wm*t) * cos(e*wc*t) + G_ijfQ*sin(f*wc*t) * V_jeQC*cos(wm*t) * sin(e*wc*t)
                // + 0.5*G_ijfI*cos(delta*wc*t) * V_jeIC*cos(wm*t) + 0.5*G_ijfQ*cos(delta*wc*t) * V_jeQC*cos(wm*t)
                // RHS has positive sign
                harmonicSpaceMatrix_omegaC_0_->block(4*deltaAbs-2, 4*e-2)[i][j] = +(-omega_*f)*Cf_[i]->block(j)[2*f+1]; // in-phase cosine translated from in-phase cosine by in-phase G
                harmonicSpaceMatrix_omegaC_0_->block(4*deltaAbs-2, 4*e  )[i][j] = +(+omega_*f)*Cf_[i]->block(j)[2*f]; // in-phase cosine translated from quadrature cosine by quadrature G

                // - G_ijfI*cos(f*wc*t) * V_jeIS*sin(wm*t) * cos(e*wc*t) - G_ijfQ*sin(f*wc*t) * V_jeQS*sin(wm*t) * sin(e*wc*t)
                // - 0.5*G_ijfI*cos(delta*wc*t) * V_jeIS*sin(wm*t) - 0.5*G_ijfQ*cos(delta*wc*t) * V_jeQS*sin(wm*t)
                // RHS has negative sign
                harmonicSpaceMatrix_omegaC_0_->block(4*deltaAbs-1,4*e-1)[i][j] = +(-omega_*f)*Cf_[i]->block(j)[2*f+1]; // in-phase sine translated from in-phase sine by in-phase G
                harmonicSpaceMatrix_omegaC_0_->block(4*deltaAbs-1,4*e+1)[i][j] = +(+omega_*f)*Cf_[i]->block(j)[2*f]; // in-phase sine translated from quadrature sine by quadrature G

                // - G_ijfI*cos(f*wc*t) * V_jeQC*cos(wm*t) * sin(e*wc*t) - G_ijfQ*sin(f*wc*t) * V_jeIC*cos(wm*t) * cos(e*wc*t)
                // + 0.5*G_ijfI*sin(delta*wc*t) * V_jeQC*cos(wm*t) - 0.5*G_ijfQ*sin(delta*wc*t) * V_jeIC*cos(wm*t)
                // RHS has positive sign
                harmonicSpaceMatrix_omegaC_0_->block(4*deltaAbs  ,4*e  )[i][j] = +sign*(-omega_*f)*Cf_[i]->block(j)[2*f+1]; // quadrature cosine translated from quadrature cosine by in-phase G
                harmonicSpaceMatrix_omegaC_0_->block(4*deltaAbs  ,4*e-2)[i][j] = -sign*(+omega_*f)*Cf_[i]->block(j)[2*f]; // quadrature cosine translated from in-phase cosine by quadrature G

                // + G_ijfI*cos(f*wc*t) * V_jeQS*sin(wm*t) * sin(e*wc*t) + G_ijfQ*sin(f*wc*t) * V_jeIS*sin(wm*t) * cos(e*wc*t)
                // - 0.5*G_ijfI*sin(delta*wc*t) * V_jeQS*sin(wm*t) + 0.5*G_ijfQ*sin(delta*wc*t) * V_jeIS*sin(wm*t)
                // RHS has negative sign
                harmonicSpaceMatrix_omegaC_0_->block(4*deltaAbs+1,4*e+1)[i][j] = +sign*(-omega_*f)*Cf_[i]->block(j)[2*f+1]; // quadrature sine translated from quadrature sine by in-phase G
                harmonicSpaceMatrix_omegaC_0_->block(4*deltaAbs+1,4*e-1)[i][j] = -sign*(+omega_*f)*Cf_[i]->block(j)[2*f]; // quadrature sine translated from in-phase sine by quadrature G
              }
            } // end of f>0 and e>0
          } // end of e loop
        } // end of f!=0
      } // end of f loop of Gf_
    } // end of column loop
  } // end of row loop

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
  // which parameters should be swapped to reuse the G matrix loop?
  harmonicSpaceMatrix_omegaC_1_ = Xyce::Linear::createBlockMatrix( numBlocks, offset, blockPattern, blockGraph.get(), baseFullGraph);
  harmonicSpaceMatrix_omegaC_1_->put( 0.0 );

  for (int i=0; i<numRows; i++)
  { // selectring a row of C Matrix
    for (int j=0; j<numRows; j++)
    { // selecting a column of C Matrix
      for (int f=0; f<=numHarms; f++)
      { // selecting a frequency point of Cf_
        // now we fill the matrix
        if (f==0)
        { // linear transformation
          // C_ij0 * [ ( V_jeQC*cos(wm*t) - V_jeQS*sin(wm*t) )*e*wc*cos(e*wc*t) - ( V_jeIC*cos(wm*t) - V_jeIS*sin(wm*t) )*e*wc*sin(e*wc*t) ]
          for (int e=1; e<numHarms; e++){
            // the first two blocks are baseband and zero, hence starting from e=1
            // + C_ij0 * V_jeQC*cos(wm*t) * e*wc*cos(e*wc*t) 
            // RHS has positive sign
            harmonicSpaceMatrix_omegaC_1_->block(4*e-2, 4*e  )[i][j] = +(e*omega_)*Cf_[i]->block(j)[0]; // in-phase cosine

            // - C_ij0 * V_jeQS*sin(wm*t) * e*wc*cos(e*wc*t)
            // RHS has negative sign
            harmonicSpaceMatrix_omegaC_1_->block(4*e-1, 4*e+1)[i][j] = +(e*omega_)*Cf_[i]->block(j)[0]; // in-phase sine

            // - C_ij0 * V_jeIC*cos(wm*t) * e*wc*sin(e*wc*t) ]
            // RHS has positive sign
            harmonicSpaceMatrix_omegaC_1_->block(4*e  , 4*e-2)[i][j] = -(e*omega_)*Cf_[i]->block(j)[0]; // quadrature cosine

            // + C_ij0 * V_jeIS*sin(wm*t) * e*wc*sin(e*wc*t) ]
            // RHS has negative sign
            harmonicSpaceMatrix_omegaC_1_->block(4*e+1, 4*e-1)[i][j] = -(e*omega_)*Cf_[i]->block(j)[0]; // quadrature sine
          }
        } else
        { // now the mixing parts (harmonic coupling)
          for (int e=1; e<=numHarms; e++)
          { // now f>0 and e>0
            // [ C_ijfI*cos(f*wc*t) - C_ijfQ*sin(f*wc*t) ] * [ - ( V_jeQC*cos(wm*t) - V_jeQS*sin(wm*t) )*cos(e*wc*t) - ( V_jeIC*cos(wm*t) - V_jeIS*sin(wm*t) )*sin(e*wc*t) ]*e*wc
            int sigma = f+e;
            // int delta = f-e; this line is just for the sake of understanding the code. Equations are with reference to delta, not deltaAbs.
            int deltaAbs = std::abs(f-e);
            int sign = f>=e ? 1 : -1;

            // case sigma
            if (sigma<=numHarms) 
            { // f+e should not be larget than numHarms, otherwise ignore it
              // in finite-harmonics space the system stil shows nonlineary and some mixing products have to be ignored

              // - C_ijfI*cos(f*wc*t) * V_jeQC*cos(wm*t) * cos(e*wc*t) + C_ijfQ*sin(f*wc*t) * V_jeIC*cos(wm*t) * sin(e*wc*t)  // all *e*wc
              // - 0.5*C_ijfI*cos(sigma*wc*t) * V_jeQC*cos(wm*t) - 0.5*C_ijfQ*cos(sigma*wc*t) * V_jeIC*cos(wm*t)              // all *e*wc
              // RHS has positive sign
              harmonicSpaceMatrix_omegaC_1_->block(4*sigma-2, 4*e  )[i][j] = -(e*omega_)*Cf_[i]->block(j)[2*f  ]; // in-phase cosine translated from quadrature cosine
              harmonicSpaceMatrix_omegaC_1_->block(4*sigma-2, 4*e-2)[i][j] = -(e*omega_)*Cf_[i]->block(j)[2*f+1]; // in-phase cosine translated from in-phase cosine

              // + C_ijfI*cos(f*wc*t) * V_jeQS*sin(wm*t) * cos(e*wc*t) - C_ijfQ*sin(f*wc*t) * V_jeIS*sin(wm*t) * sin(e*wc*t)  // all *e*wc
              // + 0.5*C_ijfI*cos(sigma*wc*t) * V_jeQS*sin(wm*t) + 0.5*C_ijfQ*cos(sigma*wc*t) * V_jeIS*sin(wm*t)              // all *e*wc
              // RHS has negative sign
              harmonicSpaceMatrix_omegaC_1_->block(4*sigma-1,4*e+1)[i][j] = -(e*omega_)*Cf_[i]->block(j)[2*f  ]; // in-phase sine translated from quadrature sine
              harmonicSpaceMatrix_omegaC_1_->block(4*sigma-1,4*e-1)[i][j] = -(e*omega_)*Cf_[i]->block(j)[2*f+1]; // in-phase sine translated from in-phase sine

              // - C_ijfI*cos(f*wc*t) * V_jeIC*cos(wm*t) * sin(e*wc*t) + C_ijfQ*sin(f*wc*t) * V_jeQC*cos(wm*t) * cos(e*wc*t)  // all *e*wc
              // - 0.5*C_ijfI*sin(sigma*wc*t) * V_jeIC*cos(wm*t) + 0.5*C_ijfQ*sin(sigma*wc*t) * V_jeQC*cos(wm*t)              // all *e*wc
              // RHS has positive sign
              harmonicSpaceMatrix_omegaC_1_->block(4*sigma  ,4*e-2)[i][j] = -(e*omega_)*Cf_[i]->block(j)[2*f  ]; // quadrature cosine translated from in-phase cosine
              harmonicSpaceMatrix_omegaC_1_->block(4*sigma  ,4*e  )[i][j] = +(e*omega_)*Cf_[i]->block(j)[2*f+1]; // quadrature cosine translated from quadrature cosine

              // + C_ijfI*cos(f*wc*t) * V_jeIS*sin(wm*t) * sin(e*wc*t) - C_ijfQ*sin(f*wc*t) * V_jeQS*sin(wm*t) * cos(e*wc*t)  // all *e*wc
              // + 0.5*C_ijfI*sin(sigma*wc*t) * V_jeIS*sin(wm*t) - 0.5*C_ijfQ*sin(sigma*wc*t) * V_jeQS*sin(wm*t)              // all *e*wc
              // RHS has negative sign
              harmonicSpaceMatrix_omegaC_1_->block(4*sigma+1,4*e-1)[i][j] = -(e*omega_)*Cf_[i]->block(j)[2*f  ]; // quadrature sine translated from in-phase sine
              harmonicSpaceMatrix_omegaC_1_->block(4*sigma+1,4*e+1)[i][j] = +(e*omega_)*Cf_[i]->block(j)[2*f+1]; // quadrature sine translated from quadrature sine
            }

            // case delta
            // now I have to deal with the index issue! if delta==0 then the indices become negative!
            if (deltaAbs==0)
            {
              // translation from harmonics to baseband
              // - C_ijfI*cos(f*wc*t) * V_jeQC*cos(wm*t) * cos(e*wc*t) + C_ijfQ*sin(f*wc*t) * V_jeIC*cos(wm*t) * sin(e*wc*t)  // all *e*wc
              // - 0.5*C_ijfI * V_jeQC*cos(wm*t) + 0.5*C_ijfQ * V_jeIC*cos(wm*t)                                              // all *e*wc
              // RHS has positive sign
              harmonicSpaceMatrix_omegaC_1_->block(0, 4*e  )[i][j] = -(e*omega_)*Cf_[i]->block(j)[2*f  ]; // baseband cosine translated from quadrature cosine
              harmonicSpaceMatrix_omegaC_1_->block(1, 4*e-2)[i][j] = +(e*omega_)*Cf_[i]->block(j)[2*f+1]; // baseband cosine translated from in-phase cosine

              // + C_ijfI*cos(f*wc*t) * V_jeQS*sin(wm*t) * cos(e*wc*t) - C_ijfQ*sin(f*wc*t) * V_jeIS*sin(wm*t) * sin(e*wc*t)  // all *e*wc
              // + 0.5*C_ijfI * V_jeQS*sin(wm*t) - 0.5*C_ijfQ * V_jeIS*sin(wm*t)                                              // all *e*wc
              // RHS has negative sign
              harmonicSpaceMatrix_omegaC_1_->block(0, 4*e+1)[i][j] = -(e*omega_)*Cf_[i]->block(j)[2*f  ]; // baseband sine translated from quadrature sine
              harmonicSpaceMatrix_omegaC_1_->block(1, 4*e-1)[i][j] = +(e*omega_)*Cf_[i]->block(j)[2*f+1]; // baseband sine translated from in-phase sine
            }
            else
            {
              // - C_ijfI*cos(f*wc*t) * V_jeQC*cos(wm*t) * cos(e*wc*t) + C_ijfQ*sin(f*wc*t) * V_jeIC*cos(wm*t) * sin(e*wc*t)  // all *e*wc
              // - 0.5*C_ijfI*cos(delta*wc*t) * V_jeQC*cos(wm*t) - 0.5*C_ijfQ*cos(delta*wc*t) * V_jeIC*cos(wm*t)              // all *e*wc
              // RHS has positive sign
              harmonicSpaceMatrix_omegaC_1_->block(4*deltaAbs-2, 4*e  )[i][j] = -(e*omega_)*Cf_[i]->block(j)[2*f  ]; // in-phase cosine translated from quadrature cosine
              harmonicSpaceMatrix_omegaC_1_->block(4*deltaAbs-2, 4*e-2)[i][j] = -(e*omega_)*Cf_[i]->block(j)[2*f+1]; // in-phase cosine translated from in-phase cosine

              // + C_ijfI*cos(f*wc*t) * V_jeQS*sin(wm*t) * cos(e*wc*t) - C_ijfQ*sin(f*wc*t) * V_jeIS*sin(wm*t) * sin(e*wc*t)  // all *e*wc
              // + 0.5*C_ijfI*cos(delta*wc*t) * V_jeQS*sin(wm*t) + 0.5*C_ijfQ*cos(delta*wc*t) * V_jeIS*sin(wm*t)              // all *e*wc
              // RHS has negative sign
              harmonicSpaceMatrix_omegaC_1_->block(4*deltaAbs-1,4*e+1)[i][j] = -(e*omega_)*Cf_[i]->block(j)[2*f  ]; // in-phase sine translated from quadrature sine
              harmonicSpaceMatrix_omegaC_1_->block(4*deltaAbs-1,4*e-1)[i][j] = -(e*omega_)*Cf_[i]->block(j)[2*f+1]; // in-phase sine translated from in-phase sine

              // - C_ijfI*cos(f*wc*t) * V_jeIC*cos(wm*t) * sin(e*wc*t) + C_ijfQ*sin(f*wc*t) * V_jeQC*cos(wm*t) * cos(e*wc*t)  // all *e*wc
              // - 0.5*C_ijfI*sin(delta*wc*t) * V_jeIC*cos(wm*t) + 0.5*C_ijfQ*sin(delta*wc*t) * V_jeQC*cos(wm*t)              // all *e*wc
              // RHS has positive sign
              harmonicSpaceMatrix_omegaC_1_->block(4*deltaAbs  ,4*e-2)[i][j] = -(e*omega_)*sign*Cf_[i]->block(j)[2*f  ]; // quadrature cosine translated from in-phase cosine
              harmonicSpaceMatrix_omegaC_1_->block(4*deltaAbs  ,4*e  )[i][j] = +(e*omega_)*sign*Cf_[i]->block(j)[2*f+1]; // quadrature cosine translated from quadrature cosine

              // + C_ijfI*cos(f*wc*t) * V_jeIS*sin(wm*t) * sin(e*wc*t) - C_ijfQ*sin(f*wc*t) * V_jeQS*sin(wm*t) * cos(e*wc*t)  // all *e*wc
              // + 0.5*C_ijfI*sin(delta*wc*t) * V_jeIS*sin(wm*t) - 0.5*C_ijfQ*sin(delta*wc*t) * V_jeQS*sin(wm*t)              // all *e*wc
              // RHS has negative sign
              harmonicSpaceMatrix_omegaC_1_->block(4*deltaAbs+1,4*e-1)[i][j] = -(e*omega_)*sign*Cf_[i]->block(j)[2*f  ]; // quadrature sine translated from in-phase sine
              harmonicSpaceMatrix_omegaC_1_->block(4*deltaAbs+1,4*e+1)[i][j] = +(e*omega_)*sign*Cf_[i]->block(j)[2*f+1]; // quadrature sine translated from quadrature sine
            }
            // end of f>0 and e>0
          } // end of e loop
        } // end of f!=0
      } // end of f loop of Gf_
    } // end of column loop
  } // end of row loop

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
  // This matrix has also similar signature to G, but is does inner-harmonic transformation.
  harmonicSpaceMatrix_C_2_ = Xyce::Linear::createBlockMatrix( numBlocks, offset, blockPattern, blockGraph.get(), baseFullGraph);
  harmonicSpaceMatrix_C_2_->put( 0.0 );

  for (int i=0; i<numRows; i++)
  { // selectring a row of C Matrix
    for (int j=0; j<numRows; j++)
    { // selecting a column of C Matrix
      for (int f=0; f<=numHarms; f++)
      { // selecting a frequency point of Cf_
        // now we fill the matrix
        // first the diagon of G matrix
        if (f==0)
        { // linear transformation
          // C_ij0 * [ - V_jS*cos(wm*t) - V_jC*sin(wm*t) + ( - V_jeIS*cos(wm*t) - V_jeIC*sin(wm*t) )*cos(e*wc*t) - ( - V_jeQS*cos(wm*t) - V_jeQC*sin(wm*t) )*sin(e*wc*t) ]
          for (int e=0; e<numHarms; e++){
            if (e==0)
            { // linear transformation at baseband
              // - C_ij0 * V_jS*cos(wm*t)
              // RHS has positive sign
              harmonicSpaceMatrix_C_2_->block(0, 1)[i][j] = -Cf_[i]->block(j)[0]; // in-phase cosine

              // - C_ij0 * V_jC*sin(wm*t)
              // RHS has negative sign
              harmonicSpaceMatrix_C_2_->block(1, 0)[i][j] = +Cf_[i]->block(j)[0]; // in-phase sine
            } else
            { // linear transformation at harmonics 
              // - C_ij0 * V_jeIS*cos(wm*t) * cos(e*wc*t)
              // RHS has positive sign
              harmonicSpaceMatrix_C_2_->block(4*e-2, 4*e-1)[i][j] = -Cf_[i]->block(j)[0]; // in-phase cosine

              // - C_ij0 * V_jeIC*sin(wm*t) * cos(e*wc*t)
              // RHS has negative sign
              harmonicSpaceMatrix_C_2_->block(4*e-1, 4*e-2)[i][j] = +Cf_[i]->block(j)[0]; // in-phase sine

              // + C_ij0 * V_jeQS*cos(wm*t) * sin(e*wc*t)
              // RHS has positive sign
              harmonicSpaceMatrix_C_2_->block(4*e  , 4*e+1)[i][j] = +Cf_[i]->block(j)[0]; // quadrature cosine

              // + C_ij0 * V_jeQC*sin(wm*t) * sin(e*wc*t)
              // RHS has negative sign
              harmonicSpaceMatrix_C_2_->block(4*e+1, 4*e  )[i][j] = -Cf_[i]->block(j)[0]; // quadrature sine
            }
          }
        } else
        { // now the mixing parts (harmonic coupling)
          for (int e=0; e<=numHarms; e++)
          { // harmonic index of v
            if (e==0) 
            { // baseband modulation
              // [ C_ijfI*cos(f*wc*t) - C_ijfQ*sin(f*wc*t) ] * [ - V_jS*cos(wm*t) - V_jC*sin(wm*t) ]

              // - C_ijfI*cos(f*wc*t) * V_jS*cos(wm*t)
              // RHS has positive sign
              harmonicSpaceMatrix_C_2_->block(4*f-2,0)[i][j] = -2*Cf_[i]->block(j)[2*f]; // in-phase cosine

              // - C_ijfI*cos(f*wc*t) * V_jC*sin(wm*t)
              // RHS has negative sign
              harmonicSpaceMatrix_C_2_->block(4*f-1,0)[i][j] = +2*Cf_[i]->block(j)[2*f]; // in-phase sine

              // + C_ijfQ*sin(f*wc*t) * V_jC*cos(wm*t)
              // RHS has positive sign
              harmonicSpaceMatrix_C_2_->block(4*f  ,0)[i][j] = +2*Cf_[i]->block(j)[2*f+1]; // quadrature cosine

              // + C_ijfQ*sin(f*wc*t) * V_jS*sin(wm*t)
              // RHS has negative sign
              harmonicSpaceMatrix_C_2_->block(4*f+1,0)[i][j] = -2*Cf_[i]->block(j)[2*f+1]; // quadrature sine
            } else 
            { // now f>0 and e>0
              // [ C_ijfI*cos(f*wc*t) - C_ijfQ*sin(f*wc*t) ] * [ ( - V_jeIS*cos(wm*t) - V_jeIC*sin(wm*t) )*cos(e*wc*t) - ( - V_jeQS*cos(wm*t) - V_jeQC*sin(wm*t) )*sin(e*wc*t) ]
              int sigma = f+e;
              // int delta = f-e; this line is just for the sake of understanding the code. Equations are with reference to delta, not deltaAbs.
              int deltaAbs = std::abs(f-e);
              int sign = f>=e ? 1 : -1;

              // case sigma
              if (sigma<=numHarms) 
              { // f+e should not be larget than numHarms, otherwise ignore it
                // in finite-harmonics space the system stil shows nonlineary and some mixing products have to be ignored

                // - C_ijfI*cos(f*wc*t) * V_jeIS*cos(wm*t) * cos(e*wc*t) - C_ijfQ*sin(f*wc*t) * V_jeQS*cos(wm*t) * sin(e*wc*t)
                // - 0.5*C_ijfI*cos(sigma*wc*t) * V_jeIS*cos(wm*t) + 0.5*C_ijfQ*cos(sigma*wc*t) * V_jeQS*cos(wm*t)
                // RHS has positive sign
                harmonicSpaceMatrix_C_2_->block(4*sigma-2, 4*e-1)[i][j] = -Cf_[i]->block(j)[2*f  ]; // in-phase cosine translated from in-phase sine by d/dt
                harmonicSpaceMatrix_C_2_->block(4*sigma-2, 4*e+1)[i][j] = +Cf_[i]->block(j)[2*f+1]; // in-phase cosine translated from quadrature sine by d/dt

                // - C_ijfI*cos(f*wc*t) * V_jeIC*sin(wm*t) * cos(e*wc*t) - C_ijfQ*sin(f*wc*t) * V_jeQC*sin(wm*t) * sin(e*wc*t)
                // - 0.5*C_ijfI*cos(sigma*wc*t) * V_jeIC*sin(wm*t) + 0.5*C_ijfQ*cos(sigma*wc*t) * V_jeQC*sin(wm*t)
                // RHS has negative sign
                harmonicSpaceMatrix_C_2_->block(4*sigma-1,4*e-2)[i][j] = +Cf_[i]->block(j)[2*f  ]; // in-phase sine translated from in-phase cosine by d/dt
                harmonicSpaceMatrix_C_2_->block(4*sigma-1,4*e  )[i][j] = -Cf_[i]->block(j)[2*f+1]; // in-phase sine translated from quadrature cosine by d/dt

                // + C_ijfI*cos(f*wc*t) * V_jeQS*cos(wm*t) * sin(e*wc*t) + C_ijfQ*sin(f*wc*t) * V_jeIS*cos(wm*t) * cos(e*wc*t)
                // + 0.5*C_ijfI*sin(sigma*wc*t) * V_jeQS*cos(wm*t) + 0.5*C_ijfQ*sin(sigma*wc*t) * V_jeIS*cos(wm*t)
                // RHS has positive sign
                harmonicSpaceMatrix_C_2_->block(4*sigma  ,4*e+1)[i][j] = +Cf_[i]->block(j)[2*f  ]; // quadrature cosine translated from quadrature sine by d/dt
                harmonicSpaceMatrix_C_2_->block(4*sigma  ,4*e-1)[i][j] = +Cf_[i]->block(j)[2*f+1]; // quadrature cosine translated from in-phase sine by d/dt

                // + C_ijfI*cos(f*wc*t) * V_jeQC*sin(wm*t) * sin(e*wc*t) + C_ijfQ*sin(f*wc*t) * V_jeIC*sin(wm*t) * cos(e*wc*t)
                // + 0.5*C_ijfI*sin(sigma*wc*t) * V_jeQC*sin(wm*t) + 0.5*C_ijfQ*sin(sigma*wc*t) * V_jeIC*sin(wm*t)
                // RHS has negative sign
                harmonicSpaceMatrix_C_2_->block(4*sigma+1,4*e  )[i][j] = -Cf_[i]->block(j)[2*f  ]; // quadrature sine translated from quadrature cosine by d/dt
                harmonicSpaceMatrix_C_2_->block(4*sigma+1,4*e-2)[i][j] = -Cf_[i]->block(j)[2*f+1]; // quadrature sine translated from in-phase cosine by d/dt
              }

              // case delta
              // now I have to deal with the index issue! if delta==0 then the indices become negative!
              if (deltaAbs==0)
              {
                // translation from harmonics to baseband
                // - C_ijfI*cos(f*wc*t) * V_jeIS*cos(wm*t) * cos(e*wc*t) - C_ijfQ*sin(f*wc*t) * V_jeQS*cos(wm*t) * sin(e*wc*t)
                // - 0.5*C_ijfI * V_jeIS*cos(wm*t) + 0.5*C_ijfQ * V_jeQS*cos(wm*t)
                // RHS has positive sign
                harmonicSpaceMatrix_C_2_->block(0, 4*e-1)[i][j] = -Cf_[i]->block(j)[2*f  ]; // in-phase cosine translated from in-phase sine by d/dt
                harmonicSpaceMatrix_C_2_->block(0, 4*e+1)[i][j] = +Cf_[i]->block(j)[2*f+1]; // in-phase cosine translated from quadrature sine by d/dt

                // - C_ijfI*cos(f*wc*t) * V_jeIC*sin(wm*t) * cos(e*wc*t) - C_ijfQ*sin(f*wc*t) * V_jeQC*sin(wm*t) * sin(e*wc*t)
                // - 0.5*C_ijfI * V_jeIC*sin(wm*t) + 0.5*C_ijfQ * V_jeQC*sin(wm*t)
                // RHS has negative sign
                harmonicSpaceMatrix_C_2_->block(1,4*e-2)[i][j] = +Cf_[i]->block(j)[2*f  ]; // in-phase sine translated from in-phase cosine by d/dt
                harmonicSpaceMatrix_C_2_->block(1,4*e  )[i][j] = -Cf_[i]->block(j)[2*f+1]; // in-phase sine translated from quadrature cosine by d/dt
              }
              else
              {
                // - C_ijfI*cos(f*wc*t) * V_jeIS*cos(wm*t) * cos(e*wc*t) - C_ijfQ*sin(f*wc*t) * V_jeQS*cos(wm*t) * sin(e*wc*t)
                // - 0.5*C_ijfI*cos(delta*wc*t) * V_jeIS*cos(wm*t) + 0.5*C_ijfQ*cos(delta*wc*t) * V_jeQS*cos(wm*t)
                // RHS has positive sign
                harmonicSpaceMatrix_C_2_->block(4*deltaAbs-2, 4*e-1)[i][j] = -Cf_[i]->block(j)[2*f  ]; // in-phase cosine translated from in-phase sine by d/dt
                harmonicSpaceMatrix_C_2_->block(4*deltaAbs-2, 4*e+1)[i][j] = +Cf_[i]->block(j)[2*f+1]; // in-phase cosine translated from quadrature sine by d/dt

                // - C_ijfI*cos(f*wc*t) * V_jeIC*sin(wm*t) * cos(e*wc*t) - C_ijfQ*sin(f*wc*t) * V_jeQC*sin(wm*t) * sin(e*wc*t)
                // - 0.5*C_ijfI*cos(delta*wc*t) * V_jeIC*sin(wm*t) + 0.5*C_ijfQ*cos(delta*wc*t) * V_jeQC*sin(wm*t)
                // RHS has negative sign
                harmonicSpaceMatrix_C_2_->block(4*deltaAbs-1,4*e-2)[i][j] = +Cf_[i]->block(j)[2*f  ]; // in-phase sine translated from in-phase cosine by d/dt
                harmonicSpaceMatrix_C_2_->block(4*deltaAbs-1,4*e  )[i][j] = -Cf_[i]->block(j)[2*f+1]; // in-phase sine translated from quadrature cosine by d/dt

                // + C_ijfI*cos(f*wc*t) * V_jeQS*cos(wm*t) * sin(e*wc*t) + C_ijfQ*sin(f*wc*t) * V_jeIS*cos(wm*t) * cos(e*wc*t)
                // + 0.5*C_ijfI*sin(delta*wc*t) * V_jeQS*cos(wm*t) + 0.5*C_ijfQ*sin(delta*wc*t) * V_jeIS*cos(wm*t)
                // RHS has positive sign
                harmonicSpaceMatrix_C_2_->block(4*deltaAbs  ,4*e+1)[i][j] = +sign*Cf_[i]->block(j)[2*f  ]; // quadrature cosine translated from quadrature sine by d/dt
                harmonicSpaceMatrix_C_2_->block(4*deltaAbs  ,4*e-1)[i][j] = +sign*Cf_[i]->block(j)[2*f+1]; // quadrature cosine translated from in-phase sine by d/dt

                // + C_ijfI*cos(f*wc*t) * V_jeQC*sin(wm*t) * sin(e*wc*t) + C_ijfQ*sin(f*wc*t) * V_jeIC*sin(wm*t) * cos(e*wc*t)
                // + 0.5*C_ijfI*sin(delta*wc*t) * V_jeQC*sin(wm*t) + 0.5*C_ijfQ*sin(delta*wc*t) * V_jeIC*sin(wm*t)
                // RHS has negative sign
                harmonicSpaceMatrix_C_2_->block(4*deltaAbs+1,4*e  )[i][j] = -sign*Cf_[i]->block(j)[2*f  ]; // quadrature sine translated from quadrature cosine by d/dt
                harmonicSpaceMatrix_C_2_->block(4*deltaAbs+1,4*e-2)[i][j] = -sign*Cf_[i]->block(j)[2*f+1]; // quadrature sine translated from in-phase cosine by d/dt

              }
            } // end of f>0 and e>0
          } // end of e loop
        } // end of f!=0
      } // end of f loop of Gf_
    } // end of column loop
  } // end of row loop



  harmonicSpaceB_->putScalar( 0.0 );
  harmonicSpaceB_->block( 0 ).update( 1.0, *bVecRealPtr);
  harmonicSpaceB_->block( 1 ).update( 1.0, *bVecImagPtr);


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
    for (int i=0; i<BlockSize; i++){
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