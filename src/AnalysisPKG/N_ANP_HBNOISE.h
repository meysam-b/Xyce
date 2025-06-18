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
// Purpose       : This is the Harmonic Balance Noise Analysis class
// Special Notes :
// Creator       : Meysam Bahmanian
// Creation Date : 5/13/2025
//-----------------------------------------------------------------------------

#ifndef Xyce_N_ANP_HBNOISE_h
#define Xyce_N_ANP_HBNOISE_h

#include <vector>
#include <map>

#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_RCP.hpp>

#include <N_ANP_fwd.h>
#include <N_LOA_fwd.h>
#include <N_LAS_fwd.h>
#include <N_TOP_fwd.h>

#include <N_ANP_AnalysisBase.h>
#include <N_ANP_RegisterAnalysis.h>
#include <N_ANP_StepEvent.h>
#include <N_UTL_DFTInterfaceDecl.hpp>
#include <N_UTL_FFTInterface.hpp>
#include <N_UTL_Listener.h>
#include <N_UTL_Op.h>
#include <N_UTL_OptionBlock.h>

namespace Xyce {
namespace Analysis {

class HBNOISE : public AnalysisBase, public Util::ListenerAutoSubscribe<StepEvent>
{
public:
  HBNOISE(
    AnalysisManager &                     analysis_manager,
    Linear::System &                      linear_system,
    Nonlinear::Manager &                  nonlinear_manager,
    Loader::Loader &                      loader,
    Device::DeviceMgr &                   device_manager,
    Topo::Topology &                      topology
    );

  virtual ~HBNOISE();

  void notify(const StepEvent &event);

  // Existing virtual functions
  virtual bool processSuccessfulDCOP();
  virtual bool processFailedDCOP();
  virtual void finalExpressionBasedSetup();
  virtual TimeIntg::TIAParams &getTIAParams();
  virtual const TimeIntg::TIAParams &getTIAParams() const;
  virtual bool doRun();
  virtual bool getDCOPFlag() const;
  virtual bool convertDataToSweepParams();

  // Add parameter setting methods
  bool setAnalysisParams(const Util::OptionBlock & paramsBlock);
  bool setLinSol(const Util::OptionBlock & OB);
  bool setDataStatements(const Util::OptionBlock & paramsBlock);

protected:
  virtual bool doInit();
  virtual bool doLoopProcess();
  virtual bool doProcessSuccessfulStep();
  virtual bool doProcessFailedStep();
  virtual bool doFinish();
  virtual bool doHandlePredictor();

private:
  Analysis::HB* getHBAnalysis();
  int setupSweepParam_();
  bool updateDataParams_(int stepNumber);
  bool updateCurrentFreq_(int stepNumber);

  // This helper function should be moved to a more appropriate place
  void setMatrixElement(Linear::Matrix& mat, int row, int col, double value);

  bool updateLinearTimeVariantSystem_C_and_G_();
  bool createHarmonicSpaceLinearSystem_();
  bool updateHarmonicSpaceMatrix_G_();
  bool updateHarmonicSpaceMatrix_omegaC_0_();
  bool updateHarmonicSpaceMatrix_omegaC_1_();
  bool updateHarmonicSpaceMatrix_C_2_();
  bool updateHarmonicSpaceFreq_(); // update the harmonic space matrix for the frequency
  void resetAdjointHBNOISELinearSystem_();

  void setupAdjointRHS_();
  bool solveAdjointHBNOISE_();
  void prepareHBNOISEOutputVectors_(
    Linear::BlockVector *           harmonicSpaceX,
    std::vector<std::vector<Xyce::Analysis::NoiseData*> > &noiseDataVecVec,
    std::vector<Xyce::Analysis::NoiseData*> &noiseDataVec,
    double &totalRelativeNoiseDens);

  void processOutputNodes ();
  
  void clearNoiseIntegrals_();
  inline void evalDeviceNoiseDensities(
    bool isBaseband,
    Xyce::Analysis::NoiseData& noiseData, 
    Linear::Vector& XIreal, 
    Linear::Vector& XIimag,
    double& totalRelativeNoiseDens);

private:
  AnalysisManager &                     analysisManager_;
  Loader::Loader &                      loader_;
  Linear::System &                      linearSystem_;
  Nonlinear::Manager &                  nonlinearManager_;
  Device::DeviceMgr &                   deviceManager_;
  Topo::Topology &                      topology_;
  OutputMgrAdapter &                    outputManagerAdapter_;
  Parallel::Manager *                   pdsMgrPtr_;
  AnalysisBase *                        currentAnalysisObject_;
  Loader::HBLoader *                    hbLoaderPtr_; /// HB loader, builder, system, and DFT
  Teuchos::RCP<Linear::HBBuilder>       hbBuilderPtr_;
  Linear::Builder *                     builderPtr_;
  Linear::System *                      hbLinearSystem_;
  
  // HBNOISE specific parameters
  bool outputNodeSingle_;              // Flag for single output node
  std::string outputNode1_;            // First output node
  std::string outputNode2_;            // Second output node
  double harmonicNumber_;              // Harmonic number for hbnoise analysis
  std::string type_;                   // Type of sweep (LIN, DEC, OCT)
  double np_;                          // Number of points
  double fOffsetStart_;                // Start offset frequency for hbnoise analysis
  double fOffsetStop_;                 // Stop offset frequency for hbnoise analysis
  double stepMult_;                    // Multiplier for frequency steps (DEC/OCT)
  double fstep_;                       // Step size for frequency (LIN)
  int pts_per_summary_;                // Points per summary
  bool dataSpecification_;             // Flag for data specification
  int hbnoiseLoopSize_;                // Size of the hbnoise analysis loop
  SweepVector hbnoiseSweepVector_;       // Vector of sweep parameters
  std::map< std::string, std::vector<std::string> > dataNamesMap_;  // Maps dataset name to parameter names
  std::map< std::string, std::vector< std::vector<double> > > dataTablesMap_;  // Maps dataset name to parameter values
  Analysis::HB *hbAnalysis_;

  double                        delFreq_;
  double                        lastFreq_;
  double                        currentFreq_;
  double                        lnFreq_;
  double                        lnLastFreq_;
  double                        delLnFreq_;

  double                        totalAMNoiseDens_;
  double                        totalPMNoiseDens_;

  // NOISE B-vectors
  Linear::Vector * bNoiseVecPtr;
  // Linear::Vector * bNoiseVecRealPtr;
  // Linear::Vector * bNoiseVecImagPtr;

  //time domain matrices
  std::vector<Teuchos::RCP<Linear::BlockVector> > Ct_;
  std::vector<Teuchos::RCP<Linear::BlockVector> > Gt_;

  //frequency domain matrices
  std::vector<Teuchos::RCP<Linear::BlockVector> > Cf_;
  std::vector<Teuchos::RCP<Linear::BlockVector> > Gf_;

  Linear::BlockMatrix *           harmonicSpaceMatrix_;
  Linear::BlockMatrix *           harmonicSpaceMatrixConstant_;
  Linear::BlockMatrix *           harmonicSpaceMatrix_G_; //conductance matrix
  Linear::BlockMatrix *           harmonicSpaceMatrix_omegaC_0_; // PART 0: dC(t)/dt
  Linear::BlockMatrix *           harmonicSpaceMatrix_omegaC_1_; // PART 1: carrier derivative 
  Linear::BlockMatrix *           harmonicSpaceMatrix_C_2_;   // PART 2: baseband/in-phase/quadrature components derivative
  Linear::BlockMatrix *           harmonicSpaceMatrix_omegamC_2_;   // PART 2: scaled by omegam
  Linear::BlockVector *           harmonicSpaceBI_; // in-phase
  Linear::BlockVector *           harmonicSpaceBQ_; // quadrature
  Linear::BlockVector *           harmonicSpaceXI_; // in-phase
  Linear::BlockVector *           harmonicSpaceXQ_; // quadrature
  Linear::BlockVector *           harmonicSpace_SavedXI_; // in-phase
  Linear::BlockVector *           harmonicSpace_SavedXQ_; // quadrature

  Linear::Solver *              blockSolverI_;
  Linear::Solver *              blockSolverQ_;
  Linear::Problem *             blockProblemI_;
  Linear::Problem *             blockProblemQ_;
  Util::OptionBlock             linSolOptionBlock_;

  std::vector<std::string> outputVarNames_;
  std::vector<int>    outputVarGIDs_;
  double outputValReal_;
  double outputValImag_;
  double outputValSqr_;
  double outputValCosPhi_;
  double outputValSinPhi_;

  double freq_;                                                 // primary frequency from HB analysis
  double omega_;                                                // primary angular frequency (2*pi*freq_) from HB analysis
  int                   numHarms_;                              // number of harmonics
  int                   size_;                                  // Problem Size: 2*harmonics+1
  double                period_;                                // Periodicity Information
  std::vector<double>                   times_;

  // hbnoise integrals are not calculated for DATA=<n> case if the
  // specified frequencies are not monotonically increasing
  bool calcNoiseIntegrals_;

  // Option blocks for parameters
  Util::OptionBlock saved_lsOB_;
  Util::OptionBlock saved_timeIntOB_;

  // noise contribution of each device summed over all noise harmonic frequencies for the in-phase output
  std::vector<Xyce::Analysis::NoiseData*> noiseDataVecI_;

  // noise contribution of each device summed over all noise harmonic frequencies for the quadrature output
  std::vector<Xyce::Analysis::NoiseData*> noiseDataVecQ_;

  // noise contribution of each device at each noise harmonic freuquency for the in-phase output
  // The element are the contributions from baseband
  // The other components are pairs of LSB and USB components around the carrier frequency
  // total number of elements is 2*numHarms+1
  // The first element of this vector is our old AC NOISE!
  std::vector< std::vector<Xyce::Analysis::NoiseData*> > noiseDataVecVecI_;

  // noise contribution of each device at each noise harmonic freuquency for the quadrature output
  // same as noiseDataVecVecI_ but for the quadrature output
  // when the desired output is baseband, there will be no quadrature noise.
  std::vector< std::vector<Xyce::Analysis::NoiseData*> > noiseDataVecVecQ_;
};

bool registerHBNOISEFactory(FactoryBlock &factory_block);

} // namespace Analysis
} // namespace Xyce

#endif // Xyce_N_ANP_HBNOISE_h 