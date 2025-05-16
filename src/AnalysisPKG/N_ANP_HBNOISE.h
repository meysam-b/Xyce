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
#include <N_ANP_StepEvent.h>
#include <N_ANP_RegisterAnalysis.h>
#include <N_UTL_DFTInterfaceDecl.hpp>
#include <N_UTL_FFTInterface.hpp>
#include <N_UTL_Listener.h>
#include <N_UTL_OptionBlock.h>

namespace Xyce {
namespace Analysis {

class HBNOISE : public AnalysisBase, public Util::ListenerAutoSubscribe<StepEvent>
{
public:
  HBNOISE(
    AnalysisManager &analysis_manager,
    Linear::System & linear_system,
    Nonlinear::Manager & nonlinear_manager,
    Loader::Loader &loader,
    Device::DeviceMgr & device_manager,
    Linear::Builder & builder,
    Topo::Topology & topology,
    IO::InitialConditionsManager & initial_conditions_manager,
    IO::RestartMgr & restart_manager);

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

  // Member variables similar to HB
  AnalysisManager &                     analysisManager_;
  Loader::Loader &                      loader_;
  Linear::System &                      linearSystem_;
  Nonlinear::Manager &                  nonlinearManager_;
  Device::DeviceMgr &                   deviceManager_;
  Linear::Builder &                     builder_;
  Topo::Topology &                      topology_;
  IO::InitialConditionsManager &        initialConditionsManager_;
  IO::RestartMgr &                      restartManager_;
  Parallel::Manager *                   pdsMgrPtr_;

  // HBNOISE specific parameters
  bool outputNodeSingle_;              // Flag for single output node
  std::string outputNode1_;            // First output node
  std::string outputNode2_;            // Second output node
  double harmonicNumber_;              // Harmonic number for hbnoise analysis
  std::string type_;                   // Type of sweep (LIN, DEC, OCT)
  double np_;                          // Number of points
  double fOffsetStart_;                // Start offset frequency for hbnoise analysis
  double fOffsetStop_;                 // Stop offset frequency for hbnoise analysis
  double fAbsStart_;                   // Start offset frequency for hbnoise analysis
  double fAbsStop_;                    // Stop offset frequency for hbnoise analysis
  double stepMult_;                    // Multiplier for frequency steps (DEC/OCT)
  double fstep_;                       // Step size for frequency (LIN)
  int pts_per_summary_;                // Points per summary
  bool dataSpecification_;             // Flag for data specification
  int hbnoiseLoopSize_;                // Size of the hbnoise analysis loop
  SweepVector hbnoiseSweepVector_;       // Vector of sweep parameters
  std::map< std::string, std::vector<std::string> > dataNamesMap_;  // Maps dataset name to parameter names
  std::map< std::string, std::vector< std::vector<double> > > dataTablesMap_;  // Maps dataset name to parameter values
  Analysis::HB *hbAnalysis_;

  // hbnoise integrals are not calculated for DATA=<n> case if the
  // specified frequencies are not monotonically increasing
  bool calcNoiseIntegrals_;

  // Option blocks for parameters
  Util::OptionBlock saved_lsOB_;
  Util::OptionBlock saved_timeIntOB_;

  double freq_;  // Store primary frequency from HB analysis

  int setupSweepParam_();
};

bool registerHBNOISEFactory(FactoryBlock &factory_block);

} // namespace Analysis
} // namespace Xyce

#endif // Xyce_N_ANP_HBNOISE_h 