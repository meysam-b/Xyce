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

#include <N_ANP_AnalysisManager.h>
#include <N_ANP_DCSweep.h>
#include <N_ANP_HBNOISE.h>
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
#include <N_LAS_BlockSystemHelpers.h>
#include <N_LAS_BlockVector.h>
#include <N_LAS_HBBuilder.h>
#include <N_LAS_HBPrecondFactory.h>
#include <N_LAS_HBSolverFactory.h>
#include <N_LAS_PrecondFactory.h>
#include <N_LAS_System.h>
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

using Teuchos::rcp;
using Teuchos::RCP;
using Teuchos::rcp_dynamic_cast;

namespace Xyce {
namespace Analysis {

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
  Linear::Builder &                     builder,
  Topo::Topology &                      topology,
  IO::InitialConditionsManager &        initial_conditions_manager,
  IO::RestartMgr &                      restart_manager)
  : AnalysisBase(analysis_manager, "HBNOISE"),
    StepEventListener(&analysis_manager),
    analysisManager_(analysis_manager),
    loader_(loader),
    linearSystem_(linear_system),
    nonlinearManager_(nonlinear_manager),
    deviceManager_(device_manager),
    builder_(builder),
    topology_(topology),
    initialConditionsManager_(initial_conditions_manager),
    restartManager_(restart_manager),
    pdsMgrPtr_(0),
    outputNodeSingle_(true),
    outputNode1_(""),
    outputNode2_(""),
    harmonicNumber_(1.0),
    type_("DEC"),
    np_(10.0),
    fOffsetStart_(1.0),
    fOffsetStop_(1.0),
    pts_per_summary_(0),
    dataSpecification_(false)
{
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
  return true;
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
  return true;
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
    noiseSweepVector_.push_back(parseSweepParams(paramsBlock.begin(), paramsBlock.end()));
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
public:
  HBNOISEFactory(
    Analysis::AnalysisManager &          analysis_manager,
    Linear::System &                     linear_system,
    Nonlinear::Manager &                 nonlinear_manager,
    Loader::Loader &                     loader,
    Device::DeviceMgr &                  device_manager,
    Linear::Builder &                    builder,
    Topo::Topology &                     topology,
    IO::InitialConditionsManager &       initial_conditions_manager,
    IO::RestartMgr &                     restart_manager)
    : HBNOISEFactoryBase(),
      analysisManager_(analysis_manager),
      linearSystem_(linear_system),
      nonlinearManager_(nonlinear_manager),
      loader_(loader),
      deviceManager_(device_manager),
      builder_(builder),
      topology_(topology),
      initialConditionsManager_(initial_conditions_manager),
      restartManager_(restart_manager)
  {}

  virtual ~HBNOISEFactory()
  {}

  HBNOISE *create() const
  {
    analysisManager_.setAnalysisMode(ANP_MODE_HBNOISE);

    HBNOISE *hbnoise = new HBNOISE(analysisManager_, linearSystem_,
                                  nonlinearManager_, loader_, deviceManager_,
                                  builder_, topology_, initialConditionsManager_,
                                  restartManager_);

    hbnoise->setAnalysisParams(hbnoiseAnalysisOptionBlock_);
    hbnoise->setLinSol(linSolOptionBlock_);

    return hbnoise;
  }

  void setHBNOISEAnalysisOptionBlock(const Util::OptionBlock &option_block)
  {
    hbnoiseAnalysisOptionBlock_ = option_block;
  }

  // Move these to public section
  AnalysisManager &                     analysisManager_;
  Linear::System &                      linearSystem_;
  Nonlinear::Manager &                  nonlinearManager_;
  Loader::Loader &                      loader_;
  Device::DeviceMgr &                   deviceManager_;
  Linear::Builder &                     builder_;
  Topo::Topology &                      topology_;
  IO::InitialConditionsManager &        initialConditionsManager_;
  IO::RestartMgr &                      restartManager_;

private:
  Util::OptionBlock     hbnoiseAnalysisOptionBlock_;
  Util::OptionBlock     linSolOptionBlock_;
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
    factory_block.builder_, factory_block.topology_,
    factory_block.initialConditionsManager_, factory_block.restartManager_);

  addAnalysisFactory(factory_block, factory);

  // Register the command parser and processor
  factory_block.optionsManager_.addCommandParser(".HBNOISE", extractHBNOISEData);
  factory_block.optionsManager_.addCommandProcessor("HBNOISE", new HBNOISEAnalysisReg(*factory));

  return true;
}

} // namespace Analysis
} // namespace Xyce 