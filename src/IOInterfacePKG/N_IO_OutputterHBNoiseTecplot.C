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

//-------------------------------------------------------------------------
//
// Purpose        : HBNOISE Tecplot output
//
// Special Notes  : 
//
// Creator        : Meysam Bahmanian
//
// Creation Date  : 6/14/2025
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_IO_OutputterHBNoiseTecplot.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_Tecplot.h>
#include <N_IO_Op.h>
#include <N_UTL_DeleteList.h>

namespace Xyce {
namespace IO {
namespace Outputter {

//-----------------------------------------------------------------------------
// Function      : NoiseTecPlot::NoiseTecPlot
// Purpose       : Constructor
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
HBNoiseTecPlot::HBNoiseTecPlot(Parallel::Machine comm, OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    outFilename_(),
    os_(0),
    index_(0),
    currentStep_(0),
    numberOfSteps_(0)
{
  if (printParameters_.defaultExtension_.empty())
    printParameters_.defaultExtension_ = "NOISE.dat";

  fixupColumns(comm, outputManager_.getOpBuilderManager(), printParameters_, opList_);
}

//-----------------------------------------------------------------------------
// Function      : HBNoiseTecPlot::~HBNoiseTecPlot
// Purpose       : Destructor
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
HBNoiseTecPlot::~HBNoiseTecPlot()
{
  outputManager_.closeFile(os_);

  deleteList(opList_.begin(), opList_.end());
}

//-----------------------------------------------------------------------------
// Function      : HBNoiseTecplot::doOutputHBNoise
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Meysam Bahmanian
// Creation Date : 6/14/2025
//-----------------------------------------------------------------------------
void HBNoiseTecPlot::doOutputHBNoise(
  Parallel::Machine   comm,
  double              frequency,
  const Linear::Vector &real_solution_vector, 
  const Linear::Vector &imaginary_solution_vector,
  double              totalAMNoiseDens, 
  double              totalPMNoiseDens, 
  const std::vector<Xyce::Analysis::NoiseData*> & noiseDataVecI,
  const std::vector<Xyce::Analysis::NoiseData*> & noiseDataVecQ,
  const std::vector<std::vector<Xyce::Analysis::NoiseData*> > & noiseDataVecVecI,
  const std::vector<std::vector<Xyce::Analysis::NoiseData*> > & noiseDataVecVecQ)
{
  double tmpTime = outputManager_.getCircuitTime();

  if (Parallel::rank(comm) == 0 && !os_)
  {
    outFilename_ = outputFilename(printParameters_.filename_, 
                                  printParameters_.defaultExtension_,
                                  printParameters_.suffix_+outputManager_.getFilenameSuffix(), 
                                  outputManager_.getNetlistFilename(),
                                  printParameters_.overrideRawFilename_,
                                  printParameters_.formatSupportsOverrideRaw_,
                                  printParameters_.dashoFilename_,
                                  printParameters_.fallback_);
    os_ = outputManager_.openFile(outFilename_);
    os_->setf(std::ios::scientific);
    os_->precision(printParameters_.streamPrecision_);
    os_->setf(std::ios::left, std::ios::adjustfield);
  }

  if (os_ && index_ == 0)
  {
    tecplotFreqHeader(*os_, currentStep_ == 0, outputManager_.getNetlistFilename(), 
        opList_, outputManager_);
  }

  std::vector<complex> result_list;
  Util::Op::OpData op_data;
  op_data.amnoise_ = totalAMNoiseDens;
  op_data.pmnoise_ = totalPMNoiseDens;
  op_data.amnoiseDataVec_ = &noiseDataVecI;
  op_data.pmnoiseDataVec_ = &noiseDataVecQ;
  getValues(comm, opList_, op_data, result_list);

  for (int i = 0; i < result_list.size(); ++i)
  {
    result_list[i] = complex(filter(result_list[i].real(), printParameters_.filter_), 0.0);

    if (os_)
    {
      printValue(*os_, printParameters_.table_.columnList_[i], printParameters_.delimiter_, i, result_list[i].real());
    }
  }

  if (os_)
  {
    *os_ << std::endl;
  }

  ++index_;
}

//-----------------------------------------------------------------------------
// Function      : HBNoiseTecPlot::doFinishOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Meysam Bahmanian
// Creation Date : 6/14/2025
//-----------------------------------------------------------------------------
void HBNoiseTecPlot::doFinishOutput()
{
  if (os_)
  {
    if (numberOfSteps_ == 0)
    {
      if (outputManager_.getPrintFooter ())
      {
        (*os_) << "End of Xyce(TM) Simulation" << std::endl;
      }

      outputManager_.closeFile(os_);
      os_ = 0;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : HBNoiseTecPlot::doStartStep
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Meysam Bahmanian
// Creation Date : 6/14/2025
//-----------------------------------------------------------------------------
void HBNoiseTecPlot::doStartStep(
  int                           current_step,
  int                           number_of_steps)
{
  index_ = 0;
  currentStep_ = current_step;
  numberOfSteps_ = number_of_steps;
}

//-----------------------------------------------------------------------------
// Function      : HBNoiseTecPlot::doStartStep
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Meysam Bahmanian
// Creation Date : 6/14/2025
//-----------------------------------------------------------------------------
void HBNoiseTecPlot::doResetIndex()
{
  index_ = 0;
}

//-----------------------------------------------------------------------------
// Function      : HBNoiseTecPlot::doSteppingComplete
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Meysam Bahmanian
// Creation Date : 6/14/2025
//-----------------------------------------------------------------------------
void HBNoiseTecPlot::doSteppingComplete()
{
  // close the noise file.
  if (os_)
  {
    if ( outputManager_.getPrintFooter () )
    {
      (*os_) << "End of Xyce(TM) Parameter Sweep" << std::endl;
    }
  }

  outputManager_.closeFile(os_);
  os_ = 0;
}

} // namespace Outputter
} // namespace IO
} // namespace Xyce
