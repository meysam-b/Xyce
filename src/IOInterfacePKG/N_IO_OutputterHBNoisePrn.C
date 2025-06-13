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
// Purpose        : HBNOISE Prn output
//
// Special Notes  :
//
// Creator        : Meysam Bahmanian
//
// Creation Date  : 6/9/2025
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_IO_OutputterHBNoisePrn.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_Op.h>
#include <N_UTL_DeleteList.h>

namespace Xyce {
namespace IO {
namespace Outputter {

//-----------------------------------------------------------------------------
// Function      : HBNoisePrn::HBNoisePrn
// Purpose       : Constructor
// Special Notes :
// Scope         :
// Creator       : Meysam Bahmanian
// Creation Date : 6/9/2025
//-----------------------------------------------------------------------------
HBNoisePrn::HBNoisePrn(Parallel::Machine comm, OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    outFilename_(),
    os_(0),
    index_(0),
    currentStep_(0),
    numberOfSteps_(0)
{
  if (printParameters_.defaultExtension_.empty())
    printParameters_.defaultExtension_ = "NOISE.prn";

  fixupColumns(comm, outputManager_.getOpBuilderManager(), printParameters_, opList_);
}

//-----------------------------------------------------------------------------
// Function      : HBNoisePrn::~HBNoisePrn
// Purpose       : Destructor
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
HBNoisePrn::~HBNoisePrn()
{
  outputManager_.closeFile(os_);

  deleteList(opList_.begin(), opList_.end());
}

//-----------------------------------------------------------------------------
// Function      : HBNoisePrn::hbNoiseHeader
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Meysam Bahmanian
// Creation Date : 6/9/2025
//-----------------------------------------------------------------------------
void HBNoisePrn::hbNoiseHeader()
{
  if (os_ && currentStep_ == 0)
  {
    int column_index = 0;
    for (Table::ColumnList::const_iterator
        it = printParameters_.table_.columnList_.begin();
        it != printParameters_.table_.columnList_.end();
        ++it, ++column_index)
    {
      if (it != printParameters_.table_.columnList_.begin())
      {
        *os_ << (printParameters_.delimiter_.empty() ? " " : printParameters_.delimiter_);
      }
      printHeader(*os_, (*it));
    }

    for (Table::ColumnList::const_iterator it2 = columnList_.begin(); it2 != columnList_.end(); ++it2)
    {
      if (it2 != columnList_.begin())
      {
        *os_ << printParameters_.delimiter_;
      }
      printHeader(*os_, (*it2));
    }
    *os_ << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : HBNoisePrn::doOutputHBNoise
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Meysam Bahmanian
// Creation Date : 6/9/2025
//-----------------------------------------------------------------------------
void HBNoisePrn::doOutputHBNoise(
  Parallel::Machine   comm,
  double              frequency,
  const Linear::Vector &real_solution_vector, 
  const Linear::Vector &imaginary_solution_vector,
  double              totalAMNoiseDens, 
  double              totalPMNoiseDens, 
  const std::vector<Xyce::Analysis::NoiseData*> & noiseDataVecI,
  const std::vector<Xyce::Analysis::NoiseData*> & noiseDataVecQ,
  // I have not yet found a proper way to do a detailed noise separation for every noise source for every noise sideband
  // Noise separation is not just a feature for users, it is extremely helpful for regression testing and debugging.
  // I think best would be to have .print options to specify noise separation.
  // So these two vectors are not used in the outputter.
  const std::vector<std::vector<Xyce::Analysis::NoiseData*> > & noiseDataVecVecI,
  const std::vector<std::vector<Xyce::Analysis::NoiseData*> > & noiseDataVecVecQ)
{
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
    if (outputManager_.getPrintHeader())
    {
      printHeader(*os_, printParameters_);
    }
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
      printValue(*os_, printParameters_.table_.columnList_[i], printParameters_.delimiter_, i, result_list[i].real());
  }

  if (os_)
    *os_ << std::endl;

  ++index_;
}

//-----------------------------------------------------------------------------
// Function      : HBNoisePrn::doFinishOutput
// Purpose       : Output the footer, and close the stream if there is no
//               : .STEP loop.  This function is also called after each step,
//               : if there is a .STEP loop, but currently does nothing in 
//               : that case.
// Special Notes :
// Scope         :
// Creator       : Meysam Bahmanian
// Creation Date : 6/9/2025
//-----------------------------------------------------------------------------
void HBNoisePrn::doFinishOutput()
{
  if (os_)
  {
    if (numberOfSteps_ == 0)
    {
      if (outputManager_.getPrintFooter ())
      {
        // this end-of-simulation footer is used if there is no .STEP loop
        (*os_) << "End of Xyce(TM) Simulation" << std::endl;
      }
      outputManager_.closeFile(os_);
      os_ = 0;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : HBNoisePrn::doStartStep
// Purpose       : This function is executed at the start of each step.
// Special Notes :
// Scope         :
// Creator       : Meysam Bahmanian
// Creation Date : 6/9/2025
//-----------------------------------------------------------------------------
void HBNoisePrn::doStartStep(
  int                           current_step,
  int                           number_of_steps)
{
  index_ = 0;
  currentStep_ = current_step;
  numberOfSteps_ = number_of_steps;

  // If using Format::GNUPLOT then add two blank lines before the output for
  // steps 1, 2, ... if there is a .STEP loop.  (Note: currentStep_ goes from 
  // 0 to numberOfSteps_-1 if there is a .STEP loop.)
  //
  // If using Format::SPLOT then add a single blank line before the output for
  // steps 1, 2, ... if there is a .STEP loop.  
  if (os_)
  {
    if ( (printParameters_.addGnuplotSpacing_) && (currentStep_ > 0 ) )
    {
      *os_ << std::endl << std::endl;
    }
    else if ( (printParameters_.addSplotSpacing_) && (currentStep_ > 0 ) )
    {
      *os_ << std::endl;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : HBNoisePrn::doResetIndex
// Purpose       : Reset the value for the Index column to zero
// Special Notes :
// Scope         :
// Creator       : Meysam Bahmanian
// Creation Date : 6/9/2025
//-----------------------------------------------------------------------------
void HBNoisePrn::doResetIndex()
{
  index_ = 0;
}

//-----------------------------------------------------------------------------
// Function      : HBNoisePrn::doSteppingComplete
// Purpose       : Output footer and close the stream  when a .STEP loop 
//               : is used.
// Special Notes :
// Scope         :
// Creator       : Meysam Bahmanian
// Creation Date : 6/9/2025
//-----------------------------------------------------------------------------
void HBNoisePrn::doSteppingComplete()
{
  // close the file.
  if (os_)
  {
    // this end-of-simulation footer is used if there is a .STEP loop
    if ( outputManager_.getPrintFooter() )
    {
      (*os_) << "End of Xyce(TM) Parameter Sweep" << std::endl;
    }

    outputManager_.closeFile(os_);
    os_ = 0;
  }
}

} // namespace Outputter
} // namespace IO
} // namespace Xyce
