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
// Purpose        : Generate HBNoise output in CSV format
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

#include <N_IO_OutputterHBNoise.h>
#include <N_IO_OutputterHBNoiseCSV.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_Op.h>
#include <N_UTL_DeleteList.h>

namespace Xyce {
namespace IO {
namespace Outputter {

//-----------------------------------------------------------------------------
// Class         : HBNoiseCSV
// Purpose       : Outputter class for HBNoise output, (CSV) output
//                 format
// Special Notes :
// Creator       : Meysam Bahmanian
// Creation Date : 6/14/2025
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Function      : HBNoiseCSV::HBNoiseCSV
// Purpose       : constructor
// Special Notes :
// Scope         :
// Creator       : Meysam Bahmanian
// Creation Date : 6/14/2025
//-----------------------------------------------------------------------------
HBNoiseCSV::HBNoiseCSV(Parallel::Machine comm, OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    outFilename_(),
    os_(0),
    index_(0),
    currentStep_(0),
    numberOfSteps_(0)
{
  if (printParameters_.defaultExtension_.empty())
    printParameters_.defaultExtension_ = ".NOISE.csv";

  fixupColumns(comm, outputManager_.getOpBuilderManager(), printParameters_, opList_);
}

//-----------------------------------------------------------------------------
// Function      : HBNoiseCSV::~HBNoiseCSV
// Purpose       : destructor
// Special Notes :
// Scope         :
// Creator       : Meysam Bahmanian
// Creation Date : 6/14/2025
//-----------------------------------------------------------------------------
HBNoiseCSV::~HBNoiseCSV()
{
  outputManager_.closeFile(os_);

  deleteList(opList_.begin(), opList_.end());
}

//-----------------------------------------------------------------------------
// Function      : HBNoiseCSV::hbNoiseHeader
// Purpose       : Print out the header line in the .csv file
// Special Notes :
// Scope         :
// Creator       : Meysam Bahmanian
// Creation Date : 6/14/2025
//-----------------------------------------------------------------------------
void HBNoiseCSV::hbNoiseHeader()
{
  if (os_ && currentStep_ == 0)
  {
    printHBNoiseHeader(*os_, printParameters_);
  }
}

//-----------------------------------------------------------------------------
// Function      : HBNoiseCSV::doOutputHBNoise
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Meysam Bahmanian
// Creation Date : 6/14/2025
//-----------------------------------------------------------------------------
void HBNoiseCSV::doOutputHBNoise(
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

    printHBNoiseHeader(*os_, printParameters_);
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
// Function      : HBNoiseCSV::doFinishOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Meysam Bahmanian
// Creation Date : 6/14/2025
//-----------------------------------------------------------------------------
void HBNoiseCSV::doFinishOutput()
{
  if (os_)
  {
    if (numberOfSteps_ == 0)
    {      
      outputManager_.closeFile(os_);
      os_ = 0;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : HBNoiseCSV::doStartStep
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Meysam Bahmanian
// Creation Date : 6/14/2025
//-----------------------------------------------------------------------------
void HBNoiseCSV::doStartStep(
  int                           current_step,
  int                           number_of_steps)
{
  index_ = 0;
  currentStep_ = current_step;
  numberOfSteps_ = number_of_steps;
}

//-----------------------------------------------------------------------------
// Function      : HBNoiseCSV::doStartStep
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Meysam Bahmanian
// Creation Date : 6/14/2025
//-----------------------------------------------------------------------------
void HBNoiseCSV::doResetIndex()
{
  index_ = 0;
}

//-----------------------------------------------------------------------------
// Function      : HBNoiseCSV::doSteppingComplete
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Meysam Bahmanian
// Creation Date : 6/14/2025
//-----------------------------------------------------------------------------
void HBNoiseCSV::doSteppingComplete()
{
  // close the file.
  if (os_)
  {
    outputManager_.closeFile(os_);
    os_ = 0;
  }
}

} // namespace Outputter
} // namespace IO
} // namespace Xyce
