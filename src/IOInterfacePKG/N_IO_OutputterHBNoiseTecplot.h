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
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_OutputterHBNoiseTecplot_h
#define Xyce_N_IO_OutputterHBNoiseTecplot_h

#include <N_IO_OutputterLocal.h>

namespace Xyce {
namespace IO {
namespace Outputter {

class HBNoiseTecPlot : public Interface
{
public:
  HBNoiseTecPlot(Parallel::Machine comm, OutputMgr &output_manager, const PrintParameters &print_parameters);

  virtual ~HBNoiseTecPlot();

private:
  HBNoiseTecPlot(const HBNoiseTecPlot &);
  HBNoiseTecPlot &operator=(const HBNoiseTecPlot &);

public:

  virtual void doSetAnalysisMode(Analysis::Mode analysis_mode)
  {
    printParameters_.analysisMode_ = analysis_mode;
  }

  virtual void doFinishOutput();

  virtual void doStartStep(int current_step, int number_of_step);

  virtual void doResetIndex();

  virtual void doSteppingComplete();

  virtual void doOutputHBNoise(
    Parallel::Machine   comm,
    double              frequency,
    const Linear::Vector &real_solution_vector, 
    const Linear::Vector &imaginary_solution_vector,
    double              totalAMNoiseDens, 
    double              totalPMNoiseDens, 
    const std::vector<Xyce::Analysis::NoiseData*> & noiseDataVecI,
    const std::vector<Xyce::Analysis::NoiseData*> & noiseDataVecQ,
    const std::vector<std::vector<Xyce::Analysis::NoiseData*> > & noiseDataVecVecI,
    const std::vector<std::vector<Xyce::Analysis::NoiseData*> > & noiseDataVecVecQ);

private:
  OutputMgr &           outputManager_;
  PrintParameters       printParameters_;
  std::string           outFilename_;
  std::ostream *        os_;
  int                   index_;
  int                   currentStep_;
  int                   numberOfSteps_;

  Util::Op::OpList      opList_;
};


} // namespace Outputter
} // namespace IO
} // namespace Xyce

#endif // Xyce_N_IO_OutputterHBNoisePrn_h
