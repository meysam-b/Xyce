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
// Purpose        : HBNOISE output
//
// Special Notes  :
//
// Creator        : Meysam Bahmanian
//
// Creation Date  : 6/10/2025
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_OutputterHBNoise_h
#define Xyce_N_IO_OutputterHBNoise_h

#include <N_IO_OutputterLocal.h>

namespace Xyce {
namespace IO {
namespace Outputter {

void enableHBNoiseOutput(Parallel::Machine comm, OutputMgr &output_manager, Analysis::Mode analysis_mode);

std::ostream &printHBNoiseHeader(std::ostream &os, const PrintParameters &print_parameters);



} // namespace Outputter
} // namespace IO
} // namespace Xyce

#endif // Xyce_N_IO_OutputterHBNoise_h
