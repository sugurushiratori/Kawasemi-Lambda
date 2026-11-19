// Eigenfunction output: ASCII tables (gnuplot-friendly) and ParaView XML
// StructuredGrid files (.vts).
#pragma once

#include <string>
#include <vector>

#include "io/Config.h"
#include "post/Eigenfunction.h"

namespace lambda {

// Writes, for one mode and one m, the files
//   <prefix>_l<l>_radial.dat              u(r), v(r), p(r)            (ASCII, m-independent)
//   <prefix>_l<l>_m<m>_meridional.dat     fields in the plane phi = 0 (ASCII)
//   <prefix>_l<l>_m<m>_surface.dat        d Y and velocity on r = 1   (ASCII)
//   <prefix>_l<l>_m<m>.vts                fields in the drop          (ParaView)
//   <prefix>_l<l>_m<m>_surface.vts        surface, d Y e_r for Warp By Vector (ParaView)
// Returns the names of the files written.
std::vector<std::string> write_eigenfunction(const Eigenfunction& ef,
                                             const EigenfunctionOutput& opt);

}  // namespace lambda
