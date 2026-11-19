#pragma once

#include <ostream>
#include <vector>

#include "io/Config.h"
#include "solver/Solver.h"

namespace lambda {

void write_yaml(std::ostream& os, const Config& cfg, const std::vector<ModeResult>& res);
void write_csv(std::ostream& os, const Config& cfg, const std::vector<ModeResult>& res);
void write_spectrum(std::ostream& os, const Config& cfg, int l, const std::vector<cplx>& ev);

// Short table to stdout, and the full result to cfg.output_file if set.
void report(const Config& cfg, const std::vector<ModeResult>& res);

}  // namespace lambda
