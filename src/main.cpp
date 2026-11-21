// LAMBDA: linear eigenvalues of a viscous drop oscillating in vacuum.
//
//   lambdaDrop input.yaml
//   lambdaDrop We l [m] [Nr] [guessRe guessIm]     (We = Oh^2; m has no effect)
//   lambdaDrop --all We l [Nr]                           (all eigenvalues of degree l)
#include <cstdio>
#include <exception>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "io/Config.h"
#include "io/FieldWriter.h"
#include "io/Output.h"
#include "post/Eigenfunction.h"
#include "solver/Solver.h"

namespace {

void usage() {
    std::fputs(
        "LAMBDA " LAMBDA_VERSION
        " - eigenvalues of a free viscous drop\n"
        "usage:\n"
        "  lambdaDrop input.yaml\n"
        "  lambdaDrop We l [m] [Nr] [guessRe guessIm]   (We = Oh^2; m has no effect)\n"
        "  lambdaDrop --all We l [Nr]                         (all finite eigenvalues of degree l)\n"
        "  lambdaDrop --version\n",
        stderr);
}

bool ends_with(const std::string& s, const std::string& suf) {
    return s.size() >= suf.size() && s.compare(s.size() - suf.size(), suf.size(), suf) == 0;
}

}  // namespace

int main(int argc, char** argv) {
    std::vector<std::string> args(argv + 1, argv + argc);
    if (args.empty() || args[0] == "-h" || args[0] == "--help") {
        usage();
        return args.empty() ? 1 : 0;
    }
    if (args[0] == "--version") {
        std::printf("LAMBDA %s\n", LAMBDA_VERSION);
        return 0;
    }
    try {
        const bool yaml =
            args.size() == 1 && (ends_with(args[0], ".yaml") || ends_with(args[0], ".yml"));
        const lambda::Config cfg =
            yaml ? lambda::load_yaml(args[0]) : lambda::parse_command_line(args);

        if (cfg.method == lambda::Method::Spectrum) {
            std::ofstream f;
            if (!cfg.output_file.empty()) f.open(cfg.output_file);
            std::ostream& os = cfg.output_file.empty() ? std::cout : f;
            for (const auto& md : cfg.modes)
                lambda::write_spectrum(os, cfg, md.l, lambda::spectrum(md.l, cfg.Oh, cfg.Nr));
            return 0;
        }

        const auto res = lambda::solve(cfg);
        lambda::report(cfg, res);
        if (cfg.eigenfunction.enabled)
            for (const auto& r : res)
                for (int m : cfg.eigenfunction.m) {
                    const lambda::Eigenfunction ef(r, m);
                    for (const auto& fn : lambda::write_eigenfunction(ef, cfg.eigenfunction))
                        std::printf("wrote %s\n", fn.c_str());
                }
        for (const auto& r : res)
            if (!r.converged) return 2;
        return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "lambdaDrop: %s\n", e.what());
        return 1;
    }
}
