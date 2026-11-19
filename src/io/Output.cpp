#include "io/Output.h"

#include <cstdio>
#include <fstream>
#include <iomanip>
#include <numbers>
#include <stdexcept>

namespace lambda {

namespace {

constexpr int kDigits = 16;

// Gamma [1/s] and f [Hz] from lambda in viscous units, if physical properties are given.
bool dimensional(const Config& cfg, cplx lam, double& gamma, double& freq) {
    if (!cfg.physical) return false;
    const double tv = cfg.physical->viscous_time();
    gamma = -lam.real() / tv;
    freq = lam.imag() / (2.0 * std::numbers::pi * tv);
    return true;
}

}  // namespace

void write_yaml(std::ostream& os, const Config& cfg, const std::vector<ModeResult>& res) {
    os << std::setprecision(kDigits);
    os << "# LAMBDA " << LAMBDA_VERSION << " output.  lambda = -Gamma + i omega\n"
       << "# viscous units: time R^2/nu;  capillary units: time sqrt(rho R^3/sigma)\n";
    os << "name: " << cfg.name << "\n";
    os << "Oh: " << cfg.Oh << "\n";
    if (cfg.physical) {
        const auto& p = *cfg.physical;
        os << "physical: {density: " << p.density << ", viscosity: " << p.viscosity
           << ", surface_tension: " << p.surface_tension << ", radius: " << p.radius << "}\n";
    }
    os << "Nr: " << cfg.Nr << "\n";
    os << "results:\n";
    for (const auto& r : res) {
        const cplx lc = r.lambda * cfg.Oh;  // capillary units
        os << "  - l: " << r.l << "\n"
           << "    lambda_viscous: [" << r.lambda.real() << ", " << r.lambda.imag() << "]\n"
           << "    lambda_capillary: [" << lc.real() << ", " << lc.imag() << "]\n";
        double g, f;
        if (dimensional(cfg, r.lambda, g, f))
            os << "    Gamma_per_s: " << g << "\n"
               << "    frequency_Hz: " << f << "\n";
        os << "    backward_error: " << r.backward_error << "\n";
        if (r.resolution_change_gamma >= 0) {
            os << "    resolution_change:   # relative change with Nr -> Nr+" << cfg.check_dNr << "\n"
               << "      Gamma: " << r.resolution_change_gamma << "\n";
            if (r.resolution_change_omega >= 0)
                os << "      omega: " << r.resolution_change_omega << "\n";
        }
        os << "    converged: " << (r.converged ? "true" : "false") << "\n";
    }
}

void write_csv(std::ostream& os, const Config& cfg, const std::vector<ModeResult>& res) {
    os << std::setprecision(kDigits);
    os << "l,Oh,Nr,re_viscous,im_viscous,re_capillary,im_capillary,"
          "Gamma_per_s,frequency_Hz,backward_error,dNr_change_Gamma,dNr_change_omega,converged\n";
    for (const auto& r : res) {
        const cplx lc = r.lambda * cfg.Oh;
        os << r.l << ',' << cfg.Oh << ',' << cfg.Nr << ',' << r.lambda.real() << ','
           << r.lambda.imag() << ',' << lc.real() << ',' << lc.imag() << ',';
        double g, f;
        if (dimensional(cfg, r.lambda, g, f)) os << g << ',' << f;
        else os << ',';
        os << ',' << r.backward_error << ',';
        if (r.resolution_change_gamma >= 0) os << r.resolution_change_gamma;
        os << ',';
        if (r.resolution_change_omega >= 0) os << r.resolution_change_omega;
        os << ',' << (r.converged ? 1 : 0) << '\n';
    }
}

void write_spectrum(std::ostream& os, const Config& cfg, int l, const std::vector<cplx>& ev) {
    os << std::setprecision(kDigits);
    os << "# finite eigenvalues (viscous units), Oh = " << cfg.Oh << ", l = " << l
       << ", Nr = " << cfg.Nr << "\n# re im\n";
    for (const auto& z : ev) os << z.real() << ' ' << z.imag() << '\n';
}

void report(const Config& cfg, const std::vector<ModeResult>& res) {
    std::printf("# Oh = %.10e   Nr = %d   (viscous units)\n", cfg.Oh, cfg.Nr);
    std::printf("%-3s %24s %24s %9s", "l", "Re(lambda)", "Im(lambda)", "bwd.err");
    if (cfg.check_resolution) std::printf(" %9s %9s", "dNr:dG/G", "dNr:dw/w");
    std::printf("\n");
    for (const auto& r : res) {
        std::printf("%-3d %24.16e %24.16e %9.1e", r.l, r.lambda.real(), r.lambda.imag(),
                    r.backward_error);
        if (r.resolution_change_gamma >= 0) {
            std::printf(" %9.1e", r.resolution_change_gamma);
            if (r.resolution_change_omega >= 0) std::printf(" %9.1e", r.resolution_change_omega);
            else std::printf(" %9s", "-");
        }
        std::printf("%s\n", r.converged ? "" : "  NOT CONVERGED");
    }
    if (cfg.output_file.empty()) return;
    std::ofstream f(cfg.output_file);
    if (!f) throw std::runtime_error("cannot open " + cfg.output_file);
    if (cfg.output_format == OutputFormat::Csv) write_csv(f, cfg, res);
    else write_yaml(f, cfg, res);
}

}  // namespace lambda
