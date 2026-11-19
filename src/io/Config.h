// Input description for LAMBDA: one drop, several degrees l.
#pragma once

#include <complex>
#include <optional>
#include <string>
#include <vector>

namespace lambda {

// Physical properties in SI units.
struct Physical {
    double density = 0.0;          // rho   [kg/m^3]
    double viscosity = 0.0;        // mu    [Pa s]
    double surface_tension = 0.0;  // sigma [N/m]
    double radius = 0.0;           // R     [m]

    double ohnesorge() const;       // mu / sqrt(rho sigma R)
    double viscous_time() const;    // rho R^2 / mu          [s]
    double capillary_time() const;  // sqrt(rho R^3 / sigma) [s]
};

struct Mode {
    int l = 2;
    // Initial guess in viscous units (time R^2/nu). Empty: automatic.
    std::optional<std::complex<double>> guess;
};

enum class Method { Mode, Spectrum };

// Spatial form of the eigenmodes (see io/FieldWriter.h).
struct EigenfunctionOutput {
    bool enabled = false;
    std::string prefix;       // file name prefix; default: the case name
    std::vector<int> m{0};    // azimuthal orders, 0 <= m <= l
    bool ascii = true;        // .dat tables
    bool vtk = true;          // ParaView .vts
    int n_r = 33, n_theta = 49, n_phi = 65;  // output points (phi includes 0 and 2 pi)
};
enum class OutputFormat { Yaml, Csv };

struct Config {
    std::string name = "case";

    std::optional<Physical> physical;  // if given, Oh is computed from it
    double Oh = 0.0;                   // Ohnesorge number

    std::vector<Mode> modes;

    int Nr = 48;  // radial Chebyshev functions per variable

    Method method = Method::Mode;
    double tolerance = 1e-10;       // accept when relative backward error < tolerance
    bool check_resolution = false;  // also solve with Nr + check_dNr and report the change
    int check_dNr = 16;

    std::string output_file;  // empty: stdout only
    OutputFormat output_format = OutputFormat::Yaml;
    EigenfunctionOutput eigenfunction;

    double We() const { return Oh * Oh; }  // legacy parameter name (= Oh^2)
    void validate() const;                 // throws std::runtime_error
};

// Read a YAML input file.
Config load_yaml(const std::string& path);
Config parse_yaml_string(const std::string& text);

// Command line (compatible with dropEigen):
//   lambdaDrop We l [m] [Nr] [Nq] [guessRe guessIm]   m and Nq are accepted and ignored
//   lambdaDrop --all We l [Nr]                          all eigenvalues of degree l
Config parse_legacy_args(const std::vector<std::string>& args);

}  // namespace lambda
