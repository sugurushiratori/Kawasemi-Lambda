// Unit tests for the input layer (no solver needed).
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>

#include "io/Config.h"

using namespace lambda;

static int failures = 0;
#define CHECK(cond)                                                                   \
    do {                                                                              \
        if (!(cond)) {                                                                \
            std::fprintf(stderr, "%s:%d: CHECK failed: %s\n", __FILE__, __LINE__, #cond); \
            ++failures;                                                               \
        }                                                                             \
    } while (0)

template <class F>
static bool throws(F f) {
    try {
        f();
    } catch (const std::runtime_error&) {
        return true;
    }
    return false;
}

static bool close(double a, double b, double rtol = 1e-14) {
    return std::abs(a - b) <= rtol * std::abs(b);
}

int main() {
    {  // physical properties -> Oh (Ni at 1728 K, R = 1 mm)
        auto c = parse_yaml_string(R"(
name: Ni
drop:
  physical: {density: 7890, viscosity: 4.7e-3, surface_tension: 1.739, radius: 1.0e-3}
modes: [2, 3, 4]
resolution: {Nr: 40, check: true}
output: {file: out.csv}
)");
        CHECK(c.name == "Ni");
        CHECK(c.physical.has_value());
        CHECK(close(c.Oh, 4.7e-3 / std::sqrt(7890 * 1.739 * 1e-3)));
        CHECK(c.modes.size() == 3 && c.modes[2].l == 4);
        CHECK(c.Nr == 40 && c.check_resolution && c.check_dNr == 16);
        CHECK(c.output_format == OutputFormat::Csv);
        CHECK(close(c.physical->capillary_time() / c.physical->viscous_time(), c.Oh, 1e-13));
    }
    {  // We, list of maps with a guess, defaults
        auto c = parse_yaml_string(R"(
drop: {We: 1.0e-3}
modes:
  - 2
  - {l: 4, guess: [-22.9, 264.4]}
solver: {method: mode, tolerance: 1.0e-8}
)");
        CHECK(close(c.Oh, std::sqrt(1e-3)));
        CHECK(c.modes.size() == 2 && !c.modes[0].guess);
        CHECK(c.modes[1].guess && c.modes[1].guess->imag() == 264.4);
        CHECK(c.tolerance == 1e-8);
        CHECK(c.Nr == 48);
        CHECK(c.output_file.empty());
    }
    {  // eigenfunction output
        auto c = parse_yaml_string(R"(
name: X
drop: {Oh: 0.1}
modes: [2, 3]
output:
  eigenfunction: {m: [0, 2], format: vtk, points: {r: 11}}
)");
        CHECK(c.eigenfunction.enabled && c.eigenfunction.prefix == "X");
        CHECK(c.eigenfunction.m.size() == 2 && c.eigenfunction.m[1] == 2);
        CHECK(c.eigenfunction.vtk && !c.eigenfunction.ascii);
        CHECK(c.eigenfunction.n_r == 11 && c.eigenfunction.n_theta == 49);
    }
    CHECK(throws([] { parse_yaml_string("drop: {Oh: 0.1}\nmodes: [2]\noutput: {eigenfunction: {m: 3}}"); }));
    CHECK(throws([] { parse_yaml_string("drop: {Oh: 0.1}\nmodes: [2]\noutput: {eigenfunction: {format: png}}"); }));

    // errors
    CHECK(throws([] { parse_yaml_string("drop: {Oh: 0.1, We: 0.01}\nmodes: [2]"); }));
    CHECK(throws([] { parse_yaml_string("drop: {Oh: 0.1}\nmodes: [1]"); }));
    CHECK(throws([] { parse_yaml_string("drop: {Oh: -1}\nmodes: [2]"); }));
    CHECK(throws([] { parse_yaml_string("drop: {Oh: 0.1}\nmodes: 2"); }));
    CHECK(throws([] { parse_yaml_string("drop: {Oh: 0.1}\nmodes: [2]\nsolver: {method: x}"); }));
    CHECK(throws([] { parse_yaml_string("drop: {Oh: 0.1}\nmodes: [2]\nresolution: {Nq: 8}"); }));
    CHECK(throws([] { parse_yaml_string("drop: {Oh: [0.1}"); }));

    {  // command line compatible with dropEigen
        auto c = parse_legacy_args({"1e-3", "3", "1", "40", "36", "-12.2", "171.4"});
        CHECK(close(c.Oh, std::sqrt(1e-3)));
        CHECK(c.modes.size() == 1 && c.modes[0].l == 3);
        CHECK(c.Nr == 40);
        CHECK(c.modes[0].guess && c.modes[0].guess->real() == -12.2);
        auto a = parse_legacy_args({"--all", "1e-1", "2", "24"});
        CHECK(a.method == Method::Spectrum && a.Nr == 24);
        CHECK(throws([] { parse_legacy_args({"1e-3"}); }));
        CHECK(throws([] { parse_legacy_args({"1e-3", "2", "0", "32", "32", "1.0"}); }));
        CHECK(throws([] { parse_legacy_args({"1e-3x", "2"}); }));
    }

    if (failures) std::fprintf(stderr, "%d check(s) failed\n", failures);
    else std::puts("test_config: all checks passed");
    return failures ? 1 : 0;
}
