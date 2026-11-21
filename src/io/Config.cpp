#include "io/Config.h"

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <stdexcept>

#include <yaml-cpp/yaml.h>

namespace lambda {

double Physical::ohnesorge() const {
    return viscosity / std::sqrt(density * surface_tension * radius);
}
double Physical::viscous_time() const { return density * radius * radius / viscosity; }
double Physical::capillary_time() const {
    return std::sqrt(density * radius * radius * radius / surface_tension);
}

void Config::validate() const {
    auto fail = [](const std::string& msg) { throw std::runtime_error("input: " + msg); };
    if (!(Oh > 0.0) || !std::isfinite(Oh)) fail("Oh must be positive");
    if (physical) {
        const auto& p = *physical;
        if (!(p.density > 0 && p.viscosity > 0 && p.surface_tension > 0 && p.radius > 0))
            fail("physical properties must all be positive");
    }
    if (Nr < 8) fail("Nr must be >= 8");
    if (Nr > 400) fail("Nr must be <= 400");
    if (modes.empty()) fail("no degree l given");
    for (const auto& md : modes)
        if (md.l < 2) fail("l must be >= 2 (l = 0, 1 are not shape oscillations)");
    if (!(tolerance > 0)) fail("tolerance must be positive");
    if (eigenfunction.enabled) {
        const auto& e = eigenfunction;
        if (method != Method::Mode) fail("eigenfunction output needs solver.method: mode");
        if (!e.ascii && !e.vtk) fail("eigenfunction.format is empty");
        if (e.n_r < 2 || e.n_theta < 3 || e.n_phi < 3) fail("eigenfunction.points too few");
        for (int m : e.m) {
            if (m < 0) fail("eigenfunction.m must be >= 0 (m and -m give mirror images)");
            for (const auto& md : modes)
                if (m > md.l)
                    fail("eigenfunction.m = " + std::to_string(m) + " exceeds l = " +
                         std::to_string(md.l));
        }
    }
    if (check_dNr < 1) fail("check_dNr must be >= 1");
}

namespace {

void check_keys(const YAML::Node& n, const std::string& where,
                std::initializer_list<const char*> allowed) {
    if (!n.IsMap()) throw std::runtime_error("input: '" + where + "' must be a map");
    for (const auto& kv : n) {
        const auto key = kv.first.as<std::string>();
        if (std::none_of(allowed.begin(), allowed.end(), [&](const char* a) { return key == a; }))
            throw std::runtime_error("input: unknown key '" + key + "' in '" + where + "'");
    }
}

template <class T>
T get(const YAML::Node& n, const char* key, T def) {
    return n[key] ? n[key].as<T>() : def;
}

Config from_node(const YAML::Node& root) {
    Config c;
    check_keys(root, "top level", {"name", "drop", "modes", "resolution", "solver", "output"});
    c.name = get<std::string>(root, "name", c.name);

    // ---- drop: exactly one of physical, Oh, We ----
    const YAML::Node d = root["drop"];
    if (!d) throw std::runtime_error("input: missing 'drop'");
    check_keys(d, "drop", {"physical", "Oh", "We"});
    const int given = (d["physical"] ? 1 : 0) + (d["Oh"] ? 1 : 0) + (d["We"] ? 1 : 0);
    if (given != 1)
        throw std::runtime_error("input: give exactly one of drop.physical, drop.Oh, drop.We");
    if (const YAML::Node p = d["physical"]) {
        check_keys(p, "drop.physical", {"density", "viscosity", "surface_tension", "radius"});
        Physical ph;
        ph.density = p["density"].as<double>();
        ph.viscosity = p["viscosity"].as<double>();
        ph.surface_tension = p["surface_tension"].as<double>();
        ph.radius = p["radius"].as<double>();
        c.physical = ph;
        c.Oh = ph.ohnesorge();
    } else if (d["Oh"]) {
        c.Oh = d["Oh"].as<double>();
    } else {
        c.Oh = std::sqrt(d["We"].as<double>());
    }

    // ---- modes: [2, 3, 4]  or  [{l: 2, guess: [re, im]}, ...] ----
    const YAML::Node ms = root["modes"];
    if (!ms || !ms.IsSequence()) throw std::runtime_error("input: 'modes' must be a list");
    for (const auto& x : ms) {
        Mode md;
        if (x.IsScalar()) {
            md.l = x.as<int>();
        } else {
            check_keys(x, "modes[]", {"l", "guess"});
            md.l = x["l"].as<int>();
            if (x["guess"]) {
                auto g = x["guess"].as<std::vector<double>>();
                if (g.size() != 2) throw std::runtime_error("input: guess must be [re, im]");
                md.guess = std::complex<double>(g[0], g[1]);
            }
        }
        c.modes.push_back(md);
    }

    if (const YAML::Node r = root["resolution"]) {
        check_keys(r, "resolution", {"Nr", "check", "check_dNr"});
        c.Nr = get<int>(r, "Nr", c.Nr);
        c.check_resolution = get<bool>(r, "check", c.check_resolution);
        c.check_dNr = get<int>(r, "check_dNr", c.check_dNr);
    }

    if (const YAML::Node s = root["solver"]) {
        check_keys(s, "solver", {"method", "tolerance"});
        const auto method = get<std::string>(s, "method", "mode");
        if (method == "mode")
            c.method = Method::Mode;
        else if (method == "spectrum")
            c.method = Method::Spectrum;
        else
            throw std::runtime_error("input: solver.method must be 'mode' or 'spectrum'");
        c.tolerance = get<double>(s, "tolerance", c.tolerance);
    }

    if (const YAML::Node o = root["output"]) {
        check_keys(o, "output", {"file", "format", "eigenfunction"});
        c.output_file = get<std::string>(o, "file", "");
        auto fmt = get<std::string>(o, "format", "");
        if (fmt.empty()) {
            const auto& f = c.output_file;
            fmt = (f.size() >= 4 && f.compare(f.size() - 4, 4, ".csv") == 0) ? "csv" : "yaml";
        }
        if (fmt == "yaml")
            c.output_format = OutputFormat::Yaml;
        else if (fmt == "csv")
            c.output_format = OutputFormat::Csv;
        else
            throw std::runtime_error("input: output.format must be yaml or csv");

        if (const YAML::Node e = o["eigenfunction"]) {
            check_keys(e, "output.eigenfunction", {"prefix", "m", "format", "points"});
            auto& ef = c.eigenfunction;
            ef.enabled = true;
            ef.prefix = get<std::string>(e, "prefix", c.name);
            if (e["m"]) {
                ef.m.clear();
                if (e["m"].IsSequence())
                    for (const auto& x : e["m"]) ef.m.push_back(x.as<int>());
                else
                    ef.m.push_back(e["m"].as<int>());
            }
            if (const YAML::Node f = e["format"]) {
                ef.ascii = ef.vtk = false;
                std::vector<std::string> fs;
                if (f.IsSequence()) fs = f.as<std::vector<std::string>>();
                else fs.push_back(f.as<std::string>());
                for (const auto& x : fs) {
                    if (x == "ascii") ef.ascii = true;
                    else if (x == "vtk" || x == "paraview") ef.vtk = true;
                    else throw std::runtime_error("input: eigenfunction.format: ascii and/or vtk");
                }
            }
            if (const YAML::Node pt = e["points"]) {
                check_keys(pt, "output.eigenfunction.points", {"r", "theta", "phi"});
                ef.n_r = get<int>(pt, "r", ef.n_r);
                ef.n_theta = get<int>(pt, "theta", ef.n_theta);
                ef.n_phi = get<int>(pt, "phi", ef.n_phi);
            }
        }
    }

    c.validate();
    return c;
}

double to_number(const std::string& s) {
    std::size_t pos = 0;
    const double v = std::stod(s, &pos);
    if (pos != s.size()) throw std::runtime_error("not a number: " + s);
    return v;
}

int to_int(const std::string& s) {
    std::size_t pos = 0;
    const int v = std::stoi(s, &pos);
    if (pos != s.size()) throw std::runtime_error("not an integer: " + s);
    return v;
}

}  // namespace

Config parse_yaml_string(const std::string& text) {
    try {
        return from_node(YAML::Load(text));
    } catch (const YAML::Exception& e) {
        throw std::runtime_error(std::string("input (YAML): ") + e.what());
    }
}

Config load_yaml(const std::string& path) {
    try {
        return from_node(YAML::LoadFile(path));
    } catch (const YAML::Exception& e) {
        throw std::runtime_error(path + ": " + e.what());
    }
}

Config parse_command_line(const std::vector<std::string>& a) {
    Config c;
    try {
        if (!a.empty() && a[0] == "--all") {
            if (a.size() < 3 || a.size() > 4) throw std::runtime_error("");
            c.method = Method::Spectrum;
            c.Oh = std::sqrt(to_number(a[1]));
            c.modes.push_back({to_int(a[2]), std::nullopt});
            if (a.size() == 4) c.Nr = to_int(a[3]);
        } else {
            if (a.size() < 2 || a.size() > 6 || a.size() == 5) throw std::runtime_error("");
            c.Oh = std::sqrt(to_number(a[0]));
            Mode md;
            md.l = to_int(a[1]);
            if (a.size() > 2) {  // m: allowed for completeness; no effect for a sphere
                const int m = to_int(a[2]);
                if (std::abs(m) > md.l) throw std::runtime_error("|m| must not exceed l");
            }
            if (a.size() > 3) c.Nr = to_int(a[3]);
            if (a.size() == 6) md.guess = std::complex<double>(to_number(a[4]), to_number(a[5]));
            c.modes.push_back(md);
        }
    } catch (const std::logic_error&) {  // from stod/stoi
        throw std::runtime_error("bad number in command line");
    } catch (const std::runtime_error& e) {
        if (std::string(e.what()).empty())
            throw std::runtime_error(
                "usage: lambdaDrop We l [m] [Nr] [guessRe guessIm] | lambdaDrop --all We l [Nr]");
        throw;
    }
    c.validate();
    return c;
}

}  // namespace lambda
