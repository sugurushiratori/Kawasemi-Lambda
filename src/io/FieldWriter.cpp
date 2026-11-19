#include "io/FieldWriter.h"

#include <array>
#include <bit>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <numbers>
#include <sstream>
#include <stdexcept>

namespace lambda {

namespace {

constexpr double kPi = std::numbers::pi;
using Vec3 = std::array<double, 3>;

std::ofstream open(const std::string& name) {
    std::ofstream f(name);
    if (!f) throw std::runtime_error("cannot open " + name);
    f << std::setprecision(12);
    return f;
}

// spherical (r, theta, phi) components -> Cartesian
std::array<cplx, 3> to_cartesian(cplx ur, cplx ut, cplx uf, double th, double ph) {
    const double st = std::sin(th), ct = std::cos(th), sp = std::sin(ph), cp = std::cos(ph);
    return {ur * st * cp + ut * ct * cp - uf * sp, ur * st * sp + ut * ct * sp + uf * cp,
            ur * ct - ut * st};
}

Vec3 point(double r, double th, double ph) {
    return {r * std::sin(th) * std::cos(ph), r * std::sin(th) * std::sin(ph), r * std::cos(th)};
}

std::string header(const Eigenfunction& ef, const char* what, bool with_m = true) {
    std::ostringstream h;
    h << std::setprecision(16);
    h << "# LAMBDA " << LAMBDA_VERSION << " eigenfunction: " << what << "\n"
      << "# l = " << ef.l();
    if (with_m) h << ", m = " << ef.m();
    h << ", lambda = " << ef.lambda().real() << " + "
      << ef.lambda().imag() << " i (viscous units)\n"
      << "# complex amplitudes q; field at time t = Re[q exp(lambda t)]\n"
      << "# normalised so that the surface is r = 1 + Re[d Y exp(lambda t)], d = 1, max|Y| = 1\n"
      << "# units: length R, velocity nu/R, pressure rho nu^2/R^2\n";
    return h.str();
}

// ---- ParaView XML StructuredGrid, binary (base64) Float32 ----
std::string base64(const unsigned char* p, std::size_t n) {
    static const char* tbl = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    std::string out;
    out.reserve(4 * ((n + 2) / 3));
    for (std::size_t i = 0; i < n; i += 3) {
        const unsigned v = (unsigned(p[i]) << 16) | ((i + 1 < n ? unsigned(p[i + 1]) : 0u) << 8) |
                           (i + 2 < n ? unsigned(p[i + 2]) : 0u);
        out += tbl[(v >> 18) & 63];
        out += tbl[(v >> 12) & 63];
        out += (i + 1 < n) ? tbl[(v >> 6) & 63] : '=';
        out += (i + 2 < n) ? tbl[v & 63] : '=';
    }
    return out;
}

// VTK "binary" inline format: base64 of [UInt32 byte count][data]
void write_binary(std::ostream& f, const std::vector<double>& data) {
    std::vector<unsigned char> buf(4 + 4 * data.size());
    const std::uint32_t nbytes = static_cast<std::uint32_t>(4 * data.size());
    std::memcpy(buf.data(), &nbytes, 4);
    for (std::size_t i = 0; i < data.size(); ++i) {
        const float x = static_cast<float>(data[i]);
        std::memcpy(buf.data() + 4 + 4 * i, &x, 4);
    }
    f << base64(buf.data(), buf.size()) << '\n';
}
struct VtsArray {
    std::string name;
    int ncomp;
    std::vector<double> data;
};

void write_vts(const std::string& name, int ni, int nj, int nk, const std::vector<Vec3>& pts,
               const std::vector<VtsArray>& arrays) {
    auto f = open(name);
    const std::string ext = "0 " + std::to_string(ni - 1) + " 0 " + std::to_string(nj - 1) +
                            " 0 " + std::to_string(nk - 1);
    f << "<?xml version=\"1.0\"?>\n"
      << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\""
      << (std::endian::native == std::endian::little ? "LittleEndian" : "BigEndian")
      << "\" header_type=\"UInt32\">\n"
      << "<StructuredGrid WholeExtent=\"" << ext << "\">\n<Piece Extent=\"" << ext << "\">\n"
      << "<PointData>\n";
    for (const auto& a : arrays) {
        f << "<DataArray type=\"Float32\" Name=\"" << a.name << "\" NumberOfComponents=\""
          << a.ncomp << "\" format=\"binary\">\n";
        write_binary(f, a.data);
        f << "</DataArray>\n";
    }
    f << "</PointData>\n<Points>\n"
      << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"binary\">\n";
    std::vector<double> xyz;
    xyz.reserve(3 * pts.size());
    for (const auto& p : pts) xyz.insert(xyz.end(), p.begin(), p.end());
    write_binary(f, xyz);
    f << "</DataArray>\n</Points>\n</Piece>\n</StructuredGrid>\n</VTKFile>\n";
}

}  // namespace

std::vector<std::string> write_eigenfunction(const Eigenfunction& ef,
                                             const EigenfunctionOutput& opt) {
    std::vector<std::string> files;
    const std::string base = opt.prefix + "_l" + std::to_string(ef.l());
    const std::string bm = base + "_m" + std::to_string(ef.m());
    const int nr = opt.n_r, nt = opt.n_theta, np = opt.n_phi;
    auto rr = [&](int i) { return static_cast<double>(i) / (nr - 1); };
    auto tt = [&](int j) { return kPi * j / (nt - 1); };
    auto pp = [&](int k) { return 2.0 * kPi * k / (np - 1); };

    if (opt.ascii) {
        if (ef.m() == opt.m.front()) {  // m-independent: once per mode
            const std::string fn = base + "_radial.dat";
            auto f = open(fn);
            f << header(ef, "radial functions (independent of m)", false)
              << "# r  Re(u) Im(u)  Re(v) Im(v)  Re(p) Im(p)\n"
              << "#   u_r = u Y, u_tan = v r grad Y, p = p Y\n";
            for (int i = 0; i < nr; ++i) {
                const auto a = ef.radial(rr(i));
                f << rr(i) << ' ' << a.u.real() << ' ' << a.u.imag() << ' ' << a.v.real() << ' '
                  << a.v.imag() << ' ' << a.p.real() << ' ' << a.p.imag() << '\n';
            }
            files.push_back(fn);
        }
        {
            const std::string fn = bm + "_meridional.dat";
            auto f = open(fn);
            f << header(ef, "meridional plane y = 0 (phi = 0 for x > 0, phi = pi for x < 0)")
              << "# r alpha x z  Re,Im of: u_x u_y u_z p  u_r u_theta u_phi\n"
              << "#   alpha in [0, 2 pi) is the polar angle measured round the plane;\n"
              << "#   u_x, u_z lie in the plane, u_y is normal to it (= u_phi on phi = 0)\n"
              << "# blocks of constant alpha (gnuplot: splot ... with pm3d)\n";
            for (int j = 0; j < 2 * nt - 1; ++j) {
                const double al = kPi * j / (nt - 1);
                const double th = al <= kPi ? al : 2.0 * kPi - al;
                const double ph = al <= kPi ? 0.0 : kPi;
                for (int i = 0; i < nr; ++i) {
                    const double r = rr(i);
                    const auto q = ef.field(r, th, ph);
                    const auto c = to_cartesian(q.ur, q.ut, q.uf, th, ph);
                    const Vec3 x = point(r, th, ph);
                    f << r << ' ' << al << ' ' << x[0] << ' ' << x[2];
                    for (cplx z : {c[0], c[1], c[2], q.p, q.ur, q.ut, q.uf})
                        f << ' ' << z.real() << ' ' << z.imag();
                    f << '\n';
                }
                f << '\n';
            }
            files.push_back(fn);
        }
        {
            const std::string fn = bm + "_surface.dat";
            auto f = open(fn);
            f << header(ef, "surface r = 1")
              << "# theta phi  Re,Im of: dY  u_r u_theta u_phi p\n";
            for (int j = 0; j < nt; ++j) {
                for (int k = 0; k < np; ++k) {
                    const cplx s = ef.surface(tt(j), pp(k));
                    const auto q = ef.field(1.0, tt(j), pp(k));
                    f << tt(j) << ' ' << pp(k) << ' ' << s.real() << ' ' << s.imag() << ' '
                      << q.ur.real() << ' ' << q.ur.imag() << ' ' << q.ut.real() << ' '
                      << q.ut.imag() << ' ' << q.uf.real() << ' ' << q.uf.imag() << ' '
                      << q.p.real() << ' ' << q.p.imag() << '\n';
                }
                f << '\n';
            }
            files.push_back(fn);
        }
    }

    if (opt.vtk) {
        {  // volume: r fastest, then theta, then phi
            std::vector<Vec3> pts;
            VtsArray vre{"velocity_re", 3, {}}, vim{"velocity_im", 3, {}};
            VtsArray pre{"pressure_re", 1, {}}, pim{"pressure_im", 1, {}};
            VtsArray sre{"u_r_theta_phi_re", 3, {}}, sim{"u_r_theta_phi_im", 3, {}};
            for (int k = 0; k < np; ++k)
                for (int j = 0; j < nt; ++j)
                    for (int i = 0; i < nr; ++i) {
                        const double r = rr(i), th = tt(j), ph = pp(k);
                        pts.push_back(point(r, th, ph));
                        const auto q = ef.field(r, th, ph);
                        const auto c = to_cartesian(q.ur, q.ut, q.uf, th, ph);
                        for (const auto& z : c) {
                            vre.data.push_back(z.real());
                            vim.data.push_back(z.imag());
                        }
                        for (cplx z : {q.ur, q.ut, q.uf}) {
                            sre.data.push_back(z.real());
                            sim.data.push_back(z.imag());
                        }
                        pre.data.push_back(q.p.real());
                        pim.data.push_back(q.p.imag());
                    }
            const std::string fn = bm + ".vts";
            write_vts(fn, nr, nt, np, pts, {vre, vim, pre, pim, sre, sim});
            files.push_back(fn);
        }
        {  // surface (theta x phi); displacement as a vector for "Warp By Vector"
            std::vector<Vec3> pts;
            VtsArray dre{"displacement_re", 1, {}}, dim{"displacement_im", 1, {}};
            VtsArray wre{"displacement_vector_re", 3, {}}, wim{"displacement_vector_im", 3, {}};
            VtsArray vre{"velocity_re", 3, {}}, vim{"velocity_im", 3, {}};
            VtsArray pre{"pressure_re", 1, {}}, pim{"pressure_im", 1, {}};
            for (int k = 0; k < np; ++k)
                for (int j = 0; j < nt; ++j) {
                    const double th = tt(j), ph = pp(k);
                    const Vec3 e = point(1.0, th, ph);
                    pts.push_back(e);
                    const cplx s = ef.surface(th, ph);
                    dre.data.push_back(s.real());
                    dim.data.push_back(s.imag());
                    for (double x : e) {
                        wre.data.push_back(s.real() * x);
                        wim.data.push_back(s.imag() * x);
                    }
                    const auto q = ef.field(1.0, th, ph);
                    for (const auto& z : to_cartesian(q.ur, q.ut, q.uf, th, ph)) {
                        vre.data.push_back(z.real());
                        vim.data.push_back(z.imag());
                    }
                    pre.data.push_back(q.p.real());
                    pim.data.push_back(q.p.imag());
                }
            const std::string fn = bm + "_surface.vts";
            write_vts(fn, nt, np, 1, pts, {dre, dim, wre, wim, vre, vim, pre, pim});
            files.push_back(fn);
        }
    }
    return files;
}

}  // namespace lambda
