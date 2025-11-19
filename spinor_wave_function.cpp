/*
Dependence: libboost-math

Build shared library:
g++ -g -O3 -shared -fPIC spinor_wave_function.cpp -o libspinor_wave_function.so -lboost_math_c99

Build test executable:
g++ -g -O3 -DNEED_TEST_MAIN spinor_wave_function.cpp -o spinor_wave_function_test -lboost_math_c99
*/
#include <vector>
#include <array>
#include <cstdint>
#include <iostream>
#include <fstream>
#include <array>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <cmath>

const double MEV_AU = 268.102; // 1 MeV = 268.102 a.u.
const size_t COLUMN_COUNT = 4; 
const double EPS = 1e-13;

static double linear_inter(double x0, double y0, double x1, double y1, double x)
{
  // exception
  if (x1 - x0 < EPS) {
    return 0.5 * (y0 + y1);
  }
  return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
}



static int nE_k_to_col(int nE, double kappa){
  // 1s  -> E=1 , j = 1/2, kappa = -1 -> _index = -1 --> col = 1
  // 2s  -> E=2 , j = 1/2, kappa = -1 -> _index = 9  --> col = 2
  // 2p- -> E=2 , j = 3/2, kappa = -2 -> _index = 8  --> col = 3
  // 2p+ -> E=2 , j = 1/2, kappa = 1  -> _index = 11 --> col = 4
  int _index = (int)( kappa + (nE - 1) * 20); 
  switch (_index)
  {
    case -1: return 0;
    case 19: return 1;
    case 18: return 2;
    case 21: return 3;
    default: return -1;
  }
}

struct FastTable {
    std::vector<double> xs;
    std::vector<std::array<double,COLUMN_COUNT>> ys;
    mutable bool M_ready[COLUMN_COUNT] = {false};
    mutable std::vector<double> M_cache[COLUMN_COUNT]; 

    bool load_bin(const char* binfile)
    {
        std::ifstream in(binfile, std::ios::binary);
        if (!in.is_open()) return false;
        uint64_t nrows = 0;
        uint32_t ncols = 0;
        in.read(reinterpret_cast<char*>(&nrows), sizeof(nrows));
        in.read(reinterpret_cast<char*>(&ncols), sizeof(ncols));
        if (nrows == 0 || ncols < 2) return false; 

        xs.resize(nrows);
        ys.assign(nrows, std::array<double, COLUMN_COUNT>{});

        in.read(reinterpret_cast<char*>(xs.data()), sizeof(double) * nrows);

        for (uint32_t col = 1; col <= COLUMN_COUNT; ++col) {
            std::vector<double> colbuf(nrows);
            if (col < ncols) {
                in.read(reinterpret_cast<char*>(colbuf.data()), sizeof(double) * nrows);
            } else {
                for (uint64_t i = 0; i < nrows; ++i) {
                    colbuf[i] = 0.0;
                }
            }
            for (uint64_t i = 0; i < nrows; ++i) {
                ys[i][col-1] = colbuf[i];
            }
        }

        // if columns more than COLUMN_COUNT + 1, skip remaining columns
        if (ncols > COLUMN_COUNT + 1) {
            uint64_t remaining = (uint64_t)(ncols - (COLUMN_COUNT + 1)) * nrows;
            in.seekg(static_cast<std::streamoff>(sizeof(double) * remaining), std::ios::cur);
        }

        in.close();
        return true;
    }
    
    double value(int col, double x, int interp=1)
    {
      if (col < 0 || col >= int(COLUMN_COUNT) || xs.empty())
        return 0.0;
      if (x <= xs.front())
        return ys.front()[col];
      if (x >= xs.back())
        return ys.back()[col];

      // Binary search for the interval [xs[i], xs[i+1]] containing x
      size_t left = 0, right = xs.size() - 1;
      while (left + 1 < right)
      {
        size_t mid = left + (right - left) / 2;
        if (xs[mid] <= x)
          left = mid;
        else
          right = mid;
      }
        switch (interp){
        
        case 1: // linear
        {
        return linear_inter(xs[left], ys[left][col], xs[left+1], ys[left+1][col], x);
        }
        case 2: // cubic spline
        {
          build_M(col);
          return cubic_spline_eval(col, left, x);
        }
        default:
            break;
        }
        // If interp is not 1 or 2, return 0.0 as a safe default
        return 0.0;
    }
    void build_M(int col) const {
        if (M_ready[col]) return;

        const size_t n = xs.size();
        auto &M = M_cache[col];
        M.assign(n, 0.0);

        if (n < 3) {
            M_ready[col] = true;
            return;
        }

        std::vector<double> h(n - 1, 0.0);
        for (size_t i = 0; i + 1 < n; ++i)
            h[i] = xs[i + 1] - xs[i];

        std::vector<double> lower(n, 0.0), diag(n, 0.0), upper(n, 0.0), rhs(n, 0.0);
        diag[0] = 1.0; rhs[0] = 0.0;
        diag[n - 1] = 1.0; rhs[n - 1] = 0.0;

        for (size_t i = 1; i + 1 < n; ++i) {
            const double hi_1 = h[i - 1];
            const double hi   = h[i];

            lower[i] = hi_1;
            diag[i]  = 2.0 * (hi_1 + hi);
            upper[i] = hi;

            const double dy_i   = (std::abs(hi)   < EPS) ? 0.0 : (ys[i + 1][col] - ys[i][col]) / hi;
            const double dy_im1 = (std::abs(hi_1) < EPS) ? 0.0 : (ys[i][col]     - ys[i - 1][col]) / hi_1;
            rhs[i] = 6.0 * (dy_i - dy_im1);


            if (std::abs(diag[i]) < EPS) diag[i] = EPS;
        }

        // Thomas Algorithm
        std::vector<double> cprime(n, 0.0), dprime(n, 0.0);
        double denom = (std::abs(diag[0]) < EPS) ? EPS : diag[0];
        cprime[0] = upper[0] / denom;
        dprime[0] = rhs[0] / denom;

        for (size_t i = 1; i < n; ++i) {
            double denom_i = diag[i] - lower[i] * cprime[i - 1];
            if (std::abs(denom_i) < EPS) denom_i = EPS;
            cprime[i] = (i + 1 < n) ? (upper[i] / denom_i) : 0.0;
            dprime[i] = (rhs[i] - lower[i] * dprime[i - 1]) / denom_i;
        }

        M[n - 1] = dprime[n - 1];
        for (size_t i = n - 1; i-- > 0; ) {
            M[i] = dprime[i] - cprime[i] * M[i + 1];
        }

        M_ready[col] = true;
    }
    inline double cubic_spline_eval(int col, size_t i, double x) const {
        const double xi  = xs[i];
        const double xi1 = xs[i + 1];
        const double yi  = ys[i][col];
        const double yi1 = ys[i + 1][col];
        const double hi  = xi1 - xi;
        if (std::abs(hi) < EPS) return 0.5 * (yi + yi1);

        const auto &M = M_cache[col];
        const double Mi  = M[i];
        const double Mi1 = M[i + 1];

        // a + b t + c t^2 + d t^3
        const double t = x - xi;
        const double a = yi;
        const double b = (yi1 - yi) / hi - (2.0 * Mi + Mi1) * hi / 6.0;
        const double c = 0.5 * Mi;
        const double d = (Mi1 - Mi) / (6.0 * hi);
        return ((d * t + c) * t + b) * t + a;
    }
};

static std::complex<double> spherical_harmonic(double l, double m, double theta, double phi)
{
  return boost::math::spherical_harmonic(round(l), round(m), theta, phi);
}

static std::array<std::complex<double>, 2> spin_half_spherical_harmonic(
    double l, double j, double m, double theta, double phi)
{
  if(j > l) {
    std::complex<double> upper = sqrt((j + m) / (2.0 * j)) * spherical_harmonic(l, m - 0.5, theta, phi);
    std::complex<double> lower = sqrt((j - m) / (2.0 * j)) * spherical_harmonic(l, m + 0.5, theta, phi);
    return { upper, lower };
  } else {
    std::complex<double> upper = -sqrt((j - m + 1.0) / (2.0 * j + 2.0)) * spherical_harmonic(l, m - 0.5, theta, phi);
    std::complex<double> lower = sqrt((j + m + 1.0) / (2.0 * j + 2.0)) * spherical_harmonic(l, m + 0.5, theta, phi);
    return { upper, lower };
  }
}

static std::array<std::complex<double>, 2> spin_half_spherical_harmonic(
    double kappa, double m, double theta, double phi)
{
  double l, j = std::abs(kappa) - 0.5;
  if(kappa > 0) {
    l = kappa;  // j = l - 0.5
  } else {
    l = -1 - kappa;  // j = l + 0.5
  }
  return spin_half_spherical_harmonic(l, j, m, theta, phi);
}

static std::array<double, 3> cartesian_to_spherical(double x, double y, double z) {
    double r = std::sqrt(x * x + y * y + z * z);
    double theta = (r == 0.0) ? 0.0 : std::acos(z / r); // polar angle
    double phi = std::atan2(y, x); // azimuthal angle
    return {r, theta, phi};
}
// pauli matrices
const std::array<std::array<std::complex<double>, 2>, 2> sigma_x = {{{0, 1}, {1, 0}}};
const std::array<std::array<std::complex<double>, 2>, 2> sigma_y = {{{0, std::complex<double>(0, -1)}, {std::complex<double>(0, 1), 0}}};
const std::array<std::array<std::complex<double>, 2>, 2> sigma_z = {{{1, 0}, {0, -1}}};



extern "C" void spinor_wave_function(double u_out[8], //output U[ x_1 real, x_1 imag, x_2 real, x_2 imag, y_1 real, y_1 imag, y_2 real, y_2 imag ]
     int nE, double kappa, double m, double p[3]){ //input
    double l, j = std::abs(kappa) - 0.5;
    if(kappa > 0) {
        l = kappa;  // j = l - 0.5
    } else {
        l = -1 - kappa;  // j = l + 0.5
    }
    // Convert Cartesian to Spherical coordinates
    auto [pval, theta, phi] = cartesian_to_spherical(p[0], p[1], p[2]);
    auto y = spin_half_spherical_harmonic(kappa, m, theta, phi);
    FastTable f_table, g_table;
    f_table.load_bin("./data/C_f_DBSR.bin");
    g_table.load_bin("./data/C_g_DBSR.bin");
    double f_val = f_table.value(nE_k_to_col(nE, kappa), pval * MEV_AU, 2);
    double g_val = g_table.value(nE_k_to_col(nE, kappa), pval * MEV_AU, 2);
    // std::cout << "f_val: " << f_val << ", g_val: " << g_val << std::endl;
    // delete[] f_table.ys.data();
    // delete[] g_table.ys.data();
    
    std::complex<double> I(0.0, 1.0);
    std::complex<double> factor = 4 * M_PI * std::pow(-I, l);
    // Compute normalized p vector
    double norm_p = std::sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
    double nx = (norm_p != 0.0) ? p[0] / norm_p : 0.0;
    double ny = (norm_p != 0.0) ? p[1] / norm_p : 0.0;
    double nz = (norm_p != 0.0) ? p[2] / norm_p : 0.0;

    // Compute normalized p dot [sigma_x, sigma_y, sigma_z] as a complex 2x2 matrix
    std::array<std::array<std::complex<double>, 2>, 2> sigma_dot_norm_p;
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            sigma_dot_norm_p[i][j] = nx * sigma_x[i][j] + ny * sigma_y[i][j] + nz * sigma_z[i][j];
        }
    }

    std::complex<double> u1 = factor * g_val * y[0];
    std::complex<double> u2 = factor * g_val * y[1];
    std::complex<double> v1 = factor * f_val * (sigma_dot_norm_p[0][0] * y[0] + sigma_dot_norm_p[0][1] * y[1]);
    std::complex<double> v2 = factor * f_val * (sigma_dot_norm_p[1][0] * y[0] + sigma_dot_norm_p[1][1] * y[1]);
    u_out[0] = u1.real();
    u_out[1] = u1.imag();
    u_out[2] = u2.real();
    u_out[3] = u2.imag();
    u_out[4] = v1.real();
    u_out[5] = v1.imag();
    u_out[6] = v2.real();
    u_out[7] = v2.imag();
}


#ifdef NEED_TEST_MAIN

int main(int argc, char** argv) {
    if (argc == 7) {
        int nE = std::atoi(argv[1]);
        double kappa = std::atof(argv[2]);
        double m = std::atof(argv[3]);
        double p[3] = {std::atof(argv[4]), std::atof(argv[5]), std::atof(argv[6])};
        double u_out[8] = {0};
        spinor_wave_function(u_out, nE, kappa, m, p);
        std::cout << u_out[0] << "," << u_out[1] << "," << u_out[2] << "," << u_out[3] << ","
                  << u_out[4] << "," << u_out[5] << "," << u_out[6] << "," << u_out[7] << std::endl;
    return 0;
    }
    // Example parameters
    int nE = 2;
    double kappa = -1;
    double m = 0.5;

    std::ofstream fout("spinor_wave_function_scan.csv");
    fout << "pz,u1_real,u1_imag,u2_real,u2_imag,v1_real,v1_imag,v2_real,v2_imag\n";
    double p[3] = {0, 0, 0.0};
    double u_out[8] = {0};
    for (double pz = 0.0; pz <= 0.05; pz += 1e-4) {
      p[2] = pz;
      spinor_wave_function(u_out, nE, kappa, m, p);
      fout << pz;
      for (int i = 0; i < 8; ++i) {
        fout << "," << u_out[i];
      }
      fout << "\n";
    }
    fout.close();
    std::cout << "Scan complete. Output written to spinor_wave_function_scan.csv" << std::endl;

    return 0;
}

#endif