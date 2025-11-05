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

static double liner_inter(double x0, double y0, double x1, double y1, double x)
{
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
    case -1: return 1;
    case 19: return 2;
    case 18: return 3;
    case 21: return 4;
    default: return -1;
  }
}

struct FastTable {
    std::vector<double> xs;
    std::vector<std::array<double,COLUMN_COUNT>> ys;

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
    
    double value(int col, double x) {
        if (col < 0 || col >= int(COLUMN_COUNT) || xs.empty()) return 0.0;
        if (x <= xs.front()) return ys.front()[col];
        if (x >= xs.back()) return ys.back()[col];

        // Binary search for the interval [xs[i], xs[i+1]] containing x
        size_t left = 0, right = xs.size() - 1;
        while (left + 1 < right) {
            size_t mid = left + (right - left) / 2;
            if (xs[mid] <= x)
                left = mid;
            else
                right = mid;
        }
        return liner_inter(xs[left], ys[left][col], xs[left+1], ys[left+1][col], x);
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
    double f_val = f_table.value(nE_k_to_col(nE, kappa), pval * MEV_AU); 
    double g_val = g_table.value(nE_k_to_col(nE, kappa), pval * MEV_AU);
    delete[] f_table.ys.data();
    delete[] g_table.ys.data();
    
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

int main() {
    // Example parameters
    int nE = 2;
    double kappa = -1;
    double m = 0.5;
    double p[3] = {0.1, 0.2, 0.3};

    double u_out[8] = {0};

    spinor_wave_function(u_out, nE, kappa, m, p);

    std::cout << "spin_wave_function output:" << std::endl;
    for (int i = 0; i < 8; ++i) {
        std::cout << "u_out:[" << i << "] = " << u_out[i] << std::endl;
    }

    return 0;
}

#endif