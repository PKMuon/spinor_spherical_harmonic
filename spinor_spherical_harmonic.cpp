// Dependence: libboost-math
// g++ -g -O3 -shared -fPIC spinor_spherical_harmonic.cpp -o libspinor_spherical_harmonic.so
// g++ -g -O3 -DNEED_TEST_MAIN spinor_spherical_harmonic.cpp -o spinor_spherical_harmonic_test
#include <array>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <cmath>

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

extern "C" void spin_half_spherical_harmonic(double y_out[4], double kappa, double m, double theta, double phi)
{
  auto y = spin_half_spherical_harmonic(kappa, m, theta, phi);
  for(size_t i_out = 0, i = 0; i < y.size(); ++i) {
    y_out[i_out++] = y[i].real();
    y_out[i_out++] = y[i].imag();
  }
}

#ifdef NEED_TEST_MAIN

#include <iostream>

using namespace std::complex_literals;

int main()
{
  // Tunable parameters.
  double kappa = 2, m = 1.5, theta = M_PI / 4, phi = M_PI / 4;
  double j = std::abs(kappa) - 0.5;

  std::complex<double> A[2][2] = { 0 };  // Eq. (21a) in PhysRevC.50.2822
  std::complex<double> B[2][2] = { 0 };  // Eq. (21b) in PhysRevC.50.2822

  for(double m = -j; m <= j; m += 1.0) {
    // This is not robust, so please check to confirm that m iterates over [-j, j].
    std::cout << "j = " << j << " m = " << m << std::endl;
    auto y = spin_half_spherical_harmonic(kappa, m, theta, phi);
    auto z = spin_half_spherical_harmonic(-kappa, m, theta, phi);
    std::cout << "y = " << y[0] << ", " << y[1] << std::endl;
    A[0][0] += y[0] * std::conj(y[0]);
    A[0][1] += y[0] * std::conj(y[1]);
    A[1][0] += y[1] * std::conj(y[0]);
    A[1][1] += y[1] * std::conj(y[1]);
    B[0][0] += y[0] * std::conj(z[0]);
    B[0][1] += y[0] * std::conj(z[1]);
    B[1][0] += y[1] * std::conj(z[0]);
    B[1][1] += y[1] * std::conj(z[1]);
  }

  std::cout << "A:" << std::endl;
  for(size_t i0 = 0; i0 < 2; ++i0) {
    for(size_t i1 = 0; i1 < 2; ++i1) std::cout << " " << std::setw(30) << A[i0][i1] / ((2 * j + 1) / (8 * M_PI));
    std::cout << std::endl;
  }

  A[0][0] = A[1][1] = 1;
  A[0][1] = A[1][0] = 0;
  std::cout << "A ref:" << std::endl;
  for(size_t i0 = 0; i0 < 2; ++i0) {
    for(size_t i1 = 0; i1 < 2; ++i1) std::cout << " " << std::setw(30) << A[i0][i1];
    std::cout << std::endl;
  }

  std::cout << "B:" << std::endl;
  for(size_t i0 = 0; i0 < 2; ++i0) {
    for(size_t i1 = 0; i1 < 2; ++i1) std::cout << " " << std::setw(30) << B[i0][i1] / (-(2 * j + 1) / (8 * M_PI));
    std::cout << std::endl;
  }

  B[0][0] = cos(theta);
  B[0][1] = sin(theta) * std::exp(-1.0i * phi);
  B[1][0] = sin(theta) * std::exp(1.0i * phi);
  B[1][1] = -cos(theta);
  std::cout << "B ref:" << std::endl;
  for(size_t i0 = 0; i0 < 2; ++i0) {
    for(size_t i1 = 0; i1 < 2; ++i1) std::cout << " " << std::setw(30) << B[i0][i1];
    std::cout << std::endl;
  }
  return 0;
}

#endif /* NEED_TEST_MAIN */
