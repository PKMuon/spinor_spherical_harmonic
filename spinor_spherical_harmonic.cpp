// Dependence: libboost-math
// g++ -g -O3 -shared -fPIC spinor_spherical_harmonic.cpp -o libspinor_spherical_harmonic.so
// g++ -g -O3 -DNEED_TEST_MAIN spinor_spherical_harmonic.cpp -o spinor_spherical_harmonic_test
#include <array>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <cmath>
#include <complex>

using namespace std;
using namespace complex_literals;

namespace {  // Internal symbols only visible within this translation unit.

complex<double> spherical_harmonic(double l, double m, double theta, double phi)
{
  return boost::math::spherical_harmonic(round(l), round(m), theta, phi);
}

array<complex<double>, 2> spin_half_spherical_harmonic(double l, double j, double m, double theta, double phi)
{
  if(j > l) {
    complex<double> upper = sqrt((j + m) / (2.0 * j)) * spherical_harmonic(l, m - 0.5, theta, phi);
    complex<double> lower = sqrt((j - m) / (2.0 * j)) * spherical_harmonic(l, m + 0.5, theta, phi);
    return { upper, lower };
  } else {
    complex<double> upper = -sqrt((j - m + 1.0) / (2.0 * j + 2.0)) * spherical_harmonic(l, m - 0.5, theta, phi);
    complex<double> lower = sqrt((j + m + 1.0) / (2.0 * j + 2.0)) * spherical_harmonic(l, m + 0.5, theta, phi);
    return { upper, lower };
  }
}

array<complex<double>, 2> spin_half_spherical_harmonic(double kappa, double m, double theta, double phi)
{
  double l, j = abs(kappa) - 0.5;
  if(kappa > 0) {
    l = kappa;  // Corresponds to j = l - 0.5.
  } else {
    l = -1 - kappa;  // Corresponds to j = l + 0.5.
  }
  return spin_half_spherical_harmonic(l, j, m, theta, phi);
}

}  // namespace

// External interface.
// Uses the C ABI for compatibility with C, C++, Fortran, and Python.
extern "C" void spin_half_spherical_harmonic(double y_out[4], double kappa, double m, double theta, double phi)
{
  auto y = spin_half_spherical_harmonic(kappa, m, theta, phi);
  for(size_t i_out = 0, i = 0; i < y.size(); ++i) {
    y_out[i_out++] = y[i].real();
    y_out[i_out++] = y[i].imag();
  }
}

#ifdef NEED_TEST_MAIN

#include <iomanip>
#include <iostream>

int main()
{
  // Adjustable test parameters.
  double kappa = 2, m = 1.5, theta = M_PI / 4, phi = M_PI / 4;
  double j = abs(kappa) - 0.5;

  complex<double> A[2][2] = { 0 };  // Eq. (21a) in Phys. Rev. C 50, 2822.
  complex<double> B[2][2] = { 0 };  // Eq. (21b) in Phys. Rev. C 50, 2822.

  for(double m = -j; m <= j; m += 1.0) {
    // This loop is numerically stable since m can be represented exactly as a double.
    cout << "j = " << j << " m = " << m << endl;
    auto y = spin_half_spherical_harmonic(kappa, m, theta, phi);
    auto z = spin_half_spherical_harmonic(-kappa, m, theta, phi);
    cout << "y = " << y[0] << ", " << y[1] << endl;

    A[0][0] += y[0] * conj(y[0]);
    A[0][1] += y[0] * conj(y[1]);
    A[1][0] += y[1] * conj(y[0]);
    A[1][1] += y[1] * conj(y[1]);

    B[0][0] += y[0] * conj(z[0]);
    B[0][1] += y[0] * conj(z[1]);
    B[1][0] += y[1] * conj(z[0]);
    B[1][1] += y[1] * conj(z[1]);
  }

  cout << "A:" << endl;
  for(size_t i0 = 0; i0 < 2; ++i0) {
    for(size_t i1 = 0; i1 < 2; ++i1) cout << " " << setw(30) << A[i0][i1] / ((2 * j + 1) / (8 * M_PI));
    cout << endl;
  }

  A[0][0] = A[1][1] = 1;
  A[0][1] = A[1][0] = 0;
  cout << "A (reference):" << endl;
  for(size_t i0 = 0; i0 < 2; ++i0) {
    for(size_t i1 = 0; i1 < 2; ++i1) cout << " " << setw(30) << A[i0][i1];
    cout << endl;
  }

  cout << "B:" << endl;
  for(size_t i0 = 0; i0 < 2; ++i0) {
    for(size_t i1 = 0; i1 < 2; ++i1) cout << " " << setw(30) << B[i0][i1] / (-(2 * j + 1) / (8 * M_PI));
    cout << endl;
  }

  B[0][0] = cos(theta);
  B[0][1] = sin(theta) * exp(-1.0i * phi);
  B[1][0] = sin(theta) * exp(1.0i * phi);
  B[1][1] = -cos(theta);
  cout << "B (reference):" << endl;
  for(size_t i0 = 0; i0 < 2; ++i0) {
    for(size_t i1 = 0; i1 < 2; ++i1) cout << " " << setw(30) << B[i0][i1];
    cout << endl;
  }
  return 0;
}

#endif /* NEED_TEST_MAIN */
