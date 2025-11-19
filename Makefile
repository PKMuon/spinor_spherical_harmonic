CXX = g++
FC = gfortran
CXXFLAGS = -g -O3
SHARED_FLAGS = -shared -fPIC
BOOST_LIB = -lboost_math_c99

.PHONY: all clean

all: libspinor_spherical_harmonic.so libspinor_wave_function.so \
	spinor_spherical_harmonic_example spinor_spherical_harmonic_valid \
	spinor_wave_function_example spinor_spherical_harmonic_test spinor_wave_function_test convert

libspinor_spherical_harmonic.so: spinor_spherical_harmonic.cpp
	$(CXX) $(CXXFLAGS) $(SHARED_FLAGS) $< -o $@

libspinor_wave_function.so: spinor_wave_function.cpp
	$(CXX) $(CXXFLAGS) $(SHARED_FLAGS) $< -o $@ $(BOOST_LIB)

spinor_spherical_harmonic_example: spinor_spherical_harmonic_example.f90 spinor_spherical_harmonic.cpp
	$(FC) $(CXXFLAGS) $^ -o $@ -lstdc++

spinor_spherical_harmonic_valid: spinor_spherical_harmonic_valid.f90 spinor_spherical_harmonic.cpp
	$(FC) $(CXXFLAGS) $^ -o $@ -lstdc++

spinor_wave_function_example: spinor_wave_function_example.f90 spinor_wave_function.cpp
	$(FC) $(CXXFLAGS) $^ -o $@ -lstdc++ $(BOOST_LIB)

spinor_spherical_harmonic_test: spinor_spherical_harmonic.cpp
	$(CXX) $(CXXFLAGS) -DNEED_TEST_MAIN $< -o $@

spinor_wave_function_test: spinor_wave_function.cpp
	$(CXX) $(CXXFLAGS) -DNEED_TEST_MAIN $< -o $@ $(BOOST_LIB)

convert: convert.cpp
	$(CXX) -o $@ $<

clean:
	rm *.so spinor_spherical_harmonic_example spinor_spherical_harmonic_valid \
		spinor_wave_function_example spinor_spherical_harmonic_test spinor_wave_function_test convert