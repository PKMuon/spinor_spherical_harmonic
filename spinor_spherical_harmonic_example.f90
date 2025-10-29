! gfortran -g -O3 spinor_spherical_harmonic_example.f90 spinor_spherical_harmonic.cpp -o spinor_spherical_harmonic_example -lstdc++
program spinor_spherical_harmonic_example
    use iso_c_binding
    implicit none

    interface
        subroutine spin_half_spherical_harmonic(y_out, kappa, m, theta, phi) bind(C)
            import :: c_double
            real(c_double), intent(out) :: y_out(4)
            real(c_double), value :: kappa, m, theta, phi
        end subroutine
    end interface

    real(c_double) :: y(4)
    call spin_half_spherical_harmonic(y, 2.0_c_double, 1.5_c_double, 0.7853981633974483_c_double, 0.7853981633974483_c_double)
    print '("y = ", 4F10.6)', y
end program spinor_spherical_harmonic_example
