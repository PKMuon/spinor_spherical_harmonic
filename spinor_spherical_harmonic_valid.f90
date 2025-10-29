! gfortran -g -O3 spinor_spherical_harmonic_valid.f90 spinor_spherical_harmonic.cpp -o spinor_spherical_harmonic_valid -lstdc++
program spinor_spherical_harmonic_valid
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
    real(c_double) :: kappa, m, theta, phi
    integer :: iargc
    character(len=100) :: arg

    if (command_argument_count() /= 4) then
        print *, "Usage: ./spinor_spherical_harmonic_valid kappa m theta phi"
        stop 1
    end if
    call get_command_argument(1, arg); read(arg, *) kappa
    call get_command_argument(2, arg); read(arg, *) m
    call get_command_argument(3, arg); read(arg, *) theta
    call get_command_argument(4, arg); read(arg, *) phi

    call spin_half_spherical_harmonic(y, kappa, m, theta, phi)
    write(*,'(ES20.12,A,ES20.12,A,ES20.12,A,ES20.12)') y(1), char(9), y(2), char(9), y(3), char(9), y(4)
end program spinor_spherical_harmonic_valid
