! Build (example):
! gfortran -g -O3 spinor_wave_function_example.f90 spinor_wave_function.cpp -o spinor_wave_function_example -lstdc++ -lboost_math_c99
program spinor_wave_function_example
    use iso_c_binding
    implicit none

    interface
        subroutine spinor_wave_function(u_out, nE, kappa, m, p) bind(C)
            use iso_c_binding, only: c_double, c_int
            real(c_double), intent(out) :: u_out(8)
            integer(c_int), value :: nE
            real(c_double), value :: kappa, m
            real(c_double), intent(in) :: p(3)
        end subroutine
    end interface

    real(c_double) :: u(8)
    integer(c_int) :: nE
    real(c_double) :: kappa, m, p(3)

    ! Example parameters
    nE = 2
    kappa = -1.0_c_double
    m = 0.5_c_double
    p = (/ 0.1_c_double, 0.2_c_double, 0.3_c_double /)

    call spinor_wave_function(u, nE, kappa, m, p)

    print '("u_out = ", 8F12.9)', u
end program spinor_wave_function_example
