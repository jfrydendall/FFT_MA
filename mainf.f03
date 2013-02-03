program testFFTMA
    use fftMA
    use omp_lib
    implicit none
    integer,parameter  :: dp = selected_real_kind(15, 307)
    integer  :: nx = 120, ny = 80, nz = 1, Nrel = 100, npar = 4
    real(dp) :: cell = 1.0_dp, CorLenMax = 35.0_dp, CorLenMed = 20.0_dp, CorLenMin = 5.0_dp
    real(dp) :: theta = 0.0_dp, beta = 0.0_dp, alpha = 45.0_dp, rho = 0.8
    real(dp) :: mean = 0.0_dp, var = 1.0_dp
    character(len=20) :: name = 'Gaussian', name1 = 'Gaussian'
    type(fields)        :: rf, cosim
    ! local variables
    real(dp), allocatable, dimension(:,:,:)   :: C, C1
    real(dp), allocatable, dimension(:,:,:,:) :: Z, Z1, A, A1
    rf = fields(Nrel, npar, nx, ny, nz, cell, CorLenMax, CorLenMed, CorLenMin, theta, beta, alpha, mean, var, name)
    cosim = fields(Nrel, npar, nx, ny, nz, cell, CorLenMax, CorLenMed, CorLenMin, theta, beta, alpha, mean, var, name1)
    call rf%set()
    call rf%corr(C)
    call rf%normal(Z)
    call rf%generator(A,C,Z)
    write(100) A
    
    call cosim%set()
    call cosim%corr(C1)
    call cosim%normal(Z1)
    call cosim%generator(A1,C1,rho*Z + sqrt(1.0_dp-rho**2.0_dp)*Z1)
    write(200) A1
end program testFFTMA
