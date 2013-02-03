module fftma
    use omp_lib
    use mod_fftw3
    implicit none
    integer,parameter                                       :: dp = selected_real_kind(15, 307)
    integer, public, save                                   :: nxx, nyy, nzz, L = 2*10
    real(dp), public, save                                  :: delta, gamma, scale, std
    private
    real(dp)                                                :: pi = 4.0_dp * atan(1.0_dp) ! Class-wide private constant

    type, public                                            :: fields
        integer                                             :: Nrel, npar
        integer                                             :: Nx, Ny, Nz
        real(dp)                                            :: cell, CorLenMax, CorLenMed, CorLenMin, theta, beta, alpha
        real(dp)                                            :: mean, var
        character(len = 11)                                 :: name
    contains
        procedure :: dist => distance
        procedure :: Corr => correlation_functions
        procedure :: Generator => randomfield
        procedure :: normal => rand_par
        procedure :: set => initialize
    end type fields
    contains
    subroutine initialize(this)
        class(fields), intent(inout) :: this
        this % CorLenMax = this % CorLenMax/this % cell
        this % CorLenMed = this % CorLenMed/this % cell
        this % CorLenMin = this % CorLenMin/this % cell
        std = sqrt(this % var)
        nxx = this % Nx + L
        nyy = this % Ny + L
        nzz = this % Nz + L
        scale = 1.0_dp/(nxx * nyy * nzz)

        if (this % CorLenMin > this % CorLenMax) then
            print *, 'Warning:'
            print *, 'Choose CorLenMin<=CorLenMax because CorLenMax is the direction of maximum continuity.'
            print *, 'Use the input "ang" to change the direction of maximum continutiy'
            stop
        end if

        gamma = this % CorLenMin / this % CorLenMax ! anistropy factor < 1 (=1 for isotropy)
        delta = this % CorLenMin / this % CorLenMed
        this % theta = pi/180.0_dp * this % theta ! Transform angle into radians
        this % beta = pi/180.0_dp * this % beta ! Transform angle into radians
        this % alpha = pi/180.0_dp * this % alpha ! Transform angle into radians
    end subroutine initialize

    subroutine distance(this,r)!     result(r)
        class(fields), intent(in)                          :: this
        ! local variables
        integer :: iret, i, j, k
        real(dp) :: Xrot, Yrot, Zrot, xx, yy, zz, t1, t2
        real(dp), intent(out), allocatable, dimension(:,:,:)   :: r

        ! Geometry:

        call omp_set_num_threads(this%npar)
        allocate(r(-floor(Nxx/2.0)+1:ceiling(Nxx/2.0), -floor(Nyy/2.0)+1:ceiling(Nyy/2.0), -floor(Nzz/2.0)+1:ceiling(Nzz/2.0)))
        t1 = omp_get_wtime()
        !$omp parallel do private(i,j,k)
        do k = -floor(Nzz/2.0)+1,ceiling(Nzz/2.0)
            do j = -floor(Nyy/2.0)+1,ceiling(Nyy/2.0)
                do i = -floor(Nxx/2.0)+1,ceiling(Nxx/2.0)
                    ! Transform into rotated coordinates:
                    xx = i + this%cell/2.0_dp; yy = j + this%cell/2.0_dp; zz = k + this%cell/2.0_dp
                    Xrot = (cos(this%theta) * cos(this%alpha) - sin(this%theta) * sin(this%beta) * sin(this%alpha)) * xx + &
                    (cos(this%theta) * sin(this%alpha) + sin(this%theta) * sin(this%beta) * cos(this%alpha)) * yy + &
                    (-sin(this%theta) * cos(this%beta)) * zz

                    Yrot = (-cos(this%beta) * sin(this%alpha)) * xx + &
                    (cos(this%beta) * cos(this%alpha)) * yy + &
                    (sin(this%beta)) * zz

                    Zrot = (sin(this%theta) * cos(this%alpha) + cos(this%theta) * sin(this%beta) * sin(this%alpha)) * xx + &
                    (sin(this%theta) * sin(this%alpha) - cos(this%theta) * sin(this%beta) * cos(this%alpha)) * yy + &
                    (cos(this%theta) * cos(this%beta)) * zz

                    r(i, j, k) = sqrt((delta * Xrot)**2.0_dp + (gamma * Yrot)**2.0_dp + Zrot**2.0_dp)
                end do
            end do
        end do
        !$omp end parallel do
        t2 = omp_get_wtime()
        print *, 'Time: ', t2 - t1
        end subroutine distance

    subroutine randomfield(this,A,CF,chi)
        class(fields), intent(in)                                    :: this
        real(dp), allocatable, intent(out), dimension(:,:,:,:)       :: A
        real(dp), intent(in), optional, dimension(:,:,:)             :: CF
        real(dp), intent(in), optional, dimension(:,:,:,:)           :: chi
        ! Local Variables
        integer                                                      :: iret, i, j, k
        integer(8)                                                   :: planZ, planC, planR
        real(dp)                                                     :: t1, t2
        real(dp), allocatable, dimension(:,:,:)                      :: C, phi
        real(dp), allocatable, dimension(:,:,:,:)                    :: Z
        complex(dp), allocatable, dimension(:,:,:)                   :: fC, fZ, fR

        allocate(Z(1:this % Nrel,1:nxx, 1:nyy, 1:nzz), C(1:nxx, 1:nyy, 1:nzz), A(1:this%Nrel,1:this%nx, 1:this%ny, 1:this%nz))
        allocate(fZ(1:nxx/2+1, 1:nyy, 1:nzz), fC(1:nxx/2+1, 1:nyy, 1:nzz), fR(1:nxx/2+1, 1:nyy, 1:nzz), phi(1:nxx, 1:nyy, 1:nzz))

        call dfftw_init_threads(iret)
        if (iret == 0) then
            print *, 'Something went wrong with initialation of the multi threaded fftw - about'
            stop
        end if
        call dfftw_plan_with_nthreads(this%npar)
        call dfftw_plan_dft_r2c_3d(planZ, nxx, nyy, nzz, Z, fZ, FFTW_MEASURE)
        call dfftw_plan_dft_r2c_3d(planC, nxx, nyy, nzz, C, fc, FFTW_MEASURE)
        call dfftw_plan_dft_c2r_3d(planR, nxx, nyy, nzz, fR, phi, FFTW_MEASURE)

        if (present(CF)) then
            C = CF
        else
            call this%corr(C)
        end if

        if (present(chi)) then
            Z = chi
        else
            call this % normal(Z)
        end if
        t1 = omp_get_wtime()
        do i = 1, this % Nrel
            call dfftw_execute_dft_r2c(planC, C, fC)
            ! make the fourie transform of the white noise
            call dfftw_execute_dft_r2c(planZ, Z(i,:,:,:), fZ)
            ! calculate the square root of the fourie transformed covaraince structure and multiply with the
            ! transformed white noise sequence
            fR = std * sqrt(fc) * fZ
            ! Back transform everything
            call dfftw_execute_dft_c2r(planR, fR, phi)
            !            ! slice to the desired size and add the mean value
            A(i,:,:,:) = phi(1:this % nx, 1:this % ny, 1:this % nz) * scale + this % mean
        end do
        t2 = omp_get_wtime()
        print *, 'Time :',t2-t1
        !Deallocate to conserve memory
        deallocate(fZ, fC, phi, fR, C, Z)
        call dfftw_destroy_plan(planC)
        call dfftw_destroy_plan(planZ)
        call dfftw_destroy_plan(planR)
    end subroutine randomfield

    subroutine correlation_functions(this,C)
        class(fields), intent(in)                             :: this
        real(dp),intent(out), allocatable, dimension(:,:,:)   :: C
        ! local variables
        real(dp), allocatable, dimension(:,:,:)               :: r

        call this%dist(r)
        allocate(C(1:nxx,1:nyy,1:nzz))
        if (trim(this % name) == 'Spherical') then
            where (r > this % CorLenMin)
                r = this % CorLenMin
            end where
            C = 1.0_dp - (1.5_dp * (r/this % CorLenMin) - 0.5_dp * (r/this % CorLenMin)**3.0_dp)
        else if (trim(this % name) == 'Exponential') then
            C = exp(-3.0_dp * r/this % CorLenMin)
        else if (trim(this % name) == 'Gaussian') then
            C = exp(-3.0_dp * r**2.0_dp/this % CorLenMin**2.0_dp)
        end if
        deallocate(r)
    end subroutine correlation_functions

    subroutine rand_par(this,ru)
        use par_seed_mod
        use par_zig_mod
        class(fields), intent(in)                                    :: this
        real(dp),intent(out),allocatable,dimension(:,:,:,:)           :: ru
        real(dp) :: t1, t2
        integer :: grainsize = 32, i, j, k, l, kpar
        integer, allocatable :: seed(:)
        allocate(seed(this%npar))
        call omp_set_num_threads(this%npar)
        call set_time_seed(this%npar)
        allocate(ru(1:this % Nrel,1:nxx,1:nyy,1:nzz))
        !$omp parallel do private(i,j,k,l,kpar)
        do k = 1, nzz
            do j = 1, nyy
                do i = 1, nxx
                    do l = 1, this % Nrel
                        kpar = omp_get_thread_num()
                        ru(l, i, j, k) = par_rnor(kpar)
                    end do
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine rand_par

End Module fftma
