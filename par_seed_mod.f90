module par_seed_mod
   use omp_lib
   use par_zig_mod, only : par_zigset
   contains
   subroutine set_time_seed(npar,nseed)
    implicit none
    integer, parameter :: grainsize=32
    integer, intent(in) :: npar
    integer, intent(in),optional ::nseed
    integer,allocatable,dimension(:) :: seed,iseed
    integer :: clock,i,n=1
    real :: z
! ..
! .. Local Scalars ..
    character (10) :: time
! ..
! .. Intrinsic Functions ..
    INTRINSIC DATE_AND_TIME
! ..
! .. Executable Statements ..

!    print *,npar
    call random_seed(size=n)
    allocate(seed(0:npar-1),iseed(1:n))
    if(present(nseed))then 
       iseed = nseed
    else
       CALL SYSTEM_CLOCK(COUNT=clock)
       iseed = clock + 37 * (/ (i - 1, i = 1, n) /)
    end if
    call random_seed(PUT=iseed)
    do i=0,npar-1
        call random_number(z)
        seed(i) = int(123456789*z)
    end do
    call par_zigset(npar, seed, grainsize)
   end subroutine set_time_seed
end module par_seed_mod
