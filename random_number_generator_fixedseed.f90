! RANDOM SEED INITIATION !

SUBROUTINE init_random_seed(base_seed)

    use datatype
    implicit none
    integer(kind=I4), intent(in) :: base_seed
    integer(kind=I4) :: n, i
    integer(kind=I4), allocatable :: seed(:)

    call random_seed(size=n)
    allocate(seed(n))

    seed=base_seed+[(37*(i-1), i=1,n)]

    call random_seed(put=seed)
    deallocate(seed)

END SUBROUTINE


! GENERATE GAUSSIAN RANDOM NUMBER !

REAL(kind=DP) FUNCTION gaussian_random_number(mean,sigma)

    use datatype
    implicit none
    real(kind=DP), intent(in) :: mean,sigma
    real(kind=DP) :: pi,r1,r2

    pi=4.0d0*datan(1.0d0)
    call random_number(r1)
    call random_number(r2)
    gaussian_random_number=mean+sigma*dsqrt(-2.0d0*dlog(r1))*dcos(2.0d0*pi*r2)

ENDFUNCTION

