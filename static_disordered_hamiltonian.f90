! THIS PROGRAM CONSTRUCTS STATIC DISORDERED ELECTRONIC HAMILTONIAN FOR RANDOM POLYMER CHAIN CONFORMATION !
! AND CALCULATES ITS EIGENSTATES/EIGENVALUES !

!****************************************************************************
! Packs two integers i and j in the lower triangular format. Assumes i >= j
!****************************************************************************

INTEGER(kind=I4) FUNCTION ijpk(i,j)

    use datatype
    implicit none
    integer(kind=I4), intent(in) :: i,j

    ijpk=(i*(i-1))/2+j

ENDFUNCTION

SUBROUTINE hamiltonian(ichain,n_site,alpha,beta,sigma_alpha_static,sigma_beta_static,random_int,obc_or_pbc,folder_path,index)

    use datatype
    implicit none

    character(len=20), intent(in) :: obc_or_pbc
    character(len=1000), intent(in) :: folder_path
    integer(kind=I4), intent(in) :: ichain,n_site,random_int !,obc_or_pbc
    real(kind=DP), intent(in) :: alpha,beta,sigma_alpha_static,sigma_beta_static
    integer(kind=I4), intent(out) :: index
    
    integer(kind=I4) :: ii,jj,icount
    integer(kind=I4), external :: ijpk
    real(kind=DP), external :: gaussian_random_number 
    character(len=20) :: str
    character(len=2000) :: filename
    real(kind=DP), allocatable :: H(:,:),E(:),C(:,:)
    real(kind=DP), allocatable :: matrix_upper_triangle(:)

    integer(kind=I4) :: mevlfnd,info
    real(kind=DP) :: abstol,VL,VU,dlamch
    integer(kind=I4), allocatable :: ifail(:),iwork(:)
    real(kind=DP), allocatable :: work(:)

!****************************************************************************
!    Allocation of the matrices needed for external diagonlization routine
!****************************************************************************

    abstol=2*dlamch('S')
    allocate(work(8*n_site),iwork(5*n_site),ifail(n_site))
    work=0.0d0; iwork=0; ifail=0


!****************************************************************************
!    Allocation of the necessary matrices
!****************************************************************************

    allocate(H(1:n_site,1:n_site),E(1:n_site),C(1:n_site,1:n_site),matrix_upper_triangle(1:int((n_site+1)*n_site/2)))
    H=0.0d0
    E=0.0d0
    C=0.0d0
    matrix_upper_triangle=0.0d0

!****************************************************************************
! Setting up the Hamiltonian matrix from input parameters  !
!****************************************************************************

!****************************************************************************
! Diagnoal elements !
!****************************************************************************

    do ii=1,n_site
        if(dabs(sigma_alpha_static) > 0.0d0) then
            H(ii,ii)=gaussian_random_number(alpha,sigma_alpha_static)
         else
            H(ii,ii)=alpha
         endif
    enddo

!****************************************************************************
! Off-diagnoal elements !
!****************************************************************************

    do ii=1,n_site-1
        if(dabs(sigma_beta_static) > 0.0d0) then
            H(ii,ii+1)=gaussian_random_number(beta,sigma_beta_static)
            H(ii+1,ii)=H(ii,ii+1)
        else
            H(ii,ii+1)=beta
            H(ii+1,ii)=H(ii,ii+1)
        endif
    enddo

    if (obc_or_pbc == 'pbc') then
        if(dabs(sigma_beta_static) > 0.0d0) then
            H(n_site,1)=gaussian_random_number(beta,sigma_beta_static)
            H(1,n_site)=H(n_site,1)
        else
            H(n_site,1)=beta
            H(1,n_site)=H(n_site,1)
        endif
    endif
    
    write(str,'(I0)') ichain 
    filename=trim(folder_path) // 'static_disordered_hamiltonian_matrix_' // trim(str) //'.out'
    open(unit=99,file=filename,form='unformatted')
    do ii=1,n_site
        write(99)(H(jj,ii),jj=1,n_site)
    enddo
    close(99)

!******************************************
! See the description of DSPEVX in the LAPACK package for upper triangulation of the matrix
! https://www.netlib.org/lapack/explore-html/d5/dba/group__hpevx_ga7a43ba4fd62f3eee58ca18baf2e9594c.html
!******************************************

    matrix_upper_triangle=0.0d0
    do ii=1,n_site
        do jj=1,ii
            icount=ijpk(ii,jj)
            matrix_upper_triangle(icount)=H(ii,jj)
       enddo
    enddo

    call dspevx('V','A','U',n_site,matrix_upper_triangle,VL,VU,1,n_site,abstol,mevlfnd,E,C,n_site,work,iwork,ifail,info)

    if(info/=0)then
        write(*,*)"matrix diagonalization failed by DSPEVX; info=",info
        index=1
    endif

    filename=trim(folder_path) // 'localized_state_energy_' // trim(str) //'.out'
    open(unit=99,file=filename)
    do ii=1,n_site
        write(99,'(F16.10)')E(ii)
    enddo
    close(99)

    filename=trim(folder_path) // 'localized_state_coefficient_' // trim(str) //'.out'
    open(unit=99,file=filename)
    do ii=1,n_site
        write(99,'(*(F16.10,3x))')(C(jj,ii),jj=1,n_site)
    enddo
    close(99)

    index=0
    deallocate(H,E,C,matrix_upper_triangle); deallocate(work,iwork,ifail)

    return

ENDSUBROUTINE

