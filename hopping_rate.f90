! THIS SUBROUTINE CALCULATES HOPPING RATES BETWEEN LOCALIZED STATES OF THE STATIC DISORDERED HAMILTONIAN !
! ENERGY IS IN ELECTRON-VOLT UNIT ! 

MODULE paras

    use datatype
    real(kind=DP), allocatable :: E(:),C(:,:),coord(:),k_hopping(:,:)

END MODULE

SUBROUTINE marcus(ichain,n_site,sru_length,sigma_alpha_dynamic,sigma_beta_dynamic,Lambda,T,field,carrier,obc_or_pbc,folder_path,index)

    use datatype
    use paras
    implicit none
    integer(kind=I4), intent(in) :: ichain,n_site
    real(kind=DP), intent(in) :: sru_length,sigma_alpha_dynamic,sigma_beta_dynamic,Lambda,T,field
    character(len=20), intent(in) :: carrier,obc_or_pbc
    character(len=1000), intent(in) :: folder_path
    integer(kind=I4), intent(out) :: index

    integer(kind=I4) :: ii,jj,kk
    real(kind=DP) :: pi,hbar,kB,prefactor,kBT,field_eva,pbc_box_length,delta_E,lambda_total,rho_fcwt,coupling_element,junk
    real(kind=DP), external :: gaussian_random_number
    real(kind=DP) :: abstol,dlamch
    real(kind=DP), allocatable :: S4(:)
    character(len=20) :: str
    character(len=2000) :: filename

! Assignment of machine floating point precision

    abstol=2*dlamch('S')

! Assignment of physical constants

    pi=4.0d0*atan(1.0d0)
    hbar=6.582119569d-16          !!!   ev.s   !!!
    kB=8.617333262d-5             !!!   ev/k   !!!
    prefactor=2.0d0*pi/hbar
    kBT=kB*T
    field_eva=field*1.0d-8            !!!  conversion from ev/cm to ev/a units  !!!

    allocate(E(1:n_site),C(1:n_site,1:n_site),coord(1:n_site),k_hopping(1:n_site,1:n_site))
    k_hopping=0.0d0

! Read the input file 

    write(str,'(I0)') ichain
    filename=trim(folder_path) // 'localized_state_coefficient_' // trim(str) //'.out'
    open(unit=101,file=filename)
    do ii=1,n_site
        read(101,*)(C(jj,ii),jj=1,n_site)
    enddo
    close(101)

    filename=trim(folder_path) // 'localized_state_energy_' // trim(str) //'.out'
    open(unit=101,file=filename)
    do ii=1,n_site
        read(101,*)E(ii)
    enddo
    close(101)

    filename=trim(folder_path) // 'localized_state_ll_' // trim(str) //'.out'
    open(unit=101,file=filename)
    do ii=1,n_site
        read(101,*)coord(ii),junk
    enddo

    if(obc_or_pbc == 'pbc') then
        pbc_box_length=n_site*sru_length
        do ii=1,n_site
            if(coord(ii).ge.0.5*pbc_box_length) then
                coord(ii)=coord(ii)-pbc_box_length
            endif
        enddo
    endif

! Precompute the per-state quartic sum S4(i) = sum_k C(k,i)**4.
! The reorganization energy separates as lambda_total(i,j) = (S4(i)+S4(j))*Lambda,
! so this replaces the O(n_site) inner loop in calculate_reorg_energy with an
! O(1) lookup, reducing the pair loop from O(n_site^3) to O(n_site^2) for this term.

    allocate(S4(1:n_site))
    do ii=1,n_site
        S4(ii)=0.0d0
        do kk=1,n_site
            S4(ii)=S4(ii)+C(kk,ii)**4
        enddo
    enddo


! Calculation of the hopping rates

    do ii=1,n_site
        do jj=1,n_site
            if(ii.ne.jj) then

! delta_E: energy gap between localized states in their lowest vibrational state, in presence of field 

                call calculate_delta_E(ii,jj,obc_or_pbc,pbc_box_length,carrier,field_eva,delta_E)

! Reorganization energy for charge transfer from intial to final states
! (precomputed separable form; equivalent to the original loop over kk)

                lambda_total=(S4(ii)+S4(jj))*Lambda

! Franck-condon factor and temperature weighted DOS

                rho_fcwt=dsqrt(1.0d0/(4*pi*lambda_total*kBT))*exp(-(((delta_E+lambda_total)**2)/(4.0d0*lambda_total*kBT)))

! Calculation of the coupling element between localized states 

                coupling_element=0.0d0

                if (dabs(sigma_alpha_dynamic) > 0.0d0) then
                    do kk=1,n_site
                        coupling_element=coupling_element+((sigma_alpha_dynamic)**2)&
                        &*((C(kk,ii)*C(kk,jj))**2)
                    enddo
                endif

                if (dabs(sigma_beta_dynamic) > 0.0d0) then
                    do kk=1,n_site-1
                        coupling_element=coupling_element+((sigma_beta_dynamic)**2)&
                        &*((C(kk,ii)*C(kk+1,jj)+C(kk+1,ii)*C(kk,jj))**2)
                    enddo
                    if(obc_or_pbc == 'pbc') then
                        coupling_element=coupling_element+((sigma_beta_dynamic)**2)&
                        &*((C(n_site,ii)*C(1,jj)+C(1,ii)*C(n_site,jj))**2)
                    endif  
                endif

! Calculation of hopping rate matrix element

                if (carrier == 'e') then
                    k_hopping(ii,jj)=prefactor*coupling_element*rho_fcwt
                    if(k_hopping(ii,jj).lt.0.0d0) then 
                        write(*,*)'Error-negative hopping rate'
                        index=1
                        return
                    if(k_hopping(ii,jj).lt.(2.0d0*abstol)) k_hopping(ii,jj)=0.0d0
                    endif
                elseif (carrier == 'h') then
                    k_hopping(n_site-ii+1,n_site-jj+1)=prefactor*coupling_element*rho_fcwt
                    if(k_hopping(n_site-ii+1,n_site-jj+1).lt.0.0d0) then 
                        write(*,*)'Error-negative hopping rate'
                        index=1
                        return
                    if(k_hopping(n_site-ii+1,n_site-jj+1).lt.(2.0d0*abstol)) k_hopping(n_site-ii+1,n_site-jj+1)=0.0d0
                    endif
                endif

            endif
        enddo
    enddo 

    filename=trim(folder_path) // 'localized_state_hopping_rate_' // trim(str) //'.out'
    open(unit=99,file=filename)
    write(99,'(*(E20.10E3,3x))')((k_hopping(ii,jj),jj=1,n_site),ii=1,n_site)
    close(99)

    filename=trim(folder_path) // 'localized_state_hopping_rate_' // trim(str) //'_unformatted.out'
    open(unit=99,file=filename,form='unformatted')
    write(99)((k_hopping(ii,jj),jj=1,n_site),ii=1,n_site)
    close(99)

    deallocate(E,C,coord,k_hopping,S4); index=0

    return

END SUBROUTINE

SUBROUTINE calculate_delta_E(ii,jj,obc_or_pbc,pbc_box_length,carrier,field,delta_E)

    use datatype
    use paras
    implicit none
    integer(kind=I4) :: ii,jj
    real(kind=DP) :: pbc_box_length,charge,r_ij,delta_E
    real(kind=DP),intent(in) :: field
    character(len=20),intent(in) :: carrier,obc_or_pbc

    r_ij=coord(jj)-coord(ii)

    if(obc_or_pbc =='pbc') then
        if(r_ij.ge.0.5*pbc_box_length) then
            r_ij=r_ij-pbc_box_length
        elseif(r_ij.lt.(-0.5*pbc_box_length)) then
            r_ij=r_ij+pbc_box_length
        endif
    endif

    delta_E=0.0d0
    if(carrier =='e') then
        charge=-1.0d0
        delta_E=E(jj)-E(ii)+charge*field*r_ij
    elseif(carrier == 'h') then
        charge=1.0d0
        delta_E=-(E(jj)-E(ii)+charge*field*r_ij)
    endif

   return

END SUBROUTINE

SUBROUTINE calculate_reorg_energy(ii,jj,n_site,Lambda,lambda_total)

   use datatype
   use paras
   implicit none
   integer(kind=I4), intent(in) :: n_site
   real(kind=DP), intent(in) :: Lambda
   integer(kind=I4) :: ii,jj,kk
   real(kind=DP) :: lambda_total

   lambda_total=0.0d0
   do kk=1,n_site
      lambda_total=lambda_total+(C(kk,ii)**4+C(kk,jj)**4)*Lambda
   enddo

   return

END SUBROUTINE