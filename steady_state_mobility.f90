! THIS SUBROUTINE SOLVES THE STEADY-STATE MASTER EQUATION ITERATIVELY AND CALCULATES THE MOBILITY

! WE ARE TRYING TO FIND NON-TRIVIAL SOLUTION OF AN EQUATION OF TYPE AX=0.

! A NON-TRIVIAL SOLUTION EXISTS IF AND ONLY IF THE MATRIX HAS AT LEAST ONE FREE VARIABLE IN ROW-ECHLEON FORM,  
! I.E. AT LEAST ONE ROW OF THE MATRIX IN THE ROW-ECHLEON FORM HAS ONLY ZEROS AS ELEMENTS.

! SINCE RANK OF A IS LOWER THAN N=DIM(A), WE ASSUME THAT IT IS N-1, AND ONE OF THE VARIABLES CAN TAKE ARBITRARY VALUE. 

! WE SET THE VALUE OF THIS VARIABLE, AND CONSEQUENTLY THE SYSTEM OF N HOMOGENEOUS EQUATIONS GETS CONVERTED
! INTO A SYSTEM OF N-1 INHOMOGENEOUS EQUATION. WE ASSUME THAT THE LATTER HAS NON-ZERO DETERMINANT AND AN 
! UNIQUE SOLUTION. THIS ASSUMPTION IS FURTHER VERFIED DURING THE LU DECOMPOSITION AND THE REST OF THE 
! FREE VARIABLES ARE SET TO ZERO.

!!!  REFERENCE: http://math.bu.edu/people/mkon/ma242/L3.pdf   !!!


SUBROUTINE main(ichain,n_site,sru_length,field,carrier,obc_or_pbc,folder_path,mobility,index)

    use datatype
    implicit none
    integer(kind=I4), intent(in) :: ichain,n_site
    real(kind=DP), intent(in) :: sru_length,field
    character(len=20), intent(in) :: obc_or_pbc,carrier
    character(len=1000), intent(in) :: folder_path
    real(kind=DP), intent(out) :: mobility
    integer(kind=I4), intent(out) :: index

    real(kind=DP) :: field_eva,condition_number,tolerance,pbc_box_length,norm,particle_velocity,delta_r,junk
    real(kind=DP), allocatable :: k_hopping(:,:),rate_matrix(:,:),coeff_matrix(:,:),x(:),occ_prob(:),coord(:)
    real(kind=DP), allocatable :: A(:,:),b(:)
    integer(kind=I4), allocatable :: ipiv(:)
    integer(kind=I4) :: ii,jj,istart,iflag,ndim,info,ik,jk
    character(len=20) :: str
    character(len=2000) :: filename


    field_eva=field*1.0d-8            !!!  conversion from ev/cm to ev/a units  !!!
    condition_number=1.0d-16
    pbc_box_length=n_site*sru_length
    allocate(k_hopping(1:n_site,1:n_site),rate_matrix(1:n_site,1:n_site),x(1:n_site),occ_prob(1:n_site),coord(1:n_site))
    rate_matrix=0.0d0
    x=0.0d0

    write(str,'(I0)') ichain
    filename=trim(folder_path) // 'localized_state_hopping_rate_' // trim(str) //'_unformatted.out'
    open(unit=101,file=filename,form='unformatted')
    read(101)((k_hopping(ii,jj),jj=1,n_site),ii=1,n_site)
    close(101)

    do ii=1,n_site
        do jj=1,n_site
            rate_matrix(ii,ii)=rate_matrix(ii,ii)-k_hopping(ii,jj)
        enddo
        do jj=1,n_site
            if(jj /= ii) then
                rate_matrix(ii,jj)=k_hopping(jj,ii)
            endif
        enddo
    enddo

    tolerance=maxval(abs(rate_matrix))*condition_number

    do ii=1,n_site
        do jj=1,n_site
            if (abs(rate_matrix(ii,jj)) <= tolerance) rate_matrix(ii,jj)=0.0d0
        enddo
    enddo

    istart=1
    do while (istart <= n_site)
        iflag=0
        do ii=istart,n_site
            if (abs(rate_matrix(ii,istart)) >= tolerance) then
                iflag=1
                exit
            endif
        enddo
        if (iflag == 1) exit
        istart = istart + 1
    enddo

    ndim=n_site-istart
    if(ndim <= 0) then
        write(*,*)'error - ndim dimension is wrong'
        index=1
        return
    endif
    allocate(A(1:ndim,1:ndim),b(1:ndim),ipiv(1:ndim))

    do ii=1,ndim
        do jj=1,ndim
            A(ii,jj)=rate_matrix(ii+istart,jj+istart)
        enddo
        b(ii)=-rate_matrix(ii+istart,istart)
    enddo

    call dgesv(ndim,1,A,ndim,ipiv,b,ndim,info)

    if(info /= 0) then
        write(*,*)'Error - occupational probabilities are not successfully calculated'
        index=1
        return
    endif

    if (carrier == 'e') then
        x(istart)=1.0d0                 !!!  value of a particular variable gets pre-fixed  !!!
        do ii=1,ndim
            x(istart+ii)=b(ii)
        enddo
    elseif (carrier == 'h') then
        x(n_site+1-istart)=1.0d0        !!!  value of a particular variable gets pre-fixed  !!!
        do ii=1,ndim
            x(n_site+1-istart-ii)=b(ii)
        enddo
    endif

    norm=0.0d0
    do ii=1,n_site
        norm=norm+x(ii)
    enddo
    occ_prob=x/norm

    filename=trim(folder_path) // 'localized_state_ll_' // trim(str) //'.out'
    open(unit=101,file=filename)
    do ii=1,n_site
        read(101,*)coord(ii),junk
    enddo

    particle_velocity=0.0d0
    do ii=1,n_site-1
        do jj=ii+1,n_site
            if (carrier == 'e') then
                delta_r=coord(jj)-coord(ii)
                if(obc_or_pbc =='pbc') then
                    if(delta_r >= 0.5*pbc_box_length) then
                        delta_r=delta_r-pbc_box_length
                    elseif(delta_r < (-0.5*pbc_box_length)) then
                        delta_r=delta_r+pbc_box_length
                    endif
                endif
                particle_velocity=particle_velocity+(occ_prob(ii)*k_hopping(ii,jj)*delta_r)-(occ_prob(jj)*k_hopping(jj,ii)*delta_r)
            elseif(carrier == 'h') then
                ik=n_site+1-ii
                jk=n_site+1-jj
                delta_r=coord(jk)-coord(ik)
                if(obc_or_pbc =='pbc') then
                    if(delta_r >= 0.5*pbc_box_length) then
                        delta_r=delta_r-pbc_box_length
                    elseif(delta_r < (-0.5*pbc_box_length)) then
                        delta_r=delta_r+pbc_box_length
                    endif
                endif
                particle_velocity=particle_velocity+(occ_prob(ik)*k_hopping(ii,jj)*delta_r)-(occ_prob(jk)*k_hopping(jj,ii)*delta_r)
            endif
        enddo
    enddo

    mobility=particle_velocity/field_eva
    mobility=dabs(mobility)*1.0d-16            !!! unit conversion from angstrom^2/(v.s) to cm^2/(v.s)
    write(*,'(A8,F20.10)')'Mobility',mobility

    filename=trim(folder_path) // 'localized_state_occupation_probability_' // trim(str) //'.out'
    open(unit=99,file=filename)
    write(99,'(*(F20.10,3x))')(occ_prob(ii),ii=1,n_site)
    close(99)

    deallocate(k_hopping,rate_matrix,x,occ_prob,coord,A,b,ipiv)

    index=0

    return

END SUBROUTINE
