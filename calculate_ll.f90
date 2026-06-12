! THIS SUBROUTINE CALCULATES THE LOCALIZATION LENGTH OF THE LOCALIZED EIGENSTATES !
! RETURNS A TWO-COMPONENT VECTOR !
! THE FIRST REPRESNTS THE CENTROID POSITION, SECOND THE SPREAD OF THE WAVEFUNCTION. !

! WITHIN PBC, VERIFY WHETHER THE STATE IS CENTERED ABOUT SITE 1, MIDDLE OF THE CHAIN OR ANY OTHER SITE !


SUBROUTINE main(ichain,n_site,sru_length,obc_or_pbc,folder_path,index)

    use datatype
    implicit none
    integer(kind=I4), intent(in) :: ichain,n_site
    real(kind=DP), intent(in) :: sru_length
    character(len=20), intent(in) :: obc_or_pbc
    character(len=1000), intent(in) :: folder_path
    integer(kind=I4), intent(out) :: index

    integer(kind=I4) :: ii,jj,ksite,kflag,lsite,chain_begin,kk,mm,caseflag,nn,lflag
    integer(kind=I4), allocatable :: sorted_index(:)
    real(kind=DP), allocatable :: C(:,:),ll(:,:),coord(:),coord_unwrapped(:),charge_density(:)
    real(kind=DP) :: pbc_box_length,r,r_sqr,delta_r,tolerance,pbc_origin_shift
    character(len=20) :: str
    character(len=2000) :: filename


    allocate(C(1:n_site,1:n_site),coord(1:n_site),ll(1:n_site,1:2))

! Generate the site coordinates in accordance to OBC or PBC

    if(obc_or_pbc == 'obc') then
        coord=[(sru_length*(ii-1), ii=1,n_site)]
    elseif(obc_or_pbc == 'pbc') then
        pbc_box_length=n_site*sru_length
        coord=[(sru_length*(ii-1), ii=1,n_site)]
        do ii=int(n_site/2)+1,n_site
            coord(ii)=coord(ii)-pbc_box_length
        enddo
    endif

    write(str,'(I0)') ichain
    filename=trim(folder_path) // 'localized_state_coefficient_' // trim(str) //'.out'
    open(unit=101,file=filename)
    do ii=1,n_site
        read(101,*)(C(jj,ii),jj=1,n_site)
    enddo
    close(101)

    if(obc_or_pbc == 'obc') then

        filename=trim(folder_path) // 'localized_state_ll_' // trim(str) //'.out'
        open(unit=99,file=filename)
        ll=0.0d0

        do ii=1,n_site

            r=0.0d0
            r_sqr=0.0d0
      
            do jj=1,n_site
                r=r+coord(jj)*(C(jj,ii)**2)
            enddo
            ll(ii,1)=r

            do jj=1,n_site
                delta_r=coord(jj)-r
                r_sqr=r_sqr+(delta_r**2)*(C(jj,ii)**2)
            enddo
  
            if(r_sqr.ge.1.0d-10) then
                ll(ii,2)=2.0d0*sqrt(r_sqr)
            endif

            write(99,'(*(F16.10,3x))')(ll(ii,jj),jj=1,2)

        enddo

        close(99)
        index=0
        deallocate(C,coord,ll)

    elseif(obc_or_pbc == 'pbc') then

! Verification whether the state is localized middle of the chain or towards boundary. It has significant 
! component at chain end sites if coefficient**2 > tolerance. It will be checked over a small length of 
! the polymer chain given by ksite
!
! kflag=0: localized state is centered about zero (site 1); insignificant component at sites near n_site/2 
! kflag=1: localized state is centered about n_site/2; insignificant component at sites 1 and n_site
! kflag=2: localized state is centered about other site in the chain. Search for the site where charge 
! density is at its minimum along with smaller tail that can be selected as the new chain boundary. 
!

        filename=trim(folder_path) // 'localized_state_ll_' // trim(str) //'.out'
        open(unit=99,file=filename)

        tolerance=1.0d0/n_site
        ksite=int(0.05*n_site)
        ll=0.0d0

        do ii=1,n_site

            kflag=0

            do jj=(int(n_site/2)-ksite),(int(n_site/2)+ksite)
                if(abs(C(jj,ii)**2).ge.tolerance) then
                    kflag=1
                    exit
                endif
            enddo

            do kk=1,ksite
                if(abs(C(kk,ii)**2).ge.tolerance) then
                    kflag=2
                    exit
                endif
            enddo

            if(kflag.eq.1) then
                do kk=n_site-ksite,n_site
                    if(abs(C(kk,ii)**2).ge.tolerance) then
                        kflag=2
                        exit
                    endif
                enddo
            endif


            if(kflag.eq.0) then

                r=0.0d0
                r_sqr=0.0d0
      
                do jj=1,n_site
                    r=r+coord(jj)*(C(jj,ii)**2)
                enddo
                ll(ii,1)=r
                if(ll(ii,1).lt.0.0d0) ll(ii,1)=ll(ii,1)+pbc_box_length

                do jj=1,n_site
                    delta_r=coord(jj)-r
                    if(delta_r.ge.0.5*pbc_box_length) then
                        delta_r=delta_r-pbc_box_length
                    elseif(delta_r.lt.(-0.5*pbc_box_length)) then
                        delta_r=delta_r+pbc_box_length
                    endif
                    r_sqr=r_sqr+(delta_r**2)*(C(jj,ii)**2)
                enddo
  
                if(r_sqr.ge.1.0d-10) then
                    ll(ii,2)=2.0d0*sqrt(r_sqr)
                endif

                write(99,'(*(F16.10,3x))')(ll(ii,jj),jj=1,2)

                cycle

            endif

            if(kflag.eq.1) then

! Coordinates are unwrapped and the origin is shifted at the middle of the chain

                allocate(coord_unwrapped(1:n_site))
                r=0.0d0
                r_sqr=0.0d0
                coord_unwrapped=0.0d0
                pbc_origin_shift=0.5*pbc_box_length

                do jj=1,n_site
                    if(coord(jj).lt.0.0d0) then
                        coord_unwrapped(jj)=coord(jj)+pbc_box_length
                    else
                        coord_unwrapped(jj)=coord(jj)
                    endif
                enddo

                do jj=1,n_site
                    r=r+(coord_unwrapped(jj)-pbc_origin_shift)*(C(jj,ii)**2)
                enddo
                r=r+pbc_origin_shift
                ll(ii,1)=r
         
                do jj=1,n_site
                    delta_r=coord_unwrapped(jj)-r
                    if(delta_r.ge.0.5*pbc_box_length) then
                        delta_r=delta_r-pbc_box_length
                    elseif(delta_r.lt.(-0.5*pbc_box_length)) then
                        delta_r=delta_r+pbc_box_length
                    endif
                    r_sqr=r_sqr+(delta_r**2)*(C(jj,ii)**2)
                enddo

                if(r_sqr.ge.1.0d-10) then
                    ll(ii,2)=2.0d0*sqrt(r_sqr)
                endif

                write(99,'(*(F16.10,3x))')(ll(ii,jj),jj=1,2)

                deallocate(coord_unwrapped)

                cycle

            endif
        
            if(kflag.eq.2) then

! Coordinates are unwrapped and the origin is shifted to the site with highest charge density

                allocate(coord_unwrapped(1:n_site),sorted_index(1:n_site),charge_density(1:n_site))
                r=0.0d0
                r_sqr=0.0d0
                coord_unwrapped=0.0d0
                lsite=0
                sorted_index=0
                charge_density=0.0d0
                chain_begin=1

                do jj=1,n_site
                    charge_density(jj)=(abs(C(jj,ii)))**2
                    sorted_index(jj)=jj
                enddo

                call quicksort(charge_density,sorted_index,chain_begin,n_site)

                outer: do mm=1,n_site
                            kk=sorted_index(mm)
                            if(kk.ge.(int(n_site/2)+1).and.kk.le.(n_site-1)) caseflag=1
                            if(kk.ge.1.and.kk.le.(int(n_site/2)-1)) caseflag=2

                            select case(caseflag)

                                case(1)

                                    lflag=0
                                    do jj=kk-ksite,kk+ksite
                                        nn=jj
                                        if(nn.gt.n_site) nn=nn-n_site
                                        if(nn.le.0) nn=nn+n_site
                                        if(abs(C(nn,ii)**2).ge.tolerance) then
                                            lflag=1
                                            exit
                                        endif
                                    enddo

                                    if(lflag.eq.0) then
                                        lsite=kk
                                    endif

                                    if(lsite.ne.0) then
                                        pbc_origin_shift=coord(lsite)+0.5*pbc_box_length
                                        do kk=1,n_site
                                            if(kk.ge.(int(n_site/2)+1).and.kk.lt.lsite) then
                                                coord_unwrapped(kk)=coord(kk)+pbc_box_length
                                            else
                                                coord_unwrapped(kk)=coord(kk)
                                            endif
                                        enddo
                                        exit outer
                                    endif


                                case(2)

                                    lflag=0
                                    do jj=kk-ksite,kk+ksite
                                        nn=jj
                                        if(nn.gt.n_site) nn=nn-n_site
                                        if(nn.le.0) nn=nn+n_site
                                        if(abs(C(nn,ii)**2).ge.tolerance) then
                                            lflag=1
                                            exit
                                        endif
                                    enddo

                                    if(lflag.eq.0) then
                                        lsite=kk
                                    endif

                                    if(lsite.ne.0) then
                                        pbc_origin_shift=coord(lsite)-0.5*pbc_box_length
                                        do kk=1,n_site
                                            if(kk.gt.lsite.and.kk.le.int(n_site/2)) then
                                                coord_unwrapped(kk)=coord(kk)-pbc_box_length
                                            else
                                                coord_unwrapped(kk)=coord(kk)
                                            endif
                                        enddo
                                        exit outer
                                    endif

                            end select 

                enddo outer 
       
                if(lsite.eq.0) then
                    pbc_origin_shift=0.0d0
                    coord_unwrapped=coord
                endif

                do kk=1,n_site
                    r=r+(coord_unwrapped(kk)-pbc_origin_shift)*(C(kk,ii)**2)
                enddo
                r=r+pbc_origin_shift
                ll(ii,1)=r
                if(ll(ii,1).lt.0.0d0) ll(ii,1)=ll(ii,1)+pbc_box_length
         
                do kk=1,n_site
                    delta_r=coord_unwrapped(kk)-r
                    if(delta_r.ge.0.5*pbc_box_length) then
                        delta_r=delta_r-pbc_box_length
                    elseif(delta_r.lt.(-0.5*pbc_box_length)) then
                        delta_r=delta_r+pbc_box_length
                    endif
                    r_sqr=r_sqr+(delta_r**2)*(C(kk,ii)**2)
                enddo

                if(r_sqr.ge.1.0d-10) then
                    ll(ii,2)=2.0d0*sqrt(r_sqr)
                endif

                write(99,'(*(F16.10,3x))')(ll(ii,jj),jj=1,2)

                deallocate(coord_unwrapped,sorted_index,charge_density)

            endif

        enddo

        close(99)
        index=0
        deallocate(C,coord,ll)

    endif

    return

END SUBROUTINE

RECURSIVE SUBROUTINE quicksort(a,o,first,last)

    use datatype
    implicit none
    integer(kind=I4) :: first,last,i,j,itemp
    integer(kind=I4), intent(inout) :: o(*)
    real(kind=DP), intent(inout) :: a(*)
    real(kind=DP) :: pivot,temp

    if(first.ge.last) then
        return
    endif

    i=first
    j=last
    pivot=a(int((first+last)/2))
  
    do

        do while (i.le.last.and.a(i).lt.pivot)
            i=i+1
        enddo 

        do while (j.ge.first.and.a(j).gt.pivot)
            j=j-1
        enddo 

        if(i.ge.j) exit

        temp=a(i)
        a(i)=a(j)
        a(j)=temp
        itemp=o(i)
        o(i)=o(j)
        o(j)=itemp

        i=i+1
        j=j-1

    enddo 
  
    if(first.lt.(i-1)) call quicksort(a,o,first,(i-1))
    if((j+1).lt.last) call quicksort(a,o,(j+1),last)

    return

END SUBROUTINE

