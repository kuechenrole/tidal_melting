      subroutine tidpot_load
!  this program will calculate the combined equilibrium +
!  load +self-attraction tide for inclusion into my revision of fundy5sp
!  The load +self attraction tide is taken from Ray's solution
!  while the equilibrium tide is calculated
!  for specific constituents following the formulae
!  given by Schwiderski (NWSC, 1978, pages 22 and 3),
!  and cross-checked with Flather and Cartwright
!  See also my 1993 JGR paper.
!  The loadtide values are every half degree and must be
!  interpolated to the node locations
!  NB: If using only the load tide (ie no self attr), then must
!  switch sign to get total result. See atot
!
      dimension v0(8),nload(8)
      dimension xlat(418000),xlon(418000),amp(8),                       &
     &          ph(8),xk(8),xh(8),hn(8),amphload(720,361,8,2),          &
     &          xlatload(720,361), xlonload(720,361),                   &
     &          ampeq(418000,8),pheq(418000,8)
      complex a1,a2,a3,a4,a5,a6,aint,aeq,atot(418000,8),eye
      character*80 fname
!  The precise astronomical argument, v0, values chosen here are
!  for 000GMT August 1, 1984. This corresponds to the approximate
!  time when current observations were taken at sites W01 and VJ3
!  near Cape Flattery.
!     data v0/49.47,48.92,140.21,219.79,269.26,268.71,0.00,1259.58/
!  these v0 values are set to zero
      data v0/8*0./
      eye=cmplx(0.,1.)
      fac=3.1415926536/180.
      twpi=2.*3.1415926536
      open(unit=3,file='run_tidpot_load_sal.inp',status='old',          &
     &     form='formatted')
      read(3,'(a)') fname
      open(unit=7,file=fname,status='old',form='formatted')
      read(3,'(a)') fname
      open(unit=8,file=fname,status='unknown',form='formatted')
      read(3,'(a)') fname
      open(unit=10,file=fname,status='unknown',form='formatted')
      read(3,*) npts
      do i=1,8
        read(3,*) xk(i),xh(i),hn(i)
      end do
!  we assume the constituents are q1,o1,p1,k1,n2,m2,s2,k2,
!  in that order read in the load + self =attr tide files
!  longitudes are expected to be east
        do 20 i=1,8
          read(3,'(a)') fname
          open(unit=9,file=fname,status='old',form='formatted')
!  amplitudes: skip 1st 7 lines
          do k=1,7
            read(9,'(a)') fname
          end do
          do k=1,361
            do l=1,720
              xlatload(l,k)=(k-1)*0.5
              xlonload(l,k)=(l-1)*0.5
            end do ! l
            read(9,200) (amphload(l,k,i,1),l=1,720)
!  convert amplitudes from mm to m, and phases from degrees to
!  radians
            do l=1,720
              amphload(l,k,i,1)=amphload(l,k,i,1)*.001
            end do
          end do
 200      format(11f7.2)
!  phases: skip 1st 7 lines
          do k=1,7
            read(9,'(a)') fname
          end do
          do k=1,361
            read(9,200) (amphload(l,k,i,2),l=1,720)
            do l=1,720
              amphload(l,k,i,2)=amphload(l,k,i,2)*fac
            end do
          end do
          close(unit=9)
  20    continue
!       nload=361*720
!
!  read in the lat/lons where we want the potential & load/sal values
!  ngh file
!  read(7,*)
!  read(7,*)
        read(7,*) fname
        do i=1,npts
          read(7,*) xlon(i),xlat(i)
        end do
!
!        Hedstrom & Curchitser bering ROMS grid
!        for the arctic grid file we must add 360 to the longitudes
!        xlon(i)=xlon(i)+360.
!        input longitudes must be east
       read(7,'(a)') fname
       do i=1,770
         kk=542*(i-1)
         read(7,*) (xlat(kk+j),j=1,542)
       end do
       read(7,'(a)') fname
       read(7,'(a)') fname
       do i=1,770
         kk=542*(i-1)
         read(7,*) (xlon(kk+j),j=1,542)
       end do
        close(unit=7)
!
        do k=1,8
          do i=1,npts
            arg=xlat(i)*fac
            carg=cos(arg)
            s2arg=sin(2.*arg)
!  calculate the 4 surrounding lat/lons to use for the
!  loadtide interpolation
!  (l1,k1) and (l2,k2) sw & ne points
            j1=int((xlat(i)+90.)/0.5)+1
            j2=j1+1
            l1=int(xlon(i)/0.5)+1
            if(xlon(i).lt.359.5) then
              l2=l1+1
            else
              l2=1
            end if
            if (i.eq.1) then
              write(6,*) xlat(i),xlon(i),l1,l2,j1,j2
            end if
!        interpolate for the load tides
!        note that their phases are lags        (Ray, 1998, Fig 1)
            a1=amphload(l1,j1,k,1)*cexp(-eye*amphload(l1,j1,k,2))
            a2=amphload(l2,j1,k,1)*cexp(-eye*amphload(l2,j1,k,2))
            a3=amphload(l1,j2,k,1)*cexp(-eye*amphload(l1,j2,k,2))
            a4=amphload(l2,j2,k,1)*cexp(-eye*amphload(l2,j2,k,2))
            theta=((xlat(i)+90.)-0.5*(j1-1))/.5
            phi=(xlon(i)-0.5*(l1-1))/.5
            a5=phi*a2+(1.-phi)*a1
            a6=phi*a4+(1.-phi)*a3
            aint=theta*a6+(1.-theta)*a5
!        Roy's programs assume the tidal argument is
!        exp(-i(omega*t+vo-g)) so we have make our
!        phases negative for that program
!        But fundy5 assumes exp(i(omega*t+vo-g))
            if(k.le.4) then
              amp(k)=hn(k)*s2arg
!        tide3d
!        ph(k)=-(xlon(i)+v0(k))*fac
!        fundy5
              ph(k)=(xlon(i)+v0(k))*fac
            else
              amp(k)=hn(k)*carg*carg
!        tide3d
!        ph(k)=-(2.*xlon(i)+v0(k))*fac
!        fundy5
              ph(k)=(2.*xlon(i)+v0(k))*fac
            end if
            aeq=(1.+xk(k)-xh(k))*amp(k)*cexp(eye*ph(k))
!        as we are calculating elevation-load-potential
!        we simply add here
!        see Ray (1998) equn 6
!        fundy5
            atot(i,k)=aeq+aint
!        tide3d: switch SAL phases to be consistent with tidpot
!        atot(i,k)=aeq+conjg(aint)
!        for load tide only we subtract because Ray's loadtide
!        values already include a 180 shift from the ocean tides
!        atot(i,k)=aeq-aint
            if(k.eq.6.and.i.eq.1) then
              write(6,32) aeq,aint,atot(i,k)
32            format(' aeq,aint,atot for m2 at node 1 =',6f10.5)
            end if
!        may also want to create a file of zero values
!        for nontidal calculations
!        atot(i,k)=0.
!        or we may just want the tidal potential and earth tides
!        atot(i,k)=aeq
!
            ampeq(i,k)=cabs(atot(i,k))
            x=real(atot(i,k))
            y=aimag(atot(i,k))
            if(ampeq(i,k).gt.0.0001) then
              pheq(i,k)=atan2(y,x)/fac
            else
              pheq(i,k)=0.
            end if
            if(pheq(i,k).lt.0.) pheq(i,k)=pheq(i,k)+360.
          end do
        end do
        do i=1,npts
!         write(6,*) i
!         write(8,12) i,(atot(i,k),k=1,8)
!         write(8,12) i,xlon(i),xlat(i),(ampeq(i,k),pheq(i,k),k=1,8)
! 12      format(i5,2f10.4/5x,4(f10.5,f10.2)/5x,4(f10.5,f10.2))
          write(8,*) (ampeq(i,k),pheq(i,k),k=1,8)
        end do
        write(10,30) (ampeq(i,6),i=1,npts)
        write(10,31) (pheq(i,6),i=1,npts)
  30    format(15f8.5)
  31    format(15f8.2)
        return
        end
