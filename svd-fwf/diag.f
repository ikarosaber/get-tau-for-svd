
        parameter (nx=134 ,ny=61  )
        parameter (m1=4425,m2=4425,np =1 )
	parameter (imt = nx, jmt = ny, km = 1 , mmax = 12)

        parameter (iyst = 1979, iyend = 2008)
	parameter (rmiss = -1.e34)

	dimension mask(nx,ny)

	dimension sst(nx,ny),u(nx,ny),odata(nx,ny),w(nx*ny)
	dimension neign(np),lsi(nx,ny),  index0(nx,ny)

        common /eigen/ist1(m1),ist2(m2),arr(m1,20,np) 
     &,         brr(m2,20,np),stdv(20,np),sd,ud(np)

! zhang
        data neign/5/

	real   data(imt,jmt,km), datam(imt,jmt,km,mmax)
	real   ssta(imt,jmt)
	character fin*68, fout*68, fclim*68
    
 	fout = 'recon_sflx_anom.1979-2008.dta'

	irec_len1= imt*jmt*1*1

!       read(9,88) lsi
 88     format(60i1) 

      open(3,file='data/mask-sst.dat',   
     &       form='formatted')
      read(3,*) ((mask(i,j),i=1,nx),j=1,ny)
      close(3)

!----------------------------
      do i=1,nx
        do j=1,ny
        lsi(i,j)=int(mask(i,j))
        enddo
        enddo


 
        do  i=1, nx
        do  j=1, ny
        index0(i,j) = lsi(i,j)     
        end do
        end do

      open(15,file='data/ersst_anom.1979-2008.mgrid.dat',  
     &     form='unformatted',access='direct',recl=imt*jmt)


!       read in SVD statistics
!
c  	open(96,file='fort.86-sflx',form='unformatted')
   	open(96,file='fort.87',form='unformatted')
        read(96) ist1,ist2,arr,brr,stdv,sd,ud  ! ann case
        close(96)

	write(*,*) ist1(1),ist2(1),arr(1,1,1),ud(1)

        write(*,*) sd
        write(*,*) ud

	do mode=1, 20
	write(6,567) mode, (stdv(mode,np0), np0=1,np )
	end do
  567   format(1x,i2,1x,6(f9.2,1x))
!
!-------------------------------------------------------------
  
	open(20, file=fout, form='unformatted',  
     &        access='direct', recl=irec_len1)
  
  
	do iy=iyst, iyend            ! 1950- 1997
!	write(fout(16:21),'(i6.6)') iy

	print*,'anomaly iy=',iy, fout
	
	do m = 1, mmax 

	  nrec1 = m + (iy-iyst)*12
	  read(15,rec=nrec1) ((odata(i,j),i=1,imt),j=1,jmt)

	  kk    = m + (iy-iyst)*12
  	  write(*,*) iy,m,nrec1,kk

        do 10 i=1,nx
        do 10 j=1,ny
10      sst(i,j)=odata(i,j)
!c
!c       read in number of SVD modes used to construct each atmospheric variable
!c
        do 20 l=1,1                    ! tox
!c       do 20 l=2,2                    ! toy

        ne=neign(l)
!c
!c       call statistical atmosphere
!c
        call atmos(sst,u,ne,l)

	do j = 1, jmt
	do i = 1, imt
	ssta(i,j) = rmiss
	enddo
	enddo

	do j = 1, jmt
	do i = 1, imt
	if(index0(i,j) .eq. 1 ) then
	  ssta(i,j) = u(i,j)
	endif
	enddo
	enddo

!c	write(20,rec=m) ssta
 	write(20,rec=kk) ssta
        write(6,289) iy,m,(u(i,15),i=15,20)

20      continue 

	enddo
	enddo
	close(20)
!c
 289    format(1x,i4,1x,i2,1x,6(f9.4,1x) )
        end
!c
!c..............................................................................
!c
!c       atmospheric model subroutine
!c
!c..............................................................................
!c
	subroutine atmos(sst,u,neign,n)
        parameter(nx=134,ny=61)
        parameter (m1=4425,m2=4425,np =1,m3=nx*ny )
        dimension sst(nx,ny),u(nx,ny)
        dimension w(m3),s(m1),out(m2),af(40)
        common /eigen/ist1(m1),ist2(m2),arr(m1,20,np) 
     &,         brr(m2,20,np),stdv(20,np),sd,ud(np)

        ind=0
        do 10 j=1,ny
        do 10 i=1,nx
        ind=ind+1
        w(ind)=sst(i,j)
10      continue

        do i=1,m1
        s(i)=w(ist1(i))/sd
        enddo

        do 20 k=1,neign
        work=0.0
        do 22 i=1,m1
 22     work=work+arr(i,k,n)*s(i)
        af(k)=work/stdv(k,n)/stdv(k,n)
 20     continue

        do 30 i=1,m2
        work=0.0
        do 40 l=1,neign
         work=work+af(l)*brr(i,l,n) 
40      continue
        out(i)=work
30      continue
    
        do i=1,nx*ny
!c
!c       w(i)= -1.e9 for land points.
!c       w(i) can be set to zero in the coupled model
!c
        w(i)= 0.0    
!c       w(i)= -1.e7 
        enddo
        do i=1,m2
        w(ist2(i))=out(i)*ud(n)
        enddo
        ind=0
        do 150 j=1,ny
        do 150 i=1,nx
        ind=ind+1
150     u(i,j)=w(ind)
	return
	end
