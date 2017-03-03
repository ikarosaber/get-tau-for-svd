        
	
	parameter (iyst = 1982, iyend = 2014, mmax0 = 5 )   ! zhang: need to be changed
	parameter (nx=134 ,ny=61 )
        parameter (m1=4425,m2=4425,np =2 )
	parameter (imt = nx, jmt = ny, km = 1 , mmax = 12)
	parameter (rmiss = -1.e34)

	dimension sst(nx,ny),u(nx,ny),odata(nx,ny),ssta(imt,jmt)
	dimension w(nx*ny),neign(np)
	dimension lsi(nx,ny),  index0(nx,ny)

! zhang
        common /eigen/ist1(m1),ist2(m2),arr(m1,20,np)
     &,         brr(m2,20,np),stdv(20,np),sd,ud(np)
	common /eigen12/arr12(m1,20,np,12),brr12(m2,20,np,12)
     &                 ,  stdv12(20,np,12),sd12(12),ud12(np,12)

! zhang
!       data neign/3,3/
        data neign/5,5/
!       data neign/10,10/
!       data neign/16,17/

	character fin*68, fout*68,fout1*68,fout2*68, fclim*68
    
 	fout1= 'taux_anom.jan1982-may2014.svdmonth.dta'
 	fout2= 'tauy_anom.jan1982-may2014.svdmonth.dta'

	open(15,file='ssta.nov1981-may2014.kl.dta', 
     &  form='unformatted',access='direct',recl=imt*jmt*1 )

	irec_len1= imt*jmt*1*1

	read(9,88) lsi
 88     format(60i1)

        do  i=1, nx
	do  j=1, ny
	index0(i,j) = lsi(i,j)
	end do
	end do

       ij0=0
       do  j=1, ny
       do  i=1, nx
       if(lsi(i,j) .eq. 1) then
         ij0=ij0+1
	 ist1(ij0)=i+(j-1)*nx
       endif
       enddo
       enddo

       ic=0
       ij0=0
       do  j=1, ny
       do  i=1, nx
       if(index0(i,j) .eq. 1) then
         ic=ic+1
	 ij0=ij0+1
	 ist2(ij0)=i+(j-1)*nx
       endif
       enddo
       enddo



! read in monthly svd statistics

	open(97,file='fort.97-month',form='unformatted')
	read(97) arr12,brr12,stdv12,sd12,ud12
	close(97)
	
	write(*,*)"sd12", sd12
        write(*,*) sd
        write(*,*) ud

	do mode=1, 20
	write(6,567) mode, (stdv(mode,np0), np0=1,np )
	end do
  567   format(1x,i2,1x,6(f9.2,1x))

	open(21, file=fout1, form='unformatted',
     &        access='direct', recl=irec_len1)
	open(22, file=fout2, form='unformatted',
     &        access='direct', recl=irec_len1)
 
c-------------------------------------------------------------
  
	do iy=iyst, iyend 

	  mmax00=mmax

	  if(iy.eq.2014) mmax00=mmax0    

	do m = 1, mmax00

	  nrec1 = m + (iy-1982)*12   + 2     ! begining from nov 1981
	  read(15,rec=nrec1) ((odata(i,j),i=1,imt),j=1,jmt)
	  kk=m+(iy-iyst)*12
 	  write(*,*) iy,m,nrec1,kk

        do 10 i=1,nx
        do 10 j=1,ny
10      sst(i,j)=odata(i,j)

!       read in number of SVD modes used to construct each atmospheric variable

        do 30 l=1,2                    

        ne=neign(l)

	do mode=1,20
          do ii=1, m1
          arr(ii,mode,l)=arr12(ii,mode,l,m)
  	  end do
          do ii=1, m2
          brr(ii,mode,l)=brr12(ii,mode,l,m)
	  end do
  
	  stdv(mode,l)  =stdv12 ( mode,l,m)
	  sd            =sd12   (        m)
	  ud  (     l)  =ud12   (      l,m)
	end do                              !  mode=1,20

        call atmos(sst,u,ne,l)

	do j = 1, jmt
	do i = 1, imt
	ssta(i,j) = 0.0
!	ssta(i,j) = rmiss
	enddo
	enddo

	do j = 1, jmt
	do i = 1, imt
	if(index0(i,j) .eq. 1 ) then
	  ssta(i,j) = u(i,j)
	endif
	enddo
	enddo

 	if(l.eq.1) write(21,rec=kk) ssta
 	if(l.eq.2) write(22,rec=kk) ssta

!       write(6,289) iy,m,(u(i,15),i=15,20)

 30      continue 

	enddo                   ! m
	enddo                   ! iy
	close(21)
	close(22)
 
 289    format(1x,i4,1x,i2,1x,6(f9.4,1x) )
        end
!----------------------------------
!       atmospheric model subroutine

	subroutine atmos(sst,u,neign,n)
        parameter (nx=134 ,ny=61  )
        parameter (m1=4425,m2=4425,np =2 ,m3=nx*ny )
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
22      work=work+arr(i,k,n)*s(i)
20      af(k)=work/stdv(k,n)/stdv(k,n)

        do 30 i=1,m2
        work=0.0
        do 40 l=1,neign
         work=work+af(l)*brr(i,l,n) 
40      continue
        out(i)=work
30      continue
    
        do i=1,nx*ny

!       w(i)= -1.e9 for land points.
!       w(i) can be set to zero in the coupled model
        w(i)= 0.0    
!       w(i)= -1.e34 
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
