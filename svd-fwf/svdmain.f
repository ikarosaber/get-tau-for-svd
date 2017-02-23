C THIS PROGRAM IS FOR CALCULATING SVD BETWEEN TWO FIELDS

        PROGRAM SVD2
c 
c grid point: NL for SST; NR for Tox

c      time: 1979 - 2008; 30 years= 360 point
      parameter (imonth=360, mmax = 12)

	PARAMETER(NK=20,NL=4425,NR=4425,NT=360, nx=134, ny=61)
	parameter (m1=4425,m2=4425)
	parameter (imt= 134, jmt=61 )

	dimension sst1(nx,ny),u(nx,ny),tmp(nx,ny)
	dimension ww(nx*ny)

        dimension    lsi(nx,ny),     lsi0(nx,ny)
        dimension index1(nx,ny),   index2(nx,ny)

        DIMENSION sst(NT,NL),tox(NT,NR),rgr(nl,nk)
        DIMENSION test1(nx,ny),test2(nx,ny)

	DIMENSION COVAZS20(NL,NK),SIGVET20(NR,NK)

	DIMENSION COVAZS(NL,NR),SIGVET2(NR,NR)
	DIMENSION   CORR(NL,NR),VARL(NL),VARR(NR)
	DIMENSION SCF(NK),W(NR)
        DIMENSION AK(NT,NK),BK(NT,NK)
        DIMENSION CSCF(NK)

c zhang
        DIMENSION AKnorm(NT,NK),BKnorm(NT,NK),sdt1(nk),sdt2(nk)
	DIMENSION  varmode(nk), sdvmode(nk) 	
	dimension para1(NK)
	DIMENSION  ZAVE(NR),DEVIZ(NR)	

c zhang
        REAL HEL(NL,NK),HER(NR,NK),AAK(NK),BBK(NR)

	REAL AVZ(NL),AVS(NR),R(NK)
        real mask(imt,jmt)
        integer im

      character cdir*4, fname*16

      parameter (imt0=nx, jmt0=ny)

      real    odata(imt0,jmt0)


c..........................................................
      common/plt/tt(nx,ny,nk),qq(nx,ny,nk),tst(nk),tsq(nk)
      DIMENSION toxtox(nx,ny,nk)

      character*4 varnam(4)
      character*16 mname(4)
      character*20 nvert


      open(unit=6,file='svd.out',form='formatted')
c////////////////////////////////////////////////////////////////////////

      open(3,file='data/mask-sst.dat',
     &       form='formatted')
      read(3,*) ((mask(i,j),i=1,nx),j=1,ny)
      close(3)

c----------------------------
      do i=1,nx
        do j=1,ny
        lsi(i,j)=int(mask(i,j))
        enddo
        enddo



 789    format(1x,17i3)
         do j=1, jmt
	 k=jmt+1-j
	 write(6,238) k,(lsi(i,k),i=1,40)
	 end do
 238        format(1x,i2,1x,40i2  )  

       ic=0
       do  i=1, nx
       do  j=1, ny
       index2(i,j)=lsi(i,j)

       if(lsi(i,j) .ne. 0) then
         ic=ic+1
       endif
       enddo
       enddo
       write(*,*)"sst grid point=", ic

 
c      stop
c
c       read in Tdan anomaly
c
      open(20,file='data/ersst_anom.1979-2008.mgrid.dat',
     &     form='unformatted',access='direct',recl=imt0*jmt0)

      do im = 1, imonth      ! 1979-2008

	 read(20,rec=im) ((odata(i,j),i=1,imt0),j=1,jmt0)
         
	 kk=im
	 write(*,*) im, odata(17,15), kk

         ic=0
         do i=1,nx
         do j=1,ny 
         if(lsi(i,j).eq.1 ) then
	   ic=ic+1
 	   sst(kk,ic)=odata(i,j)
         end if
         end do
         end do

       end do

      close (20)

	do j=1,jmt,3
	jj0=jmt+1-j
	write(6,139)jj0,(odata(i,jj0),i=1,9)
	end do
 139    format(1x,i2,1x,9(f4.2,1x))


c-------------------------------------------------------------


c       read in fwf anomaly data:                          should be e-p

      open(20,file='data/fwf_anom.1979-2008.kl.dat',       !!!! p-e
     &     form='unformatted',access='direct',recl=imt0*jmt0*1 )
!      open(20,file='data/Salt_Argo.monthly.200501-201009.mgrid.txt',
!     &     form='formatted')

      do im = 1, imonth      ! 

	 read(20,rec=im) ((odata(i,j),i=1,imt0),j=1,jmt0)
!         read(20,1001) ((odata(i,j),i=1,imt0),j=1,jmt0)
         
	 kk=im
	 write(*,*) im, odata(17,15), kk

         ic=0
         do i=1,nx
         do j=1,ny 
         if(index2(i,j).eq.1 ) then
	   ic=ic+1
!	   tox(kk,ic)=  odata(i,j)   !    p-e
 	   tox(kk,ic)=- odata(i,j)   !!!  e-p
         end if
         end do
         end do

       end do

      close (20)
	do j=1,jmt
	jj0=jmt+1-j
	write(6,139)jj0,(odata(i,jj0),i=1,9)
	end do


c ----------------------------------------


c remove the time mean at every spatial point 

        do 9 i=1, nl
        aa=0.
        do 10 n=1, nt
        aa=aa+sst(n,i)/float(nt)
10      continue
        do 12 n=1, nt
        sst(n,i)=sst(n,i)-aa
12      continue
9       continue

        do 8 i=1,nr
        aa=0.
        do 6 n=1, nt
        aa=aa+tox(n,i)/float(nt)
6       continue
        do 7 n=1, nt
        tox(n,i)=tox(n,i)-aa
7       continue

8       continue

 234    format(1x,i2,1x,34(i2,1x))
 237    format(//1x)
 235    format(1x,i2,1x,17(f3.0,1x))

c	stop

c/////////////////////////////////////////////////////////////////////////

c   zhang: normalized !    there are problems when normalizing !
c    zhang-sdv of hh  0.2666154

	do j=1,nr
	DEVIZ(j)=0.0
	enddo

c       CALL NORMAL(sst,NL,NT,ZAVE,DEVIZ)
	aa=0.0
        do j=1,nl
        aa=aa+DEVIZ(j)/float(nl)
	enddo
	write(*,*)"zhang-sdv of sst", aa

c  spatially averaged standard deviation of sst is
c  zhang-sdv of tdan  0.7872732  

	do j=1,nr
	DEVIZ(j)=0.0
	enddo

c       CALL NORMAL(tox,NR,NT,ZAVE,DEVIZ)
	aa=0.0
        do j=1,nr
        aa=aa+DEVIZ(j)/float(nr)
	enddo
	write(*,*)"zhang-sdv of tox", aa
c  spatially averaged standard deviation of tox is

c	stop
 891    continue

!        sstsdv=0.4980432
!        ssttox=0.1585198
!       sstsdv=5.855863
!       ssttox=0.5136144

        sstsdv=0.5550922
        ssttox=2.139869
! zhang-sdv of sst  0.5550922
! zhang-sdv of tox   2.139869     Zheng Fei prepared data

        do i=1, nl
        do n=1, nt
        sst(n,i)=sst(n,i)/sstsdv
	enddo
	enddo

        do n=1,nt
        do i=1,nr
        tox(n,i)=tox(n,i)/ssttox
	enddo
	write(*,*) nt, n, sst(n,200), tox(n,200)
	enddo


C GET COVARIANCE MATRIX	

        CALL COVARI(sst, tox,NL,NR,NT,COVAZS,AVZ,AVS) 

 	WRITE(6,*) 'COVARIANCE FINISHED'


 
 
C CALCULATE SVD, GET SINGULAR VECTER ONE (LEFT VECTER,SIGVET1),SINGULAR
C VECTER TWO (RIGHT VECTER ,SIGVET2) 
 
  	CALL SVDCMP(COVAZS,NL,NR,NL,NR,W,SIGVET2)
	do mode=1, nk
	do n0  =1, nl
        COVAZS20(n0,mode)=COVAZS(n0,mode)
c       print *, "zhengf test",COVAZS20(n0,mode)
	end do
	do n0  =1,nr
	SIGVET20(n0,mode)=SIGVET2(n0,mode)
c       print *, "zhengf test",SIGVET20(n0,mode)
	end do
	end do
 
  
	WRITE(6,*) 'the first 20 singular values,  W(I)='
	WRITE(6,*) (W(I),I=1,20)
 	WRITE(6,*)'SVD FINISHED'
11      format (1x,5f10.5)
 
C CALCULATE SCF AND CSCF(SQUARED COVARIANCE FRACTION AND ACUMULATED SCF)
  
	CALL SCF1(W,NR,NK,SCF,CSCF)
  
	write(6,*)
	WRITE(6,*) 'the first 10 squared covariance fraction,SCF='
	WRITE(6,11) (SCF(I),I=1,NK)
	WRITE(6,*) 'the first 10 accumulated squared covariance fraction,CSCF='
	WRITE(6,11) (CSCF(I),I=1,NK)

C CALCULATE THE COEFFICIENTS (time series of the two fields) AK AND BK 

c	CALL COEFFI(COVAZS ,NL,NR,NT,NK,AK,sst)
c   	CALL COEFFI(SIGVET2,NR,NR,NT,NK,BK,tox)
	CALL COEFFI(COVAZS20,NL,NR,NT,NK,AK,sst)
   	CALL COEFFI(SIGVET20,NR,NR,NT,NK,BK,tox)

c zhang
	do j=1, nk

c (1)
        aa=0.
        do n=1, nt
        aa=aa+AK(n,j)/float(nt)
        end do

        do  n=1, nt
        AKnorm(n,j)=AK(n,j)-aa
        end do

	bb=0.0
	do  n=1, nt
	bb=bb + AKnorm(n,j)*AKnorm(n,j)
	end do
	bb=bb/FLOAT(nt-1)
	sdt1(j)=SQRT(bb)

        do  n=1, nt
        AKnorm(n,j)=AK(n,j)/sdt1(j)
        end do

        end do

c (2)
	do j=1, nk

        aa=0.
        do n=1, nt
        aa=aa+BK(n,j)/float(nt)
        end do

        do  n=1, nt
        BKnorm(n,j)=BK(n,j)-aa
        end do

	bb=0.0
	do  n=1, nt
	bb=bb + BKnorm(n,j)*BKnorm(n,j)
	end do
	bb=bb/FLOAT(nt-1)
	sdt2(j)=SQRT(bb)

        do  n=1, nt
        BKnorm(n,j)=BK(n,j)/sdt2(j)
        end do

        end do

	do j=1,nk
	write(*,*) j, sdt1(j), sdt2(j)
	end do

c ----------------------------

	DO  J=1,NK
  	    CALL CORREL(AK(1,J),BK(1,J),NT,R(J))
        END DO	
	write(6,*)
	WRITE(6,*)' corr. coef. of the 1st ten pairs'
	WRITE(6,11) R



C CALCULATE THE HETEROGENEOUS CORRELATION FIELDS

	CALL HOMO(BK,sst,NL,NT,NK,HEL,AAK,BBK)
	CALL HOMO(AK,tox,NR,NT,NK,HER,AAK,BBK)

	DO 20 J=1,NK
	   CALL CORREL(AK(1,J),BK(1,J),NT,R(J))
20      continue
	write(6,*)
	WRITE(6,*)' corr. coef. of the 1st ten pairs'
	WRITE(6,11) R

	CALL HOMO(AK,sst,NL,NT,NK,HEL,AAK,BBK)
	CALL HOMO(BK,tox,NR,NT,NK,HER,AAK,BBK)

	DO J=1,NK
	   CALL CORREL(AK(1,J),BK(1,J),NT,R(J))
	end do
	write(6,*)
	WRITE(6,*)' corr. coef. of the 1st ten pairs'
	WRITE(6,11) R

c  	stop
c1890   continue




c calculation of regression maps.....

        do 200 n=1, nk

c first get the variance of sst time series 

        cc=0.
        do m=1, nt
        cc=cc+bk(m,n)/float(nt)
        end do

        aa=0.
        do 201 m=1, nt
        aa=aa+(bk(m,n)-cc)*(bk(m,n)-cc)/float(nt-1)
201     continue
	varmode(n)=aa
	sdvmode(n)=sqrt(aa)

        do 202 i=1, nl
        rgr(i,n)=w(n)*covazs(i,n)/aa
202     continue

200     continue

	write(*,*) "zhang-w(n),sdvmode(n),varmode(n)"
        do n=1,nk
	write(6,342)w(n),sdvmode(n),varmode(n),w(n)/varmode(n)
        end do
 342    format(1x, 4(f8.2,1x) )


c calculate the normalized root-mean-squared covariance

c 	 CALL VAR(sst,NL,NT,VARL,AVZ)
c        CALL VAR(tox,NR,NT,VARR,AVS)

c       CALL NORMAL(sst,NL,NT,ZAVE,DEVIZ)
c       CALL NORMAL(tox,NR,NT,ZAVE,DEVIZ)
c       CALL COVARI(sst,tox,NL,NR,NT,CORR,AVZ,AVS)
 
c       CALL ABCOVA(W,VARL,VARR,CORR,NL,NR,NK,PARA1,PARA2)
        write(6,*)
c       WRITE(6,*)'absolute covariance explained by these fields:'
c       WRITE(6,66) para1
66      format(5f12.4)

c////////////////////////////////////////////////////////////////////////
c  now hook up the variables to sumry for plotting......


c now cycle through the times and convert to full 2-d grid...

       do  l=1, nk                     ! nk=20
       do  i=1, nx
       do  j=1, ny
         tt    (i,j,l)=2000.
         qq    (i,j,l)=2000.
	 toxtox(i,j,l)=2000.
       enddo
       enddo
       enddo

       do 101 l=1, nk                     ! nk=20
       ic=0
       do 102 i=1, nx
       do 103 j=1, ny
       if(index2(i,j) .eq. 1) then
         ic=ic+1
         qq    (i,j,l)=rgr    (ic,l)
         toxtox(i,j,l)=sigvet2(ic,l)*sdt2(l)
       endif
103    continue
102    continue
101    continue

       do  l=1, nk                     ! nk=20
       ic=0
       do 106 i=1, nx
       do 105 j=1, ny
       if(lsi(i,j) .eq. 1) then
         ic=ic+1
	 tt    (i,j,l)=COVAZS(ic,l)*sdt1(l)
       endif
105    continue
106    continue
       end do

c      write(243) COVAZS,SIGVET2,W,sdt1,sdt2,ak,bk
       write(243) COVAZS20,SIGVET20,W,sdt1,sdt2,ak,bk
       
c      write(236) tt
c      write(237) toxtox

 	STOP
	END
