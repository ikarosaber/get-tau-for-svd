c grid point: NL for SST; NR for Tox

	PARAMETER(NK=20,NL=4425,NR=4425 ,nr2=nr )
	PARAMETER(NT=30 , nx=134, ny=61 )
	parameter (ny1=61, m1=4425,m2=4425,np =1 )

        dimension    lsi(nx,ny),index0(nx,ny)
	dimension sst1(nx,ny1),u(nx,ny)
	dimension mask(nx,ny)

	DIMENSION COVAZS (NL,NR),SIGVET2 (NR,NR),w (nr)
	DIMENSION COVAZS2(NL,NR),SIGVET22(NR,NR),w2(nr)

        DIMENSION COVAZS20(NL,NK),SIGVET20(NR,NK)   

	dimension ww(nx*ny), sdt1(nk),sdt2(nk)
        DIMENSION AK(NT,NK),BK(NT,NK)

c zhang
	DIMENSION AKnorm(NT,NK),BKnorm(NT,NK)

      parameter (imt0=134, jmt0=61)

c..........................................................
      DIMENSION toxtox(nx,ny,nk)
      dimension zc(nx,ny,nk)
      dimension change(nx,ny),arr0(m1,20,np),brr0(m2,20,np)

      dimension sst_d(12),tox_d(12),toy_d(12)

      dimension neign(np)
      common /eigen/ist1(m1),ist2(m2),arr(m1,20,np)
     &,         brr(m2,20,np),stdv(20,np),sd,ud(np)
      
      common /eigen12/arr12(m1,20,np,12),brr12(m2,20,np,12)
     &               ,  stdv12(20,np,12),sd12(12),ud12(np,12)

c       read(89,88) lsi
 88     format(60i1)	 

      open(3,file='data/mask-sst.dat',
     &       form='formatted')
      read(3,*) ((mask(i,j),i=1,nx),j=1,ny)
      close(3)

        do i=1,nx
        do j=1,ny
        lsi(i,j)=int(mask(i,j))
        enddo
        enddo

       ic=0
       do  i=1, nx
       do  j=1, ny
       index0(i,j)=lsi(i,j)
       if(lsi(i,j) .ne. 0) then
         ic=ic+1
       endif
       enddo
       enddo
       write(*,*)"sst grid point=", ic

       ic=0
       do  i=1, nx
       do  j=1, ny
       if(index0(i,j) .ne. 0) then
         ic=ic+1
       endif
       enddo
       enddo
       write(*,*)"tox grid point=", ic

c zhang: unit from N/M/M => Dyn/cm/cm if necessary !!!!!

       ud_sst  =0.5017650
       ud_tox  =0.1839231  
       ud_toy  =0.1352975  

       ij0=0
       do  j=1, ny
       do  i=1, nx
       if(lsi(i,j) .eq. 1) then
         ij0=ij0+1
	 ist1(ij0)=i+(j-1)*nx
       endif
       enddo
       enddo
c      write(*,*) ij0, (ist1(ij),ij=1,10)

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
       write(*,*)"tox grid point=", ic
       write(*,*) (ist2(ij),ij=1,10)
 
       read(978) sst_d,tox_d
c      read(978) sst_d,tox_d,toy_d
       write(*,*)"zhang-sst_d", sst_d
       write(*,*)"zhang-tox_d", tox_d

c ----------------------------------------------
       do 789 month = 1, 12

	sd   =sst_d(month)
	ud(1)=tox_d(month)

c	ud(2)=toy_d(month)
c    tox & toy:  1, 2

        nread0=540+month
 	read(nread0) COVAZS20,SIGVET20,W2,sdt1,sdt2,ak,bk     

	do n=1,nk
	write(6,*)"toy",month, w2(n),sdt1(n),sdt2(n)
	end do

	do mode=1,nk

         stdv(mode,1)=sdt1(mode)
 	 do m0  =1,m1
	 arr(m0,mode,1)=COVAZS20(m0,mode)*sdt1(mode)
	 end do
	 do m0  =1,m2
	 brr(m0,mode,1)=SIGVET20(m0,mode)*sdt2(mode)
	 end do

c        stdv(mode,2)=sdt1(mode)
c	 do m0  =1,m1
c	 arr(m0,mode,2)=COVAZS20(m0,mode)*sdt1(mode)
c	 end do
c	 do m0  =1,m2
c	 brr(m0,mode,2)=SIGVET20(m0 + m2,mode)*sdt2(mode)
c	 end do

	end do

  234  continue


c---------------------------
       do nvar=1, np
       do mode=1, nk

       ij0=0
       do  i=1, nx
       do  j=1, ny
       if(lsi(i,j) .eq. 1) then
         ij0=ij0+1
         change(i,j)=arr(ij0,mode,nvar)
       endif
       enddo
       enddo

       ij0=0
       do  j=1, ny1
       do  i=1, nx
       if(lsi(i,j) .eq. 1) then
         ij0=ij0+1
         arr0(ij0,mode,nvar)=change(i,j)
       endif
       enddo
       enddo

       ij0=0
       do  i=1, nx
       do  j=1, ny
       if(index0(i,j) .eq. 1) then
         ij0=ij0+1
         change(i,j)=brr(ij0,mode,nvar)
       endif
       enddo
       enddo

       ij0=0
       do  j=1, ny
       do  i=1, nx
       if(index0(i,j) .eq. 1) then
         ij0=ij0+1
         brr0(ij0,mode,nvar)=change(i,j)
       endif
       enddo
       enddo

       enddo
       enddo
c---------------------------

       do np0=1, np
       do k  =1, nk
	 do ij=1, m1
	 arr12(ij,k,np0,month)=arr0(ij,k,np0)
	 end do

	 do ij=1, m2
	 brr12(ij,k,np0,month)=brr0(ij,k,np0)
	 end do

	 stdv12(k,np0,month)  =stdv(   k,np0)
       end do
       end do

	sd12(    month)=sd
	ud12(1  ,month)=tox_d(month) 

c	ud12(2  ,month)=toy_d(month) 

 789    continue            !   month  cycle


        write(*,*) "sd12", sd12

c	write(96) arr12,brr12,stdv12,sd12,ud12
	write(97) arr12,brr12,stdv12,sd12,ud12

 242    format(1x, 3(f8.2,1x) )
  243    format(1x,i3,1x,6(f8.2,1x) )
 	STOP
	END
