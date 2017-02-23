!C THIS PROGRAM IS FOR CALCULATING SVD BETWEEN TWO FIELDS
 
!c grid point: NL for SST; NR for Tox

	PARAMETER(NK=20,NL=4425,NR=4425)
 	PARAMETER(NT=360, nx=134, ny=61) 
	parameter (ny1=61,m1=4425,m2=4425,np=1)
        parameter (imt0=nx,jmt0=ny)

        dimension    lsi(nx,ny),index0(nx,ny)
	dimension sst1(nx,ny1),u(nx,ny),tmp(nx,ny)

	DIMENSION COVAZS20(NL,NK),SIGVET20(NR,NK)

	DIMENSION COVAZS(NL,NR),SIGVET2(NR,NR),w(nr)
	dimension ww(nx*ny), sdt1(nk),sdt2(nk)
        DIMENSION AK(NT,NK),BK(NT,NK), w20(NK)

      DIMENSION AKnorm(NT,NK),BKnorm(NT,NK)
      DIMENSION toxtox(nx,ny,nk)
      dimension zc(nx,ny,nk)
      dimension change(nx,ny),arr0(m1,20,np),brr0(m2,20,np)

      dimension neign(np)
      common /eigen/ist1(m1),ist2(m2),arr(m1,20,np),
     &              brr(m2,20,np),stdv(20,np),sd,ud(np)
      real mask(imt0,jmt0),temp1(imt0,jmt0,nk),temp2(imt0,jmt0,nk)
      

        open(3,file='data/mask-sst.dat',form='formatted')
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

!      ud_sst=5.855863
!      ud_tox=0.5136144

        ud_sst=0.5550922     ! freshwater flux
        ud_tox=2.139869


       ij0=0
       do  j=1, ny
       do  i=1, nx
       if(lsi(i,j) .eq. 1) then
         ij0=ij0+1
	 ist1(ij0)=i+(j-1)*nx
       endif
       enddo
       enddo
       write(*,*) ij0, (ist1(ij),ij=1,10)

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

	sd   =ud_sst
	ud(1)=ud_tox

!       goto 234
  233  continue

 	read(243) COVAZS20,SIGVET20,W,sdt1,sdt2,ak,bk   
!	read(243) COVAZS,SIGVET2,W,sdt1,sdt2,ak,bk   
        close(243)

	nbyte=1
        open(10,file='svd.akbk.1.dta',
     1         form='unformatted',access='direct',recl=2*nbyte)
        do i=1, nt
        write(10,rec=i) ak(i,1),bk(i,1)
        end do
        close(10)

	open(10,file='svd.akbk.2.dta',
     1         form='unformatted',access='direct',recl=2*nbyte)
	do i=1, nt
        write(10,rec=i) ak(i,2),bk(i,2)
	end do
	close(10)

	open(10,file='svd.w20.dta',
     1         form='unformatted',access='direct',recl=nbyte)
	do i=1,20
	w20(i)=W(i)
	end do
	do i=1,20
	write(10,rec=i) w20(i)
	end do

!	stop

	do n=1,nk
	write(6,*)"tox", w(n),sdt1(n),sdt2(n)
	end do

	do mode=1,nk
         stdv(mode,1)=sdt1(mode)

 	 do m0  =1,m1
!	 arr(m0,mode,1)=COVAZS(m0,mode)*sdt1(mode)
         arr(m0,mode,1)=COVAZS20(m0,mode)*sdt1(mode)
	 end do
	 do m0  =1,m2
!	 brr(m0,mode,1)=SIGVET2(m0,mode)*sdt2(mode)
         brr(m0,mode,1)=SIGVET20(m0,mode)*sdt2(mode)
	 end do
	end do
  234  continue


       temp1=999.0
       temp2=999.0

       do nvar=1, np
       do mode=1, nk

       ij0=0
       do  i=1, nx
       do  j=1, ny
       if(lsi(i,j) .eq. 1) then
         ij0=ij0+1
         change(i,j)=arr(ij0,mode,nvar)
         temp1(i,j,mode)=change(i,j)*sd
       endif
       enddo
       enddo

       ij0=0
       do  j=1, ny 
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
         temp2(i,j,mode)=change(i,j)*ud(1)
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

c       write(86) ist1,ist2,arr0,brr0,stdv,sd,ud
        write(87) ist1,ist2,arr0,brr0,stdv,sd,ud
!c       read (86) ist1,ist2,arr ,brr ,stdv,sd,ud

	do n=1,nk
	write(*,243) n,(stdv(n,np0), np0=1,4 )
	end do

        write(*,*) arr(1,1,1),arr(10,1,1),brr(1,1,1),brr(10,1,1)

 242    format(1x, 3(f8.2,1x) )
 243    format(1x,i3,1x,6(f8.2,1x) )

        open(10,file='svd.spatial.dta',
     1         form='unformatted',access='direct',recl=nx*ny)
        do mode=1,nk
         write(10,rec=mode*2-1) ((temp1(i,j,mode),i=1,nx),j=1,ny)
         write(10,rec=mode*2)   ((temp2(i,j,mode),i=1,nx),j=1,ny) 
        enddo
        close(10)

 	STOP
	END
