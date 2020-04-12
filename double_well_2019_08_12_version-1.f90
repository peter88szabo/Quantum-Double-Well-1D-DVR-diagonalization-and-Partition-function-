!========================================================================================
      module parameters
!========================================================================================
      save
      real*8 :: rmin,rmax,dr
      real*8 :: pi,redmass,omega
      real*8 :: c1,c2,c3,c4,c5,c6,c7
      end module
!========================================================================================



      program boundstates 
      use parameters

      implicit real*8(a-h,o-z)
      parameter (ndim=5000)
      dimension Vpot(ndim),H(ndim,ndim),Energy(ndim),rdist(ndim)

      pi=4.0*atan(1.d0) 

!     Constants and unit conversion factors:

!     [Anstrom]*c1=[bohr]
      c1=1.0d0/0.5291772d0
!     [kcal/mol]*c2=[Hartree]
      c2=1.d0/627.51d0
!     [g/mol]*c3=[electron mass unit]
      c3=1838.6836605d0
!     [Hartree]*c4=[eV]
      c4=27.2114
!     [Hartree]*c5=[cm-1]
      c5=219474.d0
!     [femto-sec]*c6=[time in au]
      c6=41.341105
!     [Hartree]*c7=[kJ/mol]
      c7=2625.5d0



      write(*,*) "give ndvr:"
      read(*,*) ndvr

      rmin=-5.d0
      rmax=5.d0
      redmass=1.d0



!----------------------------------------------------------------------------------
      dr=(rmax-rmin)/float(ndvr)

      write(6,*) "====================================================="
      write(6,"(a16,f8.3,a5,f8.3/,a7,i5,/,a5,f10.5)")                   &
              "rmin and rmax =",rmin,"and",rmax,"ndvr =",ndvr,"dr =",dr
      write(6,*) 
!----------------------------------------------------------------------------------
      if(ndvr >= ndim ) stop "WARNING: Too many DVR-points!"

      do i=1,ndvr
         rr=rmin+float(i)*dr
         rdist(i)=rr
         Vpot(i)=Vdbw(rr)
         write(777,*) rr,Vpot(i)
      enddo



      call Hamiltonian_DVR(pi,ndim,ndvr,Vpot,dr,redmass,H,Energy)

!------------------------------------------------------------------------------------
!     Print of the results of the diagonalization
!------------------------------------------------------------------------------------
      do i=1,ndvr
         write(66,"(f12.4,i5,f15.8,2e15.5)") rdist(i),i-1,Energy(i),H(i,1)/sqrt(dr),H(i,2)/sqrt(dr)
      enddo


     
      do i=1,20
         beta=float(i)/2.d0
         call partition_func(ndim,ndvr/2,beta,Energy,prtf)

         write(67,"(f10.1,2e15.5,a40)")                                 &
         beta,prtf,dexp(-beta*Energy(1)),"<--beta,prtf,exp(-beta*E(1))"

         write(*,"(f10.1,2e15.5,a40)")                                 &
         beta,prtf,dexp(-beta*Energy(1)),"<--beta,prtf,exp(-beta*E(1))"

      enddo


 
      end program
       


!===================================================================================
      subroutine partition_func(ndim,nmax,beta,Energy,prtf)
!===================================================================================
      implicit none
      integer :: i,ndim,nmax
      real*8 :: beta,dum,prtf
      real*8, dimension(ndim) :: Energy

      dum=0.d0
      do i=1,nmax
         dum=dum+dexp(-beta*Energy(i))
      enddo

      prtf=dum

      end subroutine
!==================================================================================


 



!===================================================================================
      real*8 function Vdbw(x)
!===================================================================================
      use parameters
      implicit none
      real*8 x,ss

      ss=(x*x-1.d0)
      Vdbw=ss*ss

      return
      end function
!==================================================================================



!===================================================================================
      subroutine Hamiltonian_DVR(pi,ndim,ndvr,Vpot,dr,redmass,H,E)
!===================================================================================
!     This routine finds the eigenvalues for a bound 1D potential
!     
!     Hamiltononian is represented wth the Sinc-DVR basis
!     (D.T. Colbert, H.J. Miller, J. Chem. Phys, Vol. 96, pp 1982, 1992)
!      
!     After diagonalization the columns of H contain the eigenvectors and
!     E contains the corresponding eigen-energies
!===================================================================================
      implicit none
      integer ndim, ndvr
      integer i,j
      real*8 H(ndim,ndim),E(ndim),Vpot(ndim)
      real*8 redmass,dr,pi


!     Build up the Hamiltonian
      do i=1,ndvr
         do j=i,ndvr !only upper triangular is needed
            if (i .eq. j) then
                  H(i,j)= pi*pi/(6.0*redmass*dr*dr) + Vpot(i)
            else
               if (mod(i-j, 2) .eq. 0) then
                  H(i,j)=  1.0/(redmass*dr*dr*float(i-j)*float(i-j))
               else
                  H(i,j)= -1.0/(redmass*dr*dr*float(i-j)*float(i-j))
               endif
            endif

         enddo
      enddo

      write(*,*) "DVR representation of the Hamiltionian is done"
      write(*,*) "Please, wait for the diagonalization!"

!     Diagonalization of the Hamiltonian
      call seigen(ndvr,H,ndim,E)

      return
      end
!==================================================================================





!==================================================================================
      subroutine seigen(n, H, nmax, E)
!==================================================================================
!     Lapack Eigensolver for a real symmetric matrix H
      implicit none
      integer n, nmax
      real*8 H(nmax, nmax), E(nmax)
      real*8 work(34*nmax)
      integer ierr
      call dsyev('V', 'U', n, H, nmax, E, work, 34*nmax, ierr)
      end
!==================================================================================

