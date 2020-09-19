      program test_sneut5
      implicit none

! tests the neutrino loss rate routine

! ionmax  = number of isotopes in the network
! xmass   = mass fractions
! ymass   = molar fractions
! aion    = number of nucleons
! zion    = number of protons

      integer, parameter :: ionmax=3
      real*8             :: xmass(ionmax),ymass(ionmax), &
                            aion(ionmax),zion(ionmax),wion(ionmax), &
                            temp,den,abar,zbar,wbar,ye, nxcess, &
                            snu,dsnudt,dsnudd,dsnuda,dsnudz, &
                            spair,splas,sphot,sbrem,sreco

! formats
 01   format(1x,t2,a6,1pe11.4,t22,a6,1pe11.4, &
               t42,a6,1pe11.4,t62,a6,1pe11.4)
 02   format(1x,a,1pe11.4)


! set the mass fractions, z's and a's of the composition
! hydrogen
      xmass(1) = 0.75d0
      zion(1)  = 1.0d0
      aion(1)  = 1.0d0
      wion(1)  = aion(1)

! helium
      xmass(2) = 0.23d0
      zion(2)  = 2.0d0
      aion(2)  = 4.0d0
      wion(2)  = aion(2)

! carbon 12
      xmass(3) = 0.02d0
      zion(3)  = 6.0d0
      aion(3)  = 12.0d0
      wion(3)  = aion(3)


! get abar and zbar

      call azbar(xmass,aion,zion,wion,ionmax, &
                 ymass,abar,zbar,wbar,ye,nxcess)

! set the thermodynamic state
      temp = 1.0d9
      den  = 1.0d6


! get the neutrino losses
      call sneut5_aa(temp,den,abar,zbar, &
                  snu,dsnudt,dsnudd,dsnuda,dsnudz, &
                  spair,splas,sphot,sbrem,sreco)


! report the results
      write(6,01) 'temp =',temp,'den  =',den, &
                  'abar =',abar,'zbar =',zbar
      write(6,*)
      write(6,02) 'snu   =',snu, &
                  'dsnudt=',dsnudt, &
                  'dsnudd=',dsnudd, &
                  'dsnuda=',dsnuda, &
                  'dsnudz=',dsnudz
      write(6,*)
      write(6,02) 'spair   =',spair, &
                  'splas   =',splas, &
                  'sphot   =',sphot, &
                  'sbrem   =',sbrem, &
                  'sreco   =',sreco

      write(6,*)
      stop 'normal termination'
      end


      include 'azbar.f90'
      include 'sneut5.f90'
      include 'funct_fermi1.f90'
