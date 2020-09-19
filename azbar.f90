      subroutine azbar(xmass,aion,zion,wion,ionmax, &
                       ymass,abar,zbar,wbar,ye,nxcess)
      include 'implno.dek'

! this routine calculates composition variables

! input:
! mass fractions               = xmass(1:ionmax)  dimensionless
! number of nucleons           = aion(1:ionmax)   dimensionless
! charge of nucleus            = zion(1:ionmax)   dimensionless
! atomic weight or molar mass  = wion(1:ionmax)   dimensionless
! number of isotopes           = ionmax
!
! output:
! molar abundances        = ymass(1:ionmax)   dimensionless
! mean number of nucleons = abar              dimensionless
! mean nucleon charge     = zbar              dimensionless
! mean weight             = wbar              dimensionless
! electron fraction       = ye                dimensionaless
! neutron excess          = xcess             dimensionless


! declare the pass
      integer :: ionmax
      real*8  :: xmass(ionmax),aion(ionmax),zion(ionmax), &
                 wion(ionmax),ymass(ionmax),abar,zbar,wbar,ye,nxcess

! local variables
      real*8  ::  sum1

! molar abundances
      ymass(1:ionmax) = xmass(1:ionmax)/aion(1:ionmax)

! mean numner of nucleons
      abar  = 1.0d0/sum(ymass(1:ionmax))
      wbar  = abar

! mean charge
      ye  = sum(zion(1:ionmax)*ymass(1:ionmax))
      zbar  = abar * ye

! neutron excess
      nxcess = 1.0d0 - 2.0d0 * ye

      return
      end subroutine azbar




      subroutine azbar_simple(ymass,aion,zion,ionbeg,ionend, &
                              abar,zbar,ye)
      include 'implno.dek'

! this routine calculates common composition variables

! input:
! molar abundances             = ymass(1:ionmax)  dimensionless
! number of nucleons           = aion(1:ionmax)   dimensionless
! charge of nucleus            = zion(1:ionmax)   dimensionless
! start of isotopes            = ionbeg
! end of isotopes              = ionend
!
! output:
! mean number of nucleons = abar              dimensionless
! mean nucleon charge     = zbar              dimensionless
! electron fraction       = ye                dimensionless


! declare the pass
      integer :: ionbeg,ionend
      real*8  :: ymass(ionbeg:ionend),aion(ionbeg:ionend),zion(ionbeg:ionend),abar,zbar,ye

! local variables
      real*8  :: sum1

! mean number of nucleons, numerator assumes sum xmass = 1.0d0
      abar  = 1.0d0/sum(ymass(ionbeg:ionend))

! mean electron fraction and mean charge
      ye    = sum(zion(ionbeg:ionend)*ymass(ionbeg:ionend))
      zbar  = abar * ye

      return
      end subroutine azbar_simple



