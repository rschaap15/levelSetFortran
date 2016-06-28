!> \brief
!! This program generates a 1D distribution of points for strand grids.
!!
!! Comments:
!!
!! Versions:
!!   - 1.0 Katz 07/23/2010
!!
!! Additional notes:\par
!!   none
!!
!! Source code:
!!   \include Meshing/strand1dDist.F90


SUBROUTINE strand1dDist( &
     myid, &
     nmax, &
     xs, &
     nPtsPerStrand, &
     stretchRatio, &
     strandDist, &
     strandLength, &
     wallSpacing)

IMPLICIT NONE

INTEGER,INTENT(IN   ) :: myid,nmax,strandDist
INTEGER,INTENT(INOUT) :: nPtsPerStrand
REAL   ,INTENT(INOUT) :: stretchRatio,strandLength,wallSpacing
REAL   ,DIMENSION(nmax),INTENT(  OUT) :: xs
INTEGER :: npst0,func,n,k
REAL    :: xi,xf,dxi,dxf,dxfa,dxia


xi    = 0.
xf    = strandLength
dxi   = wallspacing
dxf   = 0.
npst0 = 0
IF      (strandDist == 1   ) THEN
   IF (myid == 0) WRITE(*,*)'*** using uniform strand spacing ***'
   func  = 0
   dxi   = 0.
   dxf   = 0.
   npst0 = nPtsPerStrand+1 ! add 1 to index from 0
   IF (nPtsPerStrand < 1) THEN
      IF (myid == 0) &
           WRITE(*,*)'*** for uniform spacing, must have nPtsPerStrand >= 1 ***'
      STOP
   END IF
ELSE IF (strandDist == 2 ) THEN
   IF (myid == 0) WRITE(*,*)'*** using geometric strand spacing ***'
   func = 1
   IF (nPtsPerStrand < 1) THEN
      IF (myid == 0) &
           WRITE(*,*)'*** using stretchRatio to determine nPtsPerStrand***'
      IF (stretchRatio < 1.) THEN
         IF (myid == 0) WRITE(*,*)'*** stretchRatio must be >= 1. ***'
         STOP
      END IF
   ELSE
      npst0 = nPtsPerStrand+1 ! add 1 to index from 0
   END IF
ELSE IF (strandDist == 3) THEN
   IF (myid == 0) WRITE(*,*)'*** using hyperbolic strand spacing ***'
   func = 2
   IF (nPtsPerStrand < 1) THEN
      IF (myid == 0) &
           WRITE(*,*)'*** using stretchRatio to determine nPtsPerStrand***'
      IF (stretchRatio < 1.) THEN
         IF (myid == 0) WRITE(*,*)'*** stretchRatio must be >= 1. ***'
         STOP
      END IF
   ELSE
      npst0 = nPtsPerStrand+1 ! add 1 to index from 0
   END IF
ELSE
   IF (myid == 0) &
        WRITE(*,*)'*** strandDist not recognized in subroutine strandDist ***'
   STOP
END IF

CALL sfuns(func,xi,xf,dxi,dxf,npst0,stretchRatio,nmax,dxia,dxfa,xs)

nPtsPerStrand = npst0-1


END SUBROUTINE strand1dDist
