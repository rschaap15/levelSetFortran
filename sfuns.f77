C***********************************************************************
      SUBROUTINE SFUNS(ISF,X0,X1,DX0,DX1,N,RMAX,NMAX,DX0A,DX1A,X)
C
C    Authors : William M. Chan and Pieter G. Buning
C    Version : 1.0e
C    Date    : January, 99.
C
C    This routine computes a stretched arc length distribution (1-D
C    stretching function) in the range [X0,X1]. Three stretching
C    functions are available: uniform, geometric and hyperbolic tangent.
C    The geometric stretching allows grid spacing specification at either
C    one or the other end of the domain while the hyperbolic tangent
C    stretching allows grid spacing specification at either end or both
C    ends of the domain.
C
C    For the hyperbolic tangent stretching, the discrete spacing (output
C    in the X array) is usually slightly different from the analytic
C    spacing (input from user). The discrepancy typically gets smaller
C    as the number of points increases.
C
C    When an end point spacing is specified, the analytic spacing is the
C    specified spacing. When an end point spacing is not specified, the
C    analytic spacing is computed. In both cases, the end point analytic
C    spacings are returned.
C
C WMC
C 1.0e Converted EPSIL to always used double precision arithmetic.
C      Updated GEOUNI to treat case when DX1/DX0 is less than RMAX*RMAX.
C
C    Input variables
C
C     ISF  = 0  Uniform stretching (DX0 and DX1 are disregarded)
C          = 1  Geometric stretching
C          = 2  Hyperbolic tangent stretching
C          = 3  Geometric stretching up to end spacing, then uniform
C               for the rest of the distance. Fall back to hyperbolic
C               tangent if end spacing cannot be achieved with monotonic
C               increase of spacing.
C          = 4  New
C     X0   =    Value of arc length function at end point X0
C     X1   =    Value of arc length function at end point X1
C     DX0  =    Grid spacing at end point X0
C     DX1  =    Grid spacing at end point X1
C     N    > 0  Number of grid points including ends. RMAX is disregarded.
C          <=0  The smallest N such that RMAX is not exceeded will be
C               computed. RMAX>0 is assumed, otherwise a default is used.
C     RMAX =    Max stretching ratio allowed. This is assumed to be a
C               number greater than 1. If RMAX<1, 1/RMAX will be used.
C               (input RMAX is used only if input N <= 0)
C     NMAX =    Max number of points used if input N <= 0
C               (also equals max dimension of X)
C     DX0A =    Max uniform spacing for ISF=4
C
C    Returned variables
C
C     N     =  Number of grid points including ends
C     RMAX  =  Max stretching ratio (returned if input N > 0)
C     DX0A  =  Analytic grid spacing at end point X0
C     DX1A  =  Analytic grid spacing at end point X1
C     X     =  Array returned with coordinates of points in the
C              interval [X0,X1]
C
C    Let D = X1-X0. There are 8 acceptable combinations of ranges of values
C    for the variables D,DX0,DX1. If DX0 and DX1 are both non-zero when
C    ISF=1 is selected, DX1 is disregarded, i.e. grid spacing is fixed at
C    X0 only. In all cases, X(1) = X0 and X(N) = X1.
C
C       D  DX0  DX1               Description
C    ===============================================================
C      >0   >0   =0    Grid spacing fixed at X0, floated at X1 (*)
C      >0   =0   >0    Grid spacing fixed at X1, floated at X0 (*)
C      >0   >0   >0    Grid spacing fixed at X0, fixed at X1   (*)
C      >0   =0   =0    Uniform grid spacing                    (*)
C      <0   <0   =0    Grid spacing fixed at X0, floated at X1 (**)
C      <0   =0   <0    Grid spacing fixed at X1, floated at X0 (**)
C      <0   <0   <0    Grid spacing fixed at X0, fixed at X1   (**)
C      <0   =0   =0    Uniform grid spacing                    (**)
C
C    (*)  Increasing x coordinate from X(1) to X(N)
C    (**) Decreasing x coordinate from X(1) to X(N)
C    Cases 4 and 8 for ISF=1 or 2 are equivalent to selecting ISF=0.
C
C    Parameters
C
C     RDEF = Default RMAX used if input N and RMAX are not positive.
C
#include "precis.h"
C
      PARAMETER ( RDEF=1.3 )
      LOGICAL COMPSR,DIALOG
      DIMENSION X(NMAX)
C
      DIALOG = .FALSE.
C
C    Check input
C
      CALL SFICHK(ISF,X0,X1,DX0,DX1,N,RMAX,NMAX,RDEF,COMPSR,DIALOG)
C
C    Check for degenerate cases including uniform spacing
C
      CALL SFDGEN(X0,X1,DX0,DX1,N,RMAX,DX0A,DX1A,X,DIALOG,IEXIT)
C
C    Select stretching function
C
      IF (IEXIT.EQ.0) THEN
C
        IF      (ISF.EQ.1) THEN
C
         CALL GEOMS(X0,X1,DX0,DX1,N,NMAX,RMAX,DX0A,DX1A,X,DIALOG)
C
        ELSE IF (ISF.EQ.2) THEN
C
         CALL TANHS(X0,X1,DX0,DX1,N,NMAX,RMAX,DX0A,DX1A,X,DIALOG)
C
        ELSE IF (ISF.EQ.3) THEN
C
         CALL GEOUNI(X0,X1,DX0,DX1,N,NMAX,RMAX,DX0A,DX1A,X,DIALOG)
C
        ELSE IF (ISF.EQ.4) THEN
C
         CALL HYPUNI(X0,X1,DX0,DX1,N,NMAX,RMAX,DX0A,DX1A,X,DIALOG)
C
        ENDIF
C
        IF (COMPSR) CALL STRRAT(N,RMAX,X)
C
      ENDIF
C
      RETURN
      END
C***********************************************************************
      SUBROUTINE SFICHK(ISF,X0,X1,DX0,DX1,N,RMAX,NMAX,RDEF,
     &                  COMPSR,DIALOG)
C
#include "precis.h"
C
      PARAMETER (ZERO=0.0, ONE=1.0)
      LOGICAL COMPSR,DIALOG
C
C    Check input
C
      IF ( (ISF.LT.0) .OR. (ISF.GT.4) ) THEN
       IF (DIALOG)
     &  WRITE(*,*)'Input ISF is not 0/1/2/3/4. ISF is reset to 2.'
       ISF = 2
      ENDIF
C
      IF ( (ISF.EQ.1) .AND. (DX0.NE.ZERO) .AND. (DX1.NE.ZERO) ) THEN
       IF (DIALOG) THEN
         WRITE(*,*)'Warning: Grid spacing specified at both ends for ',
     &             '         geometric stretching.'
         WRITE(*,*)'Final grid spacing will be disregarded.'
       ENDIF
       DX1 = ZERO
      ENDIF
C
      IF (ISF.EQ.3) THEN
       IF (DX0.EQ.ZERO) THEN
        IF (DIALOG)
     &   WRITE(*,*)'Warning: Start spacing must be supplied for ISF=3.'
        STOP
       ENDIF
       IF (DX1.EQ.ZERO) THEN
        IF (DIALOG) THEN
         WRITE(*,*)'Warning: End spacing must be supplied for ISF=3.'
         WRITE(*,*)'Resetting end spacing to start spacing.'
        ENDIF
        DX1 = DX0
       ENDIF
       IF (N.GT.0) THEN
        IF (DIALOG) WRITE(*,*)'Warning: Cannot specify N with ISF=3.'
        N = 0
       ENDIF
      ENDIF
C
      IF (ISF.EQ.0) THEN
       DX0 = ZERO
       DX1 = ZERO
      ENDIF
C
      IF (N.LE.0) THEN
C
       IF (RMAX.LE.ZERO) THEN
        IF (DIALOG) THEN
          WRITE(*,*)'A positive RMAX must be supplied for N<=0.'
          WRITE(*,*)'RMAX is reset to ',RDEF
        ENDIF
        RMAX = RDEF
       ENDIF
       IF (RMAX.LT.ONE) THEN
        RMAX = ONE/RMAX
        IF (DIALOG) THEN
          WRITE(*,*)'Input stretching ratio < 1. Transformed to ',RMAX
        ENDIF
       ENDIF
       IF ( (DX0.EQ.ZERO) .AND. (DX1.EQ.ZERO) ) THEN
        IF (DIALOG) THEN
           WRITE(*,*)'N must be positive for DX0=DX1=0.'
           WRITE(*,*)'N is reset to ',NMAX
        ENDIF
        N = NMAX
       ENDIF
C
      ELSE IF (N.GT.NMAX) THEN
C
       IF (DIALOG) THEN
          WRITE(*,*)'Too many points. Make NMAX >= ',N
          WRITE(*,*)'For this run, N is reset to ',NMAX
       ENDIF
       N = NMAX
C
      ENDIF
C
      IF (N.GT.0) THEN
         COMPSR = .TRUE.
      ELSE
         COMPSR = .FALSE.
      ENDIF
C
      RETURN
      END
C***********************************************************************
      SUBROUTINE SFDGEN(X0,X1,DX0,DX1,N,RMAX,DX0A,DX1A,X,DIALOG,IEXIT)
C
#include "precis.h"
C
      PARAMETER (ZERO=0.0, ONE=1.0)
      LOGICAL DIALOG
      DIMENSION X(*)
C
C    Check for degenerate cases including uniform spacing
C
      IEXIT = 0
      D     = X1 - X0
C
C    Zero total arc length.
C
      IF (D.EQ.ZERO) THEN
       IF (DIALOG) WRITE(*,*)'Zero total arc length in input. X0=X1.'
       DO 10 I=1,N
        X(I) = X0
 10    CONTINUE
       IEXIT = 1
      ENDIF
C
C    D*DX0 and D*DX1 cannot be negative.
C
      IF ( (D*DX0.LT.ZERO) .OR. (D*DX1.LT.ZERO) ) THEN
       IF (DIALOG) THEN
         WRITE(*,*)'Error: total arc length and initial/final spacing ',
     &             'cannot be of opposite signs.'
       ENDIF
       IEXIT = 2
      ENDIF
C
C    End spacings do not fit inside domain.
C
      IF ( (ABS(DX0).GT.ABS(D)) .OR. (ABS(DX1).GT.ABS(D)) ) THEN
       IF (DIALOG) THEN
         WRITE(*,*)'Error: specified end spacing is larger than total ',
     &             'arc length.'
       ENDIF
       IEXIT = 2
      ENDIF
C
C    N is 1 or 2.
C
      IF (N.EQ.1) THEN
         X(1)  = X0
         IEXIT = 1
      ELSE IF (N.EQ.2) THEN
         X(1)  = X0
         X(2)  = X1
         IEXIT = 1
      ENDIF
C
      IF (IEXIT.GT.0) THEN
       DX0A  = ZERO
       DX1A  = ZERO
      ENDIF
      IF ((IEXIT.EQ.2).AND.(DIALOG)) WRITE(*,*)'Program terminated.'
C
C    Uniform spacing.
C
      IF ( (DX0.EQ.ZERO) .AND. (DX1.EQ.ZERO) .AND. (N.GT.2) ) THEN
C
        DXU   = D/(N-1)
        X(1)  = X0
        X(N)  = X1
        DO 20 I=2,N-1
         X(I) = X(I-1) + DXU
 20     CONTINUE
        RMAX  = 1.0
        DX0A  = DXU
        DX1A  = DXU
        IEXIT = 1
C
      ENDIF
C
      RETURN
      END
C***********************************************************************
      SUBROUTINE STRRAT(N,RMAX,X)
C
#include "precis.h"
C
C    Compute max stretching ratio and return in RMAX.
C
      PARAMETER (ONE=1.0)
      DIMENSION X(N)
C
      IF (N.LT.3) THEN
C
         RMAX = 1.0
C
      ELSE
C
         R    = (X(3)-X(2))/(X(2)-X(1))
         RMAX = MAX( R, ONE/R )
         DO 10 I=3,N-1
            R    = (X(I+1)-X(I))/(X(I)-X(I-1))
            RMAX = MAX( RMAX, R, ONE/R )
 10      CONTINUE
C
      ENDIF
C
      RETURN
      END
C***********************************************************************
      SUBROUTINE MAXSP(N,DSMAX,X)
C
#include "precis.h"
C
C    Compute max spacing and return in DSMAX.
C
      DIMENSION X(N)
C
      IF (N.LT.2) THEN
C
         DSMAX = 0.0
C
      ELSE
C
         DSMAX = X(2)-X(1)
         DO 10 I=3,N
            DSMAX = MAX( DSMAX, X(I)-X(I-1) )
 10      CONTINUE
C
      ENDIF
C
      RETURN
      END
C***********************************************************************
      SUBROUTINE GEOMS(X0,X1,DX0,DX1,N,NMAX,RMAX,DX0A,DX1A,X,DIALOG)
C
#include "precis.h"
C
C  WMC (1/95)
C
C    This routine computes a geometric stretching in the range [X0,X1]
C    where the grid spacing at either X0 or X1 can be specified.
C
C    Input variables
C
C     X0   =   Value of arc length function at end point X0
C     X1   =   Value of arc length function at end point X1
C     DX0  =   Grid spacing at end point X0
C     DX1  =   Grid spacing at end point X1
C     N    > 0  Number of grid points including ends. RMAX is disregarded.
C          <=0  The smallest N such that RMAX is not exceeded will be
C               computed. RMAX>0 is assumed.
C     RMAX =   Max stretching ratio allowed
C     NMAX =   Max number of points used if input N is not positive
C              (also equals max dimension of X)
C
C    Returned variables
C
C     N     =  Number of grid points including ends,
C     DX0A  =  Analytic grid spacing at end point X0
C     DX1A  =  Analytic grid spacing at end point X1
C     X     =  Array returned with coordinates of points in the
C              interval [X0,X1]
C
C    Let D = X1-X0. There are 4 acceptable combinations of D,DX0,DX1.
C    In all 4 cases, X(1) = X0 and X(N) = X1.
C
C       D  DX0  DX1               Description
C    ===========================================================
C      >0   >0   =0    Grid spacing fixed at X0, floated at X1
C      >0   =0   >0    Grid spacing fixed at X1, floated at X0
C      <0   <0   =0    Grid spacing fixed at X0, floated at X1
C      <0   =0   <0    Grid spacing fixed at X1, floated at X0
C
      PARAMETER ( ZERO=0.0, ONE=1.0 )
      LOGICAL DIALOG
      DIMENSION X(NMAX)
C
      D    = X1 - X0
      X(1) = X0
C
      IF (DX0.NE.ZERO) THEN
C
C      Initial spacing specified
C
         CALL EPSIL(D,DX0,N,RMAX,EPS)
         IF ((N.GT.NMAX).AND.(DIALOG)) THEN
           WRITE(*,*)'Warning: no. of points needed > NMAX. Recompile.'
           STOP
         ENDIF
         R      = ONE + EPS
         X(N)   = X1
         DX     = DX0
         DO 30 I = 2,N-1
            X(I)   = X(I-1) + DX
            DX     = DX*R
   30       CONTINUE
         DX0A   = DX0
         DX1A   = X(N) - X(N-1)
C
      ELSE IF (DX1.NE.ZERO) THEN
C
C      Final spacing specified
C
         CALL EPSIL(D,DX1,N,RMAX,EPS)
         IF ((N.GT.NMAX).AND.(DIALOG)) THEN
           WRITE(*,*)'Warning: no. of points needed > NMAX. Recompile.'
           STOP
         ENDIF
         R      = ONE + EPS
         X(N)   = X1
         DX     = -DX1
         DO 40 I = N-1,2,-1
            X(I)   = X(I+1) + DX
            DX     = DX*R
   40       CONTINUE
         DX0A   = X(2) - X(1)
         DX1A   = DX1
C
      ENDIF
C
      IF (DIALOG) THEN
        WRITE(*,101) N
        WRITE(*,102) R
      ENDIF
 101  FORMAT(/,'Number of points used      = ',I4,/)
 102  FORMAT('Geometric stretching ratio = ',F11.5)
C
      RETURN
      END
C***********************************************************************
      SUBROUTINE EPSIL ( D,DX0,N,RMAX,EPS )
C
C   If the number of points N is given, this subroutine applies a
C   Newton-Raphson root-finding technique to find a value of epsilon
C   for a particular use of the geometric stretching transformation.
C   If the stretching ratio RMAX is given, the minimum number of points
C   N needed such that RMAX is not exceeded, and the corresponding
C   epsilon are returned.
C
C   This function is modified from the EPSIL function by Steger.
C   New checks are now performed for near uniform spacing (small EPS).
C   If EPS gets small, divides by a small number are avoided by
C   series expansions (WMC).
C
C   Example implementation:
C
C      CALL EPSIL(D,DX0,N,RMAX,EPS)
C      R   = 1. + EPS
C      DX  = DX0
C      X(1)= X0
C      DO 10 I= 2,N
C         X(I)= X(I-1) + DX
C         DX  = DX*R
C   10    CONTINUE
C
C     D     =   Total distance
C     DX0   =   Specified initial spacing
C     RMAX  =   Specified stretching ratio
C     N    > 0  Number of grid points including ends. RMAX is disregarded.
C          <=0  The smallest N such that RMAX is not exceeded will be
C               computed. RMAX>0 is assumed.
C     EPS   =   Returned value of epsilon
C
C     MITER  = Maximum number of iterations
C     FFRAC  = Iterative error bound
C              (fraction of smallest normalized spacing)
C     TOL    = Normalized tolerance for uniform spacing test
C     EPSMIN = Series expansion used if EPS smaller than EPSMIN
C
#include "precis.h"
C
      REAL*8 FFRAC,TOL,EPSMIN,ONE,TWTHRD
      REAL*8 DXL,DDX,ARG,FNM1,FNM2,FNM3,C0,C1,A,B,REM,R,RNM2,FDX,DEPS,F
      PARAMETER ( MITER=20, FFRAC=2.D-5 )
      PARAMETER ( TOL=1.D-6, EPSMIN=5.D-5 )
      PARAMETER ( ONE=1.D0, TWTHRD=2.D0/3.D0 )
C      PARAMETER ( MITER=20, FFRAC=2.E-5 )
C      PARAMETER ( TOL=1.E-6, EPSMIN=5.E-5 )
C      PARAMETER ( ONE=1., TWTHRD=2.D0/3.D0 )
C
      DXL    = ABS(DX0/D)
      DDX    = ONE/DXL
C
C    Compute minimum N if it is not positive
C
      IF (N.LE.0) THEN
       ARG = ONE - DDX*(ONE-RMAX)
       N   = 2 + INT( LOG(ARG)/LOG(RMAX) )
       IF (N.EQ.2) N = 3
      ENDIF
C
C    Compute EPS
C
      FNM1   = N-1
      FNM2   = N-2
      FNM3   = N-3
      C0     = 0.5*FNM1*FNM2
      C1     = DDX - FNM1
      A      = C1/C0
      B      = 1. + TWTHRD*A*FNM3
C
C    Set initial guess
C
      REM    = C1/DDX
      IF (ABS(REM).LT.TOL) THEN
         EPS    = 0.
      ELSE
         EPS    = DDX**(1./FNM2) - 1.
      ENDIF
C
C    Do Newton iterations
C
      DO 10 ITER = 1,MITER
         R      = 1. + EPS
         RNM2   = R**FNM2
         IF (ABS(EPS).LT.EPSMIN) THEN
            FDX    = C1 - C0*EPS
            DEPS   = A - B*EPS
         ELSE
            FDX    = DDX - (RNM2*R - 1.)/EPS
            ANUM   = FDX*EPS**2
            DENOM  = 1. + RNM2*(FNM2*EPS - 1.)
            DEPS   = ANUM/DENOM
         ENDIF
C
C    Beginning spacing is DXL, end spacing estimate is RNM2*DXL
C
         IF (ABS(FDX).LT.FFRAC*MIN(ONE,RNM2)) GOTO 20
         IF (ABS(DEPS).LT.TOL) GOTO 20
         EPS    = EPS + DEPS
   10    CONTINUE
C
      ITER   = MITER
      F      = FDX*DXL
      WRITE(*,*) 'EPSIL: Exceeded maximum number of iterations.'
      WRITE(*,*) '       Input : D,DX0,N    = ',D,DX0,N
      WRITE(*,*) '       Output: F,EPS,ITER = ',F,EPS,ITER
C
   20 CONTINUE
C
      RETURN
      END
C***********************************************************************
      SUBROUTINE TANHS(X0,X1,DX0,DX1,N,NMAX,RMAX,DX0A,DX1A,X,DIALOG)
C
#include "precis.h"
C
C    Refs: J.F. Thompson, Z.U.A. Warsi, and C.W. Mastin, NUMERICAL GRID
C          GENERATION, North-Holland, New York, pp. 307-308 (1985).
C
C          M. Vinokur, On One-Dimensional Stretching Functions for
C          Finite Difference Calculations, J. of Comp. Phys., Vol. 50,
C          pp. 215-234 (1983).
C
C  WMC (1/95)
C
C    This routine computes a hyperbolic tangent stretching in the range
C    [X0,X1] with either or both the initial and final spacings
C    specified. When an end point spacing is specified, the analytic
C    spacing is the specified spacing. When an end point spacing is
C    not specified, the analytic spacing is computed. In both cases,
C    the end point analytic spacings are returned.
C
C    This routine is derived from a combination of HYPTAN and VCLUST3.
C    Series expansions from VCLUST3 are used to form initial guesses
C    for Newton iterations used to invert y=sinh(x)/x and y=sin(x)/x.
C    Treatments for B<1, B>1 and B close to 1 are taken from Vinokur's
C    paper.
C
C    Input variables
C
C     X0   =   Value of arc length function at end point X0
C     X1   =   Value of arc length function at end point X1
C     DX0  =   Grid spacing at end point X0
C     DX1  =   Grid spacing at end point X1
C     N    > 0  Number of grid points including ends. RMAX is disregarded.
C          <=0  The smallest N such that RMAX is not exceeded will be
C               computed. RMAX>0 is assumed.
C     RMAX =   Max stretching ratio allowed
C     NMAX =   Max number of points used if input N is not positive
C              (also equals max dimension of X)
C
C    Returned variables
C
C     N     =  Number of grid points including ends,
C     DX0A  =  Analytic grid spacing at end point X0
C     DX1A  =  Analytic grid spacing at end point X1
C     X     =  Array returned with coordinates of points in the
C              interval [X0,X1]
C
C    Let D = X1-X0. There are 6 acceptable combinations of D,DX0,DX1.
C    In all 6 cases, X(1) = X0 and X(N) = X1.
C
C       D  DX0  DX1               Description
C    ===========================================================
C      >0   >0   =0    Grid spacing fixed at X0, floated at X1
C      >0   =0   >0    Grid spacing fixed at X1, floated at X0
C      >0   >0   >0    Grid spacing fixed at X0, fixed at X1
C      <0   <0   =0    Grid spacing fixed at X0, floated at X1
C      <0   =0   <0    Grid spacing fixed at X1, floated at X0
C      <0   <0   <0    Grid spacing fixed at X0, fixed at X1
C
C    Parameters
C
C     EPSB  = If B is within EPSB of 1.0, a special series expansion
C             that is to first order in (B-1) is used.
C
      PARAMETER ( EPSB=0.001 )
      PARAMETER ( ZERO=0.0, ONE=1.0, TWO=2.0, HALF=0.5 )
      LOGICAL DIALOG
C
      DIMENSION X(NMAX)
C
      XLEN  = X1 - X0
      DX0A  = DX0
      DX1A  = DX1
C
C    Scale grid spacings.
C
      DS0   = DX0/XLEN
      DS1   = DX1/XLEN
      DS    = ONE/(N-1)
C
C    Set parameters.
C
        IF      ( (DX0.NE.ZERO) .AND. (DX1.EQ.ZERO) ) THEN
C
         A     = ONE
         B     = DS/DS0
         U1    = ONE
         U2    = ONE
         U3    = HALF
         U4    = TWO
C
        ELSE IF ( (DX0.EQ.ZERO) .AND. (DX1.NE.ZERO) ) THEN
C
         A     = ONE
         B     = DS/DS1
         U1    = ONE
         U2    = ZERO
         U3    = HALF
         U4    = -ONE
C
        ELSE IF ( (DX0.NE.ZERO) .AND. (DX1.NE.ZERO) ) THEN
C
         A     = SQRT(DS1/DS0)
         B     = DS/SQRT(DS0*DS1)
         U1    = HALF
         U2    = ONE
         U3    = TWO
         U4    = HALF
C
        ENDIF
C
C    Generate stretching function.
C
        IF (N.GT.0) THEN
C
C        Generate stretching function with N points
C
          CALL TANHSG(X0,X1,XLEN,DX0,DX1,EPSB,DS,N,A,B,U1,U2,U3,U4,
     &                RMAXX,DX0A,DX1A,X)
C
        ELSE IF (N.LE.0) THEN
C
C        Generate stretching function with increasing no. of points
C        until max stretching ratio is less than or equal to RMAX.
C
          ROLD = 1.E+6

          DO 30 NT=4,NMAX
C
           DS     = ONE/(NT-1)
           IF      ( (DX0.NE.ZERO) .AND. (DX1.EQ.ZERO) ) THEN
            B     = DS/DS0
           ELSE IF ( (DX0.EQ.ZERO) .AND. (DX1.NE.ZERO) ) THEN
            B     = DS/DS1
           ELSE IF ( (DX0.NE.ZERO) .AND. (DX1.NE.ZERO) ) THEN
            B     = DS/SQRT(DS0*DS1)
           ENDIF
C
           CALL TANHSG(X0,X1,XLEN,DX0,DX1,EPSB,DS,NT,A,B,U1,U2,U3,U4,
     &                 RMAXX,DX0A,DX1A,X)
C
C          Exit loop if max stretching ratio is below RMAX
C
            IF (RMAXX.LE.RMAX) THEN
             N = NT
             IF (DIALOG) WRITE(*,101) N
             GO TO 50
            ENDIF
C
C          Exit loop if max stretching ratio starts increasing (for
C          two iterations in a row)
C
            IF (NT.GT.4) THEN
             IF (RMAXX.GT.ROLD .AND. ROLD.GT.ROLDER) THEN
              N      = NT-1
              DS     = ONE/(N-1)
              IF      ( (DX0.NE.ZERO) .AND. (DX1.EQ.ZERO) ) THEN
               B     = DS/DS0
              ELSE IF ( (DX0.EQ.ZERO) .AND. (DX1.NE.ZERO) ) THEN
               B     = DS/DS1
              ELSE IF ( (DX0.NE.ZERO) .AND. (DX1.NE.ZERO) ) THEN
               B     = DS/SQRT(DS0*DS1)
              ENDIF
              CALL TANHSG(X0,X1,XLEN,DX0,DX1,EPSB,DS,N,A,B,U1,U2,U3,U4,
     &                    RMAXX,DX0A,DX1A,X)
              IF (DIALOG) THEN
                WRITE(*,102) ROLD
                WRITE(*,101) N
              ENDIF
              GO TO 50
             ENDIF
            ENDIF
            ROLDER = ROLD
            ROLD = RMAXX
C
 30       CONTINUE
C
          N = NMAX
          IF (DIALOG) THEN
            WRITE(*,103) NMAX
            WRITE(*,104) RMAXX
          ENDIF
C
 50       CONTINUE
C
        ENDIF
C
C    Check for crossovers.
C
      DO 70 I=2,N
         DSI    = (X(I)-X(I-1))/XLEN
         IF ((DSI.LE.0.0).AND.(DIALOG)) THEN
            WRITE(*,*) 'Warning: zero or negative spacing starting at ',
     &                 'point ',I
            GO TO 90
         ENDIF
 70   CONTINUE
C
 90   CONTINUE
C
 100  CONTINUE
C
 101  FORMAT(/,'Number of points used = ',I4,/)
 102  FORMAT('Smallest max stretching ratio attainable =',F12.5)
 103  FORMAT('Max number of points reached =',I4)
 104  FORMAT('Final max stretching ratio   =',F12.5)
C
      RETURN
      END
C***********************************************************************
      SUBROUTINE TANHSG(X0,X1,XLEN,DX0,DX1,EPSB,DS,N,A,B,U1,U2,U3,U4,
     &                  RMAXX,DX0A,DX1A,X)
C
#include "precis.h"
C
      PARAMETER ( ZERO=0.0, ONE=1.0, TWO=2.0, HALF=0.5 )
      DIMENSION X(*)
C
C    Generate hyperbolic tangent arc length distribution for N number
C    of points and report maximum stretching ratio.
C
      ONEM  = ONE - EPSB
      ONEP  = ONE + EPSB
      X(1)  = X0
      X(N)  = X1
C
        IF      (B .LE. ONEM) THEN
C
          CALL ASINN(B,DELTA)
          HDELTA = HALF*DELTA
          TNH2   = TAN(HDELTA)
          DO 30 I=2,N-1
             XII     = (I-1)*DS
             X(I)    = U1*( U2 + TAN( HDELTA*(XII/U1-U2) )/TNH2 )
 30       CONTINUE
C
        ELSE IF (B .GE. ONEP) THEN
C
          CALL ASINHN(B,DELTA)
          HDELTA = HALF*DELTA
          TNH2   = TANH(HDELTA)
          DO 40 I=2,N-1
             XII     = (I-1)*DS
             X(I)    = U1*( U2 + TANH( HDELTA*(XII/U1-U2) )/TNH2 )
 40       CONTINUE
C
        ELSE
C
          UBM = U3*(ONE-B)
          DO 50 I=2,N-1
             XII     = (I-1)*DS
             X(I)    = XII*( ONE + UBM*(XII-U4)*(XII-ONE) )
 50       CONTINUE
C
        ENDIF
C
C      Rescale coordinates
C
        AM = ONE - A
        DO 60 I=2,N-1
           X(I) = X0 + XLEN*( X(I)/(A + AM*X(I)) )
 60     CONTINUE
C
C      Compute end point analytic spacing
C
        IF ( (B .LE. ONEM) .OR. (B .GE. ONEP) ) THEN
          DXAN = XLEN*DS*HDELTA/TNH2
          IF (DX1.EQ.ZERO) DX1A = DXAN
          IF (DX0.EQ.ZERO) DX0A = DXAN
        ELSE
          IF (DX1.EQ.ZERO) DX1A = XLEN*DS*(ONE + UBM*(ONE-U4))
          IF (DX0.EQ.ZERO) DX0A = XLEN*DS*(ONE + UBM*U4)
        ENDIF
C
C     Compute max stretching ratio
C
       IMAXX = 2
       DXA   = X(2)-X(1)
       DXB   = X(3)-X(2)
       RMAXX = MAX( DXB/DXA, DXA/DXB )
       DXA   = DXB
C
       DO 80 I=3,N-1
        DXB = X(I+1)-X(I)
        R   = MAX( DXB/DXA, DXA/DXB )
        IF (R.GT.RMAXX) THEN
         RMAXX = R
         IMAXX = I
        ENDIF
        DXA = DXB
 80    CONTINUE
C
C      WRITE(*,101) RMAXX,IMAXX
C101   FORMAT('Max stretching ratio = ',F10.4,' at I = ',I4)
C
      RETURN
      END
C***********************************************************************
      SUBROUTINE ASINHN(B,DELTA)
C
#include "precis.h"
C
      PARAMETER (A1= -0.15,
     &           A2=  0.0573214285714,
     &           A3= -0.024907294878,
     &           A4=  0.0077424460899,
     &           A5= -0.0010794122691)
      PARAMETER (C0= -0.0204176930892,
     &           C1=  0.2490272170591,
     &           C2=  1.9496443322775,
     &           C3= -2.629454725241,
     &           C4=  8.5679591096315)
      PARAMETER (B1=  2.7829681178603,
     &           B2= 35.0539798452776)
      PARAMETER ( ONE=1.0, TWO=2.0, SIX=6.0 )
      PARAMETER ( MITER=20, DELMIN=5.0E-5, TOL=1.0E-6 )
C
C    Use series expansions to get initial guess for delta where
C    sinh(delta)/delta = B
C    Then use Newton iterations to converge delta. Since B is never
C    below (1+epsb) in the calling sequence, delta should remain
C    well above DELMIN.
C
      IF (B.LE.B1) THEN
         BB    = B - ONE
         DELTA = SQRT(SIX*BB)
     &           *(((((A5*BB+A4)*BB+A3)*BB+A2)*BB+A1)*BB+ONE)
      ELSE
         V     = LOG(B)
         W     = ONE/B - ONE/B2
         DELTA = V + LOG(TWO*V)
     &           *(ONE+ONE/V)+(((C4*W+C3)*W+C2)*W+C1)*W+C0
      ENDIF
C
C    Newton iterations
C
      DO 10 I=1,MITER
       IF (ABS(DELTA).GE.DELMIN) THEN
        SDELTA = SINH(DELTA)
        CDELTA = COSH(DELTA)
        F  = SDELTA/DELTA - B
        FP = (DELTA*CDELTA - SDELTA)/(DELTA*DELTA)
        DD = -F/FP
       ELSE
        WRITE(*,*)'Delta fell below DELMIN.'
        STOP
       ENDIF
       IF (ABS(F).LT.TOL) GO TO 20
       DELTA = DELTA + DD
 10   CONTINUE
C
C      WRITE(*,*)'Exceeded max number of iterations.'
C      WRITE(*,*)'DELTA=',DELTA,'   F=',F,'   DDELTA=',DD
C
 20   CONTINUE
C
      RETURN
      END
C***********************************************************************
      SUBROUTINE ASINN(B,DELTA)
C
#include "precis.h"
C
      PARAMETER (A1=  0.15,
     &           A2=  0.0573214285714,
     &           A3=  0.0489742834696,
     &           A4= -0.053337753213,
     &           A5=  0.0758451335824)
      PARAMETER (C3= -2.6449340668482,
     &           C4=  6.7947319658321,
     &           C5=-13.2055008110734,
     &           C6= 11.7260952338351)
      PARAMETER (B1=  0.2693897165164)
      PARAMETER (PI=  3.14159265358981)
      PARAMETER ( ONE=1.0, SIX=6.0 )
      PARAMETER ( MITER=20, DELMIN=5.0E-5, TOL=1.0E-6 )
C
C    Use series expansions to get initial guess for delta where
C    sin(delta)/delta = B
C    Then use Newton iterations to converge delta. Since B is never
C    above (1-epsb) in the calling sequence, delta should remain
C    above DELMIN.
C
      IF (B.LE.B1) THEN
         DELTA = PI*((((((C6*B+C5)*B+C4)*B+C3)*B+ONE)*B-ONE)*B+ONE)
      ELSE
         BB    = ONE - B
         DELTA = SQRT(SIX*BB)
     &         *(((((A5*BB+A4)*BB+A3)*BB+A2)*BB+A1)*BB+ONE)
      ENDIF
C
C    Newton iterations
C
      DO 10 I=1,MITER
       IF (ABS(DELTA).GE.DELMIN) THEN
        SDELTA = SIN(DELTA)
        CDELTA = COS(DELTA)
        F  = SDELTA/DELTA - B
        FP = (DELTA*CDELTA - SDELTA)/(DELTA*DELTA)
        DD = -F/FP
       ELSE
        WRITE(*,*)'Delta fell below DELMIN.'
        STOP
       ENDIF
       IF (ABS(F).LT.TOL) GO TO 20
       DELTA = DELTA + DD
 10   CONTINUE
C
C      WRITE(*,*)'Exceeded max number of iterations.'
C      WRITE(*,*)'DELTA=',DELTA,'   F=',F,'   DDELTA=',DD
C
 20   CONTINUE
C
      RETURN
      END
C***********************************************************************
      SUBROUTINE GEOUNI(X0,X1,DX0,DX1,N,NMAX,RMAX,DX0A,DX1A,X,DIALOG)
C
#include "precis.h"
C
C    Geometric stretching up to end spacing, then uniform for the rest
C    of the distance. Input RMAX may be lowered to make points fit within
C    total distance. Fall back to hyperbolic tangent if end spacing
C    cannot be achieved with monotonic increase of spacing.
C
C    Input variables
C
C     X0   =   Value of arc length function at end point X0
C     X1   =   Value of arc length function at end point X1
C     DX0  =   Grid spacing at end point X0
C     DX1  =   Grid spacing at end point X1
C     RMAX =   Max stretching ratio allowed
C     NMAX =   Max number of points used if input N is not positive
C              (also equals max dimension of X)
C
C    Returned variables
C
C     N     =  Number of grid points including ends
C     DX0A  =  Analytic grid spacing at end point X0
C     DX1A  =  Analytic grid spacing at end point X1
C     X     =  Array returned with coordinates of points in the
C              interval [X0,X1]
C
C    Let D = X1-X0. There are 6 acceptable combinations of D,DX0,DX1.
C    In all 6 cases, X(1) = X0 and X(N) = X1.
C
C       D  DX0  DX1               Description
C    ===========================================================
C      >0   >0   =0    Grid spacing fixed at X0, floated at X1
C      >0   =0   >0    Grid spacing fixed at X1, floated at X0
C      >0   >0   >0    Grid spacing fixed at X0, fixed at X1
C      <0   <0   =0    Grid spacing fixed at X0, floated at X1
C      <0   =0   <0    Grid spacing fixed at X1, floated at X0
C      <0   <0   <0    Grid spacing fixed at X0, fixed at X1
C
C    Parameters
C
C     EPSB  = If B is within EPSB of 1.0, a special series expansion
C             that is to first order in (B-1) is used.
C
      PARAMETER ( ZERO=0.0, ONE=1.0, TWO=2.0, HALF=0.5 )
      LOGICAL DIALOG
      DIMENSION X(NMAX)
C
      N     = 0
      RMAXO = RMAX
C
C    Check ratio of DX1 to DX0.
C
      DXRAT = DX1/DX0
      RSQ   = RMAX*RMAX
C
C    First generate geometric stretching.
C
      D    = X1 - X0
      X(1) = X0
      CALL EPSIL(D,DX0,N,RMAX,EPS)
      IF ((N.GT.NMAX).AND.(DIALOG)) THEN
        WRITE(*,*)'Warning: no. of points needed > NMAX. Recompile.'
        STOP
      ENDIF
      R      = ONE + EPS
      X(N)   = X1
      DX     = DX0
      DO I = 2,N-1
         X(I)   = X(I-1) + DX
         DX     = DX*R
      ENDDO
      DX0A   = DX0
      DX1A   = X(N) - X(N-1)
C
C    Check end spacing.
C
      IF ( (DX1A.LT.DX1) .OR. (DXRAT.LT.RSQ) ) THEN
C
C      Cannot achieve end spacing with monotonic increase of spacing,
C      or end spacing ratio is less than RMAX.
C      Go to hyperbolic tangent stretching.
C
        N    = 0
        RMAX = RMAXO
        CALL TANHS(X0,X1,DX0,DX1,N,NMAX,RMAX,DX0A,DX1A,X,DIALOG)
        CALL MAXSP(N,DSMAX,X)
        RNEW = RMAX
        NTOP = MIN( INT(D/DX0) + 1, NMAX )
C
        DO WHILE (((DSMAX.GT.DX1).OR.(RNEW.GT.RMAX)) .AND. (N.LT.NTOP))
           N  = N+1
           RX = 0.0
           CALL TANHS(X0,X1,DX0,DX1,N,NMAX,RX,DX0A,DX1A,X,DIALOG)
           CALL MAXSP(N,DSMAX,X)
           CALL STRRAT(N,RNEW,X)
        ENDDO
C
      ELSE IF (DX1A.GT.DX1) THEN
C
C      Find index of point with spacing just less than DX1.
C
        IS = 2
        DO I=N-1,2,-1
           DELX = X(I)-X(I-1)
           IF (DELX.LT.DX1) THEN
              IS = I
              GO TO 30
           ENDIF
        ENDDO
 30     CONTINUE
C
C      Determine size of uniform spacing region and number of intervals.
C
        TDUNI = X(N)-X(IS)
        NUNI  = INT(TDUNI/DX1)
C
C      Determine size of modified geometric region and new number of points
C      for geometric region.
C
        TDGEO = D - (NUNI*DX1)
        NGEO  = IS + 1
C
C      Check if total no. of points exceed NMAX.
C
        N = NGEO + NUNI
        IF ((N.GT.NMAX).AND.(DIALOG)) THEN
          WRITE(*,*)'Warning: no. of points needed > NMAX. Recompile.'
          STOP
        ENDIF
C
C      Generate uniform region.
C
        X(N) = X1
        IU   = N - NUNI
        DO I=N-1,IU,-1
           X(I) = X(I+1) - DX1
        ENDDO
        DX1A = DX1
C
C      Re-generate geometric region.
C
        RGEO = 0.0
        CALL EPSIL(TDGEO,DX0,NGEO,RGEO,EPS)
        R      = ONE + EPS
        DX     = DX0
        DO I = 2,NGEO-1
           X(I)   = X(I-1) + DX
           DX     = DX*R
        ENDDO
C
      ENDIF
C
      RETURN
      END
C***********************************************************************
      SUBROUTINE HYPUNI(X0,X1,DX0,DX1,N,NMAX,RMAX,DX0A,DX1A,X,DIALOG)
C
#include "precis.h"
C
C    Input variables
C
C     X0   =   Value of arc length function at end point X0
C     X1   =   Value of arc length function at end point X1
C     DX0  =   Grid spacing at end point X0
C     DX1  =   Grid spacing at end point X1
C     RMAX =   Max stretching ratio allowed
C     NMAX =   Max number of points used if input N is not positive
C              (also equals max dimension of X)
C     DX0A  =  Max uniform spacing
C
C    Returned variables
C
C     N     =  Number of grid points including ends
C     DX0A  =  Analytic grid spacing at end point X0
C     DX1A  =  Analytic grid spacing at end point X1
C     X     =  Array returned with coordinates of points in the
C              interval [X0,X1]
C
      PARAMETER ( ZERO=0.0, ONE=1.0, TWO=2.0, HALF=0.5 )
      PARAMETER (NM=1000)
      LOGICAL DIALOG,INCN
      DIMENSION X(NMAX)
      DIMENSION XT(NM)
C
      DS = DX0A
      D  = X1-X0
C
C    First check for degenerate cases.
C
      DXT = DX0+DX1
      IF (DXT.GT.D) THEN
         N    = 2
         X(1) = X0
         X(2) = X1
         DX0A = D
         DX1A = D
         IDEGEN = 1
      ELSE IF (DXT.EQ.D) THEN
         N    = 3
         X(1) = X0
         X(2) = X0+DX0
         X(3) = X1
         DX0A = DX0
         DX1A = DX1
         IDEGEN = 1
      ELSE
         IDEGEN = 0
      ENDIF
C
C    Now treat non-degenerate cases.
C
      IF (IDEGEN.EQ.0) THEN
C
C       Compute modified stretching ratios at 0 and 1.
C
         R00 = DS/DX0
         N0  = 1
         R0  = R00
         DO WHILE (R0.GT.RMAX)
            N0 = N0 + 1
            R0 = R00**(1.0/FLOAT(N0))
         ENDDO
C
         R10 = DS/DX1
         N1  = 1
         R1  = R10
         DO WHILE (R1.GT.RMAX)
            N1 = N1 + 1
            R1 = R10**(1.0/FLOAT(N1))
         ENDDO
C
C       Modify N0 and N1 to number of points needed.
C
         N0 = N0 + 2
         N1 = N1 + 2
C
C       Compute distance needed to reach DS from X0 and X1.
C
         D0P = DX0*( ( 1.0 - R0**(FLOAT(N0-1)) ) / ( 1.0 - R0 ) )
         D1P = DX1*( ( 1.0 - R1**(FLOAT(N1-1)) ) / ( 1.0 - R1 ) )
C
C       Test for 2 cases.
C
         DPT = D0P + D1P
C
         IF (DPT.LE.D) THEN
C
C          Extra space available after reaching DS from both ends.
C
            DEXTRA = D-DPT
C
            NSEGU = INT(DEXTRA/DS)
            DADJ = DEXTRA - DS*NSEGU
C
            IF (DADJ.NE.ZERO) THEN
C
C             Modify geometric stretching distance from both ends.
C
               DMOD = 0.5*DADJ
               D0P  = D0P+DMOD
               D1P  = D1P+DMOD
C
            ENDIF
C
C          Compute geometric stretching from both ends.
C
            DM   = 0.0
            RM   = 0.0
            CALL GEOMS(X0,X0+D0P,DX0,DM,N0,NMAX,RM,DT0,DT1,X,DIALOG)
            DM   = 0.0
            RM   = 0.0
            CALL GEOMS(X1-D1P,X1,DM,DX1,N1,NM,RM,DT0,DT1,XT,DIALOG)
C
C          Compute N and check.
C
            N = N0 + NSEGU + N1 - 1
            IF ((N.GT.NMAX).AND.(DIALOG)) THEN
            WRITE(*,*)'Warning: no. of points needed > NMAX. Recompile.'
            STOP
            ENDIF
C
C          Load into X array.
C
            DO I=N0+1,N0+NSEGU
               X(I) = X(I-1) + DS
            ENDDO
C
            II = 1
            DO I=N0+NSEGU+1,N
               II = II + 1
               X(I) = XT(II)
            ENDDO
C
            DX0A = DX0
            DX1A = DX1
C
         ELSE IF (DPT.GT.D) THEN
C
C          Ran out of space before DS is reached from one or both ends.
C          Call hyperbolic tangent to get first guess.
C
           INCN = .TRUE.
           NN   = 0
           ICALL= 0
C
           DO WHILE (INCN)
C
              CALL TANHS(X0,X1,DX0,DX1,NN,NMAX,RMAX,DX0A,DX1A,X,DIALOG)
C
C            Store NN from first call.
C
              IF (ICALL.EQ.0) THEN
                 NN02  = 2*NN
                 ICALL = 1
              ENDIF
C
C            Check if DS is exceeded.
C
              NEX = 0
              DO I=1,NN-1
                 DEL = X(I+1)-X(I)
                 IF (DEL.GT.DS) NEX = NEX + 1
              ENDDO
C
              IF ((NEX.EQ.0).OR.(NN.EQ.NMAX).OR.(NN.GE.NN02)) THEN
                 INCN = .FALSE.
              ELSE
                 NN = NN + 1
              ENDIF
C
           ENDDO
C
           N = NN
C
         ENDIF
C
      ENDIF
C
      RETURN
      END
C***********************************************************************
      SUBROUTINE GEOEND(X0,X1,DX0,DX1,DXU,RMAX,NMAX,N,X)
C
#include "precis.h"
C
C    Generate stretching function in [X0,X1] given spacing at both ends,
C    a max spacing value, and a max stretching ratio. Use hyperbolic
C    tangent and increase number of points until max spacing value,
C    and max stretching ratio are not exceeded.
C
C    Input
C
C     X0   =    Value of arc length function at end point X0
C     X1   =    Value of arc length function at end point X1
C     DX0  =    Grid spacing at end point X0
C     DX1  =    Grid spacing at end point X1
C     DXU  =    Max grid spacing
C     RMAX =    Max stretching ratio allowed. This is assumed to be a
C               number greater than 1. If RMAX<1, 1/RMAX will be used.
C               (input RMAX is used only if input N <= 0)
C     NMAX =    Max number of points allowed
C               (also equals max dimension of X)
C
C    Returned variables
C
C     N     =  Number of grid points including ends
C     X     =  Array returned with coordinates of points in the
C              interval [X0,X1]
C     ISTAT =  1   status ok
C           =  0   status not ok
C
      PARAMETER (ZERO=0.0,ONE=1.0)
      LOGICAL DIALOG, NEXTN
      DIMENSION X(NMAX)
C
      DIALOG = .FALSE.
      ISTAT = 1
      DD = X1-X0
C
C    Check RMAX.
C
      IF (RMAX.LT.ZERO) THEN
         RMAX = -RMAX
      ELSE IF (RMAX.EQ.ZERO) THEN
         RMAX = 1.1
      ENDIF
      IF ((RMAX.GT.ZERO).AND.(RMAX.LT.ONE)) THEN
         RMAX = 1.0/RMAX
      ELSE IF (RMAX.EQ.ONE) THEN
         RMAX = 1.1
      ENDIF
C
C    Check different cases based on grid spacings.
C
      DXS = DX0+DX1
      DXA = DXS+DXU
C
      IF ((DX0.GE.DD).OR.(DX1.GE.DD).OR.(DXU.GE.DD).OR.(DXS.GT.DD)) THEN
C
C        Domain too small.
C
          N    = 2
          X(1) = X0
          X(2) = X1
C
      ELSE IF (DXS.EQ.DD) THEN
C
C        Domain large enough to fit end spacings only.
C
          N    = 3
          X(1) = X0
          X(2) = X0+DX0
          X(3) = X1
C
      ELSE IF ((DXS.LT.DD).AND.(DXA.GE.DD)) THEN
C
C        Domain large enough to fit end spacings plus a small middle section.
C
          N    = 4
          X(1) = X0
          X(2) = X0+DX0
          X(3) = X1-DX1
          X(4) = X1
C
      ELSE IF ((DX0.EQ.DX1).AND.(DX1.EQ.DXU)) THEN
C
C        Uniform spacing.
C
          N     = INT(DD/DXU) + 1
          IF (N.GT.NMAX) THEN
             WRITE(*,*)'N is too large. Pass larger NMAX into SFUNS.'
             STOP 'SFUNS'
          ENDIF
          DXNEW = DD/FLOAT(N-1)
          X(1)  = X0
          X(N)  = X1
          DO 5 I=2,N-1
             X(I) = X0 + FLOAT(I-1)*DXNEW
 5        CONTINUE
C
      ELSE
C
          NEXTN = .TRUE.
          N     = 4
C
          IF ((DXU.GE.DX0).AND.(DXU.GE.DX1)) THEN
C
             DO WHILE (NEXTN)
C
                CALL TANHS(X0,X1,DX0,DX1,N,NMAX,RMAX,DX0A,DX1A,X,DIALOG)
C
C              Check grid spacing.
C
                DSMAX = X(3)-X(2)
                DO 10 I=2,N-2
                   DS = X(I+1)-X(I)
                   IF (DS.GT.DSMAX) DSMAX = DS
 10             CONTINUE
C
                IF ((DSMAX.LE.DXU).OR.(N.EQ.NMAX)) THEN
                   NEXTN = .FALSE.
                ELSE
                       N = N+1
                ENDIF
C
             ENDDO
C
          ENDIF
C
      ENDIF
C
      RETURN
      END
C***********************************************************************
      SUBROUTINE GEOEND_O(X0,X1,DX0,DX1,DXU,RMAX,NMAX,N,X)
C
#include "precis.h"
C
C    Generate stretching function in [X0,X1] given spacing at both ends
C    and a mid-section value. If domain is sufficiently long, use
C    geometric from both ends with given max stretching ratio, and
C    uniform in middle; may have to re-adjust stretching ratio such that
C    mid-section fits an integer number of the given uniform spacing.
C    Otherwise, use geometric from both ends until the grids
C    points meet somewhere in the middle.
C
C    Input
C
C     X0   =    Value of arc length function at end point X0
C     X1   =    Value of arc length function at end point X1
C     DX0  =    Grid spacing at end point X0
C     DX1  =    Grid spacing at end point X1
C     RMAX =    Max stretching ratio allowed. This is assumed to be a
C               number greater than 1. If RMAX<1, 1/RMAX will be used.
C               (input RMAX is used only if input N <= 0)
C     NMAX =    Max number of points allowed
C               (also equals max dimension of X)
C
C    Returned variables
C
C     N     =  Number of grid points including ends
C     X     =  Array returned with coordinates of points in the
C              interval [X0,X1]
C     ISTAT =  1   status ok
C           =  0   status not ok
C
      PARAMETER (ZERO=0.0,ONE=1.0)
      LOGICAL DIALOG
      DIMENSION X(NMAX)
C
      DIALOG = .FALSE.
      ISTAT = 1
      DD = X1-X0
C
C    Check RMAX.
C
      IF (RMAX.LT.ZERO) THEN
         RMAX = -RMAX
      ELSE IF (RMAX.EQ.ZERO) THEN
         RMAX = 1.1
      ENDIF
      IF ((RMAX.GT.ZERO).AND.(RMAX.LT.ONE)) THEN
         RMAX = 1.0/RMAX
      ELSE IF (RMAX.EQ.ONE) THEN
         RMAX = 1.1
      ENDIF
C
C    Check different cases based on grid spacings.
C
      DXS = DX0+DX1
      DXA = DXS+DXU
C
      IF ((DX0.GE.DD).OR.(DX1.GE.DD).OR.(DXU.GE.DD).OR.(DXS.GE.DD)) THEN
C
C        Domain too small.
C
          N    = 2
          X(1) = X0
          X(2) = X1
C
      ELSE IF ((DXS.LT.DD).AND.(DXA.GE.DD)) THEN
C
C        Domain large enough to fit end spacings only.
C
          N    = 4
          X(1) = X0
          X(2) = X0+DX0
          X(3) = X1-DX1
          X(4) = X1
C
      ELSE IF ((DX0.EQ.DX1).AND.(DX1.EQ.DXU)) THEN
C
C        Uniform spacing.
C
          N     = INT(DD/DXU) + 1
          IF (N.GT.NMAX) THEN
             WRITE(*,*)'N is too large. Pass larger NMAX into SFUNS.'
             STOP 'SFUNS'
          ENDIF
          DXNEW = DD/FLOAT(N-1)
          X(1)  = X0
          X(N)  = X1
          DO 5 I=2,N-1
             X(I) = X0 + FLOAT(I-1)*DXNEW
 5        CONTINUE
C
      ELSE
C
C       Compute length needed to reach DXU from DX0 with RMAX.
C
         IF (DX0.NE.DXU) THEN
C
            SP = DX0
            D0 = SP
            N0 = 2
            IF (DXU.GT.DX0) THEN
               DO WHILE (SP.LT.DXU)
                  SP = SP*RMAX
                  D0 = D0 + SP
                  N0 = N0 + 1
               ENDDO
            ELSE IF (DXU.LT.DX0) THEN
               DO WHILE (SP.GT.DXU)
                  SP = SP/RMAX
                  D0 = D0 + SP
                  N0 = N0 + 1
               ENDDO
            ENDIF
C
            D0 = D0 - SP
            N0 = N0 - 1
C
         ENDIF
C
C       Compute length needed to reach DXU from DX1 with RMAX.
C
         IF (DX1.NE.DXU) THEN
C
            SP = DX1
            D1 = SP
            N1 = 2
            IF (DXU.GT.DX1) THEN
               DO WHILE (SP.LT.DXU)
                  SP = SP*RMAX
                  D1 = D1 + SP
                  N1 = N1 + 1
               ENDDO
            ELSE IF (DXU.LT.DX1) THEN
               DO WHILE (SP.GT.DXU)
                  SP = SP/RMAX
                  D1 = D1 + SP
                  N1 = N1 + 1
               ENDDO
            ENDIF
C
            D1 = D1 - SP
            N1 = N1 - 1
C
         ENDIF
C
        IF ((DX0.NE.DXU).AND.(DX1.EQ.DXU)) THEN
C
            DTOT = DD - DX1
C
            IF (D0.GT.DTOT) THEN
C
C             Takes more length than is available. Fall back to hyptan.
C
               N = 0
               CALL TANHS(X0,X1,DX0,DX1,N,NMAX,RMAX,DX0A,DX1A,X,DIALOG)
C
            ELSE
C
C             It fits. Determine gap, and re-do geometric stretching.
C
               GAP   = DTOT - D0
               NUINT = INT(GAP/DXU) + 1
               GAPN  = FLOAT(NUINT)*DXU
               D0N   = DTOT - GAPN
               X1N   = X0 + D0N
               DX1N  = 0.0
               R     = 0.0
               CALL GEOMS(X0,X1N,DX0,DX1N,N0,NMAX,R,DX0A,DX1A,X,DIALOG)
               N     = N0 + NUINT + 1
C               
C             Uniform spacing region.
C
               IS = N0 + 1
               IE = IS + NUINT - 1
               DO 10 I=IS,IE
                  X(I) = X(I-1) + DXU
 10            CONTINUE
               X(N) = X1
C
            ENDIF
C
        ELSE IF ((DX0.EQ.DXU).AND.(DX1.NE.DXU)) THEN
C
            DTOT = DD - DX0
C
            IF (D1.GT.DTOT) THEN
C
C             Takes more length than is available. Fall back to hyptan.
C
               N = 0
               CALL TANHS(X0,X1,DX0,DX1,N,NMAX,RMAX,DX0A,DX1A,X,DIALOG)
C
            ELSE
C
C             It fits. Determine gap, and re-do geometric stretching.
C
               GAP   = DTOT - D1
               NUINT = INT(GAP/DXU) + 1
               GAPN  = FLOAT(NUINT)*DXU
               D1N   = DTOT - GAPN
               X0N   = X1 - D1N
               DX0N  = 0.0
               R     = 0.0
               CALL GEOMS(X0N,X1,DX0N,DX1,N1,NMAX,R,DX0A,DX1A,X,DIALOG)
               N     = N1 + NUINT + 1
C               
C             Shift X array.
C
               IS = N-N1+1
               DO 20 I=N,IS,-1
                  X(I) = X(I-NUINT-1)
 20            CONTINUE
C
C             Uniform spacing region.
C
               DO 30 I=IS-1,2,-1
                  X(I) = X(I+1) - DXU
 30            CONTINUE
               X(1) = X0
C
            ENDIF
C
        ELSE IF ((DX0.NE.DXU).AND.(DX1.NE.DXU)) THEN
C
            DTOT = D0+D1
C
            IF (DTOT.GT.DD) THEN
C
C             Takes more length than is available. Fall back to hyptan.
C
               N = 0
               CALL TANHS(X0,X1,DX0,DX1,N,NMAX,RMAX,DX0A,DX1A,X,DIALOG)
C
            ELSE
C
C             It fits. Determine gap, and re-do geometric stretching.
C
               GAP   = DD - DTOT
               NUINT = INT(GAP/DXU) + 1
               GAPN  = FLOAT(NUINT)*DXU
               DEL   = 0.5*(GAPN-GAP)
               D0N   = D0-DEL
               D1N   = D1-DEL
               X0N   = X1 - D1N
               X1N   = X0 + D0N
C
               DX0N  = 0.0
               R     = 0.0
               CALL GEOMS(X0N,X1,DX0N,DX1,N1,NMAX,R,DX0A,DX1A,X,DIALOG)
C
               N     = N0 + N1 + NUINT - 1
C               
C             Shift X array.
C
               IS = N-N1+1
               DO 40 I=N,IS,-1
                  X(I) = X(I-NUINT-N0+1)
 40            CONTINUE
C
C             Re-do geometric stretching at X0.
C
               DX1N  = 0.0
               R     = 0.0
               CALL GEOMS(X0,X1N,DX0,DX1N,N0,NMAX,R,DX0A,DX1A,X,DIALOG)
C
C             Fill uniform region.
C
               IS = N0+1
               IE = N-N1
               DO 50 I=IS,IE
                  X(I) = X(I-1) + DXU
 50            CONTINUE
C
            ENDIF
C
         ENDIF
C
      ENDIF
C
      RETURN
      END
