PROGRAM set3d

!*************************************************************************************!
!
! Convective Flow and Min/Max Level Set Solver
!
! Author: Oisin Tong
!
! Version: 3.0  
!
! Date: 23/5/2016
!
! Type: H-J WENO 5 Serial Explicit TVD R-K Level Set Solver
!
!*************************************************************************************!
USE set_subs

IMPLICIT NONE

!*************************************************************************************!
! Set Data
!*************************************************************************************!

INTEGER :: j,sUnit,nbytePhi,offset,nxs,n1,n2,n3,fN,im,jm,km,ip,jp,kp,dd
INTEGER :: iter,nx,ny,nz,qG,nn,counter,orderUp,order1,order2,reIter
REAL :: a0,a1,a2,a3,a4,t,xLo(3),xHi(3),xs,ys,yp,ypp,g,gp,x0,x1,dxs,dPlus,dMinus
REAL :: pX1,pX2,pX3,pY1,pY2,pY3,pZ1,pZ2,pZ3,gX,gY,gZ,dis,minD,k1,h1,k2,k3,k4,t5
REAL :: B1,B2,B3,C1,C2,C3,pSx,pSy,pSz,pS,pX,pY,pZ,sgn,ddx,ddy,ddz,gM,gMM,phiErrS
REAL :: y1,z1,maxX,maxY,maxZ,minX,minY,minZ,h,a,b,c,d,e,f,phiErr,gradMag2,x,y,z,dxx
REAL :: dx,cX,cY,cZ,t1,t2,t3,t4,bx,phiX,phiY,phiZ,phiXX,phiYY,phiZZ,phiXZ,phiXY,phiYZ
REAL :: pi,aa,bb,period,ax,ay,az,axy,ayz,axz,bxy,by,bz,byz,bxz,cxy,cyz,cxz,pp,CFL
REAL :: aaa,bbb,ccc,ax1,ax2,ay1,ay2,az1,az2,Dijk,Fx,Fy,Fz,rTime
REAL,ALLOCATABLE,DIMENSION(:,:,:) :: phi,phiS,phiN,gradPhiMag,phiE,phi1,phi2,phi3,phiO,S
INTEGER,ALLOCATABLE,DIMENSION(:,:,:) :: phiSB,phiNB
REAL,ALLOCATABLE,DIMENSION(:,:) :: centroid,surfX
CHARACTER(LEN=1) :: lf=char(10)
CHARACTER(LEN=1024) :: extent,origin,spacing,coffset
INTEGER*4,ALLOCATABLE,DIMENSION(:,:) :: surfElem
CHARACTER header*80,filename*80
INTEGER*2 padding
INTEGER*4 ntri,iunit,nSurfNode,k,i,n,p,kk,share,nSurfElem
REAL*4,ALLOCATABLE,DIMENSION(:,:) :: normals,triangles,nodesT
REAL,ALLOCATABLE,DIMENSION(:,:,:,:) :: gridX,grad2Phi,gradMixPhi
REAL,ALLOCATABLE,DIMENSION(:,:,:,:,:) :: gradPhi
INTEGER iargc
REAL,ALLOCATABLE,DIMENSION(:,:) :: q1S,q2S,q3S,q4S,q5S,q6S,q7S,q8S
INTEGER,DIMENSION(8) :: qq


!*************************************************************************************!
! Import STL Data
!*************************************************************************************!

! start system time
CALL cpu_time(t1)

! import .stl data
call getarg(1,filename)

PRINT*,
PRINT*, " Reading in .stl Mesh "
PRINT*,

iunit=13
OPEN(unit=iunit,file=filename,status='old',access='stream',form='unformatted')

! read .stl header info 
READ(iunit) header
READ(iunit) ntri
   
ALLOCATE(normals(3,ntri))
ALLOCATE(triangles(3,ntri*3))
ALLOCATE(surfELem(ntri,3))
 
! read .stl data
k=1
DO i = 1,ntri
   READ(iunit) normals(1,i),normals(2,i),normals(3,i)
   READ(iunit) triangles(1,k),triangles(2,k),triangles(3,k)
   READ(iunit) triangles(1,k+1),triangles(2,k+1),triangles(3,k+1)
   READ(iunit) triangles(1,k+2),triangles(2,k+2),triangles(3,k+2)
   READ(iunit) padding
  k=k+3
END DO
  
CLOSE(iunit)

ALLOCATE(nodesT(3,ntri*5))
nSurfElem = ntri

! search through data and put into surfX and surfElem style arrays
DO k = 1,ntri
  nodesT(1,k) = 1000000. 
  nodesT(2,k) = 1000000. 
  nodesT(3,k) = 1000000. 
END DO

! eliminate repeated nodes and clean up arrays
i = 1
nSurfNode = 3
k = 0;
DO n = 1,ntri
   DO p = 1,3
      share = 0 
      DO kk = 1,nSurfNode
         IF ((abs(nodesT(1,kk) - triangles(1,i)) < 1.e-13) .AND. &
             (abs(nodesT(2,kk) - triangles(2,i)) < 1.e-13) .AND. &
             (abs(nodesT(3,kk) - triangles(3,i)) < 1.e-13)) THEN
            share = kk
            EXIT
         END IF
      END DO
      IF (share > 0) THEN
         surfElem(n,p) = share
      ELSE
         k             = k+1 
         nodesT(:,k)   = triangles(:,i)
         surfElem(n,p) = k !1-based
      END IF
      i = i+1
   END DO
   nSurfNode = k 
END DO

! allocate surfX
ALLOCATE(surfX(nSurfNode,3))

! fill in surface node data 
DO k = 1,nSurfNode
   surfX(k,1) = nodesT(1,k)
   surfX(k,2) = nodesT(2,k)
   surfX(k,3) = nodesT(3,k)
END DO

! deallocate unnecessary data
DEALLOCATE(nodesT)
DEALLOCATE(triangles)
DEALLOCATE(normals)


!*************************************************************************************!
! Determine xLo and xHi
!*************************************************************************************!

! initialize
x1 = surfX(1,1);
y1 = surfX(1,2);
z1 = surfX(1,3);

maxX = x1;
maxY = y1;
maxZ = z1;

minX = x1;
minY = y1;
minZ = z1;

! find the max and min
DO n = 2,nSurfNode 
   x1 = surfX(n,1)
   y1 = surfX(n,2)
   z1 = surfX(n,3)

   IF (x1 > maxX) THEN
      maxX = x1
   END IF
   IF (y1 > maxY) THEN
     maxY = y1
   END IF
   IF (z1 > maxZ) THEN
      maxZ = z1
   END IF

   IF (x1 < minX) THEN
      minX = x1
   END IF
   IF (y1 < minY) THEN
     minY = y1
   END IF
   IF (z1 < minZ) THEN
      minZ = z1
   END IF
   
END DO

!*************************************************************************************!
! Define Cartesiang Grid Size and Allocate Phi
!*************************************************************************************!

! find the characteristic size of the object to add around
ddx = maxX-minX
ddy = maxY-minY
ddz = maxZ-minZ

! set dx
dx = 0.1

! define Cartesian grid
nx = ceiling((maxX-minX)/dx)+1;
ny = ceiling((maxY-minY)/dx)+1;
nz = ceiling((maxZ-minZ)/dx)+1;

! number of cells you want to add
dd = 10.

! adding more cells edge
nx = nx+2*dd
ny = ny+2*dd
nz = nz+2*dd

! set xLo and xHi
xLo =(/minX-dd*dx,minY-dd*dx,minZ-dd*dx/)
xHi =(/maxX+dd*dx,maxY+dd*dx,maxZ+dd*dx/)

! allocate phi
ALLOCATE(phi(0:nx,0:ny,0:nz))
phi = 1.

! allocate a grid of x,y,z points 
ALLOCATE(gridX(0:nx,0:ny,0:nz,3))
DO i = 0,nx
   DO j = 0,ny
      DO k = 0,nz
         gridX(i,j,k,1) = xLo(1) + i*dx;
         gridX(i,j,k,2) = xLo(2) + j*dx;
         gridX(i,j,k,3) = xLo(3) + k*dx;
      END DO
   END DO
END DO

!*************************************************************************************!
! Determine Inside and Outside of Surface
!*************************************************************************************!

! cut down on excess and use minimum nodes
im = floor((minX-xLo(1))/dx)-3
jm = floor((minY-xLo(2))/dx)-3
km = floor((minZ-xLo(3))/dx)-3
! cut down on excess and use maximum nodes
ip = floor((maxX-xLo(1))/dx)+3
jp = floor((maxY-xLo(2))/dx)+3
kp = floor((maxZ-xLo(3))/dx)+3

! print out grid spacing
PRINT*, " Setting Grid Size "
PRINT*, " Grid Size: nx =",nx,", ny =",ny,",nz =",nz 
PRINT*, " Grid Spacing: dx =", dx
PRINT*,
PRINT*, " Determining Inside and Outside of Geometry "
PRINT*, 

! allocate centroid
ALLOCATE(centroid(nSurfElem,4))

DO n = 1,nSurfElem
   n1 = surfElem(n,1)
   n2 = surfElem(n,2)
   n3 = surfElem(n,3)
   pX1 = surfX(n1,1)
   pY1 = surfX(n1,2)
   pZ1 = surfX(n1,3)
   pX2 = surfX(n2,1)
   pY2 = surfX(n2,2)
   pZ2 = surfX(n2,3)
   pX3 = surfX(n3,1)
   pY3 = surfX(n3,2)
   pZ3 = surfX(n3,3)
   centroid(n,1) = (pX1+pX2+pX3)/3.
   centroid(n,2) = (pY1+pY2+pY3)/3.
   centroid(n,3) = (pZ1+pZ2+pZ3)/3.
END DO

! find which nodes are inside and outside
DO i = im,ip 
   DO j = jm,jp  
      DO k = km,kp  

         ! search through all surface elements to find closest element
         minD = 100000.;
         DO n = 1,nSurfElem
            pX = centroid(n,1)
            pY = centroid(n,2)
            pZ = centroid(n,3)
            gX = gridX(i,j,k,1)
            gY = gridX(i,j,k,2)
            gZ = gridX(i,j,k,3)
            dis = sqrt((pX-gX)*(pX-gX) + (pY-gY)*(pY-gY) + (pZ-gZ)*(pZ-gZ))
            IF (dis < minD) THEN
               minD = dis;
               fN   = n;
            END IF
         END DO

      ! create three vectors from our point to the surfElem points
      n1 = surfElem(fN,1)
      n2 = surfElem(fN,2)
      n3 = surfElem(fN,3)
      A1 = surfX(n1,1) - gridX(i,j,k,1)
      A2 = surfX(n1,2) - gridX(i,j,k,2)
      A3 = surfX(n1,3) - gridX(i,j,k,3)
      B1 = surfX(n2,1) - gridX(i,j,k,1)
      B2 = surfX(n2,2) - gridX(i,j,k,2)
      B3 = surfX(n2,3) - gridX(i,j,k,3)
      C1 = surfX(n3,1) - gridX(i,j,k,1)
      C2 = surfX(n3,2) - gridX(i,j,k,2)
      C3 = surfX(n3,3) - gridX(i,j,k,3)
    
      ! cross product two of the vectors
      pSx = A2*B3-A3*B2;
      pSy = -(A1*B3-B1*A3);
      pSz = A1*B2-B1*A2;

      ! dot product last vector with your cross product result
      pS =-(pSx*C1+pSy*C2+pSz*C3);
      
      gM = 1.
      ! return sign of the dot product
      CALL phiSign(pS,sgn,dx,gM)
              
      phi(i,j,k) = sgn

      END DO
   END DO
END DO


! deallocate surface mesh data
DEALLOCATE(surfX)
!DEALLOCATE(gridX)
DEALLOCATE(surfElem)
DEALLOCATE(centroid)

CALL cpu_time(t2)
PRINT*, " Search Run Time: ",t2-t1," Seconds"
PRINT*,

!*************************************************************************************!
! Fast Marching Method
!*************************************************************************************!

! allocate arrays used for FMM
ALLOCATE(phiS(0:nx,0:ny,0:nz))
ALLOCATE(phiO(0:nx,0:ny,0:nz))
ALLOCATE(phiN(0:nx,0:ny,0:nz))
ALLOCATE(gradPhiMag(0:nx,0:ny,0:nz))
ALLOCATE(phiE(0:nx,0:ny,0:nz))
ALLOCATE(phi1(0:nx,0:ny,0:nz))
ALLOCATE(phi2(0:nx,0:ny,0:nz))
ALLOCATE(phi3(0:nx,0:ny,0:nz))
ALLOCATE(S(0:nx,0:ny,0:nz))
ALLOCATE(gradPhi(0:nx,0:ny,0:nz,3,2))

PRINT*, " Level Set Time Integration "
PRINT*, 

! set the phi sign array
phiS = phi

! number of iterations
iter = 10000 !1500 ! 10000

! normalized dx
dxx = dx/sqrt(ddx*ddx+ddy*ddy+ddz*ddz)

! time step
CFL = .1
h = CFL*dxx


CALL reinit(phi,gradPhi,gradPhiMag,nx,ny,nz,iter,dx,h)

! set original phi
phiO = phi

! print out run time
CALL cpu_time(t3)
PRINT*, " Initialization Run Time: ",t3-t1," Seconds"
PRINT*,


!*************************************************************************************!
! Paraview Output
!*************************************************************************************!

PRINT*, " Writing Out Cartesian Grid to Paraview Format "
PRINT*, 

! output to Paraview
WRITE(extent,'(3(A3,I6))')' 0 ',nx,' 0 ',ny,' 0 ',nz
WRITE(origin,'(3(F20.8,A1))')xLo(1),' ',xLo(2),' ',xLo(3),' '
WRITE(spacing,'(3(F20.8,A1))')dx,' ',dx,' ',dx,' '
nbytePhi =(nx+1)**3*24
offset = 0
WRITE(coffset,'(I16)')offset


sUnit = 11
OPEN(UNIT=sUnit,FILE='signedDistanceFunction.vti',FORM='unformatted',ACCESS='stream',STATUS='replace')
WRITE(sUnit)'<?xml version="1.0"?>'//lf
WRITE(sUnit)'<VTKFile type="ImageData" version="0.1" byte_order="LittleEndian">'//lf
WRITE(sUnit)'<ImageData WholeExtent="',TRIM(extent),'" Origin="',TRIM(origin),'" Spacing="',TRIM(spacing),'">'//lf
WRITE(sUnit)'<Piece Extent="',TRIM(extent),'">'//lf
WRITE(sUnit)'<PointData Scalars="phi">'//lf
WRITE(sUnit)'<DataArray type="Float64" Name="phi" format="appended" offset="',TRIM(coffset),'"/>'//lf
WRITE(sUnit)'</PointData>'//lf
WRITE(sUnit)'</Piece>'//lf
WRITE(sUnit)'</ImageData>'//lf
WRITE(sUnit)'<AppendedData encoding="raw">'//lf
WRITE(sUnit)'_'
WRITE(sUnit)nbytePhi,(((phi(i,j,k),i=0,nx),j=0,ny),k=0,nz)
WRITE(sUnit)lf//'</AppendedData>'//lf
WRITE(sUnit)'</VTKFile>'//lf
CLOSE(sUnit)

!*************************************************************************************!
! Determine Narrow Band
!*************************************************************************************!

ALLOCATE(phiSB(0:nx,0:ny,0:nz))
ALLOCATE(phiNB(0:nx,0:ny,0:nz))

CALL narrowBand(nx,ny,nz,dx,phi,phiNB,phiSB)

!*************************************************************************************!
! Initialize Gradients
!*************************************************************************************!


ALLOCATE(grad2Phi(0:nx,0:ny,0:nz,3))
ALLOCATE(gradMixPhi(0:nx,0:ny,0:nz,3))

gradPhi = 0.
grad2Phi = 0.
gradMixPhi = 0.
gradPhiMag = 0.

phiN = phi

orderUp = 1
order1 = 2
order2 = 2

!*************************************************************************************!
! Min/Max Flow
!*************************************************************************************!

iter = 2000 !20000
CFL = .01
h1 = CFL*dxx

DO n = 1,iter

   !****************** Explicit Third Order TVD RK Stage 1 ***********************!

   ! Calculate second derivative flow if it is in the narrow band
   DO i = 0,nx
      DO j = 0,ny
         DO k = 0,nz
            IF (phiNB(i,j,k) == 1) THEN
               CALL secondDeriv(i,j,k,nx,ny,nz,dx,phi,phiXX,phiYY,phiZZ,phiXY,phiXZ,phiYZ,order2)
               grad2Phi(i,j,k,1) = phiXX
               grad2Phi(i,j,k,2) = phiYY
               grad2Phi(i,j,k,3) = phiZZ
               gradMixPhi(i,j,k,1) = phiXY
               gradMixPhi(i,j,k,2) = phiXZ
               gradMixPhi(i,j,k,3) = phiYZ          
            END IF

         END DO
      END DO
   END DO

   ! Caluclate the min/max flow
   DO i = 0,nx
      DO j = 0,ny
         DO k = 0,nz
            IF (phiNB(i,j,k) == 1) THEN 
               CALL minMax(i,j,k,nx,ny,nz,dx,phi,grad2Phi,gradMixPhi,F,gridX)
               k1 = F  
               phi1(i,j,k) = phi(i,j,k) + h1*k1
            END IF

         END DO
      END DO
   END DO

   !****************** Explicit Third Order TVD RK Stage 2 ***********************!

   ! Calculate second derivative flow if it is in the narrow band
   DO i = 0,nx
      DO j = 0,ny
         DO k = 0,nz
            IF (phiNB(i,j,k) == 1) THEN
               CALL secondDeriv(i,j,k,nx,ny,nz,dx,phi1,phiXX,phiYY,phiZZ,phiXY,phiXZ,phiYZ,order2)
               grad2Phi(i,j,k,1) = phiXX
               grad2Phi(i,j,k,2) = phiYY
               grad2Phi(i,j,k,3) = phiZZ
               gradMixPhi(i,j,k,1) = phiXY
               gradMixPhi(i,j,k,2) = phiXZ
               gradMixPhi(i,j,k,3) = phiYZ          
            END IF

         END DO
      END DO
   END DO

   ! Caluclate the min/max flow
   DO i = 0,nx
      DO j = 0,ny
         DO k = 0,nz
            IF (phiNB(i,j,k) == 1) THEN  
               CALL minMax(i,j,k,nx,ny,nz,dx,phi1,grad2Phi,gradMixPhi,F,gridX)
               k2 = F   
               phi2(i,j,k) = 3./4.*phi(i,j,k) + 1./4.*phi1(i,j,k) + 1./4.*h1*k2
            END IF

         END DO
      END DO
   END DO

   !****************** Explicit Third Order TVD RK Stage 3 ***********************!

   ! Calculate second derivative flow if it is in the narrow band
   DO i = 0,nx
      DO j = 0,ny
         DO k = 0,nz
            IF (phiNB(i,j,k) == 1) THEN
               CALL secondDeriv(i,j,k,nx,ny,nz,dx,phi2,phiXX,phiYY,phiZZ,phiXY,phiXZ,phiYZ,order2)
               grad2Phi(i,j,k,1) = phiXX
               grad2Phi(i,j,k,2) = phiYY
               grad2Phi(i,j,k,3) = phiZZ
               gradMixPhi(i,j,k,1) = phiXY
               gradMixPhi(i,j,k,2) = phiXZ
               gradMixPhi(i,j,k,3) = phiYZ          
            END IF

         END DO
      END DO
   END DO

   ! Caluclate the min/max flow
   DO i = 0,nx
      DO j = 0,ny
         DO k = 0,nz
            IF (phiNB(i,j,k) == 1) THEN  
               CALL minMax(i,j,k,nx,ny,nz,dx,phi2,grad2Phi,gradMixPhi,F,gridX)
               k3 = F 
               phiN(i,j,k) = 1./3.*phi(i,j,k) + 2./3.*phi2(i,j,k) + 2./3.*h1*k3 
            END IF

         END DO
      END DO
   END DO

   !********************************* RMS ***************************************!

   phiErr = 0.

   ! calculate RMS
   DO i = 0,nx
      DO j = 0,ny
         DO k = 0,nz            
            phiErr = phiErr + (phi(i,j,k)-phiN(i,j,k))*(phi(i,j,k)-phiN(i,j,k))
         END DO
      END DO
   END DO
     
   ! check error
   phiErr = sqrt(phiErr/(nx*ny*nz))
   IF (phiErr < 1.E-7) THEN
      PRINT*, " Min/max time integration has reached steady state "
      EXIT
   END IF
   
   ! set phi to new value
   phi = phiN
   
   PRINT*, " Iteration: ",n," ", " RMS Error: ",phiErr
  
   ! check for a NAN
   IF (isnan(phiErr)) STOP

   CALL narrowBand(nx,ny,nz,dx,phi,phiNB,phiSB)

END DO
print*,
   

!*************************************************************************************!
! Asymptotic Error
!*************************************************************************************!

phiErr = 0.

! calculate RMS
DO i = 0,nx
   DO j = 0,ny
      DO k = 0,nz
         phiErr = phiErr + (phi(i,j,k)-phiO(i,j,k))*(phi(i,j,k)-phiO(i,j,k))
      END DO
   END DO
END DO

phiErr = sqrt(phiErr/(nx*ny*nz))

PRINT*, " Asymptotic Error: ",phiErr

!*************************************************************************************!
! Output Grad Phi Mag
!*************************************************************************************!

! grad phi
order1 = 2
DO i = 0,nx
   DO j = 0,ny
      DO k = 0,nz
         CALL firstDeriv(i,j,k,nx,ny,nz,dx,phi,phiX,phiY,phiZ,order1,gMM)
         gradPhiMag(i,j,k) = gMM
      END DO
   END DO
END DO

!*************************************************************************************!
! Paraview Output
!*************************************************************************************!

PRINT*, " Writing Out Cartesian Grid to Paraview Format "
PRINT*, 

! output to Paraview
WRITE(extent,'(3(A3,I6))')' 0 ',nx,' 0 ',ny,' 0 ',nz
WRITE(origin,'(3(F20.8,A1))')xLo(1),' ',xLo(2),' ',xLo(3),' '
WRITE(spacing,'(3(F20.8,A1))')dx,' ',dx,' ',dx,' '
nbytePhi =(nx+1)**3*24
offset = 0
WRITE(coffset,'(I16)')offset

sUnit = 11
OPEN(UNIT=sUnit,FILE='smoothedDistanceFunction.vti',FORM='unformatted',ACCESS='stream',STATUS='replace')
WRITE(sUnit)'<?xml version="1.0"?>'//lf
WRITE(sUnit)'<VTKFile type="ImageData" version="0.1" byte_order="LittleEndian">'//lf
WRITE(sUnit)'<ImageData WholeExtent="',TRIM(extent),'" Origin="',TRIM(origin),'" Spacing="',TRIM(spacing),'">'//lf
WRITE(sUnit)'<Piece Extent="',TRIM(extent),'">'//lf
WRITE(sUnit)'<PointData Scalars="phi">'//lf
WRITE(sUnit)'<DataArray type="Float64" Name="phi" format="appended" offset="',TRIM(coffset),'"/>'//lf
WRITE(sUnit)'</PointData>'//lf
WRITE(sUnit)'</Piece>'//lf
WRITE(sUnit)'</ImageData>'//lf
WRITE(sUnit)'<AppendedData encoding="raw">'//lf
WRITE(sUnit)'_'
WRITE(sUnit)nbytePhi,(((phi(i,j,k),i=0,nx),j=0,ny),k=0,nz)
WRITE(sUnit)lf//'</AppendedData>'//lf
WRITE(sUnit)'</VTKFile>'//lf
CLOSE(sUnit)

!*************************************************************************************!
! Reinitialise
!*************************************************************************************!

! number of iterations
iter = 2000

! time step
CFL = .001
h = CFL*dxx

CALL reinit(phi,gradPhi,gradPhiMag,nx,ny,nz,iter,dx,h)


!*************************************************************************************!
! Paraview Output
!*************************************************************************************!

PRINT*, " Writing Out Cartesian Grid to Paraview Format "
PRINT*, 

! output to Paraview
WRITE(extent,'(3(A3,I6))')' 0 ',nx,' 0 ',ny,' 0 ',nz
WRITE(origin,'(3(F20.8,A1))')xLo(1),' ',xLo(2),' ',xLo(3),' '
WRITE(spacing,'(3(F20.8,A1))')dx,' ',dx,' ',dx,' '
nbytePhi =(nx+1)**3*24
offset = 0
WRITE(coffset,'(I16)')offset

sUnit = 11
OPEN(UNIT=sUnit,FILE='convectedDistanceFunction.vti',FORM='unformatted',ACCESS='stream',STATUS='replace')
WRITE(sUnit)'<?xml version="1.0"?>'//lf
WRITE(sUnit)'<VTKFile type="ImageData" version="0.1" byte_order="LittleEndian">'//lf
WRITE(sUnit)'<ImageData WholeExtent="',TRIM(extent),'" Origin="',TRIM(origin),'" Spacing="',TRIM(spacing),'">'//lf
WRITE(sUnit)'<Piece Extent="',TRIM(extent),'">'//lf
WRITE(sUnit)'<PointData Scalars="phi">'//lf
WRITE(sUnit)'<DataArray type="Float64" Name="phi" format="appended" offset="',TRIM(coffset),'"/>'//lf
WRITE(sUnit)'</PointData>'//lf
WRITE(sUnit)'</Piece>'//lf
WRITE(sUnit)'</ImageData>'//lf
WRITE(sUnit)'<AppendedData encoding="raw">'//lf
WRITE(sUnit)'_'
WRITE(sUnit)nbytePhi,(((phi(i,j,k),i=0,nx),j=0,ny),k=0,nz)
WRITE(sUnit)lf//'</AppendedData>'//lf
WRITE(sUnit)'</VTKFile>'//lf
CLOSE(sUnit)

!*************************************************************************************!
! Total Run Time
!*************************************************************************************!

! print out run time
CALL cpu_time(t4)
PRINT*, " Total Run Time: ",t4-t1," Seconds"
PRINT*,

!*************************************************************************************!
!
! Program End
!
!*************************************************************************************!

END PROGRAM set3d
