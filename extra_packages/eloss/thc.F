c--
c--   THC: Transverse Holding Coil
c--
      SUBROUTINE THC(vertex, cdir, pathLength)
      IMPLICIT NONE
      SAVE

      REAL vertex(3)
      REAL cdir(3)
      REAL center(3)
      DATA center /0.0, 0.0, 0.0/
      REAL a1(3), a2(3)
      REAL pathLength

      REAL RC(22)               ! radius of cylinders
      REAL ZPOS(26)             ! z-positions of cylinders
      REAL PHIPOS(1)            ! NEW! angular data for holding coil

c--   ALL DISTANCES ARE IN CENTIMETERS
c--   RC(21): NEW! Transverse holding coil shell inner radius
c--   RC(22): NEW! Transverse holding coil shell outer radius
      DATA RC/0.762, 0.8255, 1.0795, 1.143, 2.35585, 2.45745, 2.35585,
     &     2.44475, 2.35585, 2.50063, 2.35585, 2.50063, 0.0, 0.0,
     &     0.0, 0.0, 3.1395, 3.19659, 4.4577, 5.55625,
     &     2.5006, 2.51838/
c--   ZPOS(23): NEW! Transverse holding coil z center
c--   ZPOS(24): NEW! Transverse holding coil half length of racetrack section
c--   ZPOS(25): NEW! Transverse holding coil radius of curvature of inner coil
c--   ZPOS(26): NEW! Transverse holding coil radius of curvature of outer coil
      DATA ZPOS/2.64, 2.7, 3.2, 3.264, -2.64, 3.02707, 3.02707, 5.29529,
     &     5.29529, 6.15, 6.15, 6.17, -2.64, 3.02707, 3.02707, 5.29529,
     &     15.95, 15.97, -2.64, 18.9934, 21.61, 21.635,
     &     00.000, 0.0719, 0.9055, 3.9281/
c--   PHIPOS(1): NEW! Transverse holding coil angular center
      DATA PHIPOS /0.0000/

c--   Local user-defined variables

      REAL dir(3)
      DATA dir /0., 0., 1./
      INTEGER ierrInner, ierrOuter, ierr
      REAL pathLength2, pathLength1
      REAL rCylIn, rCylOut, phiCen, zLow, zHigh
      REAL zCen, zTrack, zCirIn, zCirOut
      REAL dS, dZ
      REAL myZ, myPhi

c--   Transverse holding coil
c--   Strategy: create cylinder and 'cut' non-coil parts out of it

      rCylIn  = RC(21)
      rCylOut = RC(22)
      phiCen  = PHIPOS(1)
      zCen    = ZPOS(23)
      zTrack  = ZPOS(24)
      zCirIn  = ZPOS(25)
      zCirOut = ZPOS(26)
      
c--   Check z0 and zHigh, zLow here!
      zLow    = zCen - zTrack - zCirOut
      zHigh   = zCen + zTrack + zCirOut
      call cylPath(vertex,cdir,rCylOut,zLow,zHigh,pathLength2,ierrOuter)
      call cylPath(vertex,cdir,rCylIn, zLow,zHigh,pathLength1,ierrInner)
      if (ierrOuter.ne.0) then
         pathLength = pathLength2 - pathLength1
      else
         pathLength = 0.0
      endif
      
c--   Find position of hit on inner cylinder: myZ, myPhi
      call CylInt(vertex, cdir, center, dir, rCylIn, a1, a2, ierr)
      if (ierr .eq. 0) then
         pathLength = 0.0
         RETURN
      end if
      myZ   = a1(3)
      myPhi = atan(a1(2)/a1(1))   ! image because tan returns acute angle!

c--   Snip anything beyond outside azimuthal range
      if (rCylIn*abs(myPhi-phiCen) .gt. zCirOut)
     &   pathLength = 0.0
      
c--   Racetracks
      if ( (abs(myZ-zCen) .lt. zTrack) .and.
     &     (rCylIn*abs(myPhi-phiCen) .lt. zCirIn) ) then
         pathLength = 0.0
c--   Downstream arc
      elseif (zCen-myZ .gt. zTrack) then
         dZ = -myZ + zCen - zTrack
         dS = rCylIn * abs(myPhi-phiCen)
         if (dZ**2 + dS**2 .gt. zCirOut**2)   pathLength = 0.0
         if (dZ**2 + dS**2 .lt. zCirIn**2)    pathLength = 0.0
      elseif (myZ-zCen .gt. zTrack) then
         dZ = myZ - zCen - zTrack
         dS = rCylIn * abs(myPhi-phiCen)
         if (dZ**2 + dS**2 .gt. zCirOut**2)   pathLength = 0.0
         if (dZ**2 + dS**2 .lt. zCirIn**2)    pathLength = 0.0
      endif

c--   Symmetry implies pathlength through cylinder is symmetric
c--   MUST COUNT IF PIERCES COIL TWICE, THOUGH
c--   if (ierr .eq. 2)
c--      dist(13) = pathLength * 2.0
         
      RETURN
      end
