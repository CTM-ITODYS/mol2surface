program moltosurface
implicit none
integer :: i,j,k,nbatom,nbchar,step,fix
integer,dimension(:), allocatable :: charvalue,atfix
real*8 :: xmin,xmax,dx,ymin,ymax,dy,zmin,zmax,dz
real*8, dimension(:,:), allocatable :: coord, coordtot
character*128, dimension(:), allocatable :: charatom, typeatm
character*128,dimension(117):: atomchar
character*128 :: inputxyz, outputbas, trash
logical :: stay
real*8,parameter :: borh2ang = 0.529177249
real*8, parameter :: pi = dacos(-1.0d0)
!variable paramètre surface
integer :: transa, transb
real*8 :: a,b,c,alpha,beta,gama,rotz
!variable modifier
integer :: nbatomtot   !nombre d'atomes total de la surface
character*128, dimension(:), allocatable :: charatomtot
real*8, dimension(:,:), allocatable :: newcoord
real*8, dimension(2,2) :: matrotz
real*8 :: ax,bx,by,cx,cy,cz          !box
real*8 :: amx,bmx,bmy,cmx,cmy,cmz,omega    !matrice passage
!===========================================================================!
!=============================lecture fichier xyz===========================!
!===========================================================================!
call getarg(1,inputxyz)
inputxyz=trim(inputxyz)
open(10,file=inputxyz,status='old',action='read')
read(10,*) nbatom
allocate(charatom(nbatom),coord(nbatom,3))
read(10,*)
do i = 1, nbatom
   read(10,*) charatom(i),(coord(i,j),j=1,3)
enddo
 close(10)
!==========================================================================!
!=================lecture fichier info_surface=============================!
!==========================================================================!
 open(10,file="info_surface",status='old',action='read')
read(10,*) 
read(10,*) trash, a
read(10,*) trash, b
read(10,*) trash, cz
read(10,*) trash, alpha
read(10,*) trash, beta
read(10,*) trash, gama
read(10,*)
read(10,*) nbchar
allocate(typeatm(nbchar),charvalue(nbchar))
do i = 1, nbchar
   read(10,*) typeatm(i)
enddo
 close(10)
!=========================================================================! 
!=====================conversion degre radian=============================!
!=========================================================================!
alpha = alpha*2*pi/360
beta = beta*2*pi/360
gama = gama*2*pi/360
!==========================================================================!
!============================construction box==============================!
!==========================================================================!
ax = a
bx = b*cos(gama)
by = b*sin(gama)
c = (cz**2)/(1-cos(beta)**2-((b**2)*(cos(alpha)**2)+(bx**2)*(cos(beta)**2)-2*b*bx*cos(beta)*cos(alpha))/(by**2))
c=sqrt(c)
cx = c*cos(beta)
cy = (b*c*cos(alpha)-bx*cx)/by
cz = sqrt(c*c-cx*cx-cy*cy)
 
!==========================================================================!
!==============================translation z+1=============================!
!==========================================================================!
nbatomtot = nbatom
allocate(charatomtot(nbatomtot),coordtot(nbatomtot,3),newcoord(nbatomtot,3))
step = 0
!do i = 0,1
i=1
      do k = 1, nbatom
         step =step +1
         coordtot(step,1) = coord(k,1)
         coordtot(step,2) = coord(k,2)
         coordtot(step,3) = coord(k,3) !+ (c*i)/2
         charatomtot(step) = charatom(k)
      enddo
deallocate(charatom,coord)

!==========================================================================!
!==================conversion cart coord en frac coord=====================!
!==========================================================================!
omega = a*b*c*sqrt(1-cos(alpha)*cos(alpha)-cos(beta)*cos(beta)-cos(gama)*cos(gama)+2*cos(alpha)*cos(beta)*cos(gama))
amx = 1/a
bmx = -(cos(gama)/(a*sin(gama)))
bmy = 1/(b*sin(gama))
cmx = b*c*(cos(alpha)*cos(gama)-cos(beta))/(omega*sin(gama))
cmy = a*c*(cos(beta)*cos(gama)-cos(alpha))/(omega*sin(gama))
cmz = (a*b*sin(gama))/omega
newcoord(:,:) = 0.0
do i = 1, nbatomtot
      newcoord(i,1) = coordtot(i,1)*amx + coordtot(i,2)*bmx + coordtot(i,3)*cmx
      newcoord(i,2) = coordtot(i,2)*bmy + coordtot(i,3)*cmy 
      newcoord(i,3) = coordtot(i,3)*cmz
enddo    
!==========================================================================!
!=========================calcul de charvalue==============================!
!==========================================================================!
charvalue(:) = 0
do i=1,nbchar
   do j=1,nbatomtot
         if (charatomtot(j) == typeatm(i)) then
               charvalue(i)= charvalue(i) + 1
         endif
   enddo
enddo
!==========================================================================!
!===================ecriture du fichier output=============================!
!==========================================================================!
!******************************format**************************************!
100  format("S:\CONVERSION\XYZ\TO\POSCAR")
101  format("1.000000")
102  format(f11.6,"    0.000000    0.000000")
103  format(f12.6,f12.6,"    0.000000")
104  format(f12.6, f12.6, f12.6)
105  format("  ! ")
106  format(/,"Selective Dynamics")
107  format("Direct")
108  format(3f14.9,"  T  T  T ! ",a3)
109  format(3f14.9,"  F  F  F ! ",a3)
open(55,file="POSCAR",action='write')
write(55,100)
write(55,101)
write(55,102) ax
write(55,103) bx,by
write(55,104) cx, cy, cz
do i = 1, nbchar
   write(55,'(A3)',advance='no') typeatm(i)
enddo
write(55,*)
do i = 1, nbchar
   write(55,'(I4)',advance='no') charvalue(i)
enddo
write(55,106)
write(55,107)
do i = 1, nbchar
   do j = 1, nbatomtot
      if (typeatm(i) .eq. charatomtot(j) ) then
         step = 1
         do k = 1, fix
            if ( j .eq. atfix(k)) then
               step = 0
            endif
         enddo
         if ( step .eq. 0) then
            write(55,109) newcoord(j,1),newcoord(j,2),newcoord(j,3),charatomtot(j)
         else
            write(55,108) newcoord(j,1),newcoord(j,2),newcoord(j,3),charatomtot(j)
         endif
      endif
   enddo
enddo
 close(55)
!!*************************************************************************!
!!*****************libération de la mémoire********************************!
!!*************************************************************************!
deallocate(charatomtot,coordtot)

endprogram





































