! Codigo que calcula los autovalores y autovectores del hamiltoniano
! efectio 1D para un electron.
!
! El hamiltoniano del sistema es
!
!   Hz = -hbar^2/(2*me)*d^2/dz^2 + V_{ef}(z) + Hr*S
!
! El potencial de confinamiento efectivo es
!
!                 v1    if abs(z)<=0.5*az
!   V_{ef}(z) =  -v2*vr if 0.5*az<=abs(z)<=0.5*(az+bz)
!                 0     other cases
!
! Hr es el valor de expectacion del hamiltoniano radial del 
! problema usando como funcion de onda el estado fundamental del pozo
!
! La matriz del hamiltoniano la calculo usando como base
! del espacio los B-splines. El problema de autovalores es un
! problema generalizado
!
!     H*c = E*S*c
!
! Guardo, los autovalores del hamiltoniano en funcion del
! campo magnetico aplicado.
! Algunos comentarios sobre unidades y campo magnetico.
! La cte alpha sirve para determinar el parametro l del potencial efectivo en funcion del
! campo magnetico, la relacion es:
!                 B = alpha/l**2
! donde B esta en teslas [T] y l en nanometros [nm].
!
! Para pasar de nanometros a unidades atomicas se usa la relacion:
!                 1 UAlenght = 0.0529177210 nm
!
! Para pasar de eV a unidades atomicas se use la relacion:
!                 1 UAenergia = 27.21138564 eV

!!!!!!!!!!!!!!!!!!!!!! preprocessing variable !!!!!!!
#ifndef zmin_
#define zmin_  0._pr
#endif

#ifndef zmax_
#define zmax_ 10._pr
#endif

#ifndef lum_
#define lum_ 0
#endif

#ifndef c_
#define c_ 0._pr
#endif

#ifndef gamma_
#define gamma_ 0._pr
#endif

#ifndef kord_
#define kord_ 5
#endif

#ifndef l_interval_
#define l_interval_ 100
#endif

#ifndef intg_
#define intg_ 100
#endif

#ifndef nev_
#define nev_ 10
#endif

#ifndef V1_
#define V1_ 1._pr
#endif

#ifndef V2_
#define V2_ 1._pr
#endif

#ifndef r0_
#define r0_ 1._pr
#endif

#ifndef az_
#define az_ 1._pr
#endif

#ifndef bz_
#define bz_ 1._pr
#endif

#ifndef B_campo_i_
#define B_campo_i_ 0.1_pr
#endif

#ifndef B_campo_f_
#define B_campo_f_ 1._pr
#endif

#ifndef num_puntos_B_
#define num_puntos_B_ 10
#endif

#ifndef me_
#define me_ 1._pr
#endif

!!!!!!!!!!!!!!!!!!!!!! MODULES !!!!!!!!!!!!!!!!!!!!!!!
module precision
implicit none
integer, parameter :: pr = selected_real_kind(15,307)
end module

module carga
use precision
implicit none
integer :: kord, lum, intg, nev, num_puntos_B
real(pr) :: zmin, zmax, v1, v2, sigma, B_campo_i, B_campo_f
real(pr) :: az, bz, r0
real(pr) :: me, B_campo, delta_B, c, gamma, l_campo
real(pr) :: Tr, Ur, Vr, Hr, omega
integer :: l_interval
character(1) :: tip, B1, B2
end module  carga

module matrices
use precision
implicit none
integer :: nk, nb
real(pr), allocatable, dimension(:) :: norma
real(pr), allocatable, dimension(:,:) :: s, v01, v02, ke
end module matrices

module integracion
use precision
implicit none
integer, allocatable, dimension(:) :: k
real(pr), allocatable, dimension(:) :: t, sp !, dsp
real(pr), allocatable, dimension(:,:) :: x, w, pl
end module  integracion

!!!!!!!!!!!!!!!!!!!!!!!!!!!!   MAIN !!!!!!!!!!!!!!!!!!!!!!!!!
program Bsplines
use precision
use carga
use matrices
use integracion
implicit none
integer :: i, j, ind_B, tipo
real(pr), parameter :: alpha = 658.4092645439_pr, a0 = 0.0529177210_pr, eV = 27.21138564_pr
real(pr), parameter :: c_light = 137.035999074492_pr, ua_to_T = 1.72d3
real(pr) :: time
character(150) :: archivo1e
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

zmin = real(zmin_, pr); zmax = real(zmax_, pr) ! valores de inicio y final del intervalo de integracion

tip = 'u'; tipo = 1 ! tipo de distribucion de knots, u, e, m

lum = lum_ ! # de intervalors u en m

c = real(c_, pr) ! parametod en m dist u en [a,c] y e en (c,b]

gamma = real(gamma_, pr) ! parametro de la exponencial

kord = kord_ ! orden de los B-splines

l_interval = l_interval_ ! # de intervalos

intg = intg_ ! grado de la intregracion de cuadratura gaussiana, >=k

nev = nev_ ! # de autovalores que queremos calculara

v1 = real(V1_, pr) ! profundidad del pozo en eV

v2 = real(V2_, pr)

az = real(az_, pr); bz = real(bz_, pr); r0 = real(r0_, pr);

B_campo_i = real(B_campo_i_, pr); B_campo_f = real(B_campo_f_, pr)

num_puntos_B = num_puntos_B_

me = real(me_, pr) ! en unidades atomicas

if( intg<kord ) intg = kord

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

nk = l_interval+2*kord-1   ! # de knots
! No incluimos primer   y ultimo spline debido a psi(0)=psi(R)=0
nb = kord+l_interval-3      ! size base splines
if( nev<0 ) nev = nb

write(archivo1e, '("./resultados/1e-E_vs_B-v1_",f6.4,"eV-v2_",f6.4,"eV-az_",f6.4,"-bz_",f6.4,".dat")') v1, v2, az, bz
!###########################################################
!###########################################################
!###########################################################
open(10, file=archivo1e)
!###########################################################
!###########################################################
write(10,'(A28,f8.2,A1,f8.2,A4)') "# Intervalo de integracion:[", zmin,",", zmax, "] nm"
write(10,'(A32,x,I2,x,A33)') "# Tipo de distribucion de knots:", tipo, "(1 es uniforme, 2 es exponencial)"
write(10,'(A47,x,f8.6)') "# Cte de decaimiento en distribucion exp gamma =", gamma
write(10,'(A30,x,I2)') "# Orden de los bsplines kord =", kord
write(10,'(A19,x,I3,x,A24,x,I3)') "# Num de intervalos", l_interval, ", grado de la cuadratura", intg
write(10,'(A14,x,I3,x,A24,x,I3)') "# Num de knost", nk, ", tamaño de la base nb =", nb
write(10,'(A33,x,f8.6,x,A2)') "# Masa efectiva del electron me =", me, "UA"
write(10,'(A22,x,f6.4,x,A2)') "# Altura del pozo V1 =", v1, "eV"
write(10,'(A27,x,f6.4,x,A2)') "# Profundidad del pozo V2 =", v2, "eV"
write(10,'(A22,x,f6.4,x,A9,x,f6.4,x,A2)') "# Radios del pozo az =", az, "nm y bz =", bz, "nm"
write(10,'(A27,x,f6.2,x,A19,x,f7.2)') "# Campo inicial b_campo_i =", b_campo_i, "y final b_campo_f =", b_campo_f
write(10,'(A24)') "# Autovalores calculados"
call flush();
!###########################################################
!###########################################################
!###########################################################
close(10)

! paso el parametro del potencial a unidades atomicas
zmin = zmin/a0; zmax = zmax/a0;
gamma = gamma/a0;
az = az/a0; bz = bz/a0; r0 = r0/a0
v1 = v1/eV; v2 = v2/eV;

allocate(Sp(kord))

allocate(x(l_interval,intg), w(l_interval,intg), pl(l_interval,intg))

allocate(t(nk), k(nk))

allocate(norma(nb), s(nb,nb), v01(nb,nb), v02(nb,nb), ke(nb,nb))

call KNOTS_PESOS(kord, tip, gamma, zmin, zmax, c, l_interval, lum, intg, t, k, x, w, pl);

delta_B = (B_campo_f-B_campo_i)/real(num_puntos_B-1, pr);

call matrix_elements( )

do i = 1,nb
  do j = i+1,nb
    s(j,i) = s(i,j)
    v01(j,i) = v01(i,j)
    v02(j,i) = v02(i,j)
    ke(j,i) = ke(i,j)
  end do
end do

! OJO, ESTO SOLO VALE PARA r0 = 7nm
Tr = 0.5_pr/me*9.74494d-5
Ur = 0.5_pr*me*10677.5_pr
Vr = 0.82604193_pr

do ind_B = 1, num_puntos_B

  B_campo = B_campo_i + real(ind_B-1, pr)*delta_B;

  l_campo = sqrt(alpha/B_campo)/a0;
  omega = 0.5*B_campo/(me*c_light*ua_to_T);  !! Frecuencia de oscilacion debida al campo

  Hr = Tr + omega**2*Ur

  call sener(archivo1e);

end do

deallocate(Sp, x, w, pl)
deallocate(t, k, norma, s, v01, v02, ke)

call cpu_time(time);
write(*,*)time/60._pr

end program !termina el main, termina el programa

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine matrix_elements( )
use precision
use carga
use matrices
use integracion
implicit none
integer :: i, j, m, n, im, in
real(pr) :: bm, bn, zz

s = 0._pr ; v01 = 0._pr; v02 = 0._pr; ke = 0._pr

do i = kord, kord+l_interval-1
  do j = 1, intg
    zz = x(k(i),j)
    call bsplvb(t, kord, 1, zz, i, Sp)

    do m = 1,kord
      im = i-kord+m-1
      if( im>0.and.im<nb+1 ) then
        do n = 1,kord
          in = i-kord+n-1

          if( in>0.and.in<nb+1 )then
            s(im,in) = s(im,in) + sp(m)*sp(n)*w(k(i),j)

            if(abs(zz)<=0.5*az) v02(im,in) = v02(im,in) + w(k(i),j)*sp(m)*sp(n)

            if(0.5*az<abs(zz).and.abs(zz)<=0.5*(az+bz)) v01(im,in) = v01(im,in) + w(k(i),j)*sp(m)*sp(n)

          endif
        end do
      endif
    end do
  end do
end do

do i = kord, kord+l_interval-1
  do m = i-kord+1, i
    if(m>1.and.m<nb+2)then
      do n = m,i
        if(n<nb+2)then
          do j = 1, intg

            call bder(x(k(i),j), t, bm, bn, kord, nk, m, n, i)

            ke(m-1,n-1) = ke(m-1,n-1) + 0.5_pr*w(k(i),j)*bm*bn/me

          end do
        endif
      end do
    endif
  end do
end do


end subroutine matrix_elements

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! sener.f90 basado en ener.f90
! pero es una subroutine de exp_bs_gs.f90 para calcular el GS

! Calcula E_i de una matriz que se le pasa por mudule

subroutine sener(archivo1e)
use precision
use carga
use matrices
use integracion
implicit none
integer(4) :: m, n, dp
integer(4) :: info
real(pr), parameter :: eV = 27.21138564_pr
real(pr), allocatable, dimension(:) :: e1
real(pr), allocatable, dimension(:,:) :: mh, s_copy, auvec1
character(150) :: archivo1e

dp = nb

allocate(e1(nev), auvec1(dp,nev), mh(dp,dp), s_copy(dp,dp))

!###########################################################
!###########################################################
!###########################################################
open(12, file=archivo1e, position = 'append')
!###########################################################
!###########################################################
!###########################################################

mh = 0._pr

!!! hamiltoniano de una particula
s_copy = s;
do m = 1, dp
  do n = 1, m
    mh(m,n) = ke(m,n) + v1*v01(m,n) - vr*v2*v02(m,n) + Hr*s(m,n)
    mh(n,m) = mh(m,n)
  end do
end do

call eigen_value(dp, nev, info, mh, s_copy, e1, auvec1);
e1(:) = eV*e1(:);

!!!!!!!!!!!!!!!!!!!!!!!####################################
!!!!!!!!!!!!!!!!!!!!!!!####################################
!!!!!!!!!!!!!!!!!!!!!!!####################################
write(12,6) B_campo, (e1(m), m = 1, 25)
!!!!!!!!!!!!!!!!!!!!!!!####################################
!!!!!!!!!!!!!!!!!!!!!!!####################################
!!!!!!!!!!!!!!!!!!!!!!!####################################
info = 0

close(12)

deallocate(e1, auvec1, mh)

6 format(e22.15,1x,1000(1x,e22.15))

return
end subroutine sener

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine eigen_value(dp,nev,info,nh,s,x,avec)
use precision
!
! calcula usando LAPACK N=nev autovalores
implicit none
integer*4 dp,nev
real(pr),dimension (nev)::x
real(pr),dimension (dp,nev)::avec
real(pr),dimension (dp,dp)::nh,s
! definiciones para DSYGVX
integer(4)::itype,numau,lwork,info,il,iu
character(1)::jobz,range,uplo
real(pr)::vl,vu,abstol
real(pr),dimension(dp,dp)::z
real(pr),dimension(9*dp)::work
integer(4),dimension(5*dp)::iwork
integer(4),dimension(dp)::ifail

lwork=9*dp ! lwork=dim(work)

!!!!!!!!!!!!!!!!!!!!!!!!!
! reemplaza  tred2 y tqli del NR por dsygvx de LAPACK
! parametros para dsygvx:
itype=1  ! especifica que A x=lambda Bx
jobz='V' ! V (N) (NO) calcula autovectores
!range='V'! autovalores en el intervalo [VL,VU]
!vl=-1.01_pr ! pongo  vl menor que el menor autovalor
!vu=0._pr   ! por ahora solo interezan los autovalores negativos
uplo='U' ! la matriz es storeada superior
! N es dp
!A es nh
! lda = dp
! B = s  ! NO hay que hacer cholesky !!!
! ldb = dp
! il, iu  autov entre il y iu SI range='I'
 range='I'
il=1
iu=nev
abstol=1.d-12  ! por decir algo (?)aconsejan abstol=2*DLAMCH('S') ; ????
! M = numau# de autovalores calculados : output
! W = x , de menor a mayor
! se debe definir un array Z(LDZ,M) NO conocemos M a priori! poner N=dp
! Z devuelve eigenvectors
! ldz=dp : dimension de Z
! definir work(lwork)
lwork=9*dp
! definir integer array iwork(5*N)
! ifail=output:
!If JOBZ = 'V', then if INFO = 0, the first M elements of
!          IFAIL are zero.  If INFO > 0, then IFAIL contains the
!          indices of the eigenvectors that failed to converge.
!          If JOBZ = 'N', then IFAIL is not referenced.
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  DPOTRF or DSYEVX returned an error code:
!             <= N:  if INFO = i, DSYEVX failed to converge;
!                    i eigenvectors failed to converge.  Their indices
!                    are stored in array IFAIL.
!             > N:   if INFO = N + i, for 1 <= i <= N, then the leading
!                    minor of order i of B is not positive definite.
!                    The factorization of B could not be completed and
!                    no eigenvalues or eigenvectors were computed.


call DSYGVX( ITYPE, JOBZ, RANGE, UPLO, dp, nh, dp, s, dp,            &
     &                   VL, VU, IL, IU, ABSTOL, numau, x, Z, dp, WORK,  &
     &                   LWORK, IWORK, IFAIL, INFO )


!!!!!!!!!!!!!!!!!!

avec(1:dp,1:nev)=z(1:dp,1:nev)



        return
        end subroutine eigen_value

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  bajada de http://www.atom.physto.se/~lindroth/comp08/bget.f


       subroutine bder(rr,t,dm,dn,kord,np,indexm,indexn,left)
!     returns the value of (d/dx) Bspline(kord,index) in rr
!     The first Bspline is called spline #1 (i.e. index=1)
!     the first knot point is in t(1)
!     np= number of knot points (distinct or multiple) including
!     the ghost points: N phyical points np=N +2*(kord-1)


      implicit none
      integer i,left,kord,indexm,indexn,np
      real*8 dm,dn,t(np),Sp(kord),rr

      dm=0.d0  ; dn=0.d0
!     if rr=t(np) then the routine assumes that
!     there is kord knotpoints in the last physical point and
!     returns Bder.ne.zero if index is np-kord, or np-kord-1
      if(rr.gt.t(np)) return
      if(rr.lt.t(1)) return
!      do it=1,np
!        if(rr.ge.t(it)) left=it
!      end do

      if(abs(rr-t(np)).lt.1.d-10) then
        !if(index.lt.np-kord-1) return

        if(indexm.eq.np-kord) then
          dm=dble(kord-1)/(t(np)-t(np-kord))
        else if(indexm.eq.np-kord-1) then
          dm=-dble(kord-1)/(t(np)-t(np-kord))
        end if

       if(indexn.eq.np-kord) then
          dn=dble(kord-1)/(t(np)-t(np-kord))
        else if(indexn.eq.np-kord-1) then
          dn=-dble(kord-1)/(t(np)-t(np-kord))
        end if

        return
      end if

      call bsplvb(t,kord-1,1,rr,left,Sp)

  if(indexm-left+kord.ge.1.or.indexm-left+kord.le.kord)then
      i=indexm-left+kord
      if(i.eq.1) then
        dm=dble(kord-1)*(-Sp(i)/(t(indexm+kord)-t(indexm+1)))
      else if(i.eq.kord) then
        dm=dble(kord-1)*(Sp(i-1)/(t(indexm+kord-1)-t(indexm)))
      else
        dm=dble(kord-1)*(Sp(i-1)/(t(indexm+kord-1)-t(indexm))- &
     &  Sp(i)/(t(indexm+kord)-t(indexm+1)))
      end if
  endif
  if(indexn-left+kord.ge.1.or.indexn-left+kord.le.kord)then
      i=indexn-left+kord

      if(i.eq.1) then
        dn=dble(kord-1)*(-Sp(i)/(t(indexn+kord)-t(indexn+1)))
      else if(i.eq.kord) then
        dn=dble(kord-1)*(Sp(i-1)/(t(indexn+kord-1)-t(indexn)))
      else
        dn=dble(kord-1)*(Sp(i-1)/(t(indexn+kord-1)-t(indexn))- &
     &  Sp(i)/(t(indexn+kord)-t(indexn+1)))
      end if

  endif
      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE bsplvb(t,jhigh,index,x,left,biatx)
!      INCLUDE 'standa.typ'
!      INCLUDE 'spldim.def'
!      include 'ndim_inc'
       implicit none
      integer(4),PARAMETER::JMAX=100
      integer index,jhigh,left,i,j,jp1
      real*8 t,x,biatx,deltal,deltar,saved,term
      DIMENSION biatx(jhigh),t(left+jhigh),deltal(jmax),deltar(jmax)
      SAVE deltal,deltar
      DATA j/1/
!      write(6,*) ' jmax=',jmax
      GO TO (10,20),index
 10   j = 1
      biatx(1) = 1.d0
      IF (j .GE. jhigh) GO TO 99

 20   CONTINUE
         jp1 = j + 1
         deltar(j) = t(left+j) - x
         deltal(j) = x - t(left+1-j)
!         write(6,'(1pd12.4,2(i5,1pd14.6))')
!     :   x,left+j,t(left+j),left+1-j,t(left+1-j)
!         write(6,'(i3,1p3d12.4)') j,deltal(j),deltar(j),
!     :   abs(deltal(j)-deltar(j))
         saved = 0.d0
         DO i = 1,j
!         write(6,'(2i3,1p3d12.4)') i,j,deltal(jp1-1),deltar(i),
!     :   abs(deltal(jp1-1)-deltar(i))

             term = biatx(i)/(deltar(i) + deltal(jp1-i))
             biatx(i) = saved + deltar(i)*term
             saved = deltal(jp1-i)*term
         END DO
         biatx(jp1) = saved
         j = jp1
         IF (j .LT. jhigh) GO TO 20
 99   RETURN
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE KNOTS_PESOS(kord,tip,gamma,a,b,c,l,lum,intg,t,k,x,w,pl)

!lo único que hace esta subrutina principal es dividir los casos segun se quiera knots uniformes, exponencial, mixto
!cada una de las subrutinas que llama ademas calcula los x y los w, abscisas y pesos para la cuadratura gaussiana, y los polinomios


!  INPUT:
! kord: orden de los b-splines
! tip; character(1): 'u' ; 'e' ; 'm': dist uniforme, exp o mixta de knots
! gamma : param dist. e (no usado si tip='u')
! a
! b  ; a<b todo es calculado en el intervalo [a,b]
! c solo usado si tip='m' => a<c<b ; dist u  en [a,c]; e en (c,b]; c es a u-knot
! NOTA I : dist optima: m; con c~r_0 donde v(r) es apreciablemente no nulo
! l = # total de sub_intervalos en [a,b]
! lum: solo se usa si tip=m : # subint con dist u => l-lem con dist e
! intg = # de puntos que usamos en la integral x cuadratura en c/subintervalo


!  OUTPUT:
!  t(l+2*kord-1)  ! los (l+2*kord-1) knots (contando multiplicidades en a y b)
!  k(l+2*kord-1)  ! da el # de intervalo j para el i-esimo nodo 1<=k<=l
!  x(l,intg),w(l,intg) : posiciones y pesos de los intg puntos en c/intervalo
! pl(l,intg),w(l,intg) : polinomios de Legendre en cada punto
! NOTA II : los p_l los calcula gratis la subr gauleg de int. x cuad. del NR
! y a veces hacen falta=> version modif. de gauleg que los da como output
use precision
implicit none
integer(4)::kord,l,lum,intg
character(1)::tip
real(pr)::gamma,a,b,c
integer(4),dimension(l+2*kord-1)::k
real(pr),dimension(l+2*kord-1)::t
real(pr),dimension(l,intg)::x,w,pl
!!!!!!


if(tip=='u')then
call dist_unif(kord,a,b,l,intg,t,k,x,w,pl)
elseif(tip=='e')then
call dist_exp(kord,gamma,a,b,l,intg,t,k,x,w,pl)
elseif(tip=='m')then
call dist_mix(kord,gamma,a,b,c,l,lum,intg,t,k,x,w,pl)
else
write(*,*)'error 1 en KNOTS_PESOS :',tip,' no corresponde a una distribucion'
stop
endif

return
end    SUBROUTINE KNOTS_PESOS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE DIST_UNIF(kord,a,b,l,intg,t,k,x,w,pl)
use precision
implicit none
integer(4)::kord,l,intg
real(pr)::a,b
integer(4),dimension(l+2*kord-1)::k
real(pr),dimension(l+2*kord-1)::t
real(pr),dimension(l,intg)::x,w,pl
!!!!!!
integer(4)::i,j,nk
real(pr)::ri,rf,dr
real(pr),dimension(intg)::vx,vw,vpl

nk=l+2*kord-1   ! # de knots

dr=(b-a)/dfloat(l)

! calcula los puntos y pesos para la cuadratura
x=0._pr;w=0._pr

do i=1,l

   ri=a+dfloat(i-1)*dr
   rf=ri+dr

   call gauleg_pl(ri,rf,vx,vw,vpl,intg)

   do j=1,intg
   x(i,j)=vx(j)
   w(i,j)=vw(j)
   pl(i,j)=vpl(j)
   end do

end do

t(1)=a
k(1)=1

do i=2,kord
   t(i)=t(1)
   k(i)=1
end do

do i=kord+1,kord+l
   t(i)=t(i-1)+dr
   k(i)=k(i-1)+1
end do

do i=kord+l+1,nk
   t(i)=t(i-1)
   k(i)=k(i-1)
end do

return

end SUBROUTINE DIST_UNIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE DIST_EXP(kord,gamma,a,b,l,intg,t,k,x,w,pl)
use precision
implicit none
integer(4)::kord,l,intg
real(pr)::gamma,a,b
integer(4),dimension(l+2*kord-1)::k
real(pr),dimension(l+2*kord-1)::t
real(pr),dimension(l,intg)::x,w,pl
!!!!!!
integer(4)::i,j,nk
real(pr)::ri,rf,dr,ye
real(pr),dimension(intg)::vx,vw,vpl

nk=l+2*kord-1   ! # de knots

dr=(b-a)/(dexp(gamma)-1._pr)
ye=gamma/dfloat(l)

! calcula los puntos y pesos para la cuadratura
x=0._pr;w=0._pr

do i=1,l

   ri=a+dr*(dexp(ye*dfloat(i-1))-1._pr)
   rf=a+dr*(dexp(ye*dfloat(i))-1._pr)

   call gauleg_pl(ri,rf,vx,vw,vpl,intg)

   do j=1,intg
      x(i,j)=vx(j)
      w(i,j)=vw(j)
      pl(i,j)=vpl(j)
   end do

end do

t(1)=a
k(1)=1

do i=2,kord
   t(i)=t(1)
   k(i)=1
end do

do i=kord+1,kord+l
   t(i)=a+dr*(dexp(ye*dfloat(k(i-1)))-1._pr)
   k(i)=k(i-1)+1
end do

do i=kord+l+1,nk
t(i)=t(i-1)
k(i)=k(i-1)
end do

return
end    SUBROUTINE DIST_EXP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE DIST_MIX(kord,gamma,a,b,c,l,lum,intg,t,k,x,w,pl)
use precision
implicit none
integer(4)::kord,l,lum,intg
real(pr)::gamma,a,b,c
integer(4),dimension(l+2*kord-1)::k
real(pr),dimension(l+2*kord-1)::t
real(pr),dimension(l,intg)::x,w,pl
!!!!!!
integer(4)::i,j,le,nk
real(pr),dimension(lum,intg)::xu,wu,plu
real(pr),dimension(l-lum,intg)::xe,we,ple
integer(4),dimension(lum+2*kord-1)::ku
real(pr),dimension(lum+2*kord-1)::tu
integer(4),dimension(l-lum+2*kord-1)::ke
real(pr),dimension(l-lum+2*kord-1)::te

nk=l+2*kord-1   ! # de knots

call DIST_UNIF(kord,a,c,lum,intg,tu,ku,xu,wu,plu)

le=l-lum

call DIST_EXP(kord,gamma,c,b,le,intg,te,ke,xe,we,ple)

do i=1,lum
x(i,1:intg)=xu(i,1:intg)
w(i,1:intg)=wu(i,1:intg)
pl(i,1:intg)=plu(i,1:intg)
end do

do i=1,le
x(lum+i,1:intg)=xe(i,1:intg)
w(lum+i,1:intg)=we(i,1:intg)
pl(lum+i,1:intg)=ple(i,1:intg)
end do

do i=1,kord+lum
t(i)=tu(i)
k(i)=ku(i)
end do

j=kord+lum
do i=kord+1,le+kord+1
j=j+1
t(j)=te(i)
k(j)=ke(i)+ku(kord+lum)-1
end do

do i=j,nk
t(i)=t(i-1)
k(i)=k(i-1)
end do

return
end    SUBROUTINE DIST_MIX

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! subrutina de NR modificada
! calcula coordenadas y pesos GL y los valores de P_l(x)

      SUBROUTINE gauleg_pl(x1,x2,x,w,pl,n)
      INTEGER n
      DOUBLE PRECISION x1,x2,x(n),w(n),pl(n)
      DOUBLE PRECISION EPS
      PARAMETER (EPS=1.d-14)
      INTEGER i,j,m
      DOUBLE PRECISION p1,p2,p3,pp,xl,xm,z,z1
      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do 12 i=1,m
        z=dcos(3.1415926535897932385d0*(dfloat(i)-.25d0)/(dfloat(n)+.5d0))
1       continue
          p1=1.d0
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/dfloat(j)
11        continue
          pl(i)=p1
          pp=dfloat(n)*(z*p1-p2)/(z*z-1.d0)
          z1=z
          z=z1-p1/pp
        if(dabs(z-z1).gt.EPS)goto 1
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
12    continue
      return
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
