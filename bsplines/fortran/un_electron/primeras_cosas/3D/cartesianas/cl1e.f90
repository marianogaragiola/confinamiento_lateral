! Codigo que calculos los autovalores de un electron en un pozo con forma gaussiana.
!
! El potencial es de la forma V(z) = -V0*exp(-(z-z0)**2/(2*sigma**2))
! Los parametros del pozo son los siguiente:
! V0    ---> profundidad del pozo.
! sigma ---> ancho del pozo.
! z0    ---> es donde estan centrado el pozo.
! Estos parametros son parametros de entrada.
!
! El codigo usa B-splines para calcular el hamiltoniano y luega calcula los
! autovalores y autovectores que elija.
!
! El intervalo de integracion es [zmin, zmax].
! Ya que el potencial es simetrico respecto a z=0 voy a tomar a zmin = -zmax.
!
! El codigo va variando la variable V0, por lo tanto devuelve los autovalores como
! funcion de V0.
!!!!!!!!!!!!!!!!!!!!!! preprocessing variable !!!!!!!
#ifndef zmin_
#define zmin_  0.d0
#endif

#ifndef zmax_
#define zmax_ 10.d0
#endif

#ifndef z0_
#define z0_ 0.d0
#endif

#ifndef sigma_
#define sigma_ 1.d0
#endif

#ifndef lum_
#define lum_ 0
#endif

#ifndef c_
#define c_ 0.d0
#endif

#ifndef beta_
#define beta_ 0.d0
#endif

#ifndef kord_
#define kord_ 5
#endif

#ifndef l_
#define l_ 30
#endif

#ifndef intg_
#define intg_ 100
#endif

#ifndef nev_
#define nev_ 10
#endif

#ifndef me_
#define me_ 1.d0
#endif

#ifndef V0_i_
#define V0_i_ 5.d0
#endif

#ifndef V0_f_
#define V0_f_ 10.d0
#endif

#ifndef num_puntos_V0_
#define num_puntos_V0_ 10
#endif

!!!!!!!!!!!!!!!!!!!!!! MODULES !!!!!!!!!!!!!!!!!!!!!!!

module carga
integer :: kord, lum, intg, nev, num_puntos_V0
real(kind(0.d0)) :: zmin, zmax, V0, sigma, V0_i, V0_f, z0, c, beta, me
integer :: l
character(1) :: tip, cV01, cV02, cV03
end module  carga

module matrices
integer :: nk, nb
real(kind(0.d0)), allocatable, dimension(:,:) :: s, v01, ke, aux
end module matrices

module integracion
integer, allocatable, dimension(:) :: k
real(kind(0.d0)), allocatable, dimension(:) :: t, sp
real(kind(0.d0)), allocatable, dimension(:,:) :: x, w, pl
end module  integracion

!!!!!!!!!!!!!!!!!!!!!!!!!!!!   MAIN !!!!!!!!!!!!!!!!!!!!!!!!!
program Bsplines
use carga
use matrices
use integracion
implicit none
integer :: i, j
real(kind(0.d0)), parameter :: a0 = 0.0529177210d0, eV = 27.21138564d0
real(kind(0.d0)) :: time
character(150) :: archivo
character(1) :: z1, z2, z3, z4, l1, l2, l3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

zmin = zmin_; zmax = zmax_ ! valores de inicio y final del intervalo de integracion

z0 = z0_ ! centro del pozo

sigma = sigma_ ! ancho del pozo.

tip = 'e' ! tipo de distribucion de knots, u, e, m

lum = lum_ ! # de intervalors u en m

c = c_ ! parametod en m dist u en [a,c] y e en (c,b]

beta = beta_ ! parametro de la exponencial

kord = kord_ ! orden de los B-splines

l = l_ ! # de intervalos

intg = intg_ ! grado de la intregracion de cuadratura gaussiana, >=k

nev = nev_ ! # de autovalores que queremos calculara

me = me_; ! masa del electron en UA

V0_i = V0_i_; V0_f = V0_f_

num_puntos_V0 = num_puntos_V0_

if( intg<kord ) intg = kord

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

nk = l+2*kord-1   ! # de knots
! No incluimos primer   y ultimo spline debido a psi(0)=psi(R)=0
nb = kord+l-3      ! size base splines
if( nev<0 ) nev = nb

z1 = char(modulo(int(zmax), 10) + 48);
z2 = char(modulo(int(zmax), 100)/10 + 48);
z3 = char(modulo(int(zmax), 1000)/100 + 48);
z4 = char(modulo(int(zmax), 10000)/1000 + 48);

l1 = char(modulo(int(l), 10) + 48);
l2 = char(modulo(int(l), 100)/10 + 48);
l3 = char(modulo(int(l), 1000)/100 + 48);

!!!!! Paso todo a unidades atomicas
V0_i = V0_i/eV; V0_f = V0_f/eV;
zmin = zmin/a0; zmax = zmax/a0; z0 = z0/a0; sigma = sigma/a0;
beta = beta*a0;

!###########################################################
!###########################################################
!###########################################################
archivo = './resultados/1e-E_vs_B-zmax'//z4//z3//z2//z1//'-l'//l3//l2//l1//'.dat';
!archivo = 'prueba.dat'
open(9, file = archivo)
!###########################################################
!###########################################################
write(9,*) '# Codigo que calculos los autovalores de un electron en un pozo gaussiano.'

write(9,*) '# El potencial es de la forma V(z) = -V0*exp(-(z-z0)**2/(2*sigma**2)).'
write(9,*) '# Los parametros del pozo son los siguiente:'
write(9,*) '# V0     ---> profundidad del pozo.'
write(9,*) '# sigma  ---> ancho del pozo.'
write(9,*) '# z0     ---> centro dle pozo.'
write(9,*) '# Estos parametros son parametros de entrada.'

write(9,*) '# El codigo usa B-splines para calcular el hamiltoniano y luega calcula los '
write(9,*) '# autovalores y autovectores que elija.'

write(9,*) '# El intervalo de integracion es [zmin, zmax].'
write(9,*) '# Ya que el potencial es simetrico respecto a z=0 voy a tomar a zmin = -zmax.'

write(9,*) '# El codigo va variando la variable Vl, por lo tanto devuelve los autovalores como'
write(9,*) '# funcion de V0.'

write(9,*) '# tip: entrar el tipo de distribucion, u, e, m'
write(9,*) '# tip =', tip

write(9,*) '# lum: entrar el # de intervalos u en m'
write(9,*) '# lum =', lum

write(9,*) '# c: parametro en m dist u en [a,c] y e en (c,b]'
write(9,*) '# c =', c

write(9,*) '# beta: parametro en exponencial'
write(9,*) '# beta =', beta

write(9,*) '# kord: orden de los b-splines'
write(9,*) '# kord =', kord

write(9,*) '# l: # de intervalos'
write(9,*) '# l =', l

write(9,*) '# nb: tamaño base'
write(9,*) '# nb =', nb

write(9,*) '# intg: grado de integracion de la cuadratura >= k'
write(9,*) '# intg =', intg

write(9,*) '# nev: # de autovalores a calcular <=nb'
write(9,*) '# nev =', nev

write(9,*) '# intervalo de integracion'
write(9,*) '# zmin =', a0*zmin, 'zmax =', a0*zmax

write(9,*) '# centro del pozo:'
write(9,*) '# z0 =', a0*z0

write(9,*) '# masa del electron:'
write(9,*) '# me =', me

write(9,*) '# pozo inicial=', V0_i

write(9,*) '# pozo final=', V0_f

write(9,*) '# ancho del pozo:'
write(9,*) '# sigma =', sigma*a0

write(9,*) '# autovalores calculados'
!###########################################################
!###########################################################
!###########################################################
close(9)

allocate(Sp(kord))
allocate(x(l,intg), w(l,intg), pl(l,intg))
allocate(t(nk), k(nk))
allocate(s(nb,nb), v01(nb,nb), ke(nb,nb))

call KNOTS_PESOS(kord, tip, beta, zmin, zmax, c, l, lum, intg, t, k, x, w, pl)

s(:,:) = 0.d0; v01(:,:) = 0.d0; ke(:,:) = 0.d0;
call matrix_elements( )

do i = 1, nb
  do j = i+1, nb
    s(j,i) = s(i,j)
    v01(j,i) = v01(i,j)
    ke(j,i) = ke(i,j)
  end do
end do

call sener(archivo)

deallocate(Sp)
deallocate(x, w, pl)
deallocate(t, k)
deallocate(s, v01, ke)

call cpu_time(time)
write(*,*)time/60.d0

end !termina el main, termina el programa

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine matrix_elements( )
use carga
use matrices
use integracion
implicit none
integer :: i, j, m, n, im, in
real(kind(0.d0)) :: bm, bn, zz

s = 0.d0 ; v01 = 0.d0 ; ke = 0.d0

do i = kord, kord+l-1
  do j = 1, intg
    zz = x(k(i),j)
    call bsplvb(t, kord, 1, zz, i, Sp)

    do m = 1,kord
      im = i-kord+m-1
      if( im>0.and.im<nb+1 ) then
        do n = 1,kord
          in = i-kord+n-1

          if( in>0.and.in<nb+1 )then
            s(im,in) = s(im,in) + sp(m)*sp(n)*w(k(i),j);

            v01(im,in) = v01(im,in) + w(k(i),j)*sp(m)*sp(n)*(exp(-(zz-z0)**2/(2.d0*sigma**2)) );

          endif
        end do
      endif
    end do
  end do
end do

do i = kord, kord+l-1
  do m = i-kord+1, i
    if(m>1.and.m<nb+2)then
      do n = m,i
        if(n<nb+2)then
          do j = 1, intg

            call bder(x(k(i),j), t, bm, bn, kord, nk, m, n, i)

            ke(m-1,n-1) = ke(m-1,n-1) + 0.5d0*w(k(i),j)*bm*bn/me;

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

subroutine sener(archivo)
use carga
use matrices
use integracion
implicit none
integer(4) :: m, dp, ind_V0
integer(4) :: info
real(kind(0.d0)), parameter :: eV = 27.21138564d0
real(kind(0.d0)), dimension(nev) :: e
real(kind(0.d0)), dimension(nb**3,nev) :: auvec
real(kind(0.d0)), dimension(nb**3,nb**3) :: mh, ms, ms_copy, T_cinetica, U_potencial
real(kind(0.d0)), dimension(nb**2,nb**2) :: aux1
real(kind(0.d0)), dimension(nb**3,nb**3) :: aux2
real(kind(0.d0)) :: delta_V0
character(150) :: archivo

dp = nb**3;

!###########################################################
!###########################################################
!###########################################################
open(31, file = archivo, position = 'append')
!###########################################################
!###########################################################
!###########################################################

delta_V0 = (V0_f - V0_i)/real(num_puntos_V0, kind(0.d0));

aux1 = 0.d0; aux2 = 0.d0; 
T_cinetica = 0.d0;

call prod_kron(nb, nb, S, S, aux1);
call prod_kron(nb, nb*nb, ke, aux1, aux2);

T_cinetica = T_cinetica + aux2;

call prod_kron(nb, nb, ke, S, aux1);
call prod_kron(nb, nb*nb, S, aux1, aux2);

T_cinetica = T_cinetica + aux2;

call prod_kron(nb, nb, S, ke, aux1);
call prod_kron(nb, nb*nb, S, aux1, aux2);

T_cinetica = T_cinetica + aux2;

!!!!
U_potencial = 0.d0;
aux1 = 0.d0; aux2 = 0.d0;
call prod_kron(nb, nb, v01, v01, aux1);
call prod_kron(nb, nb*nb, v01, aux1, aux2); 

U_potencial = aux2;

!!!
ms_copy = 0.d0;
aux1 = 0.d0; aux2 = 0.d0;
call prod_kron(nb, nb, S, S, aux1);
call prod_kron(nb, nb*nb, S, aux1, aux2);

ms_copy = aux2;

do ind_V0 = 0, num_puntos_V0

  V0 = V0_i + real(ind_V0, kind(0.d0))*delta_V0;

  mh(:,:) = 0.d0; ms (:,:) = 0.d0;
  
  ms = ms_copy;
  
  mh(:,:) = T_cinetica(:,:) - V0*U_potencial(:,:);
  
  call eigen_value( dp, nev, info, mh, ms, e, auvec)
 
  e(:) = eV*e(:);
  
  !!!!!!!!!!!!!!!!!!!!!!!####################################
  !!!!!!!!!!!!!!!!!!!!!!!####################################
  !!!!!!!!!!!!!!!!!!!!!!!####################################
  write(31,6) eV*V0, (e(m), m = 1, nev)
  !!!!!!!!!!!!!!!!!!!!!!!####################################
  !!!!!!!!!!!!!!!!!!!!!!!####################################
  !!!!!!!!!!!!!!!!!!!!!!!####################################
  call flush(); 
 
end do

6 format(e22.14,1x,1000(1x,e22.14))

return
end subroutine sener

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine prod_kron(NA, NB, A, B, kron)
implicit none
integer, intent(in) :: NA, NB
real(kind(0.d0)), dimension(NA,NA), intent(in) ::A
real(kind(0.d0)), dimensiOn(NB,NB), intent(in) :: B !dummy arguments
real(kind(0.d0)), dimension(NA*NB,NA*NB) :: kron !output matrix of the kronecker product
integer :: i,j !loop counters

do i = 1, NA
  do j = 1, NA
    kron(1+NB*(i-1):NB+NB*(i-1),1+NB*(j-1):NB+NB*(j-1)) = A(i,j)*B
  end do
end do

end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine eigen_value(dp,nev,info,nh,s,x,avec)

!
! calcula usando LAPACK N=nev autovalores

implicit none
integer*4 dp,nev
real(kind(0.d0)),dimension (nev)::x
real(kind(0.d0)),dimension (dp,nev)::avec
real(kind(0.d0)),dimension (dp,dp)::nh,s
! definiciones para DSYGVX
integer(4)::itype,numau,lwork,info,il,iu
character(1)::jobz,range,uplo
real(kind(0.d0))::vl,vu,abstol
real(kind(0.d0)),dimension(dp,dp)::z
real(kind(0.d0)),dimension(9*dp)::work
integer(4),dimension(5*dp)::iwork
integer(4),dimension(dp)::ifail

lwork=9*dp ! lwork=dim(work)

!!!!!!!!!!!!!!!!!!!!!!!!!
! reemplaza  tred2 y tqli del NR por dsygvx de LAPACK
! parametros para dsygvx:
itype=1  ! especifica que A x=lambda Bx
jobz='V' ! V (N) (NO) calcula autovectores
!range='V'! autovalores en el intervalo [VL,VU]
!vl=-1.01d0 ! pongo  vl menor que el menor autovalor
!vu=0.d0   ! por ahora solo interezan los autovalores negativos
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

      GO TO (10,20),index
 10   j = 1
      biatx(1) = 1.d0
      IF (j .GE. jhigh) GO TO 99

 20   CONTINUE
         jp1 = j + 1
         deltar(j) = t(left+j) - x
         deltal(j) = x - t(left+1-j)
         saved = 0.d0
         DO i = 1,j
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

SUBROUTINE KNOTS_PESOS(kord,tip,beta,a,b,c,l,lum,intg,t,k,x,w,pl)

!lo único que hace esta subrutina principal es dividir los casos segun se quiera knots uniformes, exponencial, mixto
!cada una de las subrutinas que llama ademas calcula los x y los w, abscisas y pesos para la cuadratura gaussiana, y los polinomios


!  INPUT:
! kord: orden de los b-splines
! tip; character(1): 'u' ; 'e' ; 'm': dist uniforme, exp o mixta de knots
! beta : param dist. e (no usado si tip='u')
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

implicit none
integer(4)::kord,l,lum,intg
character(1)::tip
real(kind(0.d0))::beta,a,b,c
integer(4),dimension(l+2*kord-1)::k
real(kind(0.d0)),dimension(l+2*kord-1)::t
real(kind(0.d0)),dimension(l,intg)::x,w,pl
!!!!!!


if(tip=='u')then
call dist_unif(kord,a,b,l,intg,t,k,x,w,pl)
elseif(tip=='e')then
call dist_exp(kord,beta,a,b,l,intg,t,k,x,w,pl)
elseif(tip=='m')then
call dist_mix(kord,beta,a,b,c,l,lum,intg,t,k,x,w,pl)
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

implicit none
integer(4)::kord,l,intg
real(kind(0.d0))::a,b
integer(4),dimension(l+2*kord-1)::k
real(kind(0.d0)),dimension(l+2*kord-1)::t
real(kind(0.d0)),dimension(l,intg)::x,w,pl
!!!!!!
integer(4)::i,j,nk
real(kind(0.d0))::ri,rf,dr
real(kind(0.d0)),dimension(intg)::vx,vw,vpl

nk=l+2*kord-1   ! # de knots

dr=(b-a)/dfloat(l)

! calcula los puntos y pesos para la cuadratura
x=0.d0;w=0.d0

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

end    SUBROUTINE DIST_UNIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE DIST_EXP(kord, beta, a, b, l, intg, t, k, x, w, pl)
implicit none
integer(4) :: kord,l,intg
real(kind(0.d0)) :: beta,a,b
integer(4), dimension(l+2*kord-1) :: k
real(kind(0.d0)), dimension(l+2*kord-1) :: t
real(kind(0.d0)), dimension(l,intg) :: x,w,pl
!!!!!!
integer(4) :: i, nk
real(kind(0.d0)) :: ri, dr, length
real(kind(0.d0)), dimension(intg) :: vx,vw,vpl

nk=l+2*kord-1   ! # de knots

dr = (b-a)/real(l, kind(0.d0));

! calcula los puntos y pesos para la cuadratura
x = 0.d0; w = 0.d0;

length = (b-a)/(sign(1.d0,b)*exp(beta*abs(b))-sign(1.d0,a)*exp(beta*abs(a)));

t(1) = a
k(1) = 1

do i = 2, kord
  t(i) = t(1)
  k(i) = 1
end do

do i = kord+1 , kord+l
  ri = (i-kord)*dr + a;
  t(i) = length*sign(1.d0, ri)*exp(beta*abs(ri))
  k(i) = k(i-1)+1
end do

do i = kord+l+1, nk
  t(i)=t(i-1)
  k(i)=k(i-1)
end do

call gauleg_pl(-1.d0, 1.d0, vx, vw, vpl, intg);

do i = 1, l
  x(i, :) = 0.5d0*(t(i+kord)-t(i+kord-1))*vx(:) + 0.5d0*(t(i+kord)+t(i+kord-1))
  w(i, :) = 0.5d0*(t(i+kord)-t(i+kord-1))*vw(:);
end do

pl = 0.d0;

return
end SUBROUTINE DIST_EXP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE DIST_MIX(kord,beta,a,b,c,l,lum,intg,t,k,x,w,pl)
implicit none
integer(4)::kord,l,lum,intg
real(kind(0.d0))::beta,a,b,c
integer(4),dimension(l+2*kord-1)::k
real(kind(0.d0)),dimension(l+2*kord-1)::t
real(kind(0.d0)),dimension(l,intg)::x,w,pl
!!!!!!
integer(4)::i,j,le,nk
real(kind(0.d0)),dimension(lum,intg)::xu,wu,plu
real(kind(0.d0)),dimension(l-lum,intg)::xe,we,ple
integer(4),dimension(lum+2*kord-1)::ku
real(kind(0.d0)),dimension(lum+2*kord-1)::tu
integer(4),dimension(l-lum+2*kord-1)::ke
real(kind(0.d0)),dimension(l-lum+2*kord-1)::te

nk=l+2*kord-1   ! # de knots

call DIST_UNIF(kord,a,c,lum,intg,tu,ku,xu,wu,plu)

le=l-lum

call DIST_EXP(kord,beta,c,b,le,intg,te,ke,xe,we,ple)

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

