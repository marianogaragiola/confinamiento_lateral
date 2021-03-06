!-------------------------------------------------------------------!
!Codigo cl1e.f90 confinamiento lateral un electron:                 !
!Calcula la matriz del hamiltoniano usando como base las            !
!autofunciones del oscilador armonico con frecuencia omega.         !
!                                                                   !
!El parametro omega es un parametro variacional no lineal           !
!que uso para mejorar la aproximacion.                              !
!-------------------------------------------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Variables de preprocesador
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifndef me_
#define me_ 1.d0
#endif

#ifndef V0_
#define V0_ 1.d0
#endif

#ifndef B_campo_i_
#define B_campo_i_ 1.d0
#endif

#ifndef B_campo_f_
#define B_campo_f_ 2.d0
#endif

#ifndef N_B_
#define N_B_ 10
#endif

#ifndef sigma_
#define sigma_ 1.d0
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

include 'mpfun.f90'
include 'mpmod.f90'
include 'gaussquad.f90'

module precision
  implicit none
  integer, parameter :: pr = kind(0.d0)
end module precision

program cl1e
  use precision
  use gaussquad
  implicit none
  integer, parameter :: N_base = 100, N_cuad = 50
  integer :: N_dim = ((N_base+1)*(N_base+2))/2
  integer :: n, m, N_B, ind_B
  real(pr), parameter :: alpha = 658.4092645439_pr, a0 = 0.0529177210_pr, eV = 27.21138564_pr
  real(pr) :: me, V0, omega, sigma, B_campo_i, B_campo_f, B_campo
  real(pr) :: b, l, delta_B, lambda
  real(pr) :: x(N_cuad), w(N_cuad), Hermite(N_cuad,N_base+1)
  real(pr) :: T(N_base+1,N_base+1), V(N_base+1,N_base+1), H_e(N_base+1,N_base+1)
  real(pr) :: V_ef(N_base+1,N_base+1,N_base+1,N_base+1)
  real(pr), allocatable :: e(:), H(:,:) 

!!cargo los parametros para el problema
  me = real(me_, pr);
  V0 = real(V0_, pr);
  B_campo_i = real(B_campo_i_, pr); B_campo_f = real(B_campo_f_, pr);
  N_B = N_B_;
  sigma = real(sigma_, pr);

!!calculo los parametros dependientes
  b = 0.5_pr*(1._pr/sigma)**2;
  omega = 1._pr;!0.07_pr;

  delta_B = (B_campo_f-B_campo_i)/real(N_B-1, pr);

  lambda = 1._pr;
  V0 = V0/eV;
!!calculo los puntos para la cuadratura de Gauss-Hermite y lo polinomios en esos puntos
  call cpquad(N_cuad, 1._pr, "Hermite", w, x);


  call pol_hermite(N_cuad, N_base, omega, x, Hermite);

!!calulo la energia cinetica
  call energia_cinetica(N_base, me, T);

!!!!calculo la energia potencial
  call energia_potencial(N_base, N_cuad, omega, sigma, x, w, Hermite, V);


  allocate(H(N_dim,N_dim), e(N_dim))

  open(11, file='./resultados/auval_vs_B.dat')
!  do ind_B = 1, N_B

    B_campo = B_campo_i + real(ind_B-1,pr)*delta_B;

    l = sqrt(2._pr*alpha/B_campo);
   
!!!!calculo el hamiltoniano
    H_e(:,:) = a0**2*omega*T(:,:) - V0*V(:,:)/(l**2*b+1._pr);

    call interaccion(N_base, N_cuad, omega, l, x, w, Hermite, V_ef);

    call hamiltoniano(N_base, N_dim, lambda, H_e, V_ef, H);

!!!!calculo los autovalores y autovectores
    call autovalores(N_dim, H, e);
    e(:) = eV*e(:);

    write(11,6) B_campo, (e(n), n = 1, 30)

!  end do
  close(11)

  deallocate(H, e)

6 format(e22.14,1x,1000(1x,e22.14))

end program cl1e
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine hamiltoniano(N_base, N_dim, lambda, H_e, V_ef, H)
  use precision
  implicit none
  integer, intent(in) :: N_base, N_dim
  real(pr), intent(in) :: lambda
  real(pr), intent(in) :: H_e(N_base+1,N_base+1), V_ef(N_base+1,N_base+1,N_base+1,N_base+1)
  real(pr), intent(out) :: H(N_dim,N_dim)
  !!!!
  integer :: n1, n2, m1, m2, ind1, ind2
  real(pr) :: sqrt2, kronecker

  sqrt2 = sqrt(0.5_pr);

  H(:,:) = 0._pr;
  ind1 = 1;
  do n1 = 1, N_base+1
    do m1 = 1, n1
 
      ind2 = 1;
      do n2 = 1, N_base+1
        do m2 = 1, n2

          if(m1==n1 .and. m2==n2) then

            H(ind1,ind2) = 2._pr*H_e(n1,n2)*kronecker(n1,n2) + lambda*V_ef(n1,n1,n2,n2);

          elseif(m1==n1 .and. m2.ne.n2) then

            H(ind1,ind2) = sqrt2*(2._pr*(H_e(n1,n2)*kronecker(n1,m2)+H_e(n1,m2)*kronecker(n1,n2)) +&
                         & lambda*(V_ef(n1,n1,n2,m2) + V_ef(n1,n1,m2,n2)));

          elseif(m1.ne.n1 .and. m2==n2) then

            H(ind1,ind2) = sqrt2*(2._pr*(H_e(n1,n2)*kronecker(m1,n2)+H_e(m1,n2)*kronecker(n1,n2)) +&
                         & lambda*(V_ef(n1,m1,n2,n2) + V_ef(m1,n1,n2,n2)));

          elseif(m1.ne.n1 .and. m2.ne.n2) then

            H(ind1,ind2) = H_e(m1,m2)*kronecker(n1,n2) + H_e(n1,n2)*kronecker(m1,m2) + &
                         & H_e(m1,n2)*kronecker(n1,m2) + H_e(n1,m2)*kronecker(m1,n2) + & 
                         & 0.5_pr*lambda*(V_ef(n1,m1,n2,m2) + V_ef(n1,m1,m2,n2) + &
                                       &  V_ef(m1,n1,n2,m2) + V_ef(m1,n1,n2,m2) );

          endif

          ind2 = ind2 + 1;

        end do
      end do

      ind1 = ind1 + 1;
    end do
  end do 

end subroutine hamiltoniano
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine energia_cinetica(N_base, me, T)
  use precision
  implicit none
  integer, intent(in) :: N_base
  real(pr), intent(in) :: me
  real(pr), intent(out) :: T(N_base+1,N_base+1)
  !!!!!!!
  integer :: n, m
  real(pr) :: kronecker

  T(:,:) = 0._pr;
  do n = 0, N_base
    do m = 0, n

      T(n+1,m+1) = (real(n,pr)+0.5_pr)*kronecker(n,m)-0.5_pr*(sqrt(real(n*(n-1),pr))*kronecker(n-1,m+1) +&
                                                  & sqrt(real((n+1)*(n+2),pr))*kronecker(n+1,m-1));
      T(m+1,n+1) = T(n+1,m+1);

    end do
  end do
  T(:,:) = 0.5_pr/me*T(:,:);


end subroutine energia_cinetica
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine energia_potencial(N_base, N_cuad, omega, sigma, x, w, Hermite, V)
  use precision
  implicit none
  integer, intent(in) :: N_base, N_cuad
  real(pr), intent(in) :: omega, sigma
  real(pr), intent(in) :: x(N_cuad), w(N_cuad), Hermite(N_cuad,N_base+1)
  real(pr), intent(out) :: V(N_base+1,N_base+1)
  !!!!!!!
  integer :: i, n, m
  real(pr) :: pot(N_cuad), base(N_cuad)

  pot(:) = 0._pr;
  do i = 1, N_cuad
    pot(i) = exp(-0.5_pr/omega*(x(i)/sigma)**2);
  end do

  V(:,:) = 0._pr
  do n = 1, N_base+1
    do m = 1, n 

      if(mod(n+m,2)==0) then

      base(:) = w(:)*Hermite(:,n)*Hermite(:,m);
      
      V(n,m) = dot_product(base, pot);
      V(m,n) = V(n,m);

      end if

    end do
  end do

end subroutine energia_potencial
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine interaccion(N_base, N_cuad, omega, l, x, w, Hermite, V_ef)
  use precision
  implicit none
  integer, intent(in) :: N_base, N_cuad
  real(pr), intent(in) :: omega, l, x(N_cuad), w(N_cuad)
  real(pr), intent(in) :: Hermite(N_cuad,N_base+1)
  real(pr), intent(out) :: V_ef(N_base+1,N_base+1,N_base+1,N_base+1)
  !!!!!
  integer :: n1, n2, m1, m2
  real(pr) :: pot_int(N_cuad,N_cuad)
  real(pr) :: phi_n1n2(N_cuad), phi_m1m2(N_cuad)


  call potencial_interaccion(N_cuad, omega, l, x, pot_int);

  V_ef(:,:,:,:) = 0._pr;

  do n1 = 1, N_base+1
    do n2 = 1, N_base+1

      phi_n1n2(:) = w(:)*Hermite(:,n1)*Hermite(:,n2);

      do m1 = 1, N_base+1
        do m2 = 1, N_base+1

          phi_m1m2(:) = w(:)*Hermite(:,m1)*Hermite(:,m2);

          V_ef(n1,m1,n2,m2) = dot_product(phi_n1n2, matmul(pot_int, phi_m1m2));         

        end do
      end do
    end do
  end do
  
end subroutine interaccion
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine potencial_interaccion(N_cuad, omega, l, x, pot_int);
!  use mpmodule
  use precision
  implicit none
  integer, intent(in) :: N_cuad
  real(pr), intent(in) :: omega, l, x(N_cuad)
  real(pr), intent(out) :: pot_int(N_cuad,N_cuad)
  !!!!
  integer :: i, j
  real(pr) :: pi, cte, y(N_cuad,N_cuad)
!  type(mp_real) :: xij, V_pot
!  call mpinit

  pot_int(:,:) = 0._pr;
!  pi = 2._pr*asin(1._pr);
!  cte = sqrt(0.5_pr*pi)/l;
!
!  y(:,:) = 0._pr;
!  do i = 1, N_cuad
!    y(i,:) = abs(x(i)-x(:))/(sqrt(2._pr*omega)*l);
!  end do
!
!  V_pot = mpreal(0._pr);
!  do i = 1, N_cuad
!    do j = 1, N_cuad
!      
!      xij = mpreal(y(i,j));
!
!      V_pot = exp(xij)*(mpreal(1._pr)-mpreal(erf(y(i,j))));
!
!      pot_int(i,j) = dble(V_pot);
!    end do
!  end do
!
!  pot_int(:,:) = cte*pot_int(:,:);

end subroutine potencial_interaccion
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine autovalores(N, A, e)
  !Calculamos los autovalores de una matriz A (NxN) real y simetrica
  !usando la subroutine dsyev.f de Lapack.
  use precision
  implicit none
  integer, intent(in) :: N
  real(pr), dimension(1:N,1:N), intent(inout) :: A
  real(pr), dimension(1:N), intent(out) :: e
  !!!!! Parametros para dsyev.f !!!!
  character(1) :: JOBZ, UPLO
  integer :: LDA, INFO, LWORK
  real(pr) :: W(1:N)
  real(pr), allocatable :: WORK(:)

  JOBZ = 'V';
  UPLO = 'U';
  LDA = N;
  LWORK = 3*N-1;

  allocate(WORK(1:LWORK))

  call DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )

  if(0 .ne. INFO) then
    write(*,*) 'Error en dsyev INFO = ', INFO
    stop
  end if

  e(:) = W(:);

  deallocate(WORK)

end subroutine autovalores
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine pol_hermite(N_cuad, N_dim, alpha, x, Her)
  use precision
  use mpmodule
  implicit none
  integer, intent(in) :: N_cuad, N_dim
  real(pr), intent(in) :: x(N_cuad), alpha
  real(pr), intent(out) :: Her(N_cuad,N_dim+1)
  !!!!!!!!!!!!!!
  integer :: j, n
  type(mp_real) :: xj, y(N_dim+1), dy, d2y
  call mpinit

  do j = 1, N_cuad
    xj = mpreal(dble(sqrt(alpha)*x(j)));

    call vahepo(N_dim+1, mpreal(dble(alpha)), xj, y);

    do n = 1, N_dim+1
      Her(j,n) = dble(y(n));
    end do

  end do

end subroutine pol_hermite
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE VAHEPO(N, ALPHA, X, Y)
! *************************************************************
! *   COMPUTES THE VALUE OF THE HERMITE POLYNOMIAL OF DEGREE N
! *   IN THE GIVEN POINT X
! *   N  = DEGREE OF THE POLYNOMIAL
! *   X  = POINT IN WHICH THE COMPUTATION IS PERFORMED
! *   Y  = VALUE OF THE POLYNOMIAL IN X
! *************************************************************
use mpmodule
implicit none
integer :: N
integer :: K
type(mp_real) :: ALPHA, X, Y(N)
type(mp_real) :: PI
call mpinit

PI = mpreal(2.D0)*ASIN(mpreal(1.D0))

DO K = 0, N-1
   IF(K .EQ. 0) THEN
     Y(K+1) = SQRT(SQRT(ALPHA/PI))
   ELSE IF(K .EQ. 1) THEN
     Y(K+1) = SQRT(mpreal(2.D0))*SQRT(SQRT(ALPHA/PI))*X
   ELSE IF(K .GE. 2) THEN
     Y(K+1) = SQRT(MPREAL(2)/MPREAL(K))*X*Y(K)-SQRT(MPREAL(K-1)/MPREAL(K))*Y(K-1)
   END IF
END DO

RETURN
END SUBROUTINE VAHEPO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function kronecker(n,m)
  use precision
  implicit none
  integer, intent(in) :: n, m
  real(pr) :: kronecker

  if(-1 == n .or. -1 == m) then
    kronecker = 0._pr;
  else
    if(n.eq.m) then
      kronecker = 1._pr;
    else
      kronecker = 0._pr;
    end if
  end if

end function kronecker
