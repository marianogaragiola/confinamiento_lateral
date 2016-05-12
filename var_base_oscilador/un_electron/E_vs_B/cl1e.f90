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
  integer, parameter :: N_base = 100, N_cuad = 200
  integer :: n, m, N_B, ind_B
  real(pr), parameter :: alpha = 658.4092645439_pr, a0 = 0.0529177210_pr, eV = 27.21138564_pr
  real(pr) :: me, V0, omega, sigma, B_campo_i, B_campo_f, B_campo
  real(pr) :: b, l, delta_B
  real(pr) :: x(N_cuad), w(N_cuad), Hermite(N_cuad,N_base+1)
  real(pr) :: T(N_base+1,N_base+1), V(N_base+1,N_base+1), H(N_base+1,N_base+1)
  real(pr) :: e(N_base+1)

!!cargo los parametros para el problema
  me = real(me_, pr);
  V0 = real(V0_, pr);
  B_campo_i = real(B_campo_i_, pr); B_campo_f = real(B_campo_f_, pr);
  N_B = N_B_;
  sigma = real(sigma_, pr);

!!calculo los parametros dependientes
  b = 0.5_pr*(1._pr/sigma)**2;

  omega = 0.02;

  delta_B = (B_campo_f-B_campo_i)/real(N_B-1, pr);

  V0 = V0/eV;
!!calculo los puntos para la cuadratura de Gauss-Hermite y lo polinomios en esos puntos
  call cpquad(N_cuad, 1._pr, "Hermite", w, x);


  call pol_hermite(N_cuad, N_base, 1._pr, x, Hermite);

!!calulo la energia cinetica
  call energia_cinetica(N_base, me, T);

  open(11, file='./resultados/auval_vs_B-3.dat')
  write(11,*) '# me =', me
  write(11,*) '# V0 =', V0*eV
  write(11,*) '# B_i =', B_campo_i, 'B_f =', B_campo_f
  write(11,*) '# N_B =', N_B
  write(11,*) '# omega =', omega
  write(11,*) '# sigma =', sigma

  do ind_B = 1, N_B

    B_campo = B_campo_i + real(ind_B-1, pr)*delta_B;
    l = sqrt(2._pr*alpha/B_campo);

!!!!calculo la energia potencial
    call energia_potencial(N_base, N_cuad, omega, sigma, x, w, Hermite, V);

!!!!calculo el hamiltoniano
    H(:,:) = a0**2*omega*T(:,:) - V0*V(:,:)/(l**2*b+1._pr);

!!!!calculo los autovalores y autovectores
    call autovalores(N_base+1, H, e);
    e(:) = eV*e(:);

    write(11,6) B_campo, (e(n), n = 1, 30)

  end do
  close(11)

6 format(e22.14,1x,1000(1x,e22.14))

end program cl1e
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
