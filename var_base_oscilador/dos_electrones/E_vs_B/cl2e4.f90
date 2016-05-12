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

#ifndef lambda_
#define lambda_ 1.d0
#endif

#ifndef omega_
#define omega_ 0.02d0
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
  integer, parameter :: N_base = 29, N_cuad = 100
  integer :: N_dim = ((N_base+1)*(N_base+2))/2
  integer, allocatable :: idx(:,:), idnm(:,:), idsim(:,:)
  integer :: n, m, N_B, ind_B
  real(pr), parameter :: alpha = 658.4092645439_pr, a0 = 0.0529177210_pr, eV = 27.21138564_pr
  real(pr) :: me, V0, omega, sigma, B_campo_i, B_campo_f, B_campo
  real(pr) :: b, l, delta_B, lambda
  real(pr), allocatable :: x(:), w(:), Hermite(:,:)
  real(pr), allocatable :: T(:,:), V(:,:), H_e(:,:)
  real(pr), allocatable :: V_ef(:,:), H(:,:), e(:)
  character(1) :: l1, l2, l3, l4

!!cargo los parametros para el problema
  me = real(me_, pr);
  V0 = real(V0_, pr);
  B_campo_i = real(B_campo_i_, pr); B_campo_f = real(B_campo_f_, pr);
  N_B = N_B_;
  sigma = real(sigma_, pr);
  lambda = real(lambda_, pr);
  omega = real(omega_, pr);

!!calculo los parametros dependientes
  b = 0.5_pr*(1._pr/sigma)**2;
  l = sqrt(2._pr*alpha/B_campo);

  l1 = char(modulo(int(1000._pr*lambda), 10) + 48);
  l2 = char(modulo(int(1000._pr*lambda), 100)/10 + 48);
  l3 = char(modulo(int(1000._pr*lambda), 1000)/100 + 48);
  l4 = char(modulo(int(1000._pr*lambda), 10000)/1000 + 48);

!  omega = 0.02_pr;

  delta_B = (B_campo_f-B_campo_i)/real(N_B-1, pr);

  V0 = V0/eV;

  allocate(idx((N_base+1)**2,3), idnm(N_base+1,N_base+1), idsim(N_dim,3))
  allocate(x(N_cuad), w(N_cuad))
  allocate(Hermite(N_cuad,N_base+1))
  allocate(T(N_base+1,N_base+1), V(N_base+1,N_base+1), H_e(N_base+1,N_base+1))
  allocate(V_ef((N_base+1)**2,(N_base+1)**2))
  allocate(H(N_dim,N_dim), e(N_dim))

  call indices(N_base, N_dim, idx, idnm, idsim);

!!calculo los puntos para la cuadratura de Gauss-Hermite y lo polinomios en esos puntos
  call cpquad(N_cuad, 1._pr, "Hermite", w, x);


  call pol_hermite(N_cuad, N_base, 1._pr, x, Hermite);

!!calulo la energia cinetica
  call energia_cinetica(N_base, me, T);

!!calculo la energia potencial
  call energia_potencial(N_base, N_cuad, omega, sigma, x, w, Hermite, V);


  open(11, file='./resultados/auval_vs_B-lambda'//l4//l3//l2//l1//'.dat')
  write(11,*) '# me =', me
  write(11,*) '# V0 =', V0*eV
  write(11,*) '# B_i =', B_campo_i, 'B_f =', B_campo_f
  write(11,*) '# N_B =', N_B
  write(11,*) '# omega =', omega
  write(11,*) '# sigma =', sigma
  write(11,*) '# lambda =', lambda

  do ind_B = 1, N_B

    B_campo = B_campo_i + real(ind_B-1, pr)*delta_B;

    l = sqrt(2._pr*alpha/B_campo);

!!!!calculo el hamiltoniano
    H_e(:,:) = a0**2*omega*T(:,:) - V0*V(:,:)/(l**2*b+1._pr);

    call interaccion(N_base, N_cuad, idx, omega, l, x, w, Hermite, V_ef);

    call hamiltoniano(N_base, N_dim, lambda/a0, idnm, idsim, H_e, V_ef, H);

!!!!calculo los autovalores y autovectores
    call autovalores(N_dim, H, e);
    e(:) = eV*e(:);

    write(11,6) B_campo, (e(n), n = 1, 30)

  end do
  close(11)

6 format(e22.14,1x,1000(1x,e22.14))

  deallocate(idx, idnm, idsim)
  deallocate(x, w, Hermite)
  deallocate(T, V, H_e)
  deallocate(V_ef)
  deallocate(H, e)


end program cl1e
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine indices(N_base, N_dim, idx, idnm, idsim)
  implicit none
  integer, intent(in) :: N_base, N_dim
  integer, intent(out) :: idx((N_base+1)**2,3), idnm(N_base+1,N_base+1)
  integer, intent(out) :: idsim(N_dim,3)
  !!!!
  integer :: n, m, ind

  ind = 1;
  do n = 1, N_base+1
    do m = 1, N_base+1
      idx(ind,1) = ind; idx(ind,2) = n; idx(ind,3) = m
      idnm(n,m) = ind;
      ind = ind + 1;
    end do
  end do

  ind = 1;
  do n = 1, N_base+1
    do m = 1, n
      idsim(ind,1) = ind; idsim(ind,2) = n; idsim(ind,3) = m;
      ind = ind + 1;
    end do
  end do

end subroutine indices
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine hamiltoniano(N_base, N_dim, lambda, idnm, idsim, H_e, V_ef, H)
  use precision
  implicit none
  integer, intent(in) :: N_base, N_dim
  integer, intent(in) :: idnm(N_base+1,N_base+1), idsim(N_dim,3)
  real(pr), intent(in) :: lambda
  real(pr), intent(in) :: H_e(N_base+1,N_base+1), V_ef((N_base+1)**2,(N_base+1)**2)
  real(pr), intent(out) :: H(N_dim,N_dim)
  !!!!
  integer :: n1, n2, m1, m2, ind1, ind2
  real(pr) :: sqrt2, kronecker

  sqrt2 = sqrt(0.5_pr);

  H(:,:) = 0._pr;

  do ind2 = 1, N_dim

    n2 = idsim(ind2, 2); m2 = idsim(ind2, 3);

    do ind1 = 1, N_dim

      n1 = idsim(ind1, 2); m1 = idsim(ind1, 3);

      if(m1==n1 .and. m2==n2) then

        H(ind1,ind2) = 2._pr*H_e(n1,n2)*kronecker(n1,n2) + lambda*V_ef(idnm(n1,n1),idnm(n2,n2));

      elseif(m1==n1 .and. m2.ne.n2) then

        H(ind1,ind2) = sqrt2*(2._pr*(H_e(n1,n2)*kronecker(n1,m2)+H_e(n1,m2)*kronecker(n1,n2)) +&
                     & lambda*(V_ef(idnm(n1,n1),idnm(n2,m2)) + V_ef(idnm(n1,n1),idnm(m2,n2)) ));

      elseif(m1.ne.n1 .and. m2==n2) then

        H(ind1,ind2) = sqrt2*(2._pr*(H_e(n1,n2)*kronecker(m1,n2)+H_e(m1,n2)*kronecker(n1,n2)) +&
                     & lambda*(V_ef(idnm(n1,m1),idnm(n2,n2)) + V_ef(idnm(m1,n1),idnm(n2,n2)) ));

      elseif(m1.ne.n1 .and. m2.ne.n2) then

        H(ind1,ind2) = H_e(m1,m2)*kronecker(n1,n2) + H_e(n1,n2)*kronecker(m1,m2) + &
                     & H_e(m1,n2)*kronecker(n1,m2) + H_e(n1,m2)*kronecker(m1,n2) + &
                     & 0.5_pr*lambda*(V_ef(idnm(n1,m1),idnm(n2,m2)) + V_ef(idnm(n1,m1),idnm(m2,n2)) + &
                                   &  V_ef(idnm(m1,n1),idnm(n2,m2)) + V_ef(idnm(m1,n1),idnm(n2,m2)) );

      endif

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
  do m = 0, N_base
    do n = 0, m

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

  pot(:) = exp(-0.5_pr/omega*(x(:)/sigma)**2);

  V(:,:) = 0._pr
  do m = 1, N_base+1
    do n = m, N_base+1, 2

      base(:) = w(:)*Hermite(:,n)*Hermite(:,m);

      V(n,m) = dot_product(base, pot);
      V(m,n) = V(n,m);

    end do
  end do

end subroutine energia_potencial
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine interaccion(N_base, N_cuad, idx, omega, l, x, w, Hermite, V_ef)
  use precision
  implicit none
  integer, intent(in) :: N_base, N_cuad
  integer, intent(in) :: idx((N_base+1)**2,3)
  real(pr), intent(in) :: omega, l, x(N_cuad), w(N_cuad)
  real(pr), intent(in) :: Hermite(N_cuad,N_base+1)
  real(pr), intent(out) :: V_ef((N_base+1)**2,(N_base+1)**2)
  !!!!!
  integer :: j, n1, n2, m1, m2, ind1, ind2, N_dim
  real(pr), parameter :: eV = 27.21138564_pr, a0 = 0.0529177210_pr
  real(pr) :: pot_int(N_cuad,N_cuad)
  real(pr) :: phi_n1n2(N_cuad), phi_m1m2(N_cuad)

  N_dim = (N_base+1)**2;

  call potencial_interaccion(N_cuad, omega, l, x, pot_int);

  V_ef(:,:) = 0._pr;

  open(21,file='matriz_int.dat')
  do ind2 = 1, N_dim
    do ind1 = 1, N_dim !ind2

      n1 = idx(ind1, 2); m1 = idx(ind1, 3);
      n2 = idx(ind2, 2); m2 = idx(ind2, 3);

      phi_n1n2(:) = w(:)*Hermite(:,n1)*Hermite(:,n2);
      phi_m1m2(:) = w(:)*Hermite(:,m1)*Hermite(:,m2);

      V_ef(ind1,ind2) = dot_product(phi_n1n2, matmul(pot_int, phi_m1m2));
!      V_ef(ind2,ind1) = V_ef(ind1,ind2);

      write(21,*) n1-1, m1-1, n2-1, m2-1, a0*V_ef(ind1, ind2)

    end do
  end do
  V_ef(:,:) = V_ef(:,:)*a0;

  close(21)

  ! do n1 = 1, N_base+1
  !   do n2 = 1, N_base+1
  !
  !     phi_n1n2(:) = w(:)*Hermite(:,n1)*Hermite(:,n2);
  !
  !     do m1 = 1, N_base+1
  !       do m2 = 1, N_base+1
  !
  !         phi_m1m2(:) = w(:)*Hermite(:,m1)*Hermite(:,m2);
  !
  !         V_ef((n1-1)*(N_base+1)+m1,(n2-1)*(N_base+1)+m2) = dot_product(phi_n1n2, matmul(pot_int, phi_m1m2));
  !
  !       end do
  !     end do
  !   end do
  ! end do

end subroutine interaccion
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine potencial_interaccion(N_cuad, omega, l, x, pot_int);
  use precision
  implicit none
  integer, intent(in) :: N_cuad
  real(pr), intent(in) :: omega, l, x(N_cuad)
  real(pr), intent(out) :: pot_int(N_cuad,N_cuad)
  !!!!
  integer :: i, j
  real(pr), parameter :: eV = 27.21138564_pr, a0 = 0.0529177210_pr
  real(pr) :: pi, cte, y(N_cuad,N_cuad)

  pot_int(:,:) = 0._pr;
  pi = 2._pr*asin(1._pr);
  cte = sqrt(0.5_pr*pi)/l;

  y(:,:) = 0._pr;
  do j = 1, N_cuad
  do i = 1, N_cuad
    ! y(i,j) = (x(i)-x(j))/(sqrt(2._pr*omega)*l);
    y(i,j) = abs(x(i)-x(j))/(sqrt(2._pr*omega)*l);
  end do
  end do

  do j = 1, N_cuad
    do i = 1, N_cuad

      if( abs(y(i,j)) < 4._pr) then
        ! pot_int(i,j) = cte*exp(y(i,j)**2)*(1._pr - erf(abs(y(i,j))));
        pot_int(i,j) = cte*exp(y(i,j)**2)*(1._pr - erf(y(i,j)));
      else
        ! pot_int(i,j) = 1._pr/l*(sqrt(0.5_pr)/abs(y(i,j)) - (sqrt(0.5_pr)/abs(y(i,j)))**3);
        pot_int(i,j) = 1._pr/l*(sqrt(0.5_pr)/y(i,j) - (sqrt(0.5_pr)/(y(i,j)))**3);
      end if

      pot_int(j,i) = pot_int(i,j);

    end do
  end do

  ! pot_int(:,:) = pot_int(:,:)*a0;

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
