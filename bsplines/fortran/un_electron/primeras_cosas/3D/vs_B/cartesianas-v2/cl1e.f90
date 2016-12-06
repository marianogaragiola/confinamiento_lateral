! Programa en C que resuelve el problema de una particula
! en un pozo de potencial usando B-splines.
! Usa el metodo variacional de Rayleight-Ritz usando como base del
! espacio los B-splines, como esta no es una base ortonormal el
! problema de autovalores queda de la siguiente forma:
!
!              H |psi> = e S |psi>
!
! donde H es la matriz del hamiltoniano y S es la matriz de solapamiento
! de los B-splines.
!
! Usa la funcion gauleg() del Numerical Recipies C para calcular los
! puntos de evaluacion y los pesos para la cuadratura.
! La funcion KNOTS_PESOS() esta hecha por mi guiandome de la version de
! fortran, calcula los knots de los B-splines con una distribucion
! uniforme solamente, en el caso de querer otra distribucion es solo
! cuestion de modificar el codigo para usar la distribucion que uno quiera.
! Para el calculo de los autovalores usa la funcion dsygvx_() de lapack.
! Las funciones para evaluar los B-splines y las derivadas son bsplvb() y
! bder() respectivamente, son versiones hechas por mi a partir de las
! versiones de fortran.
!
! EL programa calcula los autovalores para un potencial de la forma:
!
! V(z) = -V0*exp(-0.5*z^2/sigma^2)/(1 + 0.5*l^2/sigma^2)
!
! o sea el pozo es una gaussiano invertido. La profundidad del pozo
! esta moduloda por la intensidad de campo magnetico a travez del
! parametro l dado por:
!
! l = sqrt(2*alpha/B)
!
! donde B es la intesidad del campo y alpha es un parametro que toma el
! valor alpha = 658.4092645439 nm^2/T
!
! Se calculan los autovalores para distintos valores del campo
! magnetico.
!!!!!!!!!!!!!!!!!!!!
#ifndef ZMIN_
#define ZMIN_  0._pr
#endif

#ifndef ZMAX_
#define ZMAX_ 10._pr
#endif

#ifndef TIPO_
#define TIPO_ 1
#endif

#ifndef BETA_
#define BETA_ 0._pr
#endif

#ifndef KORD_
#define KORD_ 5
#endif

#ifndef L_INTERVAL_
#define L_INTERVAL_ 30
#endif

#ifndef N_CUAD_
#define N_CUAD_ 100
#endif

#ifndef NEV_
#define NEV_ 10
#endif

#ifndef ME_
#define ME_ 1._pr
#endif

#ifndef SIGMA_
#define SIGMA_ 1._pr
#endif

#ifndef V0_
#define V0_ 0.5_pr
#endif

#ifndef B_CAMPO_I_
#define B_CAMPO_I_ 5._pr
#endif

#ifndef B_CAMPO_F_
#define B_CAMPO_F_ 10._pr
#endif

#ifndef NUM_PUNTOS_B_
#define NUM_PUNTOS_B_ 10
#endif
!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!
program main
  use precision, pr => dp
  use cuadratura_knots
  use matrices
  use mod_eig
  implicit none
  integer :: tipo, kord, l_interval, n_cuad, nev
  integer :: num_puntos_B
  real(pr) :: zmin, zmax, beta, me, sigma, v0, b_campo_i, b_campo_f
  real(pr), parameter :: a0 = 0.0529177210_pr, eV = 27.21138564_pr
  real(pr), parameter :: c_light = 137.035999074492_pr !, ua_to_T = 2.350517464e5_pr
  real(pr), parameter :: ua_to_T = 1.72e3_pr
  !!!!!!
  integer :: i, nk, nb, ndimh, ind_b
  integer, allocatable :: k(:)
  real(pr) :: delta_b, b_campo, omega
  real(pr), allocatable :: t(:)
  real(pr), allocatable :: x(:,:), w(:,:), pl(:,:)
  real(pr), allocatable :: s(:,:), x_mat(:,:), p_mat(:,:), v_pozo(:,:), v_campo(:,:), ke(:,:)
  complex(pr), allocatable :: h(:,:), ms(:,:), auval(:)
  character(150) :: archivo1, archivo2

  zmin = real(ZMIN_, pr); zmax = real(ZMAX_, pr);
  tipo = TIPO_;
  beta = real(BETA_, pr);
  kord = KORD_;
  l_interval = L_INTERVAL_; n_cuad = N_CUAD_; nev = NEV_;
  me = real(ME_, pr); sigma = real(SIGMA_, pr);
  v0 = real(V0_, pr);
  b_campo_i = real(B_CAMPO_I_, pr); b_campo_f = real(B_CAMPO_F_, pr);
  num_puntos_B = NUM_PUNTOS_B_;

  nk = l_interval+2*kord-1;    ! # de knots
  ! No incluimos primer   y ultimo spline debido a psi(0)=psi(R)=0
  nb = kord+l_interval-3;      ! tamaño base splines
  ndimh = nb**3;

  write(archivo1, '("./resultados/1e-E_vs_B-v0_",f6.4,"eV-real.dat")') v0
  write(archivo2, '("./resultados/1e-E_vs_B-v0_",f6.4,"eV-imag.dat")') v0

  open(10, file = archivo1)
  open(11, file = archivo2)
  write(10,'(A28,f8.2,A1,f8.2,A4)') "# Intervalo de integracion:[", zmin,",", zmax, "] nm"
  write(10,'(A32,x,I2,x,A33)') "# Tipo de distribucion de knots:", tipo, "(1 es uniforme, 2 es exponencial)"
  write(10,'(A47,x,f8.6)') "# Cte de decaimiento en distribucion exp beta =", beta
  write(10,'(A30,x,I2)') "# Orden de los bsplines kord =", kord
  write(10,'(A19,x,I3,x,A24,x,I3)') "# Num de intervalos", l_interval, ", grado de la cuadratura", n_cuad
  write(10,'(A14,x,I3,x,A24,x,I3)') "# Num de knost", nk, ", tamaño de la base nb =", nb
  write(10,'(A33,x,f8.6,x,A2)') "# Masa efectiva del electron me =", me, "UA"
  write(10,'(A24,x,f5.2,x,A2)') "# Ancho del pozo sigma =", sigma, "nm"
  write(10,'(A27,x,f6.4,x,A2)') "# Profundidad del pozo V0 =", v0, "eV"
  write(10,'(A27,x,f6.2,x,A19,x,f7.2)') "# Campo inicial b_campo_i =", b_campo_i, "y final b_campo_f =", b_campo_f
  write(10,'(A24)') "# Autovalores calculados"
  call flush();

  ! Paso a unidades atomicas todo
  v0 = v0/eV;
  zmin = zmin/a0; zmax = zmax/a0; sigma = sigma/a0;
  beta = beta*a0;

  allocate(k(nk), t(nk))
  allocate(x(l_interval,n_cuad), w(l_interval,n_cuad), pl(l_interval,n_cuad))
  allocate(s(nb,nb), x_mat(nb,nb), p_mat(nb,nb), v_pozo(nb,nb), v_campo(nb,nb), ke(nb,nb))
  allocate(h(ndimh,ndimh), ms(ndimh,ndimh))
  allocate(auval(nev))

  !! calculos los knots y los puntos de la cuadratura
  call KNOTS_PESOS(kord, tipo, beta, zmin, zmax, l_interval, n_cuad, t, k, x, w, pl)

  !! ahora paso a calcular las matrices del problema
  s(:,:) = 0._pr; v_pozo(:,:) = 0._pr; v_campo(:,:) = 0._pr; ke(:,:) = 0._pr;
  call calculo_matrices(kord, l_interval, n_cuad, nk, nb, me, sigma, t, k, x, w, s, x_mat, p_mat, v_pozo, v_campo, ke);

  delta_b = (b_campo_f-b_campo_i)/real(num_puntos_B, pr);

  ind_B = 0;
  do ind_B = 0, num_puntos_B

    b_campo = b_campo_i + delta_b*real(ind_B, pr);

    omega = 0.5_pr*b_campo/(me*c_light*ua_to_T); ! asi esta en unidades atomicas

    call hamiltoniano(nb, v0, me, omega, s, x_mat, p_mat, v_pozo, v_campo, ke, h, ms);

    call eigenvalues(ndimh, nev, h, ms, auval)

    auval = eV*auval;

!    write(*,*) b_campo, real(auval(1)), aimag(auval(1))
    write(10,6) b_campo, (real(auval(i)), i = 1, nev)
    write(11,6) b_campo, (aimag(auval(i)), i = 1, nev)
    call flush();
  end do

  deallocate(k, t)
  deallocate(x, w, pl)
  deallocate(s, x_mat, p_mat, v_pozo, v_campo, ke)
  deallocate(h, ms)
  deallocate(auval)

6 format(e22.14,1x,1000(1x,e22.14))

end program
