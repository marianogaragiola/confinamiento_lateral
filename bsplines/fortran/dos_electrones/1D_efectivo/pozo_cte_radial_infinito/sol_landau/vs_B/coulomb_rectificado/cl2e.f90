
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

#ifndef AZ_
#define AZ_ 2._pr
#endif

#ifndef BZ_
#define BZ_ 5._pr
#endif

#ifndef R0_
#define R0_ 2._pr
#endif

#ifndef RMAX_
#define RMAX_ 50._pr
#endif

#ifndef V1_
#define V1_ 0._pr
#endif

#ifndef V2_
#define V2_ 1._pr
#endif

#ifndef ETA_
#define ETA_ 1._pr
#endif

#ifndef DELTA_
#define DELTA_ 0.5_pr
#endif

#ifndef BCAMPO_I_
#define BCAMPO_I_ 1._pr
#endif

#ifndef BCAMPO_F_
#define BCAMPO_F_ 10._pr
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
  integer :: num_puntos_b
  real(pr) :: zmin, zmax, beta, me, v1, v2
  real(pr) :: az, bz, r0, rmax, delta, eta
  real(pr) :: bcampo_i, bcampo_f
  real(pr) :: vr, vrmax, omega, slope_ll
  real(pr), parameter :: a0 = 0.0529177210_pr, eV = 27.21138564_pr
  real(pr), parameter :: alpha = 658.4092645439_pr, c_light = 137.035999074492_pr
  real(pr), parameter :: ua_to_T = 1.72d3
  !!!!!!
  integer :: i, nk, nb, ndimh, ind_b
  integer, allocatable :: k(:)
  real(pr) :: delta_b, bcampo, l_campo
  real(pr), allocatable :: t(:)
  real(pr), allocatable :: x(:,:), w(:,:), pl(:,:)
  real(pr), allocatable :: s(:,:), v01(:,:), v02(:,:), ke(:,:)
  real(pr), allocatable :: V_int(:,:,:,:)
  real(pr), allocatable :: h(:,:), ms(:,:), mv(:,:), auval(:), auvec(:,:)
  real(pr), allocatable :: exp_int(:,:)
  character(150) :: file_auval, file_exp

  zmin = real(ZMIN_, pr); zmax = real(ZMAX_, pr);
  tipo = TIPO_;
  beta = real(BETA_, pr);
  kord = KORD_;
  l_interval = L_INTERVAL_; n_cuad = N_CUAD_; nev = NEV_;
  me = real(ME_, pr);
  az = real(AZ_, pr); bz = real(BZ_, pr); r0 = real(R0_, pr); rmax = real(RMAX_, pr);
  v1 = real(V1_, pr); v2 = real(V2_, pr);
  delta = real(DELTA_, pr); eta = real(ETA_, pr);
  bcampo_i = real(BCAMPO_I_, pr); bcampo_f = real(BCAMPO_F_, pr);
  num_puntos_b = NUM_PUNTOS_B_;

  nk = l_interval+2*kord-1;    ! # de knots
  ! No incluimos primer   y ultimo spline debido a psi(0)=psi(R)=0
  nb = kord+l_interval-3;      ! tamaño base splines
  ndimh = nb*(nb+1)/2;

  slope_ll = eV/(me*c_light*ua_to_T)

  write(file_auval, '("./res21032017/2e-E_vs_B-v1_",f6.4,"eV-v2_",f6.4,"eV-az_",f6.4,"-bz_",f6.4,"-r0_",f6.4,"&
  &-rmax_",f6.4,"-eta_",f8.6,"-delta_",f6.4,".dat")') v1, v2, az, bz, r0, rmax, eta, delta
  write(file_exp, '("./res21032017/2e-Int_vs_B-v1_",f6.4,"eV-v2_",f6.4,"eV-az_",f6.4,"-bz_",f6.4,"-r0_",f6.4"&
  &-rmax_",f6.4,"-eta_",f8.6,"-delta_",f6.4,".dat")') v1, v2, az, bz, r0, rmax, eta, delta
  open(10, file = file_auval)
  open(11, file = file_exp)
  write(10,'(A28,f8.2,A1,f8.2,A4)') "# Intervalo de integracion:[", zmin,",", zmax, "] nm"
  write(10,'(A32,x,I2,x,A33)') "# Tipo de distribucion de knots:", tipo, "(1 es uniforme, 2 es exponencial)"
  write(10,'(A47,x,f8.6)') "# Cte de decaimiento en distribucion exp beta =", beta
  write(10,'(A30,x,I2)') "# Orden de los bsplines kord =", kord
  write(10,'(A19,x,I3,x,A24,x,I3)') "# Num de intervalos", l_interval, ", grado de la cuadratura", n_cuad
  write(10,'(A14,x,I3,x,A24,x,I3)') "# Num de knost", nk, ", tamaño de la base nb =", nb
  write(10,'(A33,x,f8.6,x,A2)') "# Masa efectiva del electron me =", me, "UA"
  write(10,'(A22,x,f6.4,x,A2)') "# Altura del pozo V1 =", v1, "eV"
  write(10,'(A27,x,f6.4,x,A2)') "# Profundidad del pozo V2 =", v2, "eV"
  write(10,'(A22,x,f6.4,x,A9,x,f6.4,x,A2)') "# Radios del pozo az =", az, "nm y bz =", bz, "nm"
  write(10,'(A41,x,f6.4,x,A2)') "# Cte de rectificacion de Coulomb delta =", delta, "nm"
  write(10,'(A31,x,f6.4,x,A11,x,f6.4,x,A2)') "# Radios de la parte radia r0 =", r0, "nm y Rmax =", Rmax, "nm"
  write(10,'(A22,x,f6.4)') "# Carga efectiva eta =", eta
  write(10,'(A27,x,f6.2,x,A19,x,f7.2)') "# Campo inicial b_campo_i =", bcampo_i, "y final b_campo_f =", bcampo_f
  write(10,'(A24)') "# Autovalores calculados"
  call flush();

  ! Paso a unidades atomicas todo
  v1 = v1/eV; v2 = v2/eV;
  zmin = zmin/a0; zmax = zmax/a0;
  beta = beta*a0; delta = delta/a0;
  az = az/a0; bz = bz/a0; r0 = r0/a0; rmax = rmax/a0;

  allocate(k(nk), t(nk))
  allocate(x(l_interval,n_cuad), w(l_interval,n_cuad), pl(l_interval,n_cuad))
  allocate(s(nb,nb), v01(nb,nb), v02(nb,nb), ke(nb,nb))
  allocate(V_int(nb,nb,nb,nb))
  allocate(h(ndimh,ndimh), ms(ndimh,ndimh), mv(ndimh,ndimh))
  allocate(auval(nev), auvec(ndimh,nev))
  allocate(exp_int(nev,nev))

  !! calculos los knots y los puntos de la cuadratura
  call KNOTS_PESOS(kord, tipo, beta, zmin, zmax, l_interval, n_cuad, t, k, x, w, pl)

  !! ahora paso a calcular las matrices del problema
  s(:,:) = 0._pr; v01(:,:) = 0._pr; v02(:,:) = 0._pr; ke(:,:) = 0._pr;
  call calculo_matrices(kord, l_interval, n_cuad, nk, nb, me, az, bz, t, k, x, w, s, v01, v02, ke);

  call interaccion(kord, l_interval, n_cuad, nk, nb, delta, k, t, x, w, V_int)

  delta_b = (bcampo_f-bcampo_i)/real(num_puntos_b, pr);

  ind_b = 0;
  do ind_b = 0, num_puntos_b

    bcampo = bcampo_i + delta_b*real(ind_b, pr);

    l_campo = sqrt(2.0_pr*alpha/bcampo)/a0;
    omega = 0.5_pr*bcampo/(me*c_light*ua_to_T);  !! Frecuencia de oscilacion debida al campo

    vr = (1._pr - exp(-(r0/l_campo)**2))
    vrmax = (1._pr - exp(-(rmax/l_campo)**2));

    call hamiltoniano(nb, vrmax*v1, vr*v2, eta, s, v01, v02, ke, v_int, h, ms, mv);

    auval = 0._pr
    call eigenvalues(ndimh, nev, h, ms, auval, auvec)

    auval = eV*(auval + 2.0_pr*omega); !slope_ll*bcampo;

    exp_int = matmul(transpose(auvec), matmul(mv, auvec))

    write(10,6) bcampo, (auval(i), i = 1, nev)
    write(11,6) bcampo, (exp_int(i,i), i = 1, nev)
    call flush();
  end do

  deallocate(k, t)
  deallocate(x, w, pl)
  deallocate(s, v01, v02, ke)
  deallocate(h, ms)
  deallocate(auval)

  close(10)
  close(11)

6 format(e22.14,1x,1000(1x,e22.14))
! 7 format(3000(1x,e22.14))
end program
