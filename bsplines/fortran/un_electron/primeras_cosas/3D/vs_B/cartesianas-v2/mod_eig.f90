module mod_eig
contains
subroutine eigenvalues(ndim, nev, A, B, auval)
  use precision, pr => dp
  implicit none
  integer, intent(in) :: ndim, nev
  complex(pr), intent(in) :: A(ndim,ndim), B(ndim,ndim)
  complex(pr), intent(out) :: auval(nev)
  !!!
  integer :: j, jmin(1), jmax(1)
  character(1) :: JOBVL, JOBVR
  integer :: LDA, LDB, LDVL, LDVR, LWORK, INFO
  complex(pr) :: ALPHA(ndim), BETA(ndim)
  complex(pr) :: VL(ndim,ndim), VR(ndim,ndim)
  complex(pr) :: WORK(2*ndim)
  real(pr) :: RWORK(8*ndim)
  complex(pr) :: auval_vec(ndim), caux(ndim)
  real(pr) :: vec_aux(ndim), raux(ndim)


  JOBVL = 'N'; JOBVR = 'N';
  LDA = ndim; LDB = ndim; LDVL = ndim; LDVR = ndim;
  LWORK = 2*ndim;

  call zggev(JOBVL, JOBVR, ndim, A, LDA, B, LDB, ALPHA, BETA, VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO)

  ! write(*,*) INFO, real(WORK(1))
  ! stop
  !if(INFO .NE. 0)then
  !  write(*,*) INFO
  !  stop
  !end if

  do j = 1, ndim
    auval_vec(j) = ALPHA(j)/BETA(j);
    !write(13,*) ALPHA(j), BETA(j)
    vec_aux(j) = real(auval_vec(j));
  end do
  !stop

  jmax = maxloc(vec_aux);
  do j = 1, ndim
    jmin = minloc(vec_aux);
    raux(j) = vec_aux(jmin(1));
    caux(j) = auval_vec(jmin(1));
    vec_aux(jmin(1)) = 2.d0*vec_aux(jmax(1)) + 1.d0;
    auval_vec(jmin(1)) = 2.d0*auval_vec(jmax(1)) + DCMPLX(100.d0, 0.d0);
  end do
  vec_aux = raux;
  auval_vec = caux;

  do j = 1, nev
    auval(j) = auval_vec(j);
  end do

end subroutine eigenvalues
end module mod_eig
