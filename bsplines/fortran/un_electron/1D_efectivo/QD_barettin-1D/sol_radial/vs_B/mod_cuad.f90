module cuadratura_knots

contains
  SUBROUTINE KNOTS_PESOS(kord, tip, beta, a, b, l, intg, t, k, x, w, pl)
  !lo único que hace esta subrutina principal es dividir los casos segun se quiera knots uniformes, exponencial, mixto
  !cada una de las subrutinas que llama ademas calcula los x y los w, abscisas y pesos para la cuadratura gaussiana, y los polinomios
  !  INPUT:
  ! kord: orden de los b-splines
  ! tip; character(1): 'u' ; 'e' ; 'm': dist uniforme, exp o mixta de knots
  ! beta : param dist. e (no usado si tip='u')
  ! a
  ! b  ; a<b todo es calculado en el intervalo [a,b]
  ! l = # total de sub_intervalos en [a,b]
  ! intg = # de puntos que usamos en la integral x cuadratura en c/subintervalo
  !  OUTPUT:
  !  t(l+2*kord-1)  ! los (l+2*kord-1) knots (contando multiplicidades en a y b)
  !  k(l+2*kord-1)  ! da el # de intervalo j para el i-esimo nodo 1<=k<=l
  !  x(l,intg),w(l,intg) : posiciones y pesos de los intg puntos en c/intervalo
  ! pl(l,intg),w(l,intg) : polinomios de Legendre en cada punto
  ! NOTA II : los p_l los calcula gratis la subr gauleg de int. x cuad. del NR
  ! y a veces hacen falta=> version modif. de gauleg que los da como output
    use precision, pr => dp
    implicit none
    integer :: kord, l, intg, tip
    integer, dimension(l+2*kord-1) :: k
    real(pr) :: beta, a, b
    real(pr), dimension(l+2*kord-1) :: t
    real(pr), dimension(l,intg) :: x, w, pl
    !!!!!!

    if(tip==1)then
      call dist_unif(kord, a, b, l, intg, t, k, x, w, pl)
    elseif(tip==2)then
      call dist_exp(kord, beta, a, b, l, intg, t, k, x, w, pl)
    else
      write(*,*)'error 1 en KNOTS_PESOS :',tip,' no corresponde a una distribucion'
      stop
    endif

    return
  end SUBROUTINE KNOTS_PESOS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE DIST_UNIF(kord,a,b,l,intg,t,k,x,w,pl)
    use precision, pr => dp
    implicit none
    integer :: kord, l, intg
    real(pr) :: a, b
    integer, dimension(l+2*kord-1) :: k
    real(pr), dimension(l+2*kord-1) :: t
    real(pr), dimension(l,intg) :: x, w, pl
    !!!!!!
    integer :: i, j, nk
    real(pr) :: ri, rf, dr
    real(pr), dimension(intg) :: vx, vw, vpl

    nk=l+2*kord-1;   ! # de knots

    dr=(b-a)/real(l, pr);

    ! calcula los puntos y pesos para la cuadratura
    x = 0._pr; w = 0._pr;

    do i = 1, l

       ri = a+ real(i-1, pr)*dr;
       rf = ri + dr;

       call gauleg_pl(ri, rf, vx, vw, vpl, intg)

       do j = 1, intg
         x(i,j) = vx(j);
         w(i,j) = vw(j);
         pl(i,j) = vpl(j);
       end do

    end do

    t(1) = a;
    k(1) = 1;

    do i = 2, kord
       t(i) = t(1);
       k(i) = 1;
    end do

    do i = kord+1, kord+l
       t(i) = t(i-1)+dr;
       k(i) = k(i-1)+1;
    end do

    do i = kord+l+1, nk
       t(i) = t(i-1);
       k(i) = k(i-1);
    end do

    return

  end SUBROUTINE DIST_UNIF

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE DIST_EXP(kord, beta, a, b, l, intg, t, k, x, w, pl)
    use precision, pr => dp
    implicit none
    integer :: kord,l,intg
    real(pr) :: beta,a,b
    integer, dimension(l+2*kord-1) :: k
    real(pr), dimension(l+2*kord-1) :: t
    real(pr), dimension(l,intg) :: x,w,pl
    !!!!!!
    integer :: i, nk
    real(pr) :: ri, dr, length
    real(pr), dimension(intg) :: vx,vw,vpl

    nk=l+2*kord-1   ! # de knots

    dr = (b-a)/real(l, pr);

    ! calcula los puntos y pesos para la cuadratura
    x = 0._pr; w = 0._pr;

    length = (b-a)/(sign(1._pr,b)*exp(beta*abs(b))-sign(1._pr,a)*exp(beta*abs(a)));

    t(1) = a
    k(1) = 1

    do i = 2, kord
      t(i) = t(1)
      k(i) = 1
    end do

    do i = kord+1 , kord+l
      ri = (i-kord)*dr + a;
      t(i) = length*sign(1._pr, ri)*exp(beta*abs(ri))
      k(i) = k(i-1)+1
    end do

    do i = kord+l+1, nk
      t(i)=t(i-1)
      k(i)=k(i-1)
    end do

    call gauleg_pl(-1._pr, 1._pr, vx, vw, vpl, intg);

    do i = 1, l
      x(i,:) = 0.5_pr*(t(i+kord)-t(i+kord-1))*vx(:) + 0.5_pr*(t(i+kord)+t(i+kord-1))
      w(i,:) = 0.5_pr*(t(i+kord)-t(i+kord-1))*vw(:);
    end do

    pl = 0._pr;

    return
  end SUBROUTINE DIST_EXP
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! subrutina de NR modificada
  ! calcula coordenadas y pesos GL y los valores de P_l(x)

        SUBROUTINE gauleg_pl(x1,x2,x,w,pl,n)
          use precision, pr => dp
          INTEGER n
          DOUBLE PRECISION x1,x2,x(n),w(n),pl(n)
          DOUBLE PRECISION EPS
          PARAMETER (EPS=1.d-14)
          INTEGER i,j,m
          DOUBLE PRECISION p1,p2,p3,pp,xl,xm,z,z1
          m=(n+1)/2
          xm=0.5_pr*(x2+x1)
          xl=0.5_pr*(x2-x1)
          do 12 i=1,m
            z=dcos(3.1415926535897932385_pr*(dfloat(i)-.25_pr)/(dfloat(n)+.5_pr))
  1         continue
              p1=1._pr
              p2=0._pr
              do 11 j=1,n
                p3=p2
                p2=p1
                p1=((2._pr*j-1._pr)*z*p2-(j-1._pr)*p3)/dfloat(j)
  11          continue
              pl(i)=p1
              pp=dfloat(n)*(z*p1-p2)/(z*z-1._pr)
              z1=z
              z=z1-p1/pp
            if(dabs(z-z1).gt.EPS)goto 1
            x(i)=xm-xl*z
            x(n+1-i)=xm+xl*z
            w(i)=2._pr*xl/((1._pr-z*z)*pp*pp)
            w(n+1-i)=w(i)
  12      continue
          return
        END

end module cuadratura_knots
