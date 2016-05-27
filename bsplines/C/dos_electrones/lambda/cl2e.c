/* Programa en C que resuelve el problema de dos electrones interactuantes
 * en un pozo de potencial usando B-splines.
 * Usa el metodo variacional de Rayleight-Ritz usando como base del
 * espacio los B-splines, como esta no es una base ortonormal el
 * problema de autovalores queda de la siguiente forma:
 *
 *              H |psi> = e S |psi>
 *
 * donde H es la matriz del hamiltoniano y S es la matriz de solapamiento
 * de los B-splines.
 *
 * Usa la funcion gauleg() del Numerical Recipies C para calcular los
 * puntos de evaluacion y los pesos para la cuadratura.
 * La funcion KNOTS_PESOS() esta hecha por mi guiandome de la version de
 * fortran, calcula los knots de los B-splines con una distribucion
 * uniforme solamente, en el caso de querer otra distribucion es solo
 * cuestion de modificar el codigo para usar la distribucion que uno quiera.
 * Para el calculo de los autovalores usa la funcion dsygvx_() de lapack.
 * Las funciones para evaluar los B-splines y las derivadas son bsplvb() y
 * bder() respectivamente, son versiones hechas por mi a partir de las
 * versiones de fortran.
 *
 * EL potencial de confinamiento de cada electron es de la forma:
 *
 * 			V(z) = -V0*exp(-0.5*z^2/sigma^2)/(1 + 0.5*l^2/sigma^2)
 *
 * o sea el pozo es una gaussiana invertida. La profundidad del pozo
 * esta moduloda por la intensidad de campo magnetico a travez del
 * parametro l dado por:
 *
 * 			l = sqrt(2*alpha/B)
 *
 * donde B es la intesidad del campo y alpha es un parametro que toma el
 * valor alpha = 658.4092645439 nm^2/T
 *
 * La interaccion entre los electrones es una interaccion efectiva en la
 * direccion de confinamiento y esta dada por:
 * 
 * 			V_ef(z1,z2) = sqrt(0.5*pi)/l*exp(x^2)*(1 - erf(x))
 * 
 * donde x = abs(z1-z2)/(sqrt(0.5)*l) y erf(x) es la funcion error.
 * 
 * El programa calcula los autavolres del hamiltoniano:
 * 
 * 			H(1,2) = h(1) + h(2) + lambda*V_ef(1,2)
 * 
 * donde nos quedamos con los estados simetricos del sistema y lambda
 * es la carga efectiva.
 *
 * Se calculan los autovalores para distintos valores de lambda
 * en el intervalo [0, 1]
 *
*/

#include <stdio.h> // prinft
#include <stdlib.h> // malloc y free
#include <malloc.h>
#include <math.h> // posibles operaciones matematicas
#include <assert.h>
#include <omp.h> //omp_get_wtime()

// defino algunas constantes para el programa //
#define EPS 3.0e-14  // EPS precision relativa para gauleg //

// defino algunos parametros del programa //

#ifndef RMIN
#define RMIN 0.0 // R minimo donde empieza el intervalo para la integracion //
#endif

#ifndef RMAX
#define RMAX 50.0 // R maximo donde termina el intervalo para la integracion //
#endif

#ifndef L
#define L 25 // numero de intervalos en el que divido al intervalo [Rmin, Rmax] //
#endif

#ifndef KORD
#define KORD 5 // orden de los B-splines, el grado es KORD-1 //
#endif

#ifndef INTG
#define INTG 100 // grado de integracion por cuadratura //
#endif

#ifndef NEV
#define NEV 10 // numero de autovalores que vamos a calcular //
#endif

#ifndef ME
#define ME 1.0 // masa de la particula //
#endif

#ifndef V_POZO
#define V_POZO 1.0 // profundidad del pozo //
#endif

#ifndef SIGMA
#define SIGMA 20.0 // ancho del pozo en nm //
#endif

#ifndef B_CAMPO
#define B_CAMPO 50.0 // campo magnetico inicial //
#endif

#ifndef LAMBDA_I
#define LAMBDA_I 0.0 // carga efectiva //
#endif

#ifndef LAMBDA_F
#define LAMBDA_F 1.0 // campo magnetico final //
#endif

#ifndef NUM_LAMBDA
#define NUM_LAMBDA 20 // numero de puntos para calcular //
#endif

#define a0 0.0529177210
#define eV 27.21138564
#define alpha 658.4092645439
#define pi 3.1415926535897

// escribo las funciones del programa //
int dsygvx_(int *itype, char *jobz, char *range, char *	uplo,
	int *n, double *a, int *lda, double *b, int *ldb,
	double *vl, double *vu, int *il, int *iu, double *abstol,
	int *m, double *w, double *z__,	int *ldz, double *work,
	int *lwork, int *iwork, int *ifail, int *info);


int idx(unsigned int y, unsigned int x, unsigned int numcolumns){
	return y*numcolumns + x;
}

int idx2(unsigned int y1, unsigned int x1, unsigned int y2, unsigned int x2, unsigned int numcolumns) {
	return idx(y1, x1, numcolumns)*numcolumns*numcolumns + idx(x2, y2, numcolumns);
}

int cleari(unsigned int N, int * __restrict__ vec) {

	for(unsigned int i = 0; i<N; ++i)
		vec[i] = 0;

	return 0;
}

int cleard(unsigned int N, double * __restrict__ vec) {

	for(unsigned int i = 0; i<N; ++i)
		vec[i] = 0.f;

	return 0;
}


void gauleg(double x1, double x2, double x[], double w[], int n) {
/* Given the lower and upper limits of integration x1 and x2,
 * and given n, this routine returns arrays x[1..n] and w[1..n]
 * of length n, containing the abscissas and weights of the Gauss-
 * Legendre n-point quadrature formula.
*/
	int m, j, i;
	double z1, z, xm, xl, pp, p3, p2, p1;

	m = (n+1)/2;
	xm = 0.5*(x2+x1);
	xl = 0.5*(x2-x1);

	for (i = 1; i<=m; i++) {
		z = cos(3.141592654*(i-0.25)/(n+0.5));
		do {
			p1 = 1.0;
			p2 = 0.0;
			for (j = 1;j<=n;j++) {
				p3 = p2;
				p2 = p1;
				p1 = ((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
			}
			pp = n*(z*p1-p2)/(z*z-1.0);
			z1 = z;
			z = z1-p1/pp;
		} while (fabs(z-z1) > EPS);
		x[i] = xm-xl*z;
		x[n+1-i] = xm+xl*z;
		w[i] = 2.0*xl/((1.0-z*z)*pp*pp);
		w[n+1-i] = w[i];
	}
}

int KNOTS_PESOS(unsigned int nk, int * __restrict__ k,
		double * __restrict__ t, double * __restrict__ x,
		double * __restrict__ w){

	double dr, ri, rf;
	double vx[INTG+1], vw[INTG+1];
  double Rmax, Rmin;

	Rmax = RMAX/a0; Rmin = RMIN/a0;

	dr = (Rmax-Rmin)/L;

	for(unsigned int i = 0; i<L; ++i) {
		ri = Rmin + i*dr;
		rf = ri + dr;
		gauleg(ri, rf, vx, vw, INTG);

		for(unsigned int j = 0; j<INTG; j += 1 ) {
			x[idx(i, j, INTG)]  = vx[j+1] ;//  x[idx(i, j+1, INTG)] = vx[j+2];

			w[idx(i, j, INTG)]  = vw[j+1] ;//  w[idx(i, j+1, INTG)] = vw[j+2];
		}
	}

	/// en esta parte controlar bien el tema de los indices de los vectores //
	t[0] = Rmin; k[0] = 0;

	for(unsigned int i = 1; i<KORD; ++i) {
		t[i] = t[i-1];
		k[i] = k[i-1];
	}

	for(unsigned int i = KORD; i<KORD+L; ++i) {
		t[i] = t[i-1] + dr;
		k[i] = k[i-1] + 1;
	}

	for(unsigned int i = KORD+L; i<nk; ++i) {
		t[i] = t[i-1];
		k[i] = k[i-1];
	}

	return 0;
}

int bsplvb(double * __restrict__ t, unsigned int jhigh, unsigned int index,
	   double rr, int left, double * __restrict__ biatx ) {

	unsigned int j, jp1, JMAX;
	double saved, term;
	double * deltar, * deltal;

	if(1 != index) printf("index no es igual a 1");

	JMAX = 100;

	deltar = (double *) calloc(JMAX, sizeof(double));
	deltal = (double *) calloc(JMAX, sizeof(double));

	biatx[0] = 1.0;

	for(j=0; j<jhigh-1; ++j) {
		jp1 = j+1;

		deltar[j] = t[left+j+1]-rr;
		deltal[j] = rr-t[left-j];

		saved = 0.0;
		for(unsigned int i = 0; i<j+1; ++i) {
			term = biatx[i]/(deltar[i]+deltal[jp1-i-1]);
			biatx[i] = saved + deltar[i]*term;
			saved = deltal[jp1-i-1]*term;
		}

		biatx[jp1] = saved;
	}

	free(deltar);
	free(deltal);

	return 0;
}

double bder(double rr, double * __restrict__ t, unsigned int korder,
	 unsigned int np, unsigned int indexm, unsigned int left,
	 double * __restrict__ Sp, double dm ) {

	unsigned int i;

	if(t[0]<rr && rr<t[np-1]) {

		if(abs(rr-t[np-1])<1.e-10) {
			if(indexm==np-korder) {
				dm = (korder-1)/(t[np-1]-t[np-1-korder]);
			}
			else {
				dm = -(korder-1)/(t[np-1]-t[np-1-korder]);
			}
		}

		bsplvb(t, korder-1, 1, rr, left, Sp);

		if(indexm-left+korder>=1 || indexm-left+korder<=korder) {
			i = indexm-left+korder;
			if(1==i) {
				dm = (korder-1)*(-Sp[i-1]/(t[indexm+korder]-t[indexm+1]));
			}
			else if(korder==i) {
				dm = (korder-1)*(Sp[i-1-1]/(t[indexm+korder-1]-t[indexm]));
			}
			else {
				dm = (korder-1)*(Sp[i-1-1]/(t[indexm+korder-1]-t[indexm])
				    - Sp[i-1]/(t[indexm+korder]-t[indexm+1]));
			}
		}

	}

	return dm;
}

void calculo_matrices(unsigned int nk, unsigned int nb,
		      int * __restrict__ k, double * __restrict__ t,
		      double * __restrict__ x, double * __restrict__ w,
		      double * __restrict__ s, double * __restrict__ v0,
		      double * __restrict__ ke) {

	unsigned int im, in;
	double rr, sigma;
	double * Sp;

  sigma = SIGMA/a0;

	Sp = (double *) calloc(KORD, sizeof(double));
	// ojo con los limites en los for's //
	for(unsigned int i = KORD-1; i<KORD+L-1; ++i) {
		for(unsigned int j = 0; j<INTG; ++j) {
			rr = x[idx(k[i], j, INTG)];

			bsplvb(t, KORD, 1, rr, i, Sp );

			for(unsigned int m = 0; m<KORD; ++m) {
				im = i-KORD+m;
				if(im<nb) {
					for(unsigned int n = 0; n<KORD; ++n) {
						in = i-KORD+n;
						if(in<nb) {

							s[idx(im, in, nb)] += Sp[m]*Sp[n]*w[idx(k[i], j, INTG)];

							v0[idx(im, in, nb)] += Sp[m]*Sp[n]*w[idx(k[i] , j, INTG)]*exp(-0.5*rr*rr/(sigma*sigma));

						}
					}
				}
			}
		}
	}

	for(unsigned int i = KORD-1; i<KORD+L-1; ++i) {
		// ojo con los indices en esta parte //
		for(unsigned int m = i-KORD+1; m<=i; ++m) {
			if(m>0 && m<nb) {
				for(unsigned int n = m; n<=i; ++n) {
					if(n<nb) {
						for(unsigned int j = 0; j<INTG; ++j) {

							double bm = 0, bn = 0;

							rr = x[idx(k[i], j, INTG)];

							bm = bder(rr, t, KORD, nk, m, i, Sp, bm);
							bn = bder(rr, t, KORD, nk, n, i, Sp, bn);

							ke[idx(m-1, n-1, nb)] += 0.5*w[idx(k[i], j, INTG)]*bm*bn/ME;

						}
					}
				}
			}
		}
	}

	for(unsigned int i = 0; i<nb; ++i) {
		for(unsigned int j = i+1; j<nb; ++j) {
			s[idx(j, i, nb)] = s[idx(i, j, nb)];
			v0[idx(j, i, nb)] = v0[idx(i, j, nb)];
			ke[idx(j, i, nb)] = ke[idx(i, j, nb)];
		}
	}

	free(Sp);
}

void calculo_interaccion(unsigned int nb,
		 int * __restrict__ k, double * __restrict__ t, double * __restrict__  x,
		 double * __restrict__ w, double * __restrict__ v_ef ) {

	unsigned int im, in, imp;
	double rr1, rr2, ww1, ww2;
	double r12, _r12, cte, l_campo;
	double * Sp;
	double * f;

	l_campo = sqrt(2.0*alpha/B_CAMPO)/a0;
	cte = sqrt(0.5*pi)/l_campo;

	Sp = (double *) calloc(KORD, sizeof(double));
	f = (double *) calloc(nb*nb, sizeof(double));
	// ojo con los limites en los for's //
	for(unsigned int i1 = KORD-1; i1<KORD+L-1; ++i1) {
		for(unsigned int j1 = 0; j1<INTG; ++j1) {
			rr1 = x[idx(k[i1], j1, INTG)]; ww1 = w[idx(k[i1], j1, INTG)];

			cleard(nb*nb, f);
			for(unsigned int i2 = KORD-1; i2<KORD+L-1; ++i2) {
				for(unsigned int j2 = 0; j2<INTG; ++j2) {
					rr2 = x[idx(k[i2], j2, INTG)]; ww2 = w[idx(k[i2], j2, INTG)];

					r12 = sqrt(0.5)*fabs(rr1 - rr2)/l_campo;
					_r12 = 1.0/r12;

					bsplvb(t, KORD, 1, rr2, i2, Sp);

					for(unsigned int m = 0; m<KORD; ++m) {
            im = i2-KORD+m;
            if(im<nb) {
              for(unsigned int n = 0; n<KORD; ++n) {
                in = i2-KORD+n;
                if(in<nb) {

									if( r12 < 5.0 ) {
										f[idx(im, in, nb)] += Sp[m]*Sp[n]*ww2*cte*exp(r12*r12)*(1.0 - erf(r12));
									}else{
										f[idx(im, in, nb)] += Sp[m]*Sp[n]*ww2/l_campo*(sqrt(0.5)*_r12 - pow(sqrt(0.5)*_r12, 3));
									}
                }
              }
            }
          }
				}
			}

			cleard(KORD ,Sp);
			bsplvb(t, KORD, 1, rr1, i1, Sp);

			for(unsigned int m = 0; m<KORD; ++m) {
				im = i1-KORD+m;
				if(im<nb) {
					for(unsigned int mp = 0; mp<KORD; ++mp) {
						imp = i1-KORD+mp;
						if(imp<nb) {
							for(unsigned int n = 0; n<nb; ++n) {
								for(unsigned int np = 0; np<nb; ++np) {
									v_ef[idx2(im, n, imp, np, nb)] += Sp[m]*Sp[mp]*ww1*f[idx(n, np, nb)];
								}
							}
						}
					}
				}
			}
		}
	}

	free(Sp);
	free(f);

}


void eigenvalues(int n, int m, double * __restrict__ a,
		 double * __restrict__ b, double * __restrict__ w,
		 double * __restrict__ z) {

	int itype, lda, ldb, ldz;
	int il, iu, lwork, info;
	char jobz, range, uplo;
	double abstol, vl, vu;
	double * work;
	int * iwork, * ifail;

	// le doy los valores correspondientes a las distintas variables //
	itype = 1; lda = n; ldb = n; ldz = n;
	vl = 0.0; vu = 0.0;
	jobz = 'V'; range = 'I'; uplo = 'U';
	il = 1; iu = m; lwork = 8*n;

	// le doy memoria a las matrices que neceista dsygvx //
	work = (double *) malloc(lwork*sizeof(double));
	iwork = (int *) malloc(5*n*sizeof(int));
	ifail = (int *) malloc(n*sizeof(int));

	dsygvx_( &itype, &jobz, &range, &uplo, &n, a, &lda, b, &ldb,
		 &vl, &vu, &il, &iu, &abstol, &m, w, z, &ldz, work,
		 &lwork, iwork, ifail, &info);

	if(0!=info) {
		printf("Hay un error en dsyvx_, info = %i\n", info);
		assert(0==info);
	}

	free(work);
	free(iwork);
	free(ifail);
}

void hamiltoniano_autovalores(unsigned int nb, double * __restrict__ s,
			      double * __restrict__ v0, double * __restrict__ ke, double * __restrict__ v_ef,
			      FILE * archivo) {

	unsigned int n_dim, ind1, ind2;
	double * h, * auval, * auvec; //, * s_copy;
	double * h_sim, * s_sim;
	double lambda, delta_lambda, l_campo;
  double b, v_pozo;
	double raiz;

	raiz = sqrt(0.5);

	n_dim = nb*(nb+1)/2; // dimension de H simetrico //

  b = 0.5*a0*a0/(SIGMA*SIGMA);
	v_pozo = V_POZO/eV;

	// doy memoria a las matrices y vectores //
	h = (double *) calloc( nb*nb, sizeof(double));
//	s_copy = (double *) calloc( nb*nb, sizeof(double));

	h_sim = (double *) calloc( n_dim*n_dim, sizeof(double));
	s_sim = (double *) calloc( n_dim*n_dim, sizeof(double));
	auval = (double *) calloc( NEV, sizeof(double));
	auvec = (double *) calloc( NEV*n_dim, sizeof(double));

	delta_lambda = (LAMBDA_F - LAMBDA_I)/(NUM_LAMBDA);

	// abro el archivo para guardar los datos //
	fprintf(archivo, "# Autovalores calculados\n");
	fprintf(archivo, "# Lambda  auval[0]   auval[1] ....\n");


	l_campo = sqrt(2.0*alpha/B_CAMPO)/a0;

	// Hamiltoniano de un electron //
	for(unsigned int m = 0; m<nb; ++m) {
		for(unsigned int n = 0; n<nb; ++n){
			h[idx(n, m, nb)] = ke[idx(n, m, nb)] - v_pozo/(1.0 + l_campo*l_campo*b)*v0[idx(n, m, nb)];
//			s_copy[idx(n, m, nb)] = s[idx(n, m, nb)];
		}
	}

	for(unsigned int ind_lambda = 0; ind_lambda<NUM_LAMBDA; ++ind_lambda) {

		lambda = LAMBDA_I + delta_lambda*ind_lambda;

		// Hamiltoniano de dos electrones simetrico //
		ind1 = 0;
		for(unsigned int n = 0; n<nb; ++n) {
			for(unsigned int m = 0; m<=n; ++m) {

				ind2 = 0;
				for(unsigned int np = 0; np<nb; ++np) {
					for(unsigned int mp = 0; mp<=np; ++mp) {

						if(m==n && mp==np){

							h_sim[idx(ind1, ind2, n_dim)] = 2.0*s[idx(n,np,nb)]*h[idx(n,np,nb)] + lambda*v_ef[idx2(n,n,np,np,nb)];
							s_sim[idx(ind1, ind2, n_dim)] = s[idx(n,np,nb)]*s[idx(n,np,nb)];

						} else if(m!=n && mp==np ) {

							h_sim[idx(ind1, ind2, n_dim)] = raiz*( 2.0*s[idx(m,np,nb)]*h[idx(n,np,nb)] + 2.0*s[idx(n,np,nb)]*h[idx(m,np,nb)]
							                              + lambda*v_ef[idx2(n,m,np,np,nb)] + lambda*v_ef[idx2(m,n,np,np,nb)] );
							s_sim[idx(ind1, ind2, n_dim)] = 2.0*raiz*s[idx(n,np,nb)]*s[idx(m,np,nb)];

						} else if(m==n && mp!=np) {

							h_sim[idx(ind1, ind2, n_dim)] = raiz*( 2.0*s[idx(n,mp,nb)]*h[idx(n,np,nb)] + 2.0*s[idx(n,np,nb)]*h[idx(n,mp,nb)]
							                              + lambda*v_ef[idx2(n,n,np,mp,nb)] + lambda*v_ef[idx2(n,n,mp,np,nb)] );
							s_sim[idx(ind1, ind2, n_dim)] = 2.0*raiz*s[idx(n,np,nb)]*s[idx(n,mp,nb)];

						} else if(m!=n && mp!=np) {

							h_sim[idx(ind1, ind2, n_dim)] = s[idx(n,np,nb)]*h[idx(m,mp,nb)] + s[idx(n,mp,nb)]*h[idx(m,np,nb)]
							                              + s[idx(m,mp,nb)]*h[idx(n,np,nb)] + s[idx(m,np,nb)]*h[idx(n,mp,nb)]
							                              + 0.5*lambda*(v_ef[idx2(n,m,np,mp,nb)] + v_ef[idx2(n,m,mp,np,nb)]
                                                        + v_ef[idx2(m,n,np,mp,nb)] + v_ef[idx2(m,n,mp,np,nb)]);
							s_sim[idx(ind1, ind2, n_dim)] = s[idx(n,np,nb)]*s[idx(m,mp,nb)] + s[idx(n,mp,nb)]*s[idx(m,np,nb)];

						}


						ind2 = ind2 + 1;
					}
				}
				ind1 = ind1 + 1;
			}
		}

		eigenvalues( n_dim, NEV, h_sim, s_sim, auval, auvec );

		fprintf(archivo, "%.5f   ", lambda);
		for(unsigned int i = 0; i < NEV; ++i) {
			fprintf(archivo, "%.15f   ", eV*auval[i]);
		}
		fprintf(archivo, "\n");

	}

	free(h);
	free(h_sim);
	free(s_sim);
	free(auval);
	free(auvec);
}

int main(void) {

	// defino algunas variables que voy a usar //
	unsigned int nk, nb;
	int *k;
	double *t, *x, *w;
	double *s, *v0, *ke;
	double *v_ef;
	double t_in, t_fin;
	FILE * archivo;
	char name [150];
	int int_name;

	nk = L + 2*KORD - 1; // numero de knots //
	nb = nk - KORD - 2; // tamaño de la base //

	// controlo algunos parametros //
	assert(INTG>KORD);
	assert(NEV>0);

	int_name = sprintf(name, "./resultados/2e-E_vs_lambda-B%2.0fT-zmax%3.0fnm.dat", B_CAMPO, RMAX );
	int_name = int_name + 1;

	archivo = fopen(name, "w");
	// imprimo los parametros //
	fprintf(archivo, "# Rmin = %.12f y RMAX = %.12f\n", RMIN, RMAX);
	fprintf(archivo, "# Numero de intervalos l = %i\n", L);
	fprintf(archivo, "# Orden los B-splines KORD = %i\n", KORD);
	fprintf(archivo, "# Grado de integracion de la cuadratura INTG = %i\n", INTG);
	fprintf(archivo, "# Numero de autovalores NEV=%i\n", NEV);
	fprintf(archivo, "# Numero de knots nk=%i\n", nk);
	fprintf(archivo, "# Tamaño de la base nb=%i\n", nb);
	fprintf(archivo, "# Masa de la particula me = %.12f\n", ME);
  fprintf(archivo, "# Profundidad del pozo V_pozo = %.12f\n", V_POZO);
	fprintf(archivo, "# Ancho del pozo sigma = %.12f\n", SIGMA);
  fprintf(archivo, "# Intensidad de campo magnetico B_CAMPO = %.12f\n", B_CAMPO);
	fprintf(archivo, "# LAMBDA_I = %.12f, LAMBDA_F = %.12f\n", LAMBDA_I, LAMBDA_F);
	fprintf(archivo, "# Numero de puntos de B NUM_LAMBDA = %i\n", NUM_LAMBDA);

	// doy memoria a las matrices que voy a necesitar //
	k = (int *) malloc( nk*sizeof(int));
	t = (double *) malloc( nk*sizeof(double));
	x = (double *) malloc( L*INTG*sizeof(double));
	w = (double *) malloc( L*INTG*sizeof(double));
	s = (double *) malloc( nb*nb*sizeof(double));
	v0 = (double *) malloc( nb*nb*sizeof(double));
	ke = (double *) malloc( nb*nb*sizeof(double));
	v_ef = (double *) malloc( nb*nb*nb*nb*sizeof(double));

	cleari(nk, k);
	cleard(nk, t);
	cleard(L*INTG, x);
	cleard(L*INTG, w);
	cleard(nb*nb, s);
	cleard(nb*nb, v0);
	cleard(nb*nb, ke);
	cleard(nb*nb*nb*nb, v_ef);

	t_in = omp_get_wtime();
	// primero calculos los knost y los pesos para hacer la cuadratura //
	KNOTS_PESOS(nk, k, t, x, w);

	// calculo las matrices que necesito para resolver el problema //
	calculo_matrices(nk, nb, k, t, x, w, s, v0, ke);

	calculo_interaccion(nb, k, t, x, w, v_ef);

	// armo el hamiltoniano y calculo los autovalores //
	hamiltoniano_autovalores(nb, s, v0, ke, v_ef, archivo);

	t_fin = omp_get_wtime();
	printf("%i     %i     %i     %.12f\n", L, nk, nb, t_fin-t_in);

	// libero la memoria //
	fclose(archivo);
	free(v_ef);
	free(ke);
	free(v0);
	free(s);
	free(w);
	free(x);
	free(t);
	free(k);

	return 0;
}
