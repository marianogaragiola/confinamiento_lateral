import numpy as _np
from index import indices

def H_antisimetrico(H, S, V_int, eta):

    N = _np.size(H[:,0])
    N_dim = N*(N+1)/2

    idx, idnm, idsim = indices(N)

    H_sim = _np.zeros((N_dim, N_dim))
    S_sim = _np.zeros((N_dim, N_dim))

    for ind1 in range(N_dim):

        n1 = idsim[ind1, 1]
        m1 = idsim[ind1, 2]

        for ind2 in range(ind1+1):

            n2 = idsim[ind2, 1]
            m2 = idsim[ind2, 2]

            if m1==n1 and m2==n2:

                i1 = idnm[n1,n1]
                j1 = idnm[n2,n2]
                H_sim[ind1,ind2] = 2*H[n1,n2]*S[n1,n2] + eta*V_int[i1,j1]
                S_sim[ind1,ind2] = S[n1,n2]*S[n1,n2]

            elif m1==n1 and m2!=n2:

                i1 = idnm[n1,n1];
                j1 = idnm[n2,m2]; j2 = idnm[m2,n2];
                H_sim[ind1,ind2] = _np.sqrt(0.5)*(2*(H[n1,n2]*S[n1,m2] + H[n1,m2]*S[n1,n2])
                + eta*(V_int[i1,j1] + V_int[i1,j2]))

                S_sim[ind1,ind2] = _np.sqrt(0.5)*2*S[n1,n2]*S[n1,m2]

            elif m1!=n1 and m2==n2:

                i1 = idnm[n1,m1]; i2 = idnm[m1,n1];
                j1 = idnm[n2,n2];
                H_sim[ind1,ind2] = _np.sqrt(0.5)*(2*(H[m1,n2]*S[n1,n2] + H[n1,n2]*S[m1,n2])
                + eta*(V_int[i1,j1] + V_int[i2,j1]));

                S_sim[ind1,ind2] = _np.sqrt(0.5)*2*S[n1,n2]*S[m1,n2];

            else:
                i1 = idnm[n1,m1]; i2 = idnm[m1,n1];
                j1 = idnm[n2,m2]; j2 = idnm[m2,n2];
                H_sim[ind1,ind2] = H[m1,m2]*S[n1,n2] + H[n1,n2]*S[m1,m2]
                + H[m1,n2]*S[n1,m2] + H[n1,m2]*S[m1,n2]
                + 0.5*eta*(V_int[i1,j1] + V_int[i1,j2] + V_int[i2,j1] + V_int[i2,j2]);

                S_sim[ind1,ind2] = S[n1,n2]*S[m1,m2] + S[n1,m2]*S[m1,n2];
        ######### Termina el if

        H_sim[ind2,ind1] = H_sim[ind1,ind2]
        S_sim[ind2,ind1] = S_sim[ind1,ind2]
        ##### Termina el for ind2 in range(ind1+1)

    ##### Termina el for ind1 in range(N_dim)

    return (H_sim, S_sim)
