import numpy as _np

def indices(N):

    M = N*(N+1)/2

    idx = _np.linspace(0, N**2-1, N**2)

    ####
    aux = _np.array([])
    # idnm = _np.array([], dtype=int)
    for n in range(N):
        aux = _np.hstack((aux, _np.linspace(n, n, N)))
    ####

    idx = _np.vstack((idx, aux))
    idx = _np.vstack((idx, _np.tile(_np.linspace(0, N-1, N), N)))
    idx = _np.transpose(idx)

    idnm = _np.linspace(0, N*N-1, N*N)
    idnm = _np.reshape(idnm, (N, N))

    idsim = _np.empty((1, 2), dtype=int)
    for n in range(N):
        aux = _np.hstack((_np.tile(n, (n+1, 1)), _np.reshape(_np.array(range(n+1)), (n+1, 1))))
        idsim = _np.vstack((idsim, aux))

    idsim = idsim[1:M+1,:]
    idsim = _np.hstack((_np.reshape(_np.array(range(M)), (M, 1)), idsim))

    idx     = _np.array(idx, dtype=int)
    idnm    = _np.array(idnm, dtype=int)
    idsim   = _np.array(idsim, dtype=int)

    return (idx, idnm, idsim)
