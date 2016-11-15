import numpy as _np

def V_interaccion(z_nodos, w_pesos, delta, bsz):

    bsp = _np.kron(bsz, bsz)
    bsp = bsp[0:_np.size(z_nodos)**2:_np.size(z_nodos)+1,:]

    z12 = _np.tile(z_nodos, (_np.size(z_nodos), 1)) - _np.transpose(_np.tile(z_nodos, (_np.size(z_nodos), 1)))
    z12 = _np.sqrt(delta**2 + z12**2)

    V_int = _np.tile(w_pesos, (_np.size(w_pesos), 1))*_np.transpose(_np.tile(w_pesos, (_np.size(w_pesos), 1)))/z12

    U_int = _np.dot(_np.transpose(bsp), _np.dot(V_int, bsp))

    return U_int
