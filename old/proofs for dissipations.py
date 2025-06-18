import numpy as np
import sympy as sym

nu_2, nu_4 = sym.symbols('nu_2, nu_4')

def diff_xi(obj, j, i):
    return obj(j, i+1/2) - obj(j, i-1/2)

def diff_eta(obj, j, i):
    return obj(j+1/2, i) - obj(j-1/2, i)

def diff2_xi(obj, j, i):
    return obj(j, i+1) - 2*obj(j, i) + obj(j, i-1)

def diff2_eta(obj, j, i):
    return obj(j+1, i) - 2*obj(j, i) + obj(j-1, i)  

def diff3_eta(obj, j, i):
    return diff2_eta(obj, j+1/2, i) - diff2_eta(obj, j-1/2, i)
    
def diff3_xi(obj, j, i):
        # "at i+1/2"
        # TODO: cut out the function call for better speed
        return diff2_xi(obj, j, i+1/2) - diff2_xi(obj, j, i-1/2)


def t1(jc,ic):
    """input - off-index values"""
    if np.isclose(nu_2, 0.0):
        return 0.0
    else:
        return switch2_xi(jc,ic) * l(jc,ic) * lamb(jc,ic,og=(i-ic) ) * diff_xi(q,jc,ic)
def t2(jc,ic):
    if np.isclose(nu_2, 0):
        return 0.0
    else:
        return switch2_eta(jc,ic) * l(jc,ic) * lamb(jc,ic, og=(j-jc)) * diff_eta(q,jc,ic)
def t3(jc, ic):
    s4 = switch4_xi(jc,ic)
    if np.isclose(s4, 0, rtol=1e-10):
        return 0.0
    # elif np.isclose(jc, 1.5) or np.isclose(jc, np.shape(cell_data)[0]-1.5):
    #     return 0.0 # lambda SHOULD be zero here.
    else:
        return s4 * l(jc,ic) * lamb(jc, ic, og=(i-ic)) * diff3_xi(q, jc,ic)
    
def t4(jc, ic):
    s4 = switch4_eta(jc,ic)
    if np.isclose(s4, 0, rtol=1e-10):
        return 0.0
    # elif np.isclose(j, 1.5) or np.isclose(jc, np.shape(cell_data)[0]-1.5):
    #     return 0.0 # lambda SHOULD be zero here.
    else:
        return s4 * l(jc,ic) * lamb(jc, ic, og=(j-jc)) * diff3_eta(q, jc, ic)

try:
    D[:] = diff_xi(t1, j,i) + diff_eta(t2, j,i) - diff_xi(t3, j, i) - diff_eta(t4, j, i)