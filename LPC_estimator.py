import numpy as np
import audiolazy.lazy_lpc as lazy

def lpc_coefficients(x, p):
    """
    Computes the p+1 LPC coeficients of the given sample x
    :param x: signal frame whose LPC coefficients are computed
    :param p: <int> Order of the LPC
    :return: The LPC coefficients
    """
    #return levinson_durbin(x, p)
    return levinson(x,p)

def levinson(x,p):
    m = len(x)
    hamm = np.hamming(m)
    x = [x[i] * hamm[i] for i in range(m)]
    ldfit = lazy.lpc.kautocor(x, p)
    coefs = ldfit.numerator
    return coefs

def levinson_durbin(t, p):
    """
    Computes the lpc coeficients using the levinson-durbin algorithm
    :param t: <[]> vector containing the sample signal
    :param p: <int> Order of the
    :return: <[]> a list conaining p+1 coeficients
    """
    m = len(t) # Number of samples in the input vector

    # 1 - Initialize matrix (1x1):
    # Apply a Hamming window:
    hamm = np.hamming(m)
    t = [t[i]*hamm[i] for i in range(m)]

    # 1.1 - Compute the  correlation matrix values
    r = [r_l(t, i) for i in range(m+1)]  # Vector with the correlation values

    A_k = [1]
    J_k = r[0]

    r_k = []

    # 2 - Iterate
    for k in range(1, p+1):
        r_k.insert(0, r[k])
        # parcor computation (it includes the lambda):
        parcor_k = -1.0 * np.dot(r_k, A_k) / J_k
        # Update the error
        J_k *= 1 - parcor_k**2
        # Updated the coeffs
        A_new = []
        for i, a in enumerate(A_k):
            A_new.append(a + parcor_k * A_k[-1-i])
        # Add the new coefficient
        A_new.append(parcor_k)
        A_k = A_new

    return A_k


def r_l(x, l):
    """
    Computes the l correlation coeficient: R[l]
    :param x: <[]> vector containing the samples
    :param l: <int> index of the correlation coeficient
    :return: <float> the value of the correlation coeficient
    """
    r = 0.0
    for i in range(len(x)-l):
        r += x[i] * x[i+l]
    return r

def autocorr(x, t=1):
    return np.corrcoef(np.array([x[0:-t], x[t:]]))

def autocorr2(x,t=1):
    res = np.correlate(x, x, mode='full')
    return res[res.size/2:]