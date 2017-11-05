import numpy as np

def Az_filter(x, alphas,x_previous=None):
    """
    Computes the prediction errors for the given set of samples using the previous samples
    :param x: <[int]> list containing the samples
    :param alphas: <[float]> list containing the lpc coefficients (without the 1)
    :param x_previous: <[int]> list containing the previous state
    :return:
        e: <[float]> Prediction errors
        x_previous: <[int]> list containing the current sate
    """
    p = len(alphas)
    if x_previous is None:
        x_previous = [0 for i in range(p)]

    e = []
    for x_n in x:
        e.append(x_n - np.dot(x_previous, alphas))
        x_previous = [x_n] + x_previous[:-1]
    
    return e, x_previous


def Sz_filter(e, alphas, x_previous=None, vmax=None, vmin=None):
    """
    Synthesis filter
    Computes the predicts samples for the given set of errors using the previous samples
    :param e: <[float]> list containing the errors to be decoded
    :param alphas: <[float]> list containing the lpc coefficients (without the 1)
    :param x_previous: <[int]> list containing the previous state
    :return:
        x: <[float]> Predicted samples
        x_previous: <[int]> list containing the current sate
    """
    p = len(alphas)

    if x_previous is None:
        x_previous = [0 for i in range(p)]

    x = []
    for e_n in e:
        x_n = int(round(e_n + np.dot(x_previous, alphas)))
        if vmax is not None:
            if x_n > vmax:
                x_n = vmax
        if vmin is not None:
            if x_n < vmin:
                x_n = vmin

        x.append(x_n)
        x_previous = [x_n] + x_previous[:-1]

    return x, x_previous
