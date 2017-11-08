import math
import numpy as np

import LPC_estimator
import LPC_filter
import Spherical_Codebook

# VALUES USED TO LIMIT THE OUTPUT:
V_MAX = 255
V_MIN = 0

def encode_simple(x, codebook, window_size, p):
    """
    Encodes the given sequence into codeword index using a simple model estimating the prediction from the unquanitzed samples
    :param x: <[]> Sequence to be encoded
    :param codebook: Codebook object
    :param window_size: <int> size of the sequence to estimate the LPC coeffs from.
    :param p: <int> order of the lpc predictor
    :return:
        codewords_indxs: <[[(int, int),..],..]> list of list of tuples of 2 ints representing the codeword indexs
        lpcs: <[[float],]> list containing lists of floats representing the lpc coeffs for each windowed samples
    """
    num_samples = len(x)
    print num_samples
    num_w_smp = int(math.ceil(num_samples*1.0 / window_size))
    print 'num_w_smp: '+str(num_w_smp)
    lv = codebook.Lv # Dimension of the VQ
    state = None
    codewords_indxs = []
    lpcs = []
    for w in range(num_w_smp):
        cws = []
        # Windowed sample
        if w == num_w_smp - 1 and num_w_smp > int(num_samples / window_size):
            # In this case, we need more samples (last window and less samples thant the window size)
            # Addig zeros on the last windowed
            x_w_not_full = x[w * window_size:]
            num_samples_windowed = len(x_w_not_full)
            num_zeros_padding = window_size - num_samples_windowed
            zeros_padd = [0 for i in range(num_zeros_padding)]
            x_w = x_w_not_full + zeros_padd
        else:
            x_w = x[w*window_size : (w+1)*window_size]
        # Get the LPC coefs for that subsamples:
        lpc_coefs = LPC_estimator.lpc_coefficients(x_w, p)
        alphas = lpc_coefs[1:] # Exclude the first 1
        # Compute the errors
        errors, state = LPC_filter.Az_filter(x_w, alphas, state)

        # Encode the errors into the codebook indexes:
        num_cws = int(math.ceil(window_size / lv))
        for c in range(num_cws):
            # TODO: Add zero padding
            e = errors[c*lv : (c+1)*lv]
            # Get the index for the current error vector
            cw_indx, g_indx = codebook.encode(e)
            cws.append((cw_indx, g_indx))
        codewords_indxs.append(cws)
        lpcs.append(lpc_coefs)

    return codewords_indxs, lpcs

def encode_simple_errors(x, codebook, window_size, p):
    """
    Same as encode_simple but also returns the errors
    :param x: <[]> Sequence to be encoded
    :param codebook: Codebook object
    :param window_size: <int> size of the sequence to estimate the LPC coeffs from.
    :param p: <int> order of the lpc predictor
    :return:
        codewords_indxs: <[[(int, int),..],..]> list of list of tuples of 2 ints representing the codeword indexs
        lpcs: <[[float],]> list containing lists of floats representing the lpc coeffs for each windowed samples
        err: List containing the errors without encoding
        num_err: <int> Number of errors encoded in each subvector (is the same as the window size)
    """
    num_samples = len(x)
    print num_samples
    num_w_smp = int(math.ceil(num_samples*1.0 / window_size))
    print 'num_w_smp: '+str(num_w_smp)
    lv = codebook.Lv # Dimension of the VQ
    state = None
    codewords_indxs = []
    lpcs = []
    err =[]
    for w in range(num_w_smp):
        cws = []
        # Windowed sample
        if w == num_w_smp - 1 and num_w_smp > int(num_samples / window_size):
            # In this case, we need more samples (last window and less samples thant the window size)
            # Addig zeros on the last windowed
            x_w_not_full = x[w * window_size:]
            num_samples_windowed = len(x_w_not_full)
            num_zeros_padding = window_size - num_samples_windowed
            zeros_padd = [0 for i in range(num_zeros_padding)]
            x_w = x_w_not_full + zeros_padd
        else:
            x_w = x[w*window_size : (w+1)*window_size]
        # Get the LPC coefs for that subsamples:
        lpc_coefs = LPC_estimator.lpc_coefficients(x_w, p)
        alphas = lpc_coefs[1:] # Exclude the first 1
        # Compute the errors
        errors, state = LPC_filter.Az_filter(x_w, alphas, state)

        # Encode the errors into the codebook indexes:
        num_cws = int(math.ceil(window_size / lv))
        for c in range(num_cws):
            # TODO: Add zero padding
            e = errors[c*lv : (c+1)*lv]
            # Get the index for the current error vector
            cw_indx, g_indx = codebook.encode(e)
            cws.append((cw_indx, g_indx))
        codewords_indxs.append(cws)
        lpcs.append(lpc_coefs)
        err += errors
    num_err = window_size
    return codewords_indxs, lpcs, err, num_err

def encode_vq(x, codebook, window_size, p):
    """
    Encoding using the quantified VECTOR values in the predictor
    :param x: <[]> Sequence to be encoded
    :param codebook: Codebook object
    :param window_size: <int> size of the sequence to estimate the LPC coeffs from.
    :param p: <int> order of the lpc predictor
    :return:
        codewords_indxs: <[[(int, int),..],..]> list of list of tuples of 2 ints representing the codeword indexs
        lpcs: <[[float],]> list containing lists of floats representing the lpc coeffs for each windowed samples
    """
    # TODO: NOT IMPLEMENTED YET!!!
    # TODO: FIND AN ALGORITHM TO SOLVE IT
    num_samples = len(x)
    num_w_smp = int(math.ceil(num_samples*1.0 / window_size))
    lv = codebook.Lv    # Dimension of the VQ
    state = np.zeros((p, lv))   # Initialize the state as 0
    codewords_indxs = []
    lpcs = []
    x_q_n = [0 for i in lv]
    for w in range(num_w_smp):
        cws = []
        # Windowed sample
        # TODO: Add zero padding
        if w == num_w_smp-1 and num_w_smp > int(num_samples/window_size):
            # In this case, we need more samples (last window and less samples thant the window size)
            # Addig zeros on the last windowed
            x_w_not_full = x[w*window_size:]
            num_samples_windowed = len(x_w_not_full)
            num_zeros_padding = window_size - num_samples_windowed
            zeros_padd = [0 for i in range(num_zeros_padding)]
            x_w = x_w_not_full + zeros_padd
        else:
            x_w = x[w*window_size:(w+1)*window_size]

        # Get the LPC coeffs for that subsamples:
        lpc_coefs = LPC_estimator.lpc_coefficients(x_w, p)
        alphas = lpc_coefs[1:]  # Exclude the first 1
        # Compute the errors
        p = len(alphas)
        for n in range(int(math.ceil(window_size/lv))):
            x_n = np.array(x_w[n*lv:(n+1)*lv])
            x_p = np.dot(state, alphas)
            e_n = x_n - x_p
            cw_indx, g_indx = codebook.encode(e_n)
            cws.append((cw_indx, g_indx))
            e_q_n = codebook.decode(cw_indx, g_indx)
            x_q_n = e_q_n + x_p
            state = [x_q_n] + state[:-1]

        codewords_indxs.append(cws)
        lpcs.append(lpc_coefs)

    return codewords_indxs, lpcs

def encode(x, codebook, window_size, p):
    """
    Encodeing using the quantified values in the predictor
    :param x: <[]> Sequence to be encoded
    :param codebook: Codebook object
    :param window_size: <int> size of the sequence to estimate the LPC coeffs from.
    :param p: <int> order of the lpc predictor
    :return:
        codewords_indxs: <[[(int, int),..],..]> list of list of tuples of 2 ints representing the codeword indexs
        lpcs: <[[float],]> list containing lists of floats representing the lpc coeffs for each windowed samples
    """
    num_samples = len(x)
    num_w_smp = int(math.ceil(num_samples*1.0/ window_size))
    lv = codebook.Lv    # Dimension of the VQ
    state = [0 for i in range(p)]   # Initialize the state as 0
    codewords_indxs = []
    lpcs = []
    for w in range(num_w_smp):
        cws = []
        # Windowed sample
        # TODO: Add zero padding
        x_w = x[w*window_size : (w+1)*window_size]
        # Get the LPC coeffs for that subsamples:
        lpc_coefs = LPC_estimator.lpc_coefficients(x_w, p)
        alphas = lpc_coefs[1:]  # Exclude the first 1
        # Compute the errors
        p = len(alphas)
        for x_n in x_w:
            x_p = np.dot(state, alphas)
            e_n = x_n - x_p
            cw_indx, g_indx = codebook.encode(e_n)
            cws.append((cw_indx, g_indx))
            e_q_n = codebook.decode(cw_indx, g_indx)
            x_q_n = e_q_n + x_p
            state = [x_q_n] + state[:-1]

        codewords_indxs.append(cws)
        lpcs.append(lpc_coefs)

    return codewords_indxs, lpcs


def decode(codewords_indxs, lpcs, codebook):
    """
    Decoded the given indexes into a decoded sequence.
    :param codewords_indxs: <[[(cw_indx, g_indx),...],...]> List containing sublist for each window.
                            Each sublist contains pairs representing the centroid and gain index
    :param lpcs: <[[1, a1,..],...]> List containing sublist with the lpc coefficient
    :param codebook: Codebook object
    :return: <[float]> List with the decoded sequence
    """
    num_w_smp = len(lpcs)
    state = None
    x_q = []
    for w in range(num_w_smp):
        # - Unwrap the subsamples:
        lpc_coefs = lpcs[w]
        alphas = lpc_coefs[1:]
        cws = codewords_indxs[w]
        errors_q = []
        # Obtain the errors
        for c in cws:
            cw_indx = c[0]
            g_indx = c[1]
            # Reconstruct the encoded elements (indx -> codeword)
            e_q = codebook.decode(cw_indx, g_indx)
            errors_q += e_q
        # 2 - Obtain the predicted value:
        x_w_q, state = LPC_filter.Sz_filter(errors_q, alphas, state, vmax=V_MAX, vmin=V_MIN)
        x_q += x_w_q

    return x_q

def decode_errors(codewords_indxs, lpcs, codebook, err, num_err):
    """
    Decodes the given errors into the sequence

    :param codewords_indxs:  <[[(cw_indx, g_indx),...],...]> List containing sublist for each window.
                            Each sublist contains pairs representing the centroid and gain index
    :param lpcs: <[[1, a1,..],...]> List containing sublist with the lpc coefficient
    :param codebook: Codebook object
    :param err: List containing the errors without encoding
    :param num_err: <int> Number of errors encoded in each subvector (is the same as the window size)
    :return: <[float]> List containing the decoded sequence
    """
    num_w_smp = len(lpcs)
    state = None
    x_q = []
    for w in range(num_w_smp):
        errors = err[w*(num_err): (w+1)*num_err]
        # - Unwrap the subsamples:
        lpc_coefs = lpcs[w]
        alphas = lpc_coefs[1:]
        cws = codewords_indxs[w]
        errors_q = []
        # Obtain the errors
        for c in cws:
            cw_indx = c[0]
            g_indx = c[1]
            # Reconstruct the encoded elements (indx -> codeword)
            e_q = codebook.decode(cw_indx, g_indx)
            errors_q += e_q
        # 2 - Obtain the predicted value:
        x_w_q, state = LPC_filter.Sz_filter(errors, alphas, state, vmax=V_MAX, vmin=V_MIN)
        x_q += x_w_q

    return x_q

def decode_error_extraction(codewords_indxs, lpcs, codebook):
    """
    Decodes the given codewords indexes into the related codeword representing the errors
    :param codewords_indxs: <[[(cw_indx, g_indx),...],...]> List containing sublist for each window.
                            Each sublist contains pairs representing the centroid and gain index
    :param lpcs: <[[1, a1,..],...]> List containing sublist with the lpc coefficient
    :param codebook: Codebook object
    :return: <[float]> List with the decoded error sequence
    """
    num_w_smp = len(lpcs)
    state = None
    errors_q = []
    for w in range(num_w_smp):
        # - Unwrap the subsamples:
        lpc_coefs = lpcs[w]
        alphas = lpc_coefs[1:]
        cws = codewords_indxs[w]
        # Obtain the errors
        for c in cws:
            cw_indx = c[0]
            g_indx = c[1]
            # Reconstruct the encoded elements (indx -> codeword)
            e_q = codebook.decode(cw_indx, g_indx)
            errors_q += e_q

    return errors_q