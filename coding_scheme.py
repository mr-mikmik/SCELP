
import LPC_estimator
import LPC_filter
import Spherical_Codebook

def encode_simple(x, codebook, window_size, p):
    """

    :param x:
    :param codebook:
    :param window_size:
    :param p: <int> order of the lpc predictor
    :return:
        codewords_indxs: <[[(int, int),..],..]> list of list of tuples of 2 ints representing the codeword indexs
        lpcs: <[[float],]> list containing lists of floats representing the lpc coeffs for each windowed samples
    """
    num_samples = len(x)
    num_w_smp = int(num_samples / window_size) + 1
    lv = codebook.Lv # Dimension of the VQ
    state = None
    codewords_indxs = []
    lpcs = []
    for w in range(num_w_smp):
        cws = []
        # Windowed sample
        # TODO: Add zero padding
        x_w = x[w*window_size : (w+1)*window_size]
        # Get the LPC coefs for that subsamples:
        lpc_coefs = LPC_estimator.lpc_coefficients(x_w, p)
        alphas = lpc_coefs[1:] # Exclude the first 1
        # Compute the errors
        errors, state = LPC_filter.Az_filter(x_w, alphas, state)

        # Encode the errors into the codebook indexes:
        num_cws = int(num_w_smp/lv) + 1
        for c in range(num_cws):
            # TODO: Add zero padding
            e = errors[c*lv : (c+1)*lv]
            # Get the index for the current error vector
            cw_indx, g_indx = codebook.encode(e)
            cws.append((cw_indx, g_indx))
        codewords_indxs.append(cws)
        lpcs.append(lpc_coefs)

    return codewords_indxs, lpcs


def encode(x, codebook, window_size, p):
    """
    Encodeing using the quantified values in the predictor
    :param x:
    :param codebook:
    :param window_size:
    :param p: <int> order of the lpc predictor
    :return:
        codewords_indxs: <[[(int, int),..],..]> list of list of tuples of 2 ints representing the codeword indexs
        lpcs: <[[float],]> list containing lists of floats representing the lpc coeffs for each windowed samples
    """
    num_samples = len(x)
    num_w_smp = int(num_samples / window_size) + 1
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

        errors, state = LPC_filter.Az_filter(x_w, alphas, state)

        # Encode the errors into the codebook indexes:
        num_cws = int(num_w_smp/lv) + 1
        for c in range(num_cws):
            # TODO: Add zero padding
            e = errors[c*lv, (c+1)*lv]
            # Get the index for the current error vector
            cw_indx, g_indx = codebook.encode(e)
            cws.append((cw_indx, g_indx))
        codewords_indxs.append(cws)
        lpcs.append(lpc_coefs)

    return codewords_indxs, lpcs


def decode(codewords_indxs, lpcs, codebook):
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
        x_w_q, state = LPC_filter.Sz_filter(errors_q, alphas, state)
        x_q += x_w_q

    return x_q