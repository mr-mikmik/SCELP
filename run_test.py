import math
import matplotlib.pyplot as plt
import numpy as np

import Spherical_Codebook
import coding_scheme
import random

def init_gain_dic(Mr, v_max,v_min):
    d_indx = {}
    d_gain = {}
    lamb = (v_max - v_min)*1.0 / Mr
    for i in range(Mr):
        d_indx[i] = lamb*i
        d_gain[i] = (lamb*i, lamb*(i+1))
    print d_gain
    print '___'
    print d_indx
    return d_indx, d_gain

def gain_quantization(g,d_indx, d_gain):
    g_q = 0
    for k,v in d_gain.items():
        if v[0] <= g <= v[1]:
            g_q = k
            break
    return d_indx[g_q]

def simulation_random():
    # VARIABLES:
    Lv = 3
    N_sp = 10
    M_r = 2**10
    p = 10
    window_size = 30

    # Create the codebook

    codebook = Spherical_Codebook.S_Codebook(Lv, N_sp, M_r)

    print 'Peel-list::::::::::::::::::::'
    print codebook.peelist
    print 'Number of centroids: '+str(codebook.centroids_count)
    print '_______________________________________________'

    print 'X'
    x = [random.randint(0,255) for i in range(90)]
    print x



    print 'Coding frames'
    #codewords_indxs, lpcs, err, num_err = coding_scheme.encode_simple_errors(x, codebook, window_size, p)
    codewords_indxs, lpcs = coding_scheme.encode_simple(x, codebook, window_size, p)

    print 'Decoding frames'
    #x_q = coding_scheme.decode_errors(codewords_indxs, lpcs, codebook, err, num_err)
    x_q = coding_scheme.decode(codewords_indxs, lpcs, codebook)



    print 'X_q:'

    print x_q

    x = np.array(x)
    x_q = np.array(x_q)
    error_p = x - x_q
    pot_e_q = np.sum(error_p ** 2) / (1.0 * len(x))
    pot_e_p = np.sum(x ** 2) / (1.0 * len(x))
    print 'Pot_e_q: %.4f' % pot_e_q
    print 'Pot_e_p: %.4f' % pot_e_p
    SNR = 10 * math.log(pot_e_p / pot_e_q, 10)
    print 'SNR: ' + str(SNR)

    plt.figure(1)
    #plt.subplot(211)
    plt.plot(x,'r-o', label='Original sequence')
    #plt.subplot(212)
    plt.plot(x_q,'g-o',label='Synthetic sequence')
    plt.legend()
    plt.title('Simple encoding simulation')
    plt.show()

    print 'DONE :)'


def simulation_errors():
    # VARIABLES:
    Lv = 5
    N_sp = 10
    M_r = 2**6
    p = 10
    window_size = 50
    TOTAL_NUM_SAMPLES = 100
    Muniform = 2**4

    # Create the codebook

    codebook = Spherical_Codebook.S_Codebook(Lv, N_sp, M_r)
    d1,d2 = init_gain_dic(Muniform,500,0)
    print 'Peel-list::::::::::::::::::::'
    #print codebook.peelist
    print 'Number of centroids: ' + str(codebook.centroids_count)
    print math.ceil(math.log(codebook.centroids_count,2))
    print '_______________________________________________'

    print 'X'
    x = [random.randint(0, 255) for i in range(TOTAL_NUM_SAMPLES)]
    print x

    print 'Coding frames'
    # codewords_indxs, lpcs, err, num_err = coding_scheme.encode_simple_errors(x, codebook, window_size, p)
    codewords_indxs, lpcs, e, _ = coding_scheme.encode_simple_errors(x, codebook, window_size, p)

    print 'Decoding frames'
    # x_q = coding_scheme.decode_errors(codewords_indxs, lpcs, codebook, err, num_err)
    e_q = coding_scheme.decode_error_extraction(codewords_indxs, lpcs, codebook)

    e_uniform = [ gain_quantization(x_i,d1,d2) for x_i in e]

    e_q = np.array(e_q)
    e_uniform = np.array(e_uniform)
    e = np.array(e)
    error = e - e_q
    error_u = e - e_uniform
    pot_e = np.sum(error ** 2) / (1.0 * len(e))
    pot_e_u = np.sum(error_u ** 2) / (1.0 * len(e))
    pot_e_p = np.sum(e ** 2) / (1.0 * len(e))

    print 'Pot_e: %.4f' % pot_e
    print 'Pot_e_u: %.4f' % pot_e_u
    print 'Pot ep: %.4f' % pot_e_p
    SNR = 10 * math.log(pot_e_p / pot_e, 10)
    SNR_u = 10 * math.log(pot_e_p / pot_e_u, 10)
    print 'SNR: ' + str(SNR)
    print 'SNR_u: ' + str(SNR_u)

    plt.figure(1)
    # plt.subplot(211)
    plt.plot(e, 'r-o', label='errors')
    # plt.subplot(212)
    plt.plot(e_q, 'g-o', label='errors_q')
    plt.plot(e_uniform, 'b-o', label='uniform_q')
    plt.legend()
    plt.title('Simple encoding simulation')
    plt.show()

    print 'DONE :)'




simulation_random()
