import Spherical_Codebook as s
import matplotlib.pyplot as plt
import random
import numpy as np
import math

# PARAMETERS TO TUNE:
lv = 10
Nsp = 3
Mr = 2**6

# Codebook creation
sc = s.S_Codebook(lv,Nsp,Mr)


print '***********************************************'
#print sc.c_indx_to_coords
print '***********************************************'
k =  sc.get_cartesians_dic()

#print sc.peelist
print 'Number of centroids = '+ str(sc.centroids_count)
print '==========================='
print '==========================='
n_sampl = 60

x = x = [random.randint(0, 255) for i in range(n_sampl)]
x_q = []
for i in range(n_sampl/lv):
    x_w = x[i*lv: (i+1)*lv]
    cw_i, g_i = sc.encode(x_w)
    res_q = sc.decode(cw_i, g_i)
    x_q += [rr for rr in res_q]

#print 'x: ' + str(x)
print '.....'
#print 'x_q:' + str(x_q)

x = np.array(x)
x_q = np.array(x_q)
error_p = x-x_q
pot_e_q = np.sum(error_p**2)/(1.0*len(x))
pot_e_p = np.sum(x**2)/(1.0*len(x))
print 'Pot_e_q: %.4f' % pot_e_q
print 'Pot_e_p: %.4f' % pot_e_p
SNR = 10*math.log(pot_e_p/pot_e_q,10)
print 'SNR: '+str(SNR)

plt.figure(1)
plt.plot(x, 'r-o', label='original samples')
plt.plot(x_q, 'g-o', label='samples quantified')
plt.legend()
plt.title('Codebook simulation using: lv=%d, Nsp=%d, Mr=%d bits'%(lv,Nsp,math.log(Mr,2)))
plt.show()






