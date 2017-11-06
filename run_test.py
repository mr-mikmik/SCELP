
import matplotlib.pyplot as plt

import Spherical_Codebook
import coding_scheme
import random

# VARIABLES:
Lv = 3
N_sp = 10
M_r = 10
p = 10
window_size = 100

# Create the codebook

codebook = Spherical_Codebook.S_Codebook(Lv, N_sp, M_r)

print 'Peel-list::::::::::::::::::::'
print codebook.peelist
print 'Number of centroids: '+str(codebook.centroids_count)
print '_______________________________________________'

print 'X'
x = [random.randint(0,255) for i in range(100)]
print x



print 'Coding frames'
#codewords_indxs, lpcs, err, num_err = coding_scheme.encode_simple_errors(x, codebook, window_size, p)
codewords_indxs, lpcs = coding_scheme.encode_simple(x, codebook, window_size, p)

print 'Decoding frames'
#x_q = coding_scheme.decode_errors(codewords_indxs, lpcs, codebook, err, num_err)
x_q = coding_scheme.decode(codewords_indxs, lpcs, codebook)



print 'X_q:'
print x_q

plt.figure(1)
#plt.subplot(211)
plt.plot(x,'r-o', label='Original sequence')
#plt.subplot(212)
plt.plot(x_q,'g-o',label='Synthetic sequence')
plt.legend()
plt.title('Simple encoding simulation')
plt.show()

print 'DONE :)'
