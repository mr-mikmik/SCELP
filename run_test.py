import Spherical_Codebook
import coding_scheme
import random

# VARIABLES:
Lv = 3
N_sp = 6
M_r = 6
p = 10
window_size = 30

# Create the codebook

codebook = Spherical_Codebook.S_Codebook(Lv, N_sp, M_r)

print 'Peel-list::::::::::::::::::::'
print codebook.peelist
print '_______________________________________________'

print 'X'
x = [random.randint(0,255) for i in range(30)]
print x



print 'Coding frames'
codewords_indxs, lpcs = coding_scheme.encode_simple(x, codebook, window_size, p)

print 'Decoding frames'
x_q = coding_scheme.decode(codewords_indxs, lpcs, codebook)



print 'X_q:'
print x_q

print 'DONE :)'
