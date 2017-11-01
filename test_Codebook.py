import Spherical_Codebook as s

sc = s.S_Codebook(3,3,0)

c = sc.get_centroids()
print c
print '***********************************************'
print sc.c_indx_to_coords
print '***********************************************'
k =  sc.get_cartesians_dic()
print '==========================='
print '==========================='
print sc.peelist
print '........................'
print sc.preselection([1,0,0])

