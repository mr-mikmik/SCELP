import numpy
import math

class S_Codebook():
    def __init__(self, Lv, N_sp, M_r ):
        """
        Inicialize the codebook
        :param Lv: <int> Dimention of the space
        :param N_sp: <int> Number of point for each pi arch
        :param M_r: <int> Number of codewords in the gain (radius) codebook
        """
        self.Lv = Lv
        self.N_sp = N_sp
        self.M_r = M_r

        # Compute the angular distance between two centroids in a pi arch (all expcept the last one)
        self.theta = math.pi/N_sp

        # Initialize the sin and cos dictionaries *(to improve the computation efficiency)
        self.sin_d = self.init_sin_dic()
        self.cos_d = self.init_cos_dic()

        # Initialize the spherical codebook
        self.centroids = []
        self.centroids_count = 0
        self.c_indx_to_coords = {}
        self.init_Centroids(self.centroids, self.c_indx_to_coords)
        self.c_indx_to_cart = self.init_cartesians_dic(self.c_indx_to_coords)

    def init_Centroids(self, c, d, lv_i=0,previous=()):
        """
        Initialize the centroids of the shape codebook. This is done following the apple-peeling method
        It stores the result in the c parameter
        :param c: <list> centroids list where the centroids are stored
        :param d: <dic> dictionary with centroid index as a key and the dimension's coordinates as value
        :param lv_i: <int> current dimension of the apple-peeling process
        :param previous: <int> index of the previous dimension
        :return: Nothing
        """
        # Number of angles -> Lv - 1
        if lv_i < self.Lv - 2 :
            for i in range(self.N_sp):
                c.append([])
                # Call the same function to compute recursively
                c[i] = self.init_Centroids(c[i], d, lv_i+1, previous=previous+(i,))
            return c
        else:
            # The last layer of the apple peeling method
            phi_p = (previous[-1] + 0.5)*self.theta
            Nspl = self.get_Nspl(phi_p)
            coords_actual_layer = []
            for i in range(Nspl):
                coords_actual_layer.append(self.centroids_count)
                d[self.centroids_count] = previous + (i,)
                self.centroids_count += 1
            return coords_actual_layer

    def init_cartesians_dic(self, d):
        cartesians_dic = {}
        for k, coords in d.items():
            cartesians = ()
            for i, c in enumerate(coords[:-1]):
                res = 1
                for j in range(i):
                    res *= self.sin_d[coords[i]]
                res = res * self.cos_d[c]
                cartesians = cartesians + (res,)
            res = 1
            for j in range(len(coords)-1):
                res *= self.sin_d[coords[j]]
            phi_p = (coords[-1] + 0.5) * 2 * math.pi / self.get_Nspl((coords[-2] + 0.5)*self.theta)
            rn_1 = res * math.cos(phi_p)
            rn = res * math.sin(phi_p)
            cartesians = cartesians + (rn_1, rn)
            cartesians_dic[k] = cartesians
            print '==================='
            print coords
            print cartesians
        return cartesians_dic




    def init_sin_dic(self):
        d = {}
        for i in range(self.N_sp):
            d[i]=math.sin((i + 0.5)*self.theta)
        return d

    def init_cos_dic(self):
        d={}
        for i in range(self.N_sp):
            d[i] = math.cos((i + 0.5)*self.theta)
        return d

    def get_Nspl(self, phi_i):
        N_sp_l = int(2 * math.pi / self.theta * math.sin(phi_i))
        return N_sp_l

    def get_centroids(self):
        return self.centroids

    def get_coords(self, index):
        return self.c_indx_to_coords[index]

    def get_cartesians_dic(self):
        return self.c_indx_to_cart