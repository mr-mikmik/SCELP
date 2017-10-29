import numpy
import math

class S_Codebook():
    def __init__(self, Lv, N_sp, M_r ):
        """
        Inicialize the codebook
        :param lv: <int> Dimention of the space
        :param n_sp: <int> Number of point for each pi arch
        :param m_r: <int> Number of codewords in the gain (radius) codebook
        """
        self.Lv = Lv
        self.N_sp = N_sp
        self.M_r = M_r

        # Compute the angular distance between two centroids in a pi arch (all expcept the last one)
        self.theta = math.pi/N_sp

        # Compute the number of centroids in the last arch (arch of 2pi rad)


        # Initialize the spherical codebook:
        self.centroids = []
        self.centroids_count = 0
        self.init_Centroids(self.centroids)

    def init_Centroids(self, c, lv_i=0,previous=-1):
        # Number of angles -> Lv - 1
        if lv_i < self.Lv - 2 :
            for i in range(self.N_sp):
                c.append([])
            for i in range(self.N_sp):
                c[i] = self.init_Centroids(c[i], lv_i+1, previous=i)
            print c
            return c
        else:
            phi_p = (previous + 0.5)*self.theta
            Nspl = self.get_Nspl(phi_p)
            coords_actual_layer = []
            for i in range(Nspl):
                coords_actual_layer.append(self.centroids_count)
                self.centroids_count += 1
            print coords_actual_layer
            return coords_actual_layer


    def get_Nspl(self, phi_i):
        N_sp_l = int(2 * math.pi / self.theta * math.sin(phi_i))
        return N_sp_l

    def get_centroids(self):
        return self.centroids