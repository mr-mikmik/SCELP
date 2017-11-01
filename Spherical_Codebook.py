import numpy as np
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
        
        self.centroids_count = 0
        self.c_indx_to_coords = {}
        # peelist: List contaning 
        self.centroids, self.peelist = self.init_centroids(self.c_indx_to_coords)
        self.c_indx_to_cart = self.init_cartesians_dic()
        self.c_indx_to_angles = self.init_angle_dic()


    def preselection(self, d0):
        """
        Gets the candidates to be quantized given the vector d0. To do so w
        :param d0: <[floats]> vector to be quantized
        :return: <[int]> a list containing the indexes (int) referring to the candidates centroids
        """
        candidates = []
        # 1 - Normalize the vector
        d0_norm = np.linalg.norm(d0)
        c0 = d0/d0_norm
        # 2 - Convert the vector to spherical coordinates
        _, sph_c0 = self.cartesian2spherical(c0)
        # 3 - Select the candidiates making use of the peelist
        peeling_centroids = self.centroids

        def search(values,lists):
            cand = []
            # Check if we are on the last layer of the lists:
            if type(lists[0][1]) is not list:
                for i in range(len(lists)):
                    if values[0] <= lists[i][0]:
                        cand.append(values[i])
                        cand.append(values[i-1])
            else:
                # If we aren't still on the last layer:
                if values[0] < lists[0][0]:
                    cand = search(values[1:], lists[0][1])
                elif values[0] > lists[-1][0]:
                    cand = search(values[1:], lists[-1][1])
                else:
                    for i in range(len(lists)-1):
                        if lists[i][0] <= values[0] <= lists[i+1][0]:
                            cand_rec_1 = search(values[1:], lists[i][1])
                            cand_rec_2 = search(values[1:], lists[i][1])
                            cand = cand + cand_rec_1 + cand_rec_2
                            break
            return cand

        candidates = search(sph_c0, self.peelist)

        return candidates

    ## INITIALIZE FUNCTIONS:

    def init_centroids(self, d, lv_i=0, previous=()):
        """
        Initialize the centroids of the shape codebook. This is done following the apple-peeling method
        It stores the result in the c parameter
        :param d: <dic> dictionary with centroid index as a key and the dimension's coordinates as value
        :param lv_i: <int> current dimension of the apple-peeling process
        :param previous: <int> index of the previous dimension
        :return:
            centroids parameter: <list> centroids list where the centroids are stored
            peelist: <[[float,list],[float, list],...]> list formed by pairs of angles (floats) and lists
                each sublist if formed by another list of the same type
                it turs to have the following shape: [float,...[[float, int],[float, int]...]..]> 
        """
        # Number of angles -> Lv - 1
        if lv_i < self.Lv - 2 :
            peelist = []
            c = []
            for i in range(self.N_sp):
                phi_0 = (i + 0.5) * self.theta
                peelist.append([phi_0,None])
                # Call the same function to compute recursively
                c[i], peelist[i][1] = self.init_centroids(d, lv_i + 1, previous=previous + (i,))
            return c, peelist
        else:
            # The last layer of the apple peeling method
            phi_p = (previous[-1] + 0.5)*self.theta
            Nspl = self.get_Nspl(phi_p)
            coords_last_layer = []
            last_layer_peelist = []
            for i in range(Nspl):
                coords_last_layer.append(self.centroids_count)
                d[self.centroids_count] = previous + (i,)
                phi_last = (i + 0.5) * 2 * math.pi / Nspl
                last_layer_peelist.append([phi_last, self.centroids_count])
                self.centroids_count += 1
            return coords_last_layer, last_layer_peelist

    def init_cartesians_dic(self):
        """
        Function to create a dictionary with centroids indexes as keys and a list of cartesian coordenates as values
        :return: <{int:(float)}> Dictionary
            keys: Centroids indexes
            values: List containing the cartesian coordenates
        """
        # Make use of the already created dictionary containing the centroid indexes and angle indexes (coords)
        d = self.c_indx_to_coords
        cartesians_dic = {}
        for k, coords in d.items():
            cartesians = ()
            # For all the components except the last 2 we make use of the already computed values of sin and cos
            for i, c in enumerate(coords[:-1]):
                res = 1
                for j in range(i):
                    res *= self.sin_d[coords[i]]
                res = res * self.cos_d[c]
                cartesians = cartesians + (res,)
            res = 1
            # For the last 2 coordinate we need to compute the last angle (phi_p) as it is different form the rest
            for j in range(len(coords)-1):
                res *= self.sin_d[coords[j]]
            phi_p = (coords[-1] + 0.5) * 2 * math.pi / self.get_Nspl((coords[-2] + 0.5)*self.theta)
            rn_1 = res * math.cos(phi_p)
            rn = res * math.sin(phi_p)
            cartesians = cartesians + (rn_1, rn)
            cartesians_dic[k] = cartesians

        return cartesians_dic



    def init_angle_dic(self,):
        d = self.c_indx_to_coords
        angle_dic = {}
        for k, coords in d.items():
            angles = ()
            for i, c in enumerate(coords[:-1]):
                a = (c + 0.5)*self.theta
                angles += (a,)
            phi_p = (coords[-1] + 0.5) * 2 * math.pi / self.get_Nspl((coords[-2] + 0.5) * self.theta)
            angles += (phi_p,)
            angle_dic[k] = angles
        return angle_dic

    def init_sin_dic(self):
        """
        Initialize a dictionary containing the sin of all the predefined angles
        :return: <{int : float}> Dictionary with
            key: angle index
            value: sin of that angle
        """
        d = {}
        for i in range(self.N_sp):
            d[i]=math.sin((i + 0.5)*self.theta)
        return d

    def init_cos_dic(self):
        """
        Initialize a dictionary containing the cos of all the predefined angles
        :return: <{int : float}> Dictionary with
            key: angle index
            value: cos of that angle
        """
        d = {}
        for i in range(self.N_sp):
            d[i] = math.cos((i + 0.5)*self.theta)
        return d


    ## GETTERS & AUXILIAR FUNCTIONS:

    def get_Nspl(self, phi_i):
        """
        Auxiliar function to compute the number of centroids in the last layer
        (centroids equally distributed along a circle)
        :param phi_i: Value of the angle form the previous apple-peeling layer coordinate
        :return: <int> Number of centroids for that layer
        """
        N_sp_l = int(2 * math.pi / self.theta * math.sin(phi_i))
        return N_sp_l

    def get_centroids(self):
        return self.centroids

    def get_coords(self, index):
        return self.c_indx_to_coords[index]

    def get_cartesians_dic(self):
        return self.c_indx_to_cart

    def cartesian2spherical(self, c0):
        """
        Converts the vector to is spherical equivalent
        :param c0: cartesian n-dimensional vector
        :return: <float, (float)> modulus, n-1-dimensional list containing the spherical coordinates
        """
        # TODO: Debugg division by Zero and ANGLES!!!!!!!
        sph_coords = ()     # Vector containg the angles
        modulus = np.linalg.norm(c0)
        for i, c in enumerate(c0[:-2]):
            mod_i = np.linalg.norm(c0[i:])
            phi_i = math.acos(c/mod_i)
            sph_coords += (phi_i,)
        # Last angle:
        last_mod = np.linalg.norm(c0[-2:])
        # Check divison by zero case:
        if c0[-2] == 0 and last_mod == 0:
            if c0[-1] >= 0:
                phi_l = math.acos(c0[-2] / last_mod)
            else:
                phi_l = 2*math.pi-math.acos(c0[-2] / last_mod)
        else:
            if c0[-1] >= 0:
                phi_l = math.acos(c0[-2] / last_mod)
            else:
                phi_l = 2*math.pi-math.acos(c0[-2] / last_mod)
        sph_coords += (phi_l,)

        return modulus, sph_coords

