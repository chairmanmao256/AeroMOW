import numpy as np
from scipy.optimize import fsolve

'''
### Description

Use the wave interaction method to design a nozzle that is free of waves.

Written for the undergraduate course: Fluid Mechanics and Aerodynamics, Tsinghua University.

Author: Chenyu Wu
'''

# some global parameters and functions
ga = 1.4
ga1m = ga - 1.0
ga1p = ga + 1.0
gar = ga1m / ga1p

def PM_res(M, nu):
    # note that the input nu should be in radians
    res = np.sqrt(1.0/gar)*np.arctan(np.sqrt(gar*(M**2 - 1.0))) - np.arctan(np.sqrt(M**2 - 1.0)) - nu
    return res

def intersect(x1: float, y1: float, k1: float,
              x2: float, y2: float, k2: float)->tuple:
    '''
    ### Description

    Calculate the intersection point of:
    y - y1 = k1(x - x1)
    y - y2 = k2(x - x2)
    '''
    x = (y2 - y1 + x1 * k1 - x2 * k2) / (k1 - k2)
    y = (x2 - x1 + y1 / k1 - y2 / k2) / (1.0/k1 - 1.0/k2)
    return x, y


class nozzle:
    '''
    ### Description
     
    This class uses the method of wave to design the supersonic nozzle

    '''
    def __init__(self, dtheta: float, dx: float, dy: float, nwave: int, nu0 = 0.0):
        '''
        ### Description

        Initialize the instance
        '''
        self.dtheta = dtheta
        self.dx = dx
        self.dy = dy
        self.nwave = nwave

        # stores the deflection angle and nu for each cell
        self.cell_theta = np.zeros((nwave+1, nwave+1))
        self.cell_nu    = np.zeros((nwave+1, nwave+1))
        self.cell_Ma    = np.zeros((nwave+1, nwave+1))

        # stores the position of every node
        self.node_x = np.zeros((nwave+1, nwave+1))
        self.node_y = np.zeros((nwave+1, nwave+1))

        # calculate the deflection angle and nu for each cell
        for i in range(nwave+1):
            for j in range(nwave + 1):
                self.cell_theta[i,j] = i * dtheta - j * dtheta
                self.cell_nu[i,j]    = i * dtheta + j * dtheta + nu0
                self.cell_Ma[i,j]    = self.calc_Ma(self.cell_nu[i,j])

        # give the coordinate on the wall
        self.node_y[1, 0] =  dy
        self.node_y[0, 1] = -dy
        for i in range(2, nwave+1):
            self.node_x[i, 0] = (i-1) * dx
            self.node_y[i, 0] = dx * np.tan(dtheta*np.pi/180.0 * (i-1)) + self.node_y[i-1, 0]

            self.node_x[0, i] = (i-1) * dx
            self.node_y[0, i] = - dx * np.tan(dtheta*np.pi/180.0 * (i-1)) + self.node_y[0, i-1]


    def calc_Ma(self, nu: float)-> float:
        '''
        ### Description

        Calculate Ma from nu
        '''
        # first convert degree to radian
        nu = nu / 180.0 * np.pi
        M  = fsolve(PM_res, 1.1, args=(nu))[0]

        return M

    def calc_wave(self, M1: float, M2: float, theta1: float, theta2: float, type_ = 0):
        '''
        ### Description

        Determine the wave angle based on the information before and behind the wave. 
        type == 0: wave from upper side, going downward
        type == 1: wave from lower side, going upward

        The returned wave angle is also in degree.
        '''
        # Mach angle in degree
        mu1 = np.arcsin(1.0 / M1) / np.pi * 180.0
        mu2 = np.arcsin(1.0 / M2) / np.pi * 180.0

        if type_ == 0:
            phi = (theta1 + theta2 - mu1 - mu2) / 2.0
        else:
            phi = (theta1 + theta2 + mu1 + mu2) / 2.0

        return phi
    
    def calc_child_node(self, i: int, j: int):
        '''
        ### Description
        Calcualte the child node (i,j), which is aligned with:

                (i, j-1)
                        \
                         \  cell(i, j-1)
                          \
          cell(i-1, j-1)   (i, j), child
                          /
                         /  cell(i-1, j)
                        /
                (i-1, j)
        '''
        # Calculate the slope and the parent nodes's position
        phi1 = self.calc_wave(self.cell_Ma[i-1, j-1], self.cell_Ma[i, j-1],
                              self.cell_theta[i-1, j-1], self.cell_theta[i,j-1],
                              type_ = 0)
        
        phi2 = self.calc_wave(self.cell_Ma[i-1, j-1], self.cell_Ma[i-1, j],
                              self.cell_theta[i-1, j-1], self.cell_theta[i-1,j],
                              type_ = 1)
        
        k1, k2 = np.tan(phi1 / 180.0 * np.pi), np.tan(phi2 / 180.0 * np.pi)
        x1, x2 = self.node_x[i, j-1], self.node_x[i-1, j]
        y1, y2 = self.node_y[i, j-1], self.node_y[i-1, j]

        x, y = intersect(x1, y1, k1, x2, y2, k2)

        self.node_x[i,j], self.node_y[i,j] = x, y

    def march(self):
        '''
        ### Description

        March from the inlet to outlet to get all child nodes' positions
        '''
        for k in range(2, 2 * self.nwave+1):
            if k <= self.nwave+1:
                for i in range(1, k):
                    # print("Calculating child node: ({},{})".format(i, k-i))
                    self.calc_child_node(i, k-i)
            else:
                for i in range(k - self.nwave, self.nwave+1):
                    # print("Calculating child node: ({},{})".format(i, k-i))
                    self.calc_child_node(i, k-i)

    def find_parent(self, i: int, j: int)->tuple:
        '''
        ### Description

        Find the parent of a child node. The nodes are aligned as follows:

                (i, j-1)
                        \
                         \  cell(i, j-1)
                          \
          cell(i-1, j-1)   (i, j), child
                          /
                         /  cell(i-1, j)
                        /
                (i-1, j)
        '''

        return [i, j-1], [i-1, j]

    def elimination(self):
        '''
        ### Description

        Calculate the nodes' position in the wave-elimination region
        '''
        nwave = self.nwave
        dtheta = self.dtheta
        xs, ys = self.node_x[nwave, 0], self.node_y[nwave, 0]

        x_list = [xs]
        y_list = [ys]

        for i in range(1, nwave+1):
            phi = self.calc_wave(self.cell_Ma[nwave, i-1], self.cell_Ma[nwave, i],
                                 self.cell_theta[nwave, i-1], self.cell_theta[nwave, i],
                                 type_ = 1)
            
            k1= np.tan(phi / 180.0 * np.pi)
            x1, y1 = self.node_x[nwave, i], self.node_y[nwave, i]

            k2 = np.tan((nwave-i+1)*dtheta / 180.0 * np.pi)
            x2, y2 = x_list[-1], y_list[-1]

            x, y = intersect(x1, y1, k1, x2, y2, k2)

            x_list.append(x)
            y_list.append(y)

        self.comp_xu = np.array(x_list)
        self.comp_yu = np.array(y_list)
        self.comp_xl = np.array(x_list)
        self.comp_yl = -1.0*np.array(y_list)
