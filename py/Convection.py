import numpy as np
from math import sin, cos, exp, pi

class Convection:
    def __init__(self, nu0, nufreq, cx, cy, cz):
        self.nu0 = nu0
        self.nufreq = nufreq
        self.cx = cx
        self.cy = cy
        self.cz = cz

    def getExact(self, x, y, z, t):
        a = self.nu0 * t \
            - (self.nu0 / (2.*self.nufreq)) * (cos(self.nufreq*t)-1.) \
          if self.nufreq > 0. \
          else self.nu0 * t
        exact = sin(2.*pi*(x-self.cx*t)) \
              * sin(2.*pi*(y-self.cy*t)) \
              * sin(2.*pi*(z-self.cz*t)) \
              * exp(-12.*pi*pi*a)
        return exact

    def getSolution(self, gridSize, t):
        solution = np.zeros((gridSize, gridSize, gridSize))
        dx = 1./gridSize
        x = dx/2.
        for i in range(gridSize):
            y = dx/2.
            for j in range(gridSize):
                z = dx/2.
                for k in range(gridSize):
                    solution[i, j, k] = self.getExact(x, y, z, t)
                    z = z+dx
                y = y+dx
            x = x+dx
        return solution
                    
if __name__ == '__main__':
    from matplotlib import pyplot as plt
    c = Convection(0.1, 100., 1., 1., 1.)
    f = c.getSolution(128, 0.1)
    plt.imshow(f[:,:,30])
    plt.colorbar()
    plt.show()

