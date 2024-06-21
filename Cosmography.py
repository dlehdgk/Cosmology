import numpy as np
import scipy.integrate as sp


class Cosmo:
    def __init__(self, h=0.67, Om=0.3, Or=1e-5, Ol=0.63):
        self.__h = h
        self.__om = Om
        self.__or = Or
        self.__ol = Ol

    def hubble(self):
        return self.__h

    def mass(self):
        return self.__om

    def radiation(self):
        return self.__or

    def constant(self):
        return self.__ol

    def curvature(self):
        return 1 - self.Om() - self.Or() - self.Ol()

    def rate(self, z):
        E = np.sqrt(
            self.mass() * (1 + z) ** 3
            + self.radiation() * (1 + z) ** 4
            + self.constant()
            + self.curvature() * (1 + z) ** 2
        )
        return 100 * self.hubble() * E

    def age(self, z):
        return sp.quad(lambda x: 1 / (1 + x) * self.rate(x), z, np.inf)
