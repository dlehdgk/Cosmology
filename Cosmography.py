import numpy as np
import scipy.integrate as sp


class Cosmo:
    # generate initial parameters of h, mass density radiation density and lambda density
    def __init__(self, h=0.67, Om=0.3, Or=1e-5, Ol=0.63):
        self.__h = h
        self.__om = Om
        self.__or = Or
        self.__ol = Ol

    # function to output h
    def hubble(self):
        return self.__h

    # function to output mass density
    def mass(self):
        return self.__om

    # function to output radiation density
    def radiation(self):
        return self.__or

    # function to output lambda density
    def constant(self):
        return self.__ol

    # function to calculate curvature density
    def curvature(self):
        return 1 - self.mass() - self.radiation() - self.constant()

    # function to calculate the Hubble parameter H at a given redshift z.
    # If SI=True then returns value in s^-1 otherwise returns km/s/Mpc
    def rate(self, z, SI=False):
        E = np.sqrt(
            self.mass() * (1 + z) ** 3
            + self.radiation() * (1 + z) ** 4
            + self.constant()
            + self.curvature() * (1 + z) ** 2
        )
        H = 100 * self.hubble() * E
        if SI == True:
            return H / 3.0857e19
        else:
            return H

    # function to find the age of the universe at a given z default gives present age of the universe.
    # If years=True gives ages in years otherwise in seconds.
    def age(self, z=0, years=True):
        time = sp.quad(lambda x: 1 / ((1 + x) * self.rate(x, SI=True)), z, np.inf)
        if years == True:
            return time[0] / 3600 / 24 / 365.25
        else:
            return time[0]

    # defined sink function which depends on curvature
    def sink(self, x):
        if self.curvature() > 0:
            return np.sin(x)
        elif self.curvature() == 0:
            return x
        else:
            return np.sinh(x)

    # function to find the luminosity distance to an object at redshift z
    # If Mpc=True, then returns distance in Mpc else in meters
    def lum(self, z, Mpc=True):
        I = sp.quad(lambda x: 1 / self.rate(x, SI=True), 0, z)
        if self.curvature() == 0:
            dl = 3e8 * (1 + z) * I[0]
        else:
            dl = (
                3e8
                * (1 + z)
                * 1
                / np.sqrt(abs(self.curvature()))
                * self.sink(np.sqrt(abs(self.curvature())) * I[0])
            )
        if Mpc == True:
            return dl / 3.08568e22
        else:
            return dl

    # function to calculate the angular diameter distance to an object at redshift z
    def ang(self, z, Mpc=True):
        da = self.lum(z, Mpc) / ((1 + z) ** 2)
        return da
