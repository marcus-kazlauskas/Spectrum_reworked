import numpy as np
import math as mt


def _cos_theta(theta, n) -> float:
    return mt.sqrt(1 - mt.pow(mt.sin(theta / n), 2))


def _beta(theta, n, lambda_const=1, lambda_var=1) -> float:
    return mt.pi * _cos_theta(theta, n) * lambda_const / (2 * lambda_var)


def _p(theta, n) -> float:
    return n * _cos_theta(theta, n)


def _q(theta, n) -> float:
    return _cos_theta(theta, n) / n


class Matrix(object):
    """Abeles matrix method for counting number of layers for dielectric mirror"""
    def __init__(self):
        self.m_te = np.array([[1, 0],
                              [0, 1]], complex)  # Abeles matrix for transverse-electric(TE) or s-polarized ray
        self.m_tm = np.array([[1, 0],
                              [0, 1]], complex)  # Abeles matrix for transverse-magnetic(TM) or p-polarized ray

    def set(self, theta, n, lambda_const=1, lambda_var=1):
        assert abs(theta) < mt.pi / 2 and n > 0 and lambda_const > 0 and lambda_var > 0, 'Wrong parameters'
        # set Abeles matrix of layer for TE-polarization
        self.m_te[0][0] = complex(mt.cos(_beta(theta, n, lambda_const, lambda_var)), 0)
        self.m_te[0][1] = complex(0, - mt.sin(_beta(theta, n, lambda_const, lambda_var)) / _p(theta, n))
        self.m_te[1][0] = complex(0, - mt.sin(_beta(theta, n, lambda_const, lambda_var)) * _p(theta, n))
        self.m_te[1][1] = complex(mt.cos(_beta(theta, n, lambda_const, lambda_var)), 0)
        # set Abeles matrix of layer for TM-polarization
        self.m_tm[0][0] = complex(mt.cos(_beta(theta, n, lambda_const, lambda_var)), 0)
        self.m_tm[0][1] = complex(0, - mt.sin(_beta(theta, n, lambda_const, lambda_var)) / _q(theta, n))
        self.m_tm[1][0] = complex(0, - mt.sin(_beta(theta, n, lambda_const, lambda_var)) * _q(theta, n))
        self.m_tm[1][1] = complex(mt.cos(_beta(theta, n, lambda_const, lambda_var)), 0)

    def unit(self):
        self.m_te = np.array([[1, 0],
                              [0, 1]], complex)  # unit matrix for transverse-electric(TE) or s-polarized ray
        self.m_tm = np.array([[1, 0],
                              [0, 1]], complex)  # unit matrix for transverse-magnetic(TM) or p-polarized ray

    def multiply(self, matrix):
        assert isinstance(matrix, Matrix), 'Input object isn\'t a Matrix instance'
        self.m_te = self.m_te @ matrix.m_te
        self.m_tm = self.m_tm @ matrix.m_tm

    def t_te(self, theta: float, n_air: float, n_ground: float) -> complex:
        assert abs(theta) < mt.pi / 2 and n_air > 0 and n_ground > 0, 'Wrong parameters'
        # transmittance coefficient for TE-polarization
        return 2 * _p(theta, n_air) / (self.m_te[0][0] * _p(theta, n_air)
                                       + self.m_te[1][1] * _p(theta, n_ground)
                                       + self.m_te[0][1] * _p(theta, n_ground) * _p(theta, n_air)
                                       + self.m_te[1][0])

    def t_tm(self, theta: float, n_air: float, n_ground: float) -> complex:
        assert abs(theta) < mt.pi / 2 and n_air > 0 and n_ground > 0, 'Wrong parameter'
        # transmittance coefficient for TM-polarization
        return 2 * _q(theta, n_air) / (self.m_tm[0][0] * _q(theta, n_air)
                                       + self.m_tm[1][1] * _q(theta, n_ground)
                                       + self.m_tm[0][1] * _q(theta, n_ground) * _q(theta, n_air)
                                       + self.m_tm[1][0])


class Mirror(object):
    """Abeles matrix for dielectric mirror"""
    def __init__(self, theta, n_ground, n_high, n_low):
        assert abs(theta) < mt.pi / 2 and n_ground > 0 and n_high >= n_low > 0, 'Wrong parameters'
        self.theta: float = theta  # angle between ray and perpendicular to mirror's surface
        self.n_air: float = 1  # refractive index of air
        self.n_ground: float = n_ground  # refractive index of ground layer
        self.n_high: float = n_high  # refractive index of layer with high optical density
        self.n_low: float = n_low  # refractive index of layer with low optical density
        self.t_te: complex = 0  # transmittance amplitude for TE-polarization
        self.t_tm: complex = 0  # transmittance amplitude for TM-polarization
        self.tt_te: float = 0  # transmittance intensity for TE-polarization
        self.tt_tm: float = 0  # transmittance intensity for TM-polarization
        self.high = Matrix()
        self.low = Matrix()
        self.mirror = Matrix()

    def set_layers(self, lambda_const=1, lambda_var=1):
        self.high.unit()
        self.low.unit()
        self.mirror.unit()
        self.high.set(self.theta, self.n_high, lambda_const, lambda_var)
        self.low.set(self.theta, self.n_low, lambda_const, lambda_var)

    def transmittance(self):
        self.t_te = self.mirror.t_te(self.theta, self.n_air, self.n_ground)
        self.tt_te = mt.pow(abs(self.t_te), 2) * _p(self.theta, self.n_ground) / _p(self.theta, self.n_air)
        self.t_tm = self.mirror.t_tm(self.theta, self.n_air, self.n_ground)
        self.tt_tm = mt.pow(abs(self.t_tm), 2) * _q(self.theta, self.n_ground) / _q(self.theta, self.n_air)

    def first_layer(self):
        self.mirror.multiply(self.high)

    def next_layers(self):
        self.mirror.multiply(self.low)
        self.mirror.multiply(self.high)
