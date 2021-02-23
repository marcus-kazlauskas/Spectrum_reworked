import unittest
import math
import numpy as np
from abeles import Matrix
from abeles import Mirror


class TestAbelesMatrix(unittest.TestCase):
    def setUp(self) -> None:
        self.matrix = Matrix()
        self.matrix.set(math.pi / 6,
                        1.457,
                        2.4,
                        1.457)

        self.mirror = Mirror(math.pi / 6,
                             1.457,
                             2.4,
                             1.457)

        self.mirror0 = Mirror(0, 1, 1, 1)

        self.number = 10
        self.precision = 3
        self.lambda_const = 632
        self.lambda_min = 400
        self.lambda_max = 900
        self.step = 10

    def test_det(self):
        self.assertEqual(1.0, np.round(np.linalg.det(self.matrix.m_te), self.precision), 'Det(m_te) should be 1')
        self.assertEqual(1.0, np.round(np.linalg.det(self.matrix.m_tm), self.precision), 'Det(m_tm) should be 1')

    def test_void(self):
        self.mirror0.set_layers()
        self.mirror0.first_layer()
        for i in range(self.number):
            self.mirror0.next_layers()
        self.mirror0.transmittance()
        self.assertEqual(1.0, np.round(self.mirror0.tt_te, self.precision), 'Air should be transparent')
        self.assertEqual(1.0, np.round(self.mirror0.tt_tm, self.precision), 'Air should be transparent')

    def test_coefficients(self):
        lambda_var = self.lambda_min
        while lambda_var <= self.lambda_max:
            self.mirror.set_layers(self.lambda_const, lambda_var)
            self.mirror.first_layer()
            for i in range(self.number):
                self.mirror.next_layers()
            self.mirror.transmittance()
            self.assertTrue(1.0 > self.mirror.tt_te > 0.0, 'Impossible transmittance')
            self.assertTrue(1.0 > self.mirror.tt_tm > 0.0, 'Impossible transmittance')
            lambda_var += self.step


if __name__ == '__main__':
    unittest.main()
