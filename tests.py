import unittest

import numpy as np
from src.stack import Stack
from src.layer import Layer


class TestCompositesCalculator(unittest.TestCase):

    def test_example_1(self):
        Q_0 = np.array([[155.7, 3.02, 0], [3.02, 12.16, 0], [0, 0, 4.40]])*10**9
        layers = [Layer(0, 0.00015, 0, 0, 0, 0, 0, Q_0), Layer(0.5*np.pi, 0.00015, 0, 0, 0, 0, 0, Q_0),Layer(0.5*np.pi, 0.00015, 0, 0, 0, 0, 0, Q_0),Layer(0, 0.00015, 0, 0, 0, 0, 0, Q_0)]
        stack = Stack(layers, np.array([50400, 1809, 0]), np.array([0, 0, 0]))
        stack.process_layers()

        # assert that the midplane is 4xthickness/2, or 2xthickness
        self.assertEqual(stack.midplane, 0.0003)

        #Q_bar of the sheet oriented at 0 should be same as Q_0
        np.testing.assert_allclose(stack.layers[0].Q_bar, Q_0,atol=10**-1)
        #Q_bar_90 was taken from the Ch.9 solutions on the google drive
        np.testing.assert_allclose(stack.layers[1].Q_bar, np.array([[12.16,3.02,0],[3.02, 155.7,0],[0,0,4.4]])*10**9,atol=10**-1)

        ABD = stack.get_ABD()
        A = ABD[0:3, 0:3]
        B = ABD[3:6, 0:3]
        D = ABD[3:6, 3:6]
        # this is an unsatifying use of absolute tolerance,but the issue is that even though it is ultimately insignficant,
        # anything is inf bigger than 0
        np.testing.assert_allclose(B,np.array([[0,0,0],[0,0,0],[0,0,0]]),atol=10**6)
        #evaluate that A matrix is calculated correctly
        np.testing.assert_allclose(A,np.array([[50.4,1.81,0],[1.81,50.4,0],[0,0,2.64]])*10**6, rtol=10**-1, atol=1)

        stack.get_global()


if __name__ == '__main__':
    unittest.main()
