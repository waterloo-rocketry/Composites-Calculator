import unittest

import numpy as np
from src.stack import Stack
from src.layer import Layer


class TestCompositesCalculator(unittest.TestCase):

    def test_example_1(self):
        Q_0 = np.array([[155.7, 3.02, 0], [3.02, 12.16, 0], [0, 0, 4.40]])
        layers = [Layer(0, 0.00015, 0, 0, 0, 0, 0, Q_0), Layer(0.5, 0.00015, 0, 0, 0, 0, 0, Q_0),Layer(0.5, 0.00015, 0, 0, 0, 0, 0, Q_0),Layer(0, 0.00015, 0, 0, 0, 0, 0, Q_0)]
        stack = Stack(0.00030, layers, np.array([50400, 1809, 0]), np.array([0, 0, 0]))
        stack.set_Q_bars()

        #Q_bar of the sheet oriented at 0 should be same as Q_0
        self.assertTrue(np.array_equal(stack.layers[0].Q_bar, Q_0))
        print(stack.layers[0].Q_bar)
        print(stack.layers[1].Q_bar)




if __name__ == '__main__':
    unittest.main()
