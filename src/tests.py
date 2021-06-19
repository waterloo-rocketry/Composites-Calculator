import unittest
import json
import numpy as np
from stack import Stack


class TestCompositesCalculator(unittest.TestCase):
    stress_tolerance = 10**2
    Q_tolerance = 10**-1
    strain_tolerance = 10**-5

    def test_example_1(self):
        with open('data/stacks/example1.json') as f:
            stack_data = json.load(f)
        stack = Stack(stack_data)
        Q_0 = np.array([[155.7, 3.02, 0], [3.02, 12.16, 0], [0, 0, 4.40]]) * 10 ** 9

        # assert that the mid-plane is 4 x thickness / 2, or 2 x thickness
        self.assertEqual(stack.midplane, 0.0003)

        # Q_bar of the sheet oriented at 0 should be same as Q_0
        np.testing.assert_allclose(stack.layers[0].Q_bar, Q_0, atol=self.Q_tolerance)
        # Q_bar_90 was taken from the Ch.9 solutions on the google drive
        np.testing.assert_allclose(stack.layers[1].Q_bar,
                                   np.array([[12.16, 3.02, 0], [3.02, 155.7, 0], [0, 0, 4.4]]) * 10 ** 9, atol=self.Q_tolerance)

        ABD = stack.get_ABD()
        A = ABD[0:3, 0:3]
        B = ABD[3:6, 0:3]
        D = ABD[3:6, 3:6]
        np.testing.assert_allclose(B, np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]]), atol=10 ** -7)
        np.testing.assert_allclose(A, np.array([[50.358, 1.812, 0], [1.812, 50.358, 0], [0, 0, 2.64]]) * 10 ** 6,
                                   atol=10**-7)

        # evaluate strains in a global reference frame
        stack.get_strains_and_stresses()

        global_estrain = np.array([10.01 ** -3, 0, 0])
        # since k_strain is expected to be 0, the actual plystrain is the estrain
        global_plystrain = global_estrain
        # in the notes, this was 0. I understand it is overall insignificant, but I believe we should carry it through,
        # or decide on better rounding
        global_plystress_0 = np.array([1.558557 * 10 ** 8, 3.023 * 10 ** 6, 0])
        global_plystress_1 = np.array([1.217216 * 10 ** 7, 3.023 * 10 ** 6, 0])

        np.testing.assert_allclose(stack.estrain_glob, global_estrain, atol=self.strain_tolerance)
        # moments are zero, and laminate is symmetric, so we expect this to be 0
        np.testing.assert_allclose(stack.kstrain_glob, np.array([0, 0, 0]), atol=self.strain_tolerance)

        # check its correct for one of the 0 oriented and one of the 90 oriented
        np.testing.assert_allclose(stack.layers[0].global_ply_strain, global_plystrain, atol=self.strain_tolerance)
        np.testing.assert_allclose(stack.layers[1].global_ply_strain, global_plystrain, atol=self.strain_tolerance)

        np.testing.assert_allclose(stack.layers[0].global_ply_stress, global_plystress_0, atol=self.stress_tolerance)
        np.testing.assert_allclose(stack.layers[1].global_ply_stress, global_plystress_1, atol=self.stress_tolerance)

        np.testing.assert_allclose(stack.layers[0].local_ply_strain, global_plystrain, atol=self.strain_tolerance)
        np.testing.assert_allclose(stack.layers[1].local_ply_strain, np.array([0, 10.01 ** -3, 0]), atol=self.strain_tolerance)
        np.testing.assert_allclose(stack.layers[0].local_ply_stress, global_plystress_0, atol=self.stress_tolerance)
        np.testing.assert_allclose(stack.layers[1].local_ply_stress, np.array([3.023 * 10 ** 6, 1.217216 * 10 ** 7, 0]),
                                   atol=self.stress_tolerance)

    def test_example_2(self):
        with open('data/stacks/example2.json') as f:
            stack_data = json.load(f)
        stack = Stack(stack_data)
        Q_0 = np.array([[155.7, 3.02, 0], [3.02, 12.16, 0], [0, 0, 4.40]]) * 10 ** 9

        ABD = stack.get_ABD()
        A = ABD[0:3, 0:3]
        B = ABD[3:6, 0:3]
        D = ABD[3:6, 3:6]
        print(B)


if __name__ == '__main__':
    unittest.main()
