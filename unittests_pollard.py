import unittest
import random
from main import CPoint, get_order_of_point, get_order_of_group, ro_algorithm_pollard

class TestPollard(unittest.TestCase):
    def test_pollard_first(self):
        random.seed(11)
        p = 59
        a = 11
        b = 12
        points = [(0, 22), (0, 37), (5, 29), (5, 30), (7, 14), (7, 45), (8, 9), (8, 50), (10, 1), (10, 58), (11, 15),
                  (11,44), (13,13), (13, 46), (14, 14), (14, 45), (15 , 22 ), (15 , 37 ), (16 , 6 ), (16 , 53 ),
                  (19 , 0 ), (21 , 8 ), (21 , 51 ), (22 , 20 ), (22 , 39 ), (25 , 10 ), (25 , 49 ), (28 , 18 ),
                  (28 , 41 ), (30 , 5 ), (30 , 54 ), (33 , 26 ), (33 , 33 ), (35 , 5 ), (35 , 54 ), (36 , 10 ),
                  (36 , 49 ), (38 , 14 ), (38 , 45 ), (41 , 0 ), (42 , 24 ), (42 , 35 ), (44 , 22 ), (44 , 37 ),
                  (45 , 8 ), (45 , 51 ), (48 , 25 ), (48 , 34 ), (52 , 8 ), (52 , 51 ), (53 , 5 ), (53 , 54 ),
                  (54 , 3 ), (54 , 56 ), (55 , 9 ), (55 , 50 ), (57 , 10 ), (57 , 49 ), (58 , 0 )]
        all_point_count = len(points)
        point_count = 0
        for point in points:
            point_count += 1
            print(f"\n\nПриступаем к точке point={point}")
            for i in range(0, p):
                P = CPoint(*point)
                order_of_point = get_order_of_point(P, p, a)
                n = i
                Q = P.fast_multiply(n, p, a)
                pollard_result = ro_algorithm_pollard(p, a, b, 10000, P, Q)
                actual_Q = P.fast_multiply(pollard_result, p, a)
                print(f"Настоящее n={n}, а полученное={pollard_result} (порядок точки={order_of_point}); Q={Q}, Q_полученное={Q}")
                if pollard_result >= order_of_point:
                    pollard_result %= order_of_point
                if n >= order_of_point:
                    n %= order_of_point
                self.assertIsNot(pollard_result, list)
                self.assertEqual(n, pollard_result)
                self.assertEqual(Q, P.fast_multiply(pollard_result, p, a))
                #if pollard_result != n:
                #    print("ВСЕ ПЛОХО!")
                #    print(f"Настоящее n={n}, а полученное={pollard_result}")
        print(f"Просмотрено точек: {point_count} / {all_point_count}")
        self.assertEqual(point_count, all_point_count)

    def test_fast_multiply(self):
        P = CPoint(24, 22)
        n = 103
        A = 7
        actual = P.fast_multiply(n, 59, A)
        expected = CPoint(24, 37)  # by sagemath
        self.assertEqual(actual, expected)

    def test_pollard_from_linux(self):
        random.seed(11)
        p = 59
        A = 11
        B = 12

        print(f"Порядок кривой: {get_order_of_group(A, B, p)}")

        P = CPoint(7, 14)
        n = 33  # 29, 27 good
        Q = P.fast_multiply(n, p, A)
        print(P)
        my_n = ro_algorithm_pollard(p, A, B, 100, P, Q)
        print(f"Полученный n от моей реализации: {my_n}")
        order_of_point = get_order_of_point(P, p, A)
        if n > order_of_point:
            n %= order_of_point
        if my_n > order_of_point:
            my_n %= order_of_point
        self.assertEqual(my_n, n)

    def test_get_order_of_group(self):
        p = 29
        A = 1
        B = 12
        actual_order = get_order_of_group(A, B, p)
        expected_order = 23  # from sage
        self.assertEqual(23, int(actual_order))


    def test_get_order_of_point(self):
        p = 29
        A = 1
        B = 12
        P1 = CPoint(3, 19)
        P2 = CPoint(11, 7)
        p1_order = get_order_of_point(P1, p, A)
        p2_order = get_order_of_point(P2, p, A)
        expected_order = 23  # from sage
        self.assertEqual(expected_order, int(p1_order))
        self.assertEqual(expected_order, int(p2_order))






if __name__ == "__main__":
    unittest.main()