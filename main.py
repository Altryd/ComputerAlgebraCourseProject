import numpy as np
import random
import math


def factorization(n):
    answer = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            count = 1
            n //= d
            while n % d == 0:
                # answer.append(d)
                count += 1
                n //= d
            answer.append({"value": d, "count": count})
        else:
            d += 1
    if n > 1:
        answer.append({"value": n, "count": 1})
    return answer


class CPoint:
    def __init__(self, x, y):
        self.__x = x
        self.__y = y

    @property
    def x(self):
        return self.__x

    @property
    def y(self):
        return self.__y

    @x.setter
    def x(self, new_x):
        self.__x = new_x

    @y.setter
    def y(self, new_y):
        self.__y = new_y

    def __eq__(self, rhs):  # ==
        if type(rhs) is type(self):
            return self.__dict__ == rhs.__dict__
        else:
            return False

    def __str__(self):
        return f'({self.x} ; {self.y})'

    def multipy_by_two(self, p, A):
        # p - модуль Fp
        try:
            lamd = (3 * pow(self.__x, 2) + A) * pow(2 * self.__y, -1, p)
            lamd %= p
            x3 = pow(lamd, 2) - 2 * self.__x
            x3 %= p
            y3 = lamd * (self.__x - x3) - self.y
            y3 %= p
            return CPoint(x3, y3)
        except ValueError:
            return CPoint(float("inf"), float("inf"))
        except TypeError:
            return CPoint(float("inf"), float("inf"))

    def add_point(self, rhs, p, A):
        if self == rhs:
            return self.multipy_by_two(p, A)
        elif type(rhs) is not type(self):
            raise ValueError("Types for add_points are incompatible !")
        elif self.__x == float("inf") and rhs.__x != float("inf"):
            return rhs
        elif self.__x != float("inf") and rhs.__x == float("inf"):
            return self
        elif self.__x == float("inf") and rhs.__x == float("inf"):
            return CPoint(float("inf"), float("inf"))
        try:
            lamd = (rhs.y - self.y) * pow(rhs.x - self.x, -1, p)
            lamd %= p
            x3 = pow(lamd, 2) - self.__x - rhs.__x
            x3 %= p
            y3 = lamd * (self.__x - x3) - self.y
            y3 %= p
            return CPoint(x3, y3)
        except ValueError:
            return CPoint(float("inf"), float("inf"))

    def fast_multiply(self, n, p, A):
        """
        :param n: степень в которую нужно возвести
        :param p: точка
        :param A: параметр А кривой
        :return: возвращает точку возведенную в степнь
        """
        if n == 1:
            return self
        elif n == 0:
            return CPoint(float("inf"), float("inf"))
        elif n <= -1:
            positive_res = self.fast_multiply(-n, p, A)
            return CPoint(positive_res.x, -positive_res.y)
        temp = CPoint(self.x, self.y)
        Q = CPoint(float("inf"), float("inf"))
        try:
            n = int(n)
            bin_n = bin(n)
            bin_n = bin_n[2:]
            bin_n = bin_n[::-1]
            for i in range(0, len(bin_n)):
                if bin_n[i] == '1':
                    Q = Q.add_point(temp, p, A)
                temp = temp.multipy_by_two(p, A)
            return Q
        except ValueError:
            return CPoint(float("inf"), float("inf"))


def get_order_of_point(point: CPoint, p: int, A: int) -> int:
    # 1.
    Q = point.fast_multiply(p + 1, p, A)
    m = int(np.floor(np.power(p, 1 / 4)) + 1)

    # 2. list of jP
    l1 = []
    for j in range(-m, m):
        # temp = j*point
        temp = point.fast_multiply(j, p, A)
        l1.append((temp, j))

    m_2_point = point.fast_multiply(2 * m, p, A)
    # 3. list of Q+k*(2mP)
    l2 = []
    for k in range(-m, m + 1):
        rhs = m_2_point.fast_multiply(k, p, A)
        temp = Q.add_point(rhs, p, A)
        l2.append((temp, k * 2 * m))

    l1 = sorted(l1, key=lambda elem: (elem[0].x, elem[0].y))
    l2 = sorted(l2, key=lambda elem: (elem[0].x, elem[0].y))
    i = 0
    j = 0
    S = []
    # 4. intersect
    while i < len(l1) and j < len(l2):
        if l1[i][0].x <= l2[j][0].x:
            if l1[i][0] == l2[j][0]:
                S.append({"elem": l1[i][0], "j": l1[i][1], "k2m": l2[j][1]})
                break
            i += 1
            while i < len(l1) - 1 and l1[i] == l1[i - 1]:
                i += 1
        else:
            j += 1
            while j < len(l2) - 1 and l2[j] == l2[j - 1]:
                j += 1

    M = p + 1 + S[0]["k2m"] - S[0]["j"]
    # factorization
    factorization_res = factorization(M)
    while True:
        divide_happened = False
        for divider in factorization_res:
            suspicious_order = 1
            if divider['count'] > 1:
                suspicious_order *= (divider['value'] ** (divider['count'] - 1))
            other_dividers = list(filter(lambda x: (x != divider and x['count'] > 0), factorization_res))
            for oth_divider in other_dividers:
                suspicious_order *= (oth_divider['value'] ** (oth_divider['count']))
            temp_P = point.fast_multiply(suspicious_order, p, A)
            if temp_P.x == float("inf"):
                divide_happened = True
                M = suspicious_order
                divider['count'] -= 1
                if divider['count'] == 0:
                    factorization_res.remove(divider)

        if not divide_happened:
            break

    return M


def is_quadratic_residue(a, p):  # является ли квадратичным вычетом. p - модуль
    while a < 0:
        a += p
    if pow(a, int((p - 1) / 2), p) == 1:
        return True
    else:
        return False


def solve_quadratic_eq(n, p):  # алгоритм тоннели-шенкса
    # x ** 2 == n (mod p)
    # from https://ru.wikipedia.org/wiki/%D0%90%D0%BB%D0%B3%D0%BE%D1%80%D0%B8%D1%82%D0%BC_%D0%A2%D0%BE%D0%BD%D0%B5%D0%BB%D0%BB%D0%B8_%E2%80%94_%D0%A8%D0%B5%D0%BD%D0%BA%D1%81%D0%B0
    assert is_quadratic_residue(n, p)
    s = 0
    temp_p = p
    temp_p -= 1
    while temp_p % 2 == 0:
        temp_p /= 2
        s += 1
    q = int(temp_p)
    if s == 1:
        y1 = pow(n, int((p + 1) / 4), p)
        return y1, p - y1
    quadratic_residue = None
    for i in range(2, p):
        if not is_quadratic_residue(i, p):
            quadratic_residue = i
            break
    c = pow(quadratic_residue, q, p)
    R = pow(n, int((q + 1) / 2), p)
    t = pow(n, q, p)
    M = s
    while True:
        while t < 0:
            t += p
        if t % p == 1:
            return R, p - R
        needed_i = None
        for i in range(0, M):
            if pow(t, 2 ** i, p) == 1:
                needed_i = i
                break
        b = pow(c, 2 ** (M - i - 1), p)
        R = R * b % p
        b_square = pow(b, 2, p)
        t = t * b_square % p
        c = b_square
        M = needed_i


def get_point_on_line(a, b, p):
    for x in range(0, p):
        f_x = pow(x, 3, p) + a * x + b
        f_x %= p
        if not is_quadratic_residue(f_x, p):
            continue
        y1, y2 = solve_quadratic_eq(f_x, p)
        return x, y1


def get_random_point_on_line(a, b, p):
    # y ** 2 = x ** 3 + a*x + b (mod p)   - подставляем x и решаем квадратич.сравнение
    # list_x = [i for i in range(0, p)]
    # random.shuffle(list_x)
    while True:  # TODO : верхние две строчки рабочие, но ооочень долгие
        x = random.randint(0, p)
        f_x = pow(x, 3, p) + a * x + b
        f_x %= p
        if not is_quadratic_residue(f_x, p):
            continue
        y1, y2 = solve_quadratic_eq(f_x, p)
        return CPoint(x, random.choice([y1, y2]))


def get_order_of_group(a, b, p):
    # y ** 2 = x ** 3 + a*x + b (mod p)
    P1 = get_random_point_on_line(a, b, p)
    P2 = get_random_point_on_line(a, b, p)
    P3 = get_random_point_on_line(a, b, p)
    while P2 == P1:
        P2 = get_random_point_on_line(a, b, p)
    P_list = [P1, P2, P3]
    m_list = [get_order_of_point(P, p, a) for P in P_list]
    lcm = math.lcm(*m_list)  # НОК
    left_border = int(p + 1 - 2 * np.sqrt(p))
    right_border = int(p + 1 + 2 * np.sqrt(p))
    length = int(right_border - left_border + 1)
    p_range = np.linspace(left_border, right_border, length, dtype="int64")  # p+1-2sqrt(p); p+1+2sqrt(p)
    maybe_be_an_answer = p_range[p_range % lcm == 0]
    if len(maybe_be_an_answer) == 1:
        return maybe_be_an_answer[0]
    else:
        return None


def func_(P: CPoint, Q: CPoint, Ri: CPoint, ai, bi, p, A, order_of_point):
    if Ri.x % 3 == 0:
        return Ri.add_point(Q, p, A), ai % order_of_point, (bi + 1) % order_of_point
    elif Ri.x % 3 == 1:
        return Ri.multipy_by_two(p, A), (2 * ai) % order_of_point, (2 * bi) % order_of_point
    else:
        return Ri.add_point(P, p, A), (ai + 1) % order_of_point, bi % order_of_point


def ro_algorithm_pollard(p, A, B, max_iteration, P: CPoint, Q: CPoint):
    order_of_point = int(get_order_of_point(P, p, A))
    print(f"Порядок точки ро-алгоритма: {order_of_point}")
    a_i = random.randint(0, order_of_point)
    b_i = random.randint(0, order_of_point)
    a_2i = a_i
    b_2i = b_i
    first_half = P.fast_multiply(a_i, p, A)
    second_half = Q.fast_multiply(b_i, p, A)
    Xi = first_half.add_point(second_half, p, A)
    X2i = Xi
    # print(f"Изначальные значения: P={P}, Q={Q}, p={p}, N={order_of_point}, A={A}")
    # print(f"Нулевая итерация: ai={a_i}, b_i={b_i}, Xi={Xi}, a2i={a_2i}, b2i={b_2i}, X2i={X2i}")
    for i in range(0, max_iteration):
        Xi, a_i, b_i = func_(P, Q, Xi, a_i, b_i, p, A, order_of_point)

        X2i, a_2i, b_2i = func_(P, Q, X2i, a_2i, b_2i, p, A, order_of_point)
        pre_X2i = X2i
        X2i, a_2i, b_2i = func_(P, Q, X2i, a_2i, b_2i, p, A, order_of_point)

        # print(f"{i}-ая итерация: ai={a_i}, b_i={b_i}, Xi={Xi}, a2i={a_2i}, b2i={b_2i}, X2i={X2i};  pre_X2i={pre_X2i}")
        if Xi == X2i:
            u = a_i - a_2i
            v = b_2i - b_i
            d = math.gcd(v, order_of_point)
            if d == 1:
                result = u * pow(v, -1, order_of_point)
                result = result % order_of_point
                while result < 0:
                    result += order_of_point
                return result % order_of_point
            else:
                if u % d != 0:
                    return None
                # u/d = x * v/d (mod N/d)
                if math.gcd(int(v/d), int(order_of_point/d)) != 1:
                    return None
                x_first = u/d * pow(int(v/d), -1, int(order_of_point/d))
                if x_first < 0:
                    x_first += order_of_point
                x_first %= order_of_point
                x_list = [x_first]
                for k in range(1, d):
                    x_temp = x_first + k * order_of_point / d
                    x_list.append(x_temp)
                for x in x_list:
                    if P.fast_multiply(x, p, A) == Q:
                        return int(x) % order_of_point
                return x_list

    return None


if __name__ == "__main__":
    random.seed(11)
    p = 59
    a = 11
    b = 12

    print(f"Порядок кривой: {get_order_of_group(a, b, p)}")

    P = CPoint(7, 14)
    n = 25
    Q = P.fast_multiply(n, p, a)
    pollard_result = ro_algorithm_pollard(p, a, b, 10000, P, Q)
    print(f"pollard_result={pollard_result}")
    for i in range(0, p):
        P = CPoint(7, 14)
        order_of_point = get_order_of_point(P, p, a)
        n = i
        Q = P.fast_multiply(n, p, a)
        pollard_result = ro_algorithm_pollard(p, a, b, 10000, P, Q)
        # print(f"Настоящее n={n}, а полученное={pollard_result}")
        if pollard_result >= order_of_point:
            pollard_result %= order_of_point
        if n >= order_of_point:
            n %= order_of_point
        if pollard_result != n:
            print("ВСЕ ПЛОХО!")
            print(f"Настоящее n={n}, а полученное={pollard_result}")
