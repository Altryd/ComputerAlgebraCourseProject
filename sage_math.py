#!/usr/bin/env sage

from sage.all import *
import random


def get_order_of_point(point, p: int) -> int:
    # 1.
    Q = (p+1) * P
    m = int(floor(pow(p, 0.25)) + 1)

    # 2. list of jP
    l1 = []
    for j in range(-m, m):  # конструировать hashmap
        temp = j*point
        l1.append((temp, j))

    m_2_point = (2*m) * point
    # 3. list of Q+k*(2mP)
    l2 = []
    for k in range(-m, m + 1):
        rhs = k * m_2_point
        temp = Q + rhs
        l2.append((temp, k * 2 * m))  # а здесь сразу искать в хэшмапе

    l1 = sorted(l1, key=lambda elem: (elem[0][0], elem[0][1]))
    l2 = sorted(l2, key=lambda elem: (elem[0][0], elem[0][1]))
    i = 0
    j = 0
    S = []
    # 4. intersect
    while i < len(l1) and j < len(l2):
        if l1[i][0][0] <= l2[j][0][0]:  # if l1[i][0].x <= l2[j][0].x:
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
    current_value = M
    while True:
        divide_happened = False
        factorization_res = factor(current_value)
        factorization_res.sort()
        factorization_res = list(factorization_res)
        for divider in factorization_res:
            temp_value = current_value / divider[0]
            temp_P = temp_value * point
            if temp_P == E(0, 1, 0):
                divide_happened = True
                current_value = temp_value
        if not divide_happened:
            break
    return current_value


def func_(P, Q, Ri, ai, bi, order_of_point):
    # здесь нужно делать другую функцию с какой-то нормальной
    # селекторной функцией (или хотя бы разбивать на больше частей)
    if int(Ri[0]) % 3 == 0:
        return Ri + Q, ai % order_of_point, (bi + 1) % order_of_point
    elif int(Ri[0]) % 3 == 1:
        return 2 * Ri, (2 * ai) % order_of_point, (2 * bi) % order_of_point
    else:
        return Ri + P, (ai + 1) % order_of_point, bi % order_of_point


def ro_algorithm_pollard(p, A, B, max_iteration, P, Q):
    order_of_point = int(P.order())
    # order_of_point = int(get_order_of_group(A, B, p))
    print(f"Порядок точки ро-алгоритма через sage: {order_of_point}")
    my_order = int(get_order_of_point(P, p))
    print(f"Порядок точки ро-алгоритма через мой алгоритм: {my_order}")
    if int(my_order) != int(order_of_point):
        print("[ALERT!] порядки не сошлись!")
    a_i = random.randint(0, order_of_point)
    b_i = random.randint(0, order_of_point)
    a_2i = a_i
    b_2i = b_i
    Xi = a_i * P + b_i * Q
    X2i = Xi
    # print(f"Изначальные значения: P={P}, Q={Q}, p={p}, N={order_of_point}, A={A}")
    # print(f"Нулевая итерация: ai={a_i}, b_i={b_i}, Xi={Xi}, a2i={a_2i}, b2i={b_2i}, X2i={X2i}")
    for i in range(0, max_iteration):
        Xi, a_i, b_i = func_(P, Q, Xi, a_i, b_i, order_of_point)

        X2i, a_2i, b_2i = func_(P, Q, X2i, a_2i, b_2i, order_of_point)
        pre_X2i = X2i
        X2i, a_2i, b_2i = func_(P, Q, X2i, a_2i, b_2i, order_of_point)

        if i % (10 ** 6) == 0:
            print(f"{i}-ая итерация: ai={a_i}, b_i={b_i}, Xi={Xi}, a2i={a_2i}, b2i={b_2i}, X2i={X2i};  pre_X2i={pre_X2i}")
        if Xi == X2i:
            u = a_i - a_2i
            v = b_2i - b_i
            d = gcd(v, order_of_point)
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
                if gcd(int(v/d), int(order_of_point/d)) != 1:
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
                    if x * P == Q:
                        return int(x) % order_of_point
                return x_list

    return None



"""    
    another example:
    p = 1208925819614629174706189
    A = 3
    B = 42

    E = EllipticCurve(GF(p), [A, B])
    P = E(1179072205265749114237101 , 945874690885348578458594 )
    Q = E(192293917757484076566111 , 887425469491771128820518 )
    # sage_n = P.discrete_log(Q)
    # print(sage_n)
    my_n = ro_algorithm_pollard(p, A, B, 10000000000000000000000000000000000000, P, Q)
    print(my_n)
    if int(sage_n) != my_n:
        print("bad")
    else:
        print("good")
"""
if __name__ == "__main__":
    random.seed(55)
    order_curve = E.order()
    print(f"Порядок кривой: {order_curve}")
    all_good = True
    bad_situations = []
    all_points = E.points()
    len_all_points = len(all_points)
    for i in range(0, p):
        P = E[random.randint(1, len_all_points-1)]
        n = i  # 29, 27 good
        Q = n * P
        print(f"\n\nЗаходим с n={i} и точкой P={P}, Q={Q}")
        # print(P)
        my_n = ro_algorithm_pollard(p, A, B, 10000, P, Q)
        print(f"Полученный n от моей реализации: {my_n}")
        sage_n = P.discrete_log(Q)
        print(f"Полученный n через sage: {sage_n}")

        if int(sage_n) !=  my_n:
            bad_situations.append(f"не совпало, у sage n={sage_n}, у меня n={my_n}")
            all_good = False
    print("all good" if all_good else "bad..")
    if not all_good:
        for situation in bad_situations:
            print(situation)
