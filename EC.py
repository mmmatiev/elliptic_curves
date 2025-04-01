POINT_INFINITY = None


def mod_min_abs(x, p):
    """
    Приводит число x к промежутку [-p//2 .. p//2] по модулю p.
    """
    r = x % p
    if r > p // 2:
        r -= p
    return r


def mod_exp(a, k, n):
    """
    Возведение a в степень k по модулю n с использованием алгоритма быстрого возведения в степень.
    """
    b = 1
    if k == 0:
        return b

    A = a % n

    if (k & 1) == 1:
        b = A

    k >>= 1
    while k > 0:
        A = (A * A) % n
        if (k & 1) == 1:
            b = (b * A) % n
        k >>= 1

    return b


def mod_inv_euler(x, p):
    """
    Вычисляет обратный элемент x^(-1) по модулю p через x^(p-2) mod p.
    """
    if x % p == 0:
        raise ValueError("Нет обратного элемента для 0 по модулю p")
    return mod_exp(x, p - 2, p)


def mod_sqrts(value, p):
    """
    Ищет все y, такие что y^2 = value (mod p). Возвращает список корней.
    """
    roots = []
    value_norm = mod_min_abs(value, p)
    for raw_y in range(p):
        y = mod_min_abs(raw_y, p)
        if mod_min_abs(y * y, p) == value_norm:
            roots.append(y)
    return list(set(roots))


def find_points_on_curve(p, a, b):
    """
    Находит все конечные точки на кривой y^2 = x^3 + a*x + b (mod p).
    """
    points = []
    for raw_x in range(p):
        x_ = mod_min_abs(raw_x, p)
        rhs = x_ ** 3 + a * x_ + b
        rhs = mod_min_abs(rhs, p)
        y_candidates = mod_sqrts(rhs, p)
        for y_ in y_candidates:
            points.append((x_, y_))
    return points


def point_add(P, Q, p, a):
    """
    Сложение двух точек P и Q на эллиптической кривой.
    """
    if P is None:
        return Q
    if Q is None:
        return P

    x1, y1 = P
    x2, y2 = Q

    if x1 == x2 and mod_min_abs(y1 + y2, p) == 0:
        return POINT_INFINITY

    if x1 != x2:
        lam_num = mod_min_abs(y2 - y1, p)
        lam_den = mod_min_abs(x2 - x1, p)
        lam = mod_min_abs(lam_num * mod_inv_euler(lam_den, p), p)
    else:
        lam_num = mod_min_abs(3 * x1 * x1 + a, p)
        lam_den = mod_min_abs(2 * y1, p)
        lam = mod_min_abs(lam_num * mod_inv_euler(lam_den, p), p)

    x3 = lam * lam
    x3 = mod_min_abs(x3, p)
    x3 = mod_min_abs(x3 - x1, p)
    x3 = mod_min_abs(x3 - x2, p)

    y3 = x1 - x3
    y3 = mod_min_abs(y3, p)
    y3 = mod_min_abs(lam * y3, p)
    y3 = mod_min_abs(y3 - y1, p)

    return (x3, y3)


def point_mul(k, P, p, a):
    """
    Умножение точки на число k.
    """
    if P is None or k == 0:
        return POINT_INFINITY
    result = POINT_INFINITY
    base = P
    while k > 0:
        # Проверка на нечётность
        if k % 2 == 1:
            result = point_add(result, base, p, a)
        base = point_add(base, base, p, a)
        k //= 2  # Делим k на 2
    return result



def point_order(P, p, a):
    """
    Находит порядок точки P.
    """
    if P is None:
        return 1
    order = 1
    current = P
    while current is not None:
        current = point_add(current, P, p, a)
        order += 1
    return order


def generate_subgroup(P, p, a):
    """
    Генерация всех точек на основе одного генератора.
    """
    subgroup = []
    current = P
    while current is not None:
        subgroup.append(current)
        current = point_add(current, P, p, a)
    subgroup.append(POINT_INFINITY)
    return subgroup


def is_prime(n):
    """ Простая проверка на простоту. """
    if n < 2:
        return False
    if n == 2 or n == 3:
        return True
    if n % 2 == 0:
        return False
    r = int(n ** 0.5)
    for i in range(3, r + 1, 2):
        if n % i == 0:
            return False
    return True


if __name__ == "__main__":
    p = 4999
    a = -2
    b = 1

    # 1) Строим все точки на кривой
    curve_points = find_points_on_curve(p, a, b)
    all_points = [POINT_INFINITY] + curve_points

    print(f"Все точки на кривой: {curve_points}")
    print(f"Порядок группы: {len(all_points)}\n")

    # 2) Умножение точки
    P = (122, 2200)
    k = 4
    kP = point_mul(k, P, p, a)
    print(f"{k} * {P} = {kP}")

    P = (122, 2200)
    k = 1000
    kP = point_mul(k, P, p, a)
    print(f"{k} * {P} = {kP}")

    # 3) Проверка подгрупп простого порядка
    prime_subgroups = []
    visited = set()
    for pt in all_points:
        if pt is not None:
            ord_pt = point_order(pt, p, a)
            if is_prime(ord_pt):
                H = generate_subgroup(pt, p, a)
                Hf = frozenset(H)
                if Hf not in visited:
                    visited.add(Hf)
                    prime_subgroups.append((pt, H))

    print("\nПодгруппы простого порядка:")
    for generator_pt, subgroup_pts in prime_subgroups:
        print(f"Точка-генератор: {generator_pt}, порядок = {len(subgroup_pts)}")
        print(f"Подгруппа: {subgroup_pts}\n")
