PRIMITIVE_POLY = 0b10011
ALPHA = 0b0010


# ============================================================
# Арифметика поля GF(2^4)
# ============================================================

def gf_add(a, b):
    return a ^ b


def gf_mul(a, b):
    result = 0

    while b:
        if b & 1:
            result ^= a

        b >>= 1
        a <<= 1

        if a & 0b10000:
            a ^= PRIMITIVE_POLY

    return result & 0b1111


def gf_pow(a, power):
    result = 1

    for _ in range(power):
        result = gf_mul(result, a)

    return result


def alpha_power(power):
    return gf_pow(ALPHA, power % 15)


def gf_to_alpha(value):
    if value == 0:
        return "0"

    for power in range(15):
        if alpha_power(power) == value:
            if power == 0:
                return "1"
            if power == 1:
                return "alpha"
            return f"alpha^{power}"

    return str(value)


# ============================================================
# Многочлены над GF(2^4)
# ============================================================

def poly_mul(a, b):
    result = [0] * (len(a) + len(b) - 1)

    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            result[i + j] ^= gf_mul(ai, bj)

    return trim_poly(result)


def trim_poly(poly):
    poly = poly[:]

    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()

    return poly


def poly_eval(poly, value):
    result = 0
    power = 1

    for coefficient in poly:
        result ^= gf_mul(coefficient, power)
        power = gf_mul(power, value)

    return result


def poly_to_string(poly):
    terms = []

    for power in range(len(poly) - 1, -1, -1):
        coefficient = poly[power]

        if coefficient == 0:
            continue

        coefficient_string = gf_to_alpha(coefficient)

        if power == 0:
            terms.append(coefficient_string)
        elif power == 1:
            if coefficient == 1:
                terms.append("x")
            else:
                terms.append(f"{coefficient_string}*x")
        else:
            if coefficient == 1:
                terms.append(f"x^{power}")
            else:
                terms.append(f"{coefficient_string}*x^{power}")

    if not terms:
        return "0"

    return " + ".join(terms)


def build_rs_generator(first_root, root_count):
    """
    Строит многочлен:

        product (x + alpha^i)

    для i = first_root, ..., first_root + root_count - 1.

    В характеристике 2 это то же самое, что product (x - alpha^i).
    """
    generator = [1]

    for power in range(first_root, first_root + root_count):
        root = alpha_power(power)
        generator = poly_mul(generator, [root, 1])

    return generator


# ============================================================
# Задачи
# ============================================================

def solve_task_1():
    print("=" * 80)
    print("ЗАДАЧА 1")
    print("=" * 80)

    n = 15
    d = 5
    root_count = d - 1

    generator = build_rs_generator(1, root_count)
    k = n - (len(generator) - 1)

    print()
    print("Порождающий многочлен для RS-кода с d = 5:")
    print("g(x) =", poly_to_string(generator))
    print("deg g =", len(generator) - 1)
    print("n =", n)
    print("k =", k)
    print(f"Параметры: [{n}, {k}, {d}]")


def solve_task_2():
    print()
    print("=" * 80)
    print("ЗАДАЧА 2")
    print("=" * 80)

    n = 15
    k = 9
    t = 3
    d = 2 * t + 1

    generator = build_rs_generator(1, d - 1)

    print()
    print("Порождающий многочлен для RS-кода (15,9):")
    print("g2(x) =", poly_to_string(generator))
    print("deg g2 =", len(generator) - 1)
    print(f"Параметры: [{n}, {k}, {d}]")

    # v(x) =
    # alpha
    # + alpha^3 x
    # + alpha^10 x^2
    # + alpha^12 x^3
    # + alpha^13 x^4
    # + alpha^3 x^5
    # + alpha^9 x^6
    # + alpha^5 x^7
    # + x^8
    # + alpha^10 x^9
    # + alpha^8 x^10
    # + alpha x^11
    # + alpha x^12
    # + alpha^13 x^13
    # + alpha^2 x^14

    received_exponents = [
        1, 3, 10, 12, 13,
        3, 9, 5, 0, 10,
        8, 1, 1, 13, 2
    ]

    received = [
        alpha_power(exponent)
        for exponent in received_exponents
    ]

    print()
    print("Принятое слово v(x):")

    for index, coefficient in enumerate(received):
        print(f"коэффициент при x^{index}: {gf_to_alpha(coefficient)}")

    syndromes = [
        poly_eval(received, alpha_power(power))
        for power in range(1, 7)
    ]

    print()
    print("Синдромы:")

    for index, syndrome in enumerate(syndromes, start=1):
        print(f"S{index} = {gf_to_alpha(syndrome)}")

    has_errors = any(syndrome != 0 for syndrome in syndromes)

    print()
    print("Есть ошибки:", has_errors)

    if not has_errors:
        error_vector = [0] * n

        print()
        print("Все синдромы равны нулю.")
        print("Следовательно, принятое слово уже является кодовым.")
        print("Вектор ошибок:")

        print(error_vector)


def main():
    solve_task_1()
    solve_task_2()


if __name__ == "__main__":
    main()
