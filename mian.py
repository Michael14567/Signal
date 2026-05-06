from collections import Counter
from itertools import product


# ============================================================
# Общие функции
# ============================================================

def bits_to_string(bits):
    return ''.join(str(bit) for bit in bits)


def hamming_weight(bits):
    return sum(bits)


def gf2_vector_matrix_mul(vector, matrix):
    """
    Умножение вектора-строки на матрицу над GF(2).
    """
    result = []

    for column in range(len(matrix[0])):
        value = 0

        for row in range(len(vector)):
            value ^= vector[row] & matrix[row][column]

        result.append(value)

    return tuple(result)


def gf2_rank(matrix):
    """
    Ранг матрицы над GF(2).
    """
    matrix = [row[:] for row in matrix]

    if not matrix:
        return 0

    row_count = len(matrix)
    column_count = len(matrix[0])

    rank = 0

    for column in range(column_count):
        pivot = None

        for row in range(rank, row_count):
            if matrix[row][column] == 1:
                pivot = row
                break

        if pivot is None:
            continue

        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]

        for row in range(row_count):
            if row != rank and matrix[row][column] == 1:
                matrix[row] = [
                    a ^ b
                    for a, b in zip(matrix[row], matrix[rank])
                ]

        rank += 1

    return rank


def nullspace_gf2(matrix):
    """
    Находит базис решений системы matrix * x^T = 0 над GF(2).
    Возвращает строки проверочной матрицы.
    """
    matrix = [row[:] for row in matrix]
    row_count = len(matrix)
    column_count = len(matrix[0])

    pivot_columns = []
    pivot_row = 0

    for column in range(column_count):
        pivot = None

        for row in range(pivot_row, row_count):
            if matrix[row][column] == 1:
                pivot = row
                break

        if pivot is None:
            continue

        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]

        for row in range(row_count):
            if row != pivot_row and matrix[row][column] == 1:
                matrix[row] = [
                    a ^ b
                    for a, b in zip(matrix[row], matrix[pivot_row])
                ]

        pivot_columns.append(column)
        pivot_row += 1

        if pivot_row == row_count:
            break

    free_columns = [
        column
        for column in range(column_count)
        if column not in pivot_columns
    ]

    basis = []

    for free_column in free_columns:
        vector = [0] * column_count
        vector[free_column] = 1

        for row, pivot_column in enumerate(pivot_columns):
            vector[pivot_column] = matrix[row][free_column]

        basis.append(vector)

    return basis


# ============================================================
# Задача 1. Полиномы над GF(2)
# ============================================================

def poly_degree(poly):
    if poly == 0:
        return -1
    return poly.bit_length() - 1


def poly_mul(a, b):
    """
    Умножение многочленов над GF(2).
    Бит i соответствует коэффициенту при x^i.
    """
    result = 0
    shift = 0

    while b:
        if b & 1:
            result ^= a << shift

        b >>= 1
        shift += 1

    return result


def poly_divmod(dividend, divisor):
    if divisor == 0:
        raise ZeroDivisionError("division by zero polynomial")

    quotient = 0
    remainder = dividend
    divisor_degree = poly_degree(divisor)

    while remainder != 0 and poly_degree(remainder) >= divisor_degree:
        shift = poly_degree(remainder) - divisor_degree
        quotient ^= 1 << shift
        remainder ^= divisor << shift

    return quotient, remainder


def poly_is_divisible(dividend, divisor):
    _, remainder = poly_divmod(dividend, divisor)
    return remainder == 0


def poly_to_string(poly):
    if poly == 0:
        return "0"

    parts = []

    for power in range(poly_degree(poly), -1, -1):
        if ((poly >> power) & 1) == 0:
            continue

        if power == 0:
            parts.append("1")
        elif power == 1:
            parts.append("x")
        else:
            parts.append(f"x^{power}")

    return " + ".join(parts)


def monic_polynomials_of_degree(degree):
    start = 1 << degree

    for tail in range(1 << degree):
        yield start | tail


def is_irreducible(poly):
    degree = poly_degree(poly)

    if degree <= 0:
        return False

    for divisor_degree in range(1, degree // 2 + 1):
        for divisor in monic_polynomials_of_degree(divisor_degree):
            if poly_is_divisible(poly, divisor):
                return False

    return True


def factor_over_gf2(poly):
    factors = []
    current = poly

    for divisor_degree in range(1, poly_degree(poly) + 1):
        for divisor in monic_polynomials_of_degree(divisor_degree):
            if not is_irreducible(divisor):
                continue

            while current != 1 and poly_is_divisible(current, divisor):
                factors.append(divisor)
                current, _ = poly_divmod(current, divisor)

    return factors


def int_to_codeword(value, length):
    return tuple((value >> i) & 1 for i in range(length))


def generate_cyclic_code(generator, n):
    r = poly_degree(generator)
    k = n - r

    codewords = []

    for message in range(2 ** k):
        codeword_poly = poly_mul(message, generator)
        codewords.append(int_to_codeword(codeword_poly, n))

    return sorted(set(codewords))


def minimum_distance_linear_code(codewords):
    nonzero_weights = [
        hamming_weight(word)
        for word in codewords
        if hamming_weight(word) != 0
    ]

    return min(nonzero_weights)


def classify_cyclic_code(n, k, d_min):
    if (n, k, d_min) == (7, 6, 2):
        return "код с проверкой на чётность"

    if (n, k, d_min) == (7, 4, 3):
        return "код Хэмминга"

    if (n, k, d_min) == (7, 3, 4):
        return "симплексный код"

    if (n, k, d_min) == (7, 1, 7):
        return "код с повторением"

    return "циклический код"


def solve_task_1():
    print("=" * 80)
    print("ЗАДАЧА 1")
    print("=" * 80)

    n = 7
    x7_minus_1 = (1 << 7) | 1

    print()
    print("Разложение x^7 - 1 над GF(2):")
    print("x^7 - 1 =", poly_to_string(x7_minus_1))

    factors = factor_over_gf2(x7_minus_1)

    print()
    print("Неприводимые множители:")

    for index, factor in enumerate(factors, start=1):
        print(f"f{index}(x) = {poly_to_string(factor)}")

    print()
    print("Проверка произведения множителей:")

    product_value = 1

    for factor in factors:
        product_value = poly_mul(product_value, factor)

    print("Произведение =", poly_to_string(product_value))
    print("Верно:", product_value == x7_minus_1)

    generators = []

    for mask in range(1, 2 ** len(factors) - 1):
        generator = 1

        for index, factor in enumerate(factors):
            if (mask >> index) & 1:
                generator = poly_mul(generator, factor)

        generators.append(generator)

    generators = sorted(generators, key=lambda p: (poly_degree(p), p))

    print()
    print("Нетривиальные порождающие многочлены и характеристики кодов:")

    for index, generator in enumerate(generators, start=1):
        codewords = generate_cyclic_code(generator, n)

        r = poly_degree(generator)
        k = n - r
        d_min = minimum_distance_linear_code(codewords)
        code_type = classify_cyclic_code(n, k, d_min)
        weights = Counter(hamming_weight(word) for word in codewords)

        print()
        print("-" * 80)
        print(f"Код {index}")
        print("g(x) =", poly_to_string(generator))
        print(f"deg g = {r}")
        print(f"Параметры: [{n}, {k}, {d_min}]")
        print("Тип:", code_type)

        print("Распределение весов:")
        for weight in sorted(weights):
            print(f"  A_{weight} = {weights[weight]}")


# ============================================================
# Задача 2. Решётки
# ============================================================

def generate_linear_code(generator_matrix):
    k = len(generator_matrix)
    messages = list(product([0, 1], repeat=k))
    codewords = [
        gf2_vector_matrix_mul(message, generator_matrix)
        for message in messages
    ]

    return messages, sorted(set(codewords))


def build_trellis_by_codewords(codewords):
    n = len(codewords[0])

    state_names = []
    prefix_to_state = []

    for level in range(n + 1):
        groups = {}
        prefix_map = {}

        prefixes = sorted(set(word[:level] for word in codewords))

        for prefix in prefixes:
            futures = tuple(
                sorted(
                    word[level:]
                    for word in codewords
                    if word[:level] == prefix
                )
            )

            if futures not in groups:
                groups[futures] = f"S{level}_{len(groups)}"

            prefix_map[prefix] = groups[futures]

        state_names.append(sorted(groups.values()))
        prefix_to_state.append(prefix_map)

    edges = []

    for level in range(n):
        edge_set = set()

        for word in codewords:
            source_prefix = word[:level]
            target_prefix = word[:level + 1]

            source = prefix_to_state[level][source_prefix]
            target = prefix_to_state[level + 1][target_prefix]
            label = word[level]

            edge_set.add((source, target, label))

        edges.append(sorted(edge_set))

    return state_names, edges


def row_span(row):
    nonzero = [
        index
        for index, value in enumerate(row)
        if value == 1
    ]

    return min(nonzero), max(nonzero)


def active_rows_after_position(generator_matrix, level):
    active = []

    for index, row in enumerate(generator_matrix):
        start, end = row_span(row)

        if start < level <= end:
            active.append(index)

    return active


def build_trellis_by_generator_matrix(generator_matrix):
    messages, codewords = generate_linear_code(generator_matrix)
    n = len(generator_matrix[0])

    states_by_level = []
    state_maps = []

    for level in range(n + 1):
        active = active_rows_after_position(generator_matrix, level)

        states = sorted(set(
            tuple(message[row] for row in active)
            for message in messages
        ))

        state_map = {
            state: f"S{level}_{index}"
            for index, state in enumerate(states)
        }

        states_by_level.append([
            f"{state_map[state]}={state}"
            for state in states
        ])
        state_maps.append((active, state_map))

    edges = []

    for level in range(n):
        edge_set = set()

        active_source, source_map = state_maps[level]
        active_target, target_map = state_maps[level + 1]

        for message, codeword in zip(messages, codewords):
            source_state = tuple(message[row] for row in active_source)
            target_state = tuple(message[row] for row in active_target)

            source = source_map[source_state]
            target = target_map[target_state]
            label = codeword[level]

            edge_set.add((source, target, label))

        edges.append(sorted(edge_set))

    return states_by_level, edges


def xor_vectors(a, b):
    return tuple(x ^ y for x, y in zip(a, b))


def multiply_bit_vector(bit, vector):
    if bit == 0:
        return tuple(0 for _ in vector)

    return tuple(vector)


def build_trellis_by_parity_check_matrix(H):
    r = len(H)
    n = len(H[0])

    columns = [
        tuple(H[row][column] for row in range(r))
        for column in range(n)
    ]

    zero_state = tuple(0 for _ in range(r))

    forward = [set() for _ in range(n + 1)]
    forward[0].add(zero_state)

    for level in range(n):
        for state in forward[level]:
            for bit in [0, 1]:
                next_state = xor_vectors(
                    state,
                    multiply_bit_vector(bit, columns[level])
                )
                forward[level + 1].add(next_state)

    backward = [set() for _ in range(n + 1)]
    backward[n].add(zero_state)

    for level in range(n - 1, -1, -1):
        for next_state in backward[level + 1]:
            for bit in [0, 1]:
                previous_state = xor_vectors(
                    next_state,
                    multiply_bit_vector(bit, columns[level])
                )
                backward[level].add(previous_state)

    valid_states = [
        sorted(forward[level] & backward[level])
        for level in range(n + 1)
    ]

    state_maps = []

    for level, states in enumerate(valid_states):
        state_map = {
            state: f"S{level}_{index}"
            for index, state in enumerate(states)
        }
        state_maps.append(state_map)

    states_by_level = []

    for level, states in enumerate(valid_states):
        states_by_level.append([
            f"{state_maps[level][state]}={state}"
            for state in states
        ])

    edges = []

    for level in range(n):
        edge_set = set()

        for state in valid_states[level]:
            for bit in [0, 1]:
                next_state = xor_vectors(
                    state,
                    multiply_bit_vector(bit, columns[level])
                )

                if next_state in state_maps[level + 1]:
                    source = state_maps[level][state]
                    target = state_maps[level + 1][next_state]
                    edge_set.add((source, target, bit))

        edges.append(sorted(edge_set))

    return states_by_level, edges


def print_codewords(codewords):
    for word in codewords:
        print(bits_to_string(word), "weight =", hamming_weight(word))


def print_trellis(title, states_by_level, edges_by_level):
    print()
    print(title)
    print("-" * len(title))

    for level, states in enumerate(states_by_level):
        print(f"Ярус {level}: {', '.join(states)}")

    print()
    print("Переходы:")

    for level, edges in enumerate(edges_by_level):
        print(f"Ярус {level} -> {level + 1}:")

        for source, target, label in edges:
            print(f"  {source} --{label}--> {target}")


def solve_task_2():
    print()
    print("=" * 80)
    print("ЗАДАЧА 2")
    print("=" * 80)

    G = [
        [1, 0, 1, 1, 0, 1],
        [0, 0, 1, 1, 1, 0],
        [0, 1, 0, 1, 1, 1],
    ]

    print()
    print("Порождающая матрица G:")

    for row in G:
        print(row)

    messages, codewords = generate_linear_code(G)

    n = len(G[0])
    k = gf2_rank(G)
    d_min = minimum_distance_linear_code(codewords)

    print()
    print(f"Параметры кода: [{n}, {k}, {d_min}]")
    print(f"Количество кодовых слов: {len(codewords)}")

    print()
    print("Кодовые слова:")
    print_codewords(codewords)

    H = nullspace_gf2(G)

    print()
    print("Проверочная матрица H:")

    for row in H:
        print(row)

    states_1, edges_1 = build_trellis_by_codewords(codewords)
    states_2, edges_2 = build_trellis_by_generator_matrix(G)
    states_3, edges_3 = build_trellis_by_parity_check_matrix(H)

    print_trellis("Решётка 1: по кодовым словам", states_1, edges_1)
    print_trellis("Решётка 2: по порождающей матрице", states_2, edges_2)
    print_trellis("Решётка 3: по проверочной матрице", states_3, edges_3)


def main():
    solve_task_1()
    solve_task_2()


if __name__ == "__main__":
    main()
