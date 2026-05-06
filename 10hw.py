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


def poly_to_bits(poly, length):
    return [(poly >> index) & 1 for index in range(length)]


def build_generator_matrix(generator_poly, n, k):
    generator_bits = poly_to_bits(generator_poly, n)
    matrix = []

    for shift in range(k):
        row = [0] * shift + generator_bits[:n - shift]
        matrix.append(row)

    return matrix


def rref_gf2(matrix):
    """
    Приведение матрицы к ступенчатому виду над GF(2).
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

    return matrix, pivot_columns


def nullspace_gf2(matrix):
    """
    Находит базис пространства решений системы matrix * x^T = 0 над GF(2).
    """
    rref_matrix, pivot_columns = rref_gf2(matrix)
    column_count = len(matrix[0])

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
            vector[pivot_column] = rref_matrix[row][free_column]

        basis.append(vector)

    return basis


def transpose(matrix):
    return [list(row) for row in zip(*matrix)]


def matrix_mul_gf2(left, right):
    """
    Умножение матриц над GF(2).
    """
    row_count = len(left)
    inner_count = len(left[0])
    column_count = len(right[0])

    result = [
        [0 for _ in range(column_count)]
        for _ in range(row_count)
    ]

    for i in range(row_count):
        for j in range(column_count):
            value = 0

            for k in range(inner_count):
                value ^= left[i][k] & right[k][j]

            result[i][j] = value

    return result


def is_zero_matrix(matrix):
    return all(
        value == 0
        for row in matrix
        for value in row
    )


def print_matrix(matrix):
    for row in matrix:
        print(" ".join(str(value) for value in row))


def main():
    n = 15

    # M1(x) = x^4 + x + 1
    M1 = 0b10011

    # M3(x) = x^4 + x^3 + x^2 + x + 1
    M3 = 0b11111

    generator_poly = poly_mul(M1, M3)
    k = n - poly_degree(generator_poly)

    print("=" * 80)
    print("Построение BCH-кода по M1(x) и M3(x)")
    print("=" * 80)

    print()
    print("M1(x) =", poly_to_string(M1))
    print("M3(x) =", poly_to_string(M3))
    print("g(x)  =", poly_to_string(generator_poly))
    print("deg g =", poly_degree(generator_poly))
    print("n =", n)
    print("k =", k)

    G = build_generator_matrix(generator_poly, n, k)

    print()
    print("Порождающая матрица G:")
    print_matrix(G)

    H = nullspace_gf2(G)

    print()
    print("Проверочная матрица H:")
    print_matrix(H)

    product = matrix_mul_gf2(G, transpose(H))

    print()
    print("Произведение G * H^T:")
    print_matrix(product)

    print()
    print("G * H^T = 0:", is_zero_matrix(product))


if __name__ == "__main__":
    main()
