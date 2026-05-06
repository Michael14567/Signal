from itertools import product, combinations
from collections import Counter


def gf2_mat_vec_mul(u, G):
    """
    Умножение вектора-строки u на матрицу G над GF(2).
    """
    n = len(G[0])
    result = []

    for j in range(n):
        value = 0
        for i in range(len(u)):
            value ^= u[i] & G[i][j]
        result.append(value)

    return result


def gf2_mat_mul(A, B):
    """
    Умножение матриц A и B над GF(2).
    """
    rows = len(A)
    cols = len(B[0])
    inner = len(B)

    C = [[0 for _ in range(cols)] for _ in range(rows)]

    for i in range(rows):
        for j in range(cols):
            value = 0
            for k in range(inner):
                value ^= A[i][k] & B[k][j]
            C[i][j] = value

    return C


def transpose(A):
    """
    Транспонирование матрицы.
    """
    return [list(row) for row in zip(*A)]


def weight(v):
    """
    Вес Хэмминга.
    """
    return sum(v)


def hamming_distance(a, b):
    """
    Расстояние Хэмминга между двумя векторами.
    """
    return sum(x != y for x, y in zip(a, b))


def bits_to_str(v):
    return ''.join(map(str, v))


def main():
    # Порождающая матрица кода Хэмминга (7,4)
    G_hamming = [
        [1, 0, 0, 0, 0, 1, 1],
        [0, 1, 0, 0, 1, 0, 1],
        [0, 0, 1, 0, 1, 1, 0],
        [0, 0, 0, 1, 1, 1, 1],
    ]

    # Проверочная матрица кода Хэмминга.
    # Она же является порождающей матрицей дуального кода.
    H = [
        [0, 1, 1, 1, 1, 0, 0],
        [1, 0, 1, 1, 0, 1, 0],
        [1, 1, 0, 1, 0, 0, 1],
    ]

    # 1. Проверяем ортогональность G * H^T = 0
    GHt = gf2_mat_mul(G_hamming, transpose(H))

    print("Проверка G * H^T:")
    for row in GHt:
        print(row)

    is_orthogonal = all(all(x == 0 for x in row) for row in GHt)

    print()
    print("Ортогональность выполнена:", is_orthogonal)

    # 2. Строим все кодовые слова дуального кода
    messages = list(product([0, 1], repeat=3))
    dual_codewords = [gf2_mat_vec_mul(u, H) for u in messages]

    print()
    print("Кодовые слова дуального кода:")

    for u, c in zip(messages, dual_codewords):
        print(f"u = {bits_to_str(u)} -> c = {bits_to_str(c)}, вес = {weight(c)}")

    # 3. Считаем параметры кода
    n = len(dual_codewords[0])
    k = 3
    M = len(dual_codewords)

    weights = [weight(c) for c in dual_codewords]
    weight_distribution = Counter(weights)

    nonzero_weights = [w for w in weights if w != 0]
    d_min = min(nonzero_weights)

    print()
    print("Параметры дуального кода:")
    print("n =", n)
    print("k =", k)
    print("M =", M)
    print("d_min =", d_min)

    # 4. Распределение весов
    print()
    print("Распределение весов:")

    for w in sorted(weight_distribution):
        print(f"A_{w} = {weight_distribution[w]}")

    # 5. Расстояния между всеми парами кодовых слов
    distances = [
        hamming_distance(c1, c2)
        for c1, c2 in combinations(dual_codewords, 2)
    ]

    distance_distribution = Counter(distances)

    print()
    print("Распределение расстояний между различными кодовыми словами:")

    for d in sorted(distance_distribution):
        print(f"расстояние {d}: {distance_distribution[d]} пар")

    # 6. Проверяем, что код является симплекс-кодом
    is_simplex = (
        n == 7
        and k == 3
        and M == 8
        and d_min == 4
        and weight_distribution[0] == 1
        and weight_distribution[4] == 7
    )

    print()
    print("Дуальный код к коду Хэмминга является симплекс-кодом:", is_simplex)


if __name__ == "__main__":
    main()
