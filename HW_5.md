# Задача 1

## Условие

Проверьте, что код, дуальный коду Хэмминга, действительно является симплекс-кодом.  
Для доказательства этого факта напишите программу, которая будет вычислять все необходимые величины.

---

## 1. Код Хэмминга \((7,4)\)

Рассмотрим двоичный код Хэмминга с параметрами

\[
(7,4,3).
\]

Для него можно взять проверочную матрицу

\[
H =
\begin{pmatrix}
0 & 1 & 1 & 1 & 1 & 0 & 0 \\
1 & 0 & 1 & 1 & 0 & 1 & 0 \\
1 & 1 & 0 & 1 & 0 & 0 & 1
\end{pmatrix}.
\]

Порождающая матрица кода Хэмминга может быть выбрана так:

\[
G =
\begin{pmatrix}
1 & 0 & 0 & 0 & 0 & 1 & 1 \\
0 & 1 & 0 & 0 & 1 & 0 & 1 \\
0 & 0 & 1 & 0 & 1 & 1 & 0 \\
0 & 0 & 0 & 1 & 1 & 1 & 1
\end{pmatrix}.
\]

Матрицы \(G\) и \(H\) должны удовлетворять условию

\[
GH^T = 0.
\]

---

## 2. Дуальный код

Дуальный код \(C^\perp\) к коду \(C\) определяется как

\[
C^\perp =
\{x \in \mathbb{F}_2^7 : xc^T = 0 \text{ для всех } c \in C\}.
\]

Если \(H\) — проверочная матрица кода Хэмминга, то строки матрицы \(H\)
порождают дуальный код.

Поэтому порождающая матрица дуального кода:

\[
G_{C^\perp} = H.
\]

Так как матрица \(H\) имеет 3 строки и 7 столбцов, дуальный код имеет параметры

\[
n = 7,\qquad k = 3.
\]

Количество кодовых слов:

\[
M = 2^k = 2^3 = 8.
\]

---

## 3. Симплекс-код

Двоичный симплекс-код имеет параметры

\[
(2^m - 1,\; m,\; 2^{m-1}).
\]

При \(m = 3\):

\[
(2^3 - 1,\; 3,\; 2^{3-1})
=
(7,3,4).
\]

Значит, чтобы проверить, что дуальный код к коду Хэмминга является симплекс-кодом, нужно получить параметры

\[
(7,3,4).
\]

То есть нужно проверить:

1. длина кода равна \(7\);
2. размерность равна \(3\);
3. количество кодовых слов равно \(8\);
4. все ненулевые кодовые слова имеют вес \(4\);
5. минимальное расстояние равно \(4\).

---

## 4. Программа проверки

```python
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


# Порождающая матрица кода Хэмминга (7,4)
G_hamming = [
    [1, 0, 0, 0, 0, 1, 1],
    [0, 1, 0, 0, 1, 0, 1],
    [0, 0, 1, 0, 1, 1, 0],
    [0, 0, 0, 1, 1, 1, 1],
]

# Проверочная матрица кода Хэмминга
# Она же порождающая матрица дуального кода
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
```

---

## 5. Результат работы программы

Программа выводит:

```text
Проверка G * H^T:
[0, 0, 0]
[0, 0, 0]
[0, 0, 0]
[0, 0, 0]

Ортогональность выполнена: True

Кодовые слова дуального кода:
u = 000 -> c = 0000000, вес = 0
u = 001 -> c = 1101001, вес = 4
u = 010 -> c = 1011010, вес = 4
u = 011 -> c = 0110011, вес = 4
u = 100 -> c = 0111100, вес = 4
u = 101 -> c = 1010101, вес = 4
u = 110 -> c = 1100110, вес = 4
u = 111 -> c = 0001111, вес = 4

Параметры дуального кода:
n = 7
k = 3
M = 8
d_min = 4

Распределение весов:
A_0 = 1
A_4 = 7

Распределение расстояний между различными кодовыми словами:
расстояние 4: 28 пар

Дуальный код к коду Хэмминга является симплекс-кодом: True
```

---

## 6. Вывод

Дуальный код к коду Хэмминга имеет параметры

\[
(7,3,4).
\]

У него:

\[
n = 7,
\qquad
k = 3,
\qquad
M = 8,
\qquad
d_{\min} = 4.
\]

Все ненулевые кодовые слова имеют вес \(4\):

\[
A_0 = 1,\qquad A_4 = 7.
\]

Также расстояние между любыми двумя различными кодовыми словами равно \(4\).

Следовательно, дуальный код к коду Хэмминга действительно является двоичным симплекс-кодом.
