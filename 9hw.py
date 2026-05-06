from collections import Counter
from math import gcd


# ============================================================
# Задача 1
# ============================================================

def element_name(power):
    if power == 0:
        return "1 = x^0"
    if power == 1:
        return "x"
    return f"x^{power}"


def order_of_power(group_order, power):
    return group_order // gcd(group_order, power)


def solve_task_1():
    print("=" * 80)
    print("ЗАДАЧА 1")
    print("=" * 80)

    q = 5 ** 2
    group_order = q - 1

    print()
    print(f"Поле GF(5^2) содержит {q} элементов.")
    print(f"Мультипликативная группа GF(25)* содержит {group_order} элемента.")
    print("Считаем, что x — примитивный элемент.")
    print()

    print("Порядки элементов:")
    print("0: порядок не определён")

    orders = []

    for power in range(group_order):
        order = order_of_power(group_order, power)
        orders.append(order)
        print(f"{element_name(power)}: порядок {order}")

    distribution = Counter(orders)

    print()
    print("Распределение по порядкам:")

    for order in sorted(distribution):
        print(f"порядок {order}: {distribution[order]} элементов")


# ============================================================
# Задача 2
# ============================================================

def next_lfsr_state(state):
    """
    Регистр для h(x) = 1 + x + x^4.

    Из h(x) получаем:
        x^4 = x + 1.

    Поэтому новый бит:
        s_{t+4} = s_{t+1} + s_t mod 2.
    """
    new_bit = state[0] ^ state[1]
    return state[1:] + (new_bit,)


def generate_lfsr_states(initial_state):
    states = [initial_state]
    current = initial_state

    while True:
        current = next_lfsr_state(current)
        states.append(current)

        if current == initial_state:
            break

    return states


def bits_to_string(bits):
    return ''.join(str(bit) for bit in bits)


def solve_task_2():
    print()
    print("=" * 80)
    print("ЗАДАЧА 2")
    print("=" * 80)

    initial_state = (0, 0, 1, 1)
    states = generate_lfsr_states(initial_state)

    print()
    print("Проверочный полином:")
    print("h(x) = 1 + x + x^4")
    print()
    print("Рекуррентное соотношение:")
    print("s_{t+4} = s_{t+1} + s_t mod 2")
    print()
    print("Начальное состояние:", bits_to_string(initial_state))
    print()

    print("Смена состояний регистра:")

    for step, state in enumerate(states):
        print(f"{step:2d}: {bits_to_string(state)}")

    period = len(states) - 1
    max_period = 2 ** 4 - 1

    print()
    print("Период регистра:", period)
    print("Максимально возможный период:", max_period)
    print("Код максимальной длины:", period == max_period)

    n = 15
    k = 4
    d_min = 2 ** (k - 1)

    print()
    print("Параметры кода:")
    print(f"n = {n}")
    print(f"k = {k}")
    print(f"d_min = {d_min}")
    print(f"Код: [{n}, {k}, {d_min}]")


def main():
    solve_task_1()
    solve_task_2()


if __name__ == "__main__":
    main()
