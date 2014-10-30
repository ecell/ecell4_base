import math


def fraction_Sp(E1, E2, K):

    if E1 == E2:
        return 0.5

    num = E1 - E2 - (E1 + E2) * K + \
        math.sqrt((E1 - E2) ** 2 + 2 * K * (E1 - E2) ** 2 + \
                      ((E1 + E2) * K) ** 2)
    den = 2 * (E1 - E2)

    return  num / den

def fraction_S(E1, E2, K):

    return 1.0 - fraction_Sp(E1, E2, K)



