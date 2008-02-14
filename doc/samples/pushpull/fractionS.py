import math


def fraction_S( E1, E2, K ):

    if E1 == E2:
        return 0.5

    num = E1 - E2 - ( E1 + E2 ) * K + \
        math.sqrt( (E1 - E2) ** 2 + 2 * K * ( E1 - E2 ) ** 2 + \
                       ( ( E1 + E2 ) * K ) ** 2 )
    den = 2 * ( E1 - E2 )

    return 1.0 - num / den



