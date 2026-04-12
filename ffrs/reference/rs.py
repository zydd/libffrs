from ffrs.reference import P


def locator(GF, w, err_pos):
    poly = P(GF, [1])
    for pos in err_pos:
        x = -w ** pos
        poly = poly * P(GF, [1, x])
    return poly
