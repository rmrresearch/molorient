from decimal import Decimal, getcontext


def pi_as_decimal():
    """
    Compute pi to current precision.
    """

    getcontext().prec += 2
    lasts, t, s, n, na, d, da = 0, Decimal(3), 3, 1, 0, 0, 24
    while s != lasts:
        lasts = s
        n, na = n+na, na+8
        d, da = d+da, da+32
        t = (t * n) / d
        s += t
    getcontext().prec -= 2
    return +s        


def arcsin_series(z):
    """
    Maclaurin series of arcsin.
    """

    getcontext().prec += 2

    if z == Decimal('1'):
        return pi_as_decimal() / 2 
    
    elif z == Decimal('-1'):
        return -pi_as_decimal() / 2

    if abs(z) < Decimal('0.9'):
        i, lasts, s, num, coeff = Decimal('0'), Decimal('0'), Decimal(z), Decimal(z), Decimal('1')
        tol = Decimal(10) ** -(Decimal(getcontext().prec) -2)
        while True:
            lasts = s
            i += 1
            num *= z * z
            coeff *= ((2*i - 1)**2 / ((2*i) * (2*i + 1)))
            s += num * coeff
            if abs(s - lasts) < tol:
                break
        result = +s
    
    elif z > 0:
        y = (1 - z**2).sqrt()
        result = (pi_as_decimal() / 2) - arcsin_series(y)
    
    else:
        y = (1 - z**2).sqrt()
        result = arcsin_series(y) - (pi_as_decimal() / 2)
    
    getcontext().prec -= 2
    return result


def arccos_series(z):
    """
    Maclauring series of arccos around 0 with the equation
    arccos(z) = pi/2 - arcsin(z).
    """

    if z == Decimal('1'):
        return Decimal('0')
    
    elif z == Decimal('-1'):
        return pi_as_decimal() 

    getcontext().prec += 2

    #Pi with 30 significant figures
    pi = pi_as_decimal()

    result = pi / 2 - arcsin_series(z)
    getcontext().prec -= 2
    return +result


def cos_series(y):
    """
    Maclaurin series of cos.
    """

    getcontext().prec += 2
    i, lasts, s, fact, num, sign = Decimal('0'), Decimal('0'), Decimal('1'), Decimal('1'), Decimal('1'), Decimal('1')
    tol = Decimal(10) ** -(Decimal(getcontext().prec) -2)
    while True:
        lasts = s
        i += 2
        fact *= i * (i-1)
        num *= y * y
        sign *= -1
        s += num / fact * sign
        if abs(s - lasts) < tol:
            break

    getcontext().prec -= 2
    return +s


def sin_series(x):
    """
    Maclaurin series of sin.
    """

    getcontext().prec += 2

    tol = Decimal(10) ** -(Decimal(getcontext().prec) -2)
    i, lasts, s, fact, num, sign = 1, 0, x, 1, x, 1
    while s != lasts:
        lasts = s
        i += 2
        fact *= i * (i-1)
        num *= x * x
        sign *= -1
        s += num / fact * sign
        if abs(s - lasts) < tol:
            break

    getcontext().prec -= 2

    return +s


def arctan_series(x):
    """
    Maclaurin series of arctan.
    """

    if x == Decimal('1'):
        return (pi_as_decimal() / 4)
    
    if x == Decimal('-1'):
        return (- pi_as_decimal() / 4)
    
    getcontext().prec += 2

    arg = x / (1 + x**2).sqrt()

    result = arcsin_series(arg)
    getcontext().prec -= 2

    return result


def arctan2(y, x):
    """
    Arctan Maclaurin series with two arguments. Returns angle from
    -pi to pi instead of -pi/2 to pi/2.
    """

    if x > 0:
        return arctan_series(y / x)
    
    elif x < 0 and y >= 0:
        return arctan_series(y / x) + pi_as_decimal()
    
    elif x < 0 and y < 0:
        return arctan_series(y / x) - pi_as_decimal()
    
    elif x == 0 and y > 0:
        return pi_as_decimal() / 2
    
    elif x == 0 and y < 0:
        return -pi_as_decimal() / 2