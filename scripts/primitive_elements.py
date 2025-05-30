def is_primitive(g, p):
    # p is prime, so order of multiplicative group is p-1
    order = p - 1
    # Find all prime divisors of order
    divisors = set()
    n = order
    d = 2
    while d * d <= n:
        if n % d == 0:
            divisors.add(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        divisors.add(n)
    # Check if g^((p-1)//q) != 1 for all prime divisors q
    for q in divisors:
        if pow(g, order // q, p) == 1:
            return False
    assert pow(g, order, p) == 1
    return True

def find_primitive_elements(p, count=10):
    found = []
    for g in range(2, p):
        if is_primitive(g, p):
            found.append(g)
            if len(found) >= count:
                break
    return found

if __name__ == "__main__":
    import sys
    p = int(sys.argv[1])
    primitives = find_primitive_elements(p, count=10)
    print("First 10 primitive elements:")
    print(primitives)
