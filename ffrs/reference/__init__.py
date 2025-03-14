import itertools


class F:
    def __init__(self, p, x):
        self.p = p
        if (type(x) == F):
            self.x = x.x % p
        elif type(x) == int:
            self.x = x % p
        else:
            raise RuntimeError(f'F: unsupported type: {type(x)}')

    def __add__(self, rhs):
        assert type(self) == type(rhs), f'incompatible types: {type(self)} and {type(rhs)}'
        v = (self.x + rhs.x) % self.p
        return F(self.p, v)

    def __eq__(self, rhs):
        assert type(self) == type(rhs), f'incompatible types: {type(self)} and {type(rhs)}'
        return self.x == rhs.x

    def __sub__(self, rhs):
        assert type(self) == type(rhs), f'incompatible types: {type(self)} and {type(rhs)}'
        v = (self.x - rhs.x) % self.p
        return F(self.p, v)

    def __neg__(self):
        return F(self.p, self.p - self.x)

    def __mul__(self, rhs):
        assert type(self) == type(rhs), f'incompatible types: {type(self)} and {type(rhs)}'
        v = (self.x * rhs.x) % self.p
        return F(self.p, v)

    def _mul(self, rhs):
        return self.__mul__(rhs)

    def __floordiv__(self, rhs):
        assert type(self) == type(rhs), f'incompatible types: {type(self)} and {type(rhs)}'
        return self.__mul__(rhs.inv())

    def __repr__(self):
        return f'F{self.p}[{self.x}]'

    def __int__(self):
        return self.x

    def inv(self):
        if self.x == 0:
            raise ZeroDivisionError

        v, newv = 0, 1
        r, newr = self.p, self.x

        while newr != 0:
            quotient = r // newr
            v, newv = (newv, v - quotient * newv)
            r, newr = (newr, r - quotient * newr)

        assert not (r > 1), self

        if v < 0:
            v += self.p

        return F(self.p, v)

    def inc(self):
        v = (self.x + 1) % self.p
        return F(self.p, v)

    def characteristic(self):
        return self.p


class P:
    def __init__(self, N, x):
        self.N = N

        if type(x) == int:
            self.x = GF.poly_from_int(N.p, x)
        else:
            self.x = [N(a) for a in x]

        self._norm()

    def __add__(self, rhs):
        assert type(self) == type(rhs), f'incompatible types: {type(self)} and {type(rhs)}'

        v = (a + b for a, b in itertools.zip_longest(self.x, rhs.x, fillvalue=self.N(0)))
        return P(self.N, v)

    def __eq__(self, rhs):
        assert type(self) == type(rhs), f'incompatible types: {type(self)} and {type(rhs)}'
        return (len(self.x) == len(rhs.x) and
            all(a == b for a, b in zip(self.x, rhs.x)))

    def __sub__(self, rhs):
        assert type(self) == type(rhs), f'incompatible types: {type(self)} and {type(rhs)}'

        v = (a - b for a, b in itertools.zip_longest(self.x, rhs.x, fillvalue=self.N(0)))
        return P(self.N, v)

    def __neg__(self):
        return P(self.N, (-a for a in self.x))

    def __mul__(self, rhs):
        assert type(self) == type(rhs), f'incompatible types: {type(self)} and {type(rhs)}'
        v = [self.N(0)] * (len(self.x) + len(rhs.x) - 1)

        for i in range(len(self.x)):
            for j in range(len(rhs.x)):
                v[i + j] += self.x[i] * rhs.x[j]

        return P(self.N, v)

    def __divmod__(self, rhs):
        assert type(self) == type(rhs), f'incompatible types: {type(self)} and {type(rhs)}'
        q, r = self._ex_synth_div(self.x, rhs.x)
        return P(self.N, q), P(self.N, r)

    def __mod__(self, rhs):
        _, r = divmod(self, rhs)
        return r

    def __floordiv__(self, rhs):
        d, _ = divmod(self, rhs)
        return d

    def __iter__(self):
        for a in self.x:
            yield a

    def __hash__(self):
        hash = 0
        for i in self.x:
            # WARNING: won't work if coefficients can be >= 2**16
            hash += int(i) << 16
        return hash

    def deg(self):
        assert len(self.x) > 0
        return len(self.x) - 1

    def eval(self, x):
        v = self.N(0)
        x0 = self.N(1)
        for a in self.x:
            v += a * x0
            x0 *= x
        return v

    def __getitem__(self, i):
        return self.x[i]

    def __setitem__(self, i, v):
        self.x[i] = v

    def scale(self, n):
        n = self.N(n)
        return P(self.N, [a * n for a in self.x])

    def deriv(self):
        ret = P(self.N, self.x)

        for i in range(len(ret.x) - 1):
            ret.x[i] = self.N(0)
            for _ in range(i + 1):
                ret.x[i] += ret.x[i + 1]

        ret.x.pop(-1)
        while int(ret.x[-1]) == 0:
            ret.x.pop(-1)

        return ret

    def __repr__(self):
        def enc_exp(n):
            exp = '⁰¹²³⁴⁵⁶⁷⁸⁹'
            r = ''
            while n:
                r += exp[n % 10]
                n //= 10
            return r[::-1]

        if len(self.x) == 0:
            return '0'

        s = ''
        for i, a in enumerate(self.x[::-1]):
            a = int(a)
            if a == 0:
                continue

            s += ' '

            if a < 0:
                s += '-'
                if i > 0:
                    s += ' '
            elif i > 0:
                s += '+ '

            a = abs(a)
            exp = len(self.x) - i - 1

            if a != 1 or exp == 0:
                s += str(a)

            if exp > 0:
                s += '𝑥'

            if exp > 1:
                s += enc_exp(exp)

        return s[1:]

    @staticmethod
    def _ex_synth_div(dividend, divisor):
        out = list(dividend)
        normalizer = divisor[-1]
        # assert int(divisor[-1]) == 1

        for i in range(len(dividend) -1, len(divisor)-1 -1, -1):
            out[i] //= normalizer # work with non-monic polynomials

            coef = out[i]
            if int(coef) != 0:
                for j in range(1, len(divisor)): # skip the first coefficient of the divisor
                    out[i - j] -= divisor[len(divisor)-1- j] * coef

        sep = len(divisor) - 1
        return out[sep:], out[:sep] # quotient, remainder

    def _norm(self):
        tailing = len(self.x)
        while tailing > 0:
            if int(self.x[tailing-1]) != 0:
                break
            tailing -= 1

        self.x = self.x[:tailing]


class Fp(P):
    def __init__(self, GF, x):
        self.GF = GF

        if type(x) == int:
            super().__init__(lambda x: F(GF.p, x), GF.poly_from_int(GF.p, x))
            _, self.x = self._ex_synth_div(self.x, GF.poly.x)
        elif type(x) == P:
            self.x = x.x
            self.N = x.N
        elif type(x) == Fp:
            assert x.GF == self.GF
            self.x = x.x
            self.N = x.N
        else:
            raise RuntimeError(f'Fp: unexpected type: {type(x)}')

        self.x_i = GF.poly_to_int(GF.p, self)

    def __eq__(self, rhs):
        x = rhs if type(rhs) is int else rhs.x_i
        return self.x_i == x

    def __add__(self, rhs):
        return Fp(self.GF, super().__add__(rhs) % self.GF.poly)

    def __sub__(self, rhs):
        return Fp(self.GF, super().__sub__(rhs) % self.GF.poly)

    def __neg__(self):
        return Fp(self.GF, super().__neg__() % self.GF.poly)

    def __mul__(self, rhs):
        if int(self) == 0 or int(rhs) == 0:
            return Fp(self.GF, 0)

        return Fp(self.GF, self._exp(self.log() + rhs.log()))

    def _mul(self, rhs):
        return Fp(self.GF, super().__mul__(rhs) % self.GF.poly)

    def _exp(self, x):
        return self.GF.exp_table[int(x) % (len(self.GF.exp_table) - 1)]

    def exp(self):
        return self.GF.exp_table[self.x_i]

    def log(self):
        return self.GF.log_table[self.x_i]

    def pow(self, power):
        return Fp(self.GF, self._exp(self.log() * power))

    def inv(self):
        if self.x_i == 0:
            raise ZeroDivisionError

        return Fp(self.GF, self.GF.exp_table[(len(self.GF.exp_table) - 1) - self.log()])

    def __floordiv__(self, rhs):
        if int(rhs) == 0:
            raise ZeroDivisionError
        if int(self) == 0:
            return Fp(self.GF, 0)

        return Fp(self.GF, self._exp(self.log() + (len(self.GF.exp_table) - 1) - rhs.log()))

    def __repr__(self):
        a = (f'{int(x):x}' for x in reversed(self.x))
        x = ('' if self.GF.p < 16 else ' ').join(a)
        return f'F{self.GF.p}[{x}]'

    def __repr__(self):
        a = (f'{int(x):x}' for x in reversed(self.x))
        x = ('' if self.GF.p < 16 else ' ').join(a)
        return f'F{self.GF.p}[{x}]'

    def __int__(self):
        return self.x_i

    def __hash__(self):
        return self.x_i


class GF:
    def __init__(self, p, k=1, a=None, poly=None):
        if a is None:
            a, poly = next(GF.irr_polynomials(p, k, p))

        if poly is None:
            a, poly = next(GF.irr_polynomials(p, k, a))

        if type(poly) == list:
            poly = P(lambda x: F(p, x), poly)
        elif type(poly) == P:
            poly = poly
        elif type(poly) == int:
            poly = self.poly_from_int(p, poly)
            poly = P(lambda x: F(p, x), poly)
        else:
            raise RuntimeError(f'unexpected type: {type(poly)}')

        self.p = p # prime
        self.k = k # power
        self.poly = poly # irreducible polynomial
        self.a = self(a) # primitive element

        self._init_exp_log_tables()

    def gen(self, i):
        return Fp(self, self.exp_table[i % (self.p ** self.k - 1)])

    def __call__(self, x):
        return Fp(self, x) if self.k > 1 else F(self.p, x)

    def __repr__(self):
        return f'GF{self.p ** self.k}({self.a}, {self.poly_to_int(self.p, self.poly)})'

    def minimal_polynomials(self):
        for i in range(1, self.p ** self.k - 1):
            a = self.gen(i)
            for j in range(self.p ** self.k + 1, self.p ** (self.k + 1) - 1):
                m = P(self, map(int, GF.poly_from_int(self.p, j)))
                mod = m.eval(a)
                if int(mod) == 0:
                    yield i, m

    def _init_exp_log_tables(self):
        field_charac = self.p ** self.k
        exp = [0] * field_charac
        log = [0] * field_charac

        a = self.a
        x = self(1)
        for i in range(field_charac):
            exp[i] = int(x)
            log[int(x)] = i
            x = x._mul(a)

        self.exp_table = exp
        self.log_table = log

    @staticmethod
    def poly_from_int(base, x):
        poly = []
        while x:
            poly.append(F(base, x))
            x //= base
        return poly

    @staticmethod
    def poly_to_int(base, poly):
        ret = 0
        p = 1

        for a in poly:
            ret += int(a) * p
            p *= base
        return ret

    @staticmethod
    def is_primitive(p, k, poly, prim):
        seen = [False] * (p ** k)

        x = P(lambda x: F(p, x), [1])

        for _ in range(p ** k - 1):
            x = (x * prim) % poly

            n = GF.poly_to_int(p, x)
            if seen[n]:
                return False
            else:
                seen[n] = True

        return True

    @staticmethod
    def primitives(p, k, poly):
        PN = lambda n: P(lambda x: F(p, x), GF.poly_from_int(p, n))
        for i in range(1, p ** k):
            prim = PN(i)

            if GF.is_primitive(p, k, poly, prim):
                yield GF.poly_to_int(p, prim)

    @staticmethod
    def irr_polynomials(p, k, primitive=None):
        PN = lambda n: P(lambda x: F(p, x), GF.poly_from_int(p, n))
        field_charac = p ** k

        if primitive is None:
            for poly in range(field_charac, 2 * field_charac):
                poly = PN(poly)

                try:
                    yield next(GF.primitives(p, k, poly)), GF.poly_to_int(p, poly)
                except StopIteration:
                    pass
        else:
            if type(primitive) == int:
                primitive = PN(primitive)

            for poly in range(field_charac, 2 * field_charac):
                poly = PN(poly)

                if GF.is_primitive(p, k, poly, primitive):
                    yield GF.poly_to_int(p, primitive), GF.poly_to_int(p, poly)

    @staticmethod
    def primitive_roots_of_unity(GF):
        res = dict()
        field_elements = GF.p ** GF.k

        for i in range(2, field_elements):
            for j in range(2, field_elements):
                if GF(i).pow(j) == 1:
                    if j not in res:
                        res[j] = []
                    res[j].append(i)
                    break
        return res
