def _vec_sub(vec, i, j, factor):
    if type(vec[i]) is list:
        vec[i] = [x - y * factor for x, y in zip(vec[i], vec[j])]
    else:
        vec[i] -= vec[j] * factor


def _vec_mul(vec, i, v):
    if type(vec[i]) is list:
        vec[i] = [x * v for x in vec[i]]
    else:
        vec[i] *= v

def _clone_mat(mat):
    if type(mat) is list:
        return [_clone_mat(x) for x in mat]
    else:
        return mat

def gaussian_elim(mat, vec):
    assert len(mat) == len(mat[0]) == len(vec)

    N = len(mat)
    mat = _clone_mat(mat)
    vec = _clone_mat(vec)

    for k in range(len(mat)):
        row_max = k
        while mat[row_max][k] == 0 and row_max < N-1:
            row_max += 1

        if k != row_max:
            t = mat[k]
            mat[k] = mat[row_max]
            mat[row_max] = t

        factor = mat[k][k].inv()
        _vec_mul(vec, k, factor)
        for col in range(k, N):
            mat[k][col] *= factor

        for row in range(k + 1, N):
            factor = mat[row][k]
            _vec_sub(vec, row, k, factor)
            for col in range(k, N):
                mat[row][col] -= mat[k][col] * factor

        for row in range(k-1, -1, -1):
            factor = mat[row][k]
            _vec_sub(vec, row, k, factor)
            for col in range(k, N):
                mat[row][col] -= mat[k][col] * factor

    return vec


def identity(GF, n):
    return [[GF(int(i == j)) for i in range(n)] for j in range(n)]


def inverse(GF, mat):
    id = identity(GF, len(mat))
    return gaussian_elim(mat, id)


def transpose(mat):
    return [[mat[i][j] for i in range(len(mat))] for j in range(len(mat[0]))]

def matmul(GF, a, b):
    assert len(a[0]) == len(b)
    vector = type(b[0]) is not list
    if vector:
        b = transpose([b])

    res = [[GF(0) for _ in range(len(b[0]))] for _ in range(len(a))]
    for i in range(len(a)):
        for j in range(len(b[0])):
            for k in range(len(b)):
                res[i][j] += a[i][k] * b[k][j]

    return transpose(res)[0] if vector else res


def print_mat(mat):
    assert(len(mat) >= 2)

    print(f"[{[int(x) for x in mat[0]]},")
    for row in mat[1:-1]:
        print(f" {[int(x) for x in row]},")
    print(f" {[int(x) for x in mat[-1]]}]")


def vec_add(a, b):
    assert len(a) == len(b)
    return [x + y for x, y in zip(a, b)]


def circulant(c):
    n = len(c)
    res = [[None for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            res[i][j] = c[(i - j) % n]
    return res
