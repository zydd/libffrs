Message polynomial:

$$
m(\alpha^j) = \sum_{i=0}^{k-1} m_i \alpha^{ji}
$$

---


$$
\begin{align*}
    k &= \text{message length}\\
    l &= \text{codeword length}\\
    n &= k + l\\
    \nu &= \text{errors in message+codeword}\\
\end{align*}
$$

$$
0 \leq j \leq l-1
$$

---


$$
\mathbf V =
\begin{bmatrix}
    \alpha^{0\cdot 0}&\alpha ^{0\cdot 1}&\cdots &\alpha ^{0\cdot (l-1)}\\
    \alpha^{1\cdot 0}&\alpha^{1\cdot 1}&\cdots &\alpha ^{1\cdot (l-1)}\\
    \vdots &\vdots &\ddots &\vdots\\
    \alpha^{(l-1)\cdot 0}&\alpha ^{(l-1)\cdot 1}&\cdots &\alpha^{(l-1)\cdot (l-1)}
\end{bmatrix}
$$


$$
\Alpha =
\begin{bmatrix}
    \alpha^{0\cdot k}\\
    \alpha^{1\cdot k}\\
    \vdots\\
    \alpha^{(l-1)\cdot k}
\end{bmatrix}

\qquad

\Alpha^{\odot-1} =
\begin{bmatrix}
    \alpha^{-0\cdot k}\\
    \alpha^{-1\cdot k}\\
    \vdots\\
    \alpha^{-(l-1)\cdot k}
\end{bmatrix}
$$

Message polynomial evaluated at $l$ points:
$$
M =
\begin{bmatrix}
    m(\alpha^0)\\
    m(\alpha^1)\\
    \vdots\\
    m(\alpha^{l-1})
\end{bmatrix}
$$

Codeword:
$$
C = \mathbf V^{-1} \cdot (\Alpha^{\odot-1} \odot M)
$$

---

Complexity analysis:

- $M$: $O(n + l\log l)$ using partial FFT algorithm (https://arxiv.org/pdf/2008.12559)

- Product with $\mathbf V^{-1}$: inverse NTT, $O(l\log l)$

    > Period of $\alpha^x$ needs to be $\geq n$ for the error locator polynomial to work, that means $\alpha^k \neq 1$.

    > Need to check if FFT algorithm is still able to compute this product when $\alpha$ is not a $k$-th primary root of unity.

    > Backup plan is to chose a primary root of unity with period $\geq n$ and use partial FFT algorithm.

- Hadamard product: element-wise multiplication, $O(l)$

---


Transmitted message:

$$
\bar M = M + S_m =
\begin{bmatrix}
    m(\alpha^0) + s_m(\alpha^0)\\
    m(\alpha^1) + s_m(\alpha^1)\\
    \vdots\\
    m(\alpha^{l-1}) + s_m(\alpha^{l-1})
\end{bmatrix}
$$

Transmitted codeword:
$$
\bar C = C + E_c
$$


Errors:

$$
E_m =
\begin{bmatrix}
    e_{0}\\
    e_{1}\\
    \vdots\\
    e_{k-1}
\end{bmatrix}

\qquad

E_c =
\begin{bmatrix}
    e_{k}\\
    e_{k+1}\\
    \vdots\\
    e_{n-1}
\end{bmatrix}
$$

---


Decompose syndromes on contributions from transmitted message and transmitted codeword:

$$
\begin{align*}
    s(\alpha^j)
    &= \sum_{i=0}^{n-1} e_i \alpha^{ji}
    \\&= \sum_{i=0}^{k-1} e_i \alpha^{ji} + \sum_{i=k}^{n-1} e_i \alpha^{ji}
    \\&= s_m(\alpha^j) + s_c(\alpha^j)
\end{align*}
$$


$$
S = S_m + S_c =
\begin{bmatrix}
    s_m(\alpha^0)\\
    s_m(\alpha^1)\\
    \vdots\\
    s_m(\alpha^{l-1})\\
\end{bmatrix}
+
\begin{bmatrix}
    s_c(\alpha^0)\\
    s_c(\alpha^1)\\
    \vdots\\
    s_c(\alpha^{l-1})
\end{bmatrix}
=
\begin{bmatrix}
    s(\alpha^0)\\
    s(\alpha^1)\\
    \vdots\\
    s(\alpha^{l-1})
\end{bmatrix}
$$


---

Recovering syndromes from transmitted message:

$$
\begin{align*}
\bar M + \Alpha \odot (\mathbf V \cdot \bar C)

&= \bar M + \Alpha \odot (\mathbf V \cdot C) + \Alpha \odot (\mathbf V \cdot E_c)

\\ &= \bar M + \Alpha \odot (\mathbf V \cdot \mathbf V^{-1} \cdot (\Alpha^{\odot-1} \odot M)) + \Alpha \odot (\mathbf V \cdot E_c)

\\ &= \bar M + \Alpha \odot \Alpha^{\odot-1} \odot M + \Alpha \odot (\mathbf V \cdot E_c)

\\ &= \bar M + M + \Alpha \odot (\mathbf V \cdot E_c)

\\ &= S_m + \Alpha \odot (\mathbf V \cdot E_c)

\\ &=
    \begin{bmatrix}
        s_m(\alpha^0)\\
        s_m(\alpha^1)\\
        \vdots\\
        s_m(\alpha^{l-1})\\
    \end{bmatrix}
    +
    \begin{bmatrix}
        \alpha^{0\cdot k}\\
        \alpha^{1\cdot k}\\
        \vdots\\
        \alpha^{(l-1)\cdot k}
    \end{bmatrix}

    \odot

    \begin{bmatrix}
        \alpha^{0\cdot 0}&\alpha ^{0\cdot 1}&\cdots &\alpha ^{0\cdot (l-1)}\\
        \alpha^{1\cdot 0}&\alpha^{1\cdot 1}&\cdots &\alpha ^{1\cdot (l-1)}\\\vdots &\vdots &\ddots &\vdots\\
        \alpha^{(l-1)\cdot 0}&\alpha ^{(l-1)\cdot 1}&\cdots &\alpha^{(l-1)\cdot (l-1)}
    \end{bmatrix}

    \begin{bmatrix}
        e_{k}\\
        e_{k+1}\\
        \vdots\\
        e_{n-1}
    \end{bmatrix}

\\ &=
    \begin{bmatrix}
        s_m(\alpha^0)\\
        s_m(\alpha^1)\\
        \vdots\\
        s_m(\alpha^{l-1})
    \end{bmatrix}
    +
    \begin{bmatrix}
        \alpha^{0\cdot k}\\
        \alpha^{1\cdot k}\\
        \vdots\\
        \alpha^{(l-1)\cdot k}
    \end{bmatrix}

    \odot

    \begin{bmatrix}
        \sum_{i=k}^{n-1} e_{i} \alpha^{0 \cdot (i-k)}\\
        \sum_{i=k}^{n-1} e_{i} \alpha^{1 \cdot (i-k)}\\
        \vdots\\
        \sum_{i=k}^{n-1} e_{i} \alpha^{(l-1) \cdot (i-k)}
    \end{bmatrix}

\\ &=
    \begin{bmatrix}
        s_m(\alpha^0)\\
        s_m(\alpha^1)\\
        \vdots\\
        s_m(\alpha^{l-1})
    \end{bmatrix}
    +
    \begin{bmatrix}
        \sum_{i=k}^{n-1} e_{i} \alpha^{0 \cdot i}\\
        \sum_{i=k}^{n-1} e_{i} \alpha^{1 \cdot i}\\
        \vdots\\
        \sum_{i=k}^{n-1} e_{i} \alpha^{(l-1) \cdot i}
    \end{bmatrix}

\\ &=
    \begin{bmatrix}
        s_m(\alpha^0)\\
        s_m(\alpha^1)\\
        \vdots\\
        s_m(\alpha^{l-1})
    \end{bmatrix}
    +
    \begin{bmatrix}
        s_c(\alpha^0)\\
        s_c(\alpha^1)\\
        \vdots\\
        s_c(\alpha^{l-1})
    \end{bmatrix}

\\ &= S_m + S_c
\\ &= S
\end{align*}
$$

