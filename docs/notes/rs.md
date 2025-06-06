Message polynomial:

$$
m(\omega^j) = \sum_{i=0}^{k-1} m_i \omega^{ji}
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
    \omega^{0\cdot 0}&\omega ^{0\cdot 1}&\cdots &\omega ^{0\cdot (l-1)}\\
    \omega^{1\cdot 0}&\omega^{1\cdot 1}&\cdots &\omega ^{1\cdot (l-1)}\\
    \vdots &\vdots &\ddots &\vdots\\
    \omega^{(l-1)\cdot 0}&\omega ^{(l-1)\cdot 1}&\cdots &\omega^{(l-1)\cdot (l-1)}
\end{bmatrix}
$$


$$
\omega =
\begin{bmatrix}
    \omega^{0\cdot k}\\
    \omega^{1\cdot k}\\
    \vdots\\
    \omega^{(l-1)\cdot k}
\end{bmatrix}

\qquad

\omega^{\odot-1} =
\begin{bmatrix}
    \omega^{-0\cdot k}\\
    \omega^{-1\cdot k}\\
    \vdots\\
    \omega^{-(l-1)\cdot k}
\end{bmatrix}
$$

Message polynomial evaluated at $l$ points:
$$
M =
\begin{bmatrix}
    m(\omega^0)\\
    m(\omega^1)\\
    \vdots\\
    m(\omega^{l-1})
\end{bmatrix}
$$

Codeword:
$$
C = \mathbf V^{-1} \cdot (\omega^{\odot-1} \odot M)
$$

---

Complexity analysis:

- $M$: $O(n + l\log l)$ using partial FFT algorithm (https://arxiv.org/pdf/2008.12559)

- Product with $\mathbf V^{-1}$: inverse NTT, $O(l\log l)$



- Hadamard product: element-wise multiplication, $O(l)$

---


Transmitted message:

$$
\bar M = M + S_m =
\begin{bmatrix}
    m(\omega^0) + s_m(\omega^0)\\
    m(\omega^1) + s_m(\omega^1)\\
    \vdots\\
    m(\omega^{l-1}) + s_m(\omega^{l-1})
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
    s(\omega^j)
    &= \sum_{i=0}^{n-1} e_i \omega^{ji}
    \\&= \sum_{i=0}^{k-1} e_i \omega^{ji} + \sum_{i=k}^{n-1} e_i \omega^{ji}
    \\&= s_m(\omega^j) + s_c(\omega^j)
\end{align*}
$$


$$
S = S_m + S_c =
\begin{bmatrix}
    s_m(\omega^0)\\
    s_m(\omega^1)\\
    \vdots\\
    s_m(\omega^{l-1})\\
\end{bmatrix}
+
\begin{bmatrix}
    s_c(\omega^0)\\
    s_c(\omega^1)\\
    \vdots\\
    s_c(\omega^{l-1})
\end{bmatrix}
=
\begin{bmatrix}
    s(\omega^0)\\
    s(\omega^1)\\
    \vdots\\
    s(\omega^{l-1})
\end{bmatrix}
$$


---

Recovering syndromes from transmitted message:

$$
\begin{align*}
\bar M + \omega \odot (\mathbf V \cdot \bar C)

&= \bar M + \omega \odot (\mathbf V \cdot C) + \omega \odot (\mathbf V \cdot E_c)

\\ &= \bar M + \omega \odot (\mathbf V \cdot \mathbf V^{-1} \cdot (\omega^{\odot-1} \odot M)) + \omega \odot (\mathbf V \cdot E_c)

\\ &= \bar M + \omega \odot \omega^{\odot-1} \odot M + \omega \odot (\mathbf V \cdot E_c)

\\ &= \bar M + M + \omega \odot (\mathbf V \cdot E_c)

\\ &= S_m + \omega \odot (\mathbf V \cdot E_c)

\\ &=
    \begin{bmatrix}
        s_m(\omega^0)\\
        s_m(\omega^1)\\
        \vdots\\
        s_m(\omega^{l-1})\\
    \end{bmatrix}
    +
    \begin{bmatrix}
        \omega^{0\cdot k}\\
        \omega^{1\cdot k}\\
        \vdots\\
        \omega^{(l-1)\cdot k}
    \end{bmatrix}

    \odot

    \begin{bmatrix}
        \omega^{0\cdot 0}&\omega ^{0\cdot 1}&\cdots &\omega ^{0\cdot (l-1)}\\
        \omega^{1\cdot 0}&\omega^{1\cdot 1}&\cdots &\omega ^{1\cdot (l-1)}\\\vdots &\vdots &\ddots &\vdots\\
        \omega^{(l-1)\cdot 0}&\omega ^{(l-1)\cdot 1}&\cdots &\omega^{(l-1)\cdot (l-1)}
    \end{bmatrix}

    \begin{bmatrix}
        e_{k}\\
        e_{k+1}\\
        \vdots\\
        e_{n-1}
    \end{bmatrix}

\\ &=
    \begin{bmatrix}
        s_m(\omega^0)\\
        s_m(\omega^1)\\
        \vdots\\
        s_m(\omega^{l-1})
    \end{bmatrix}
    +
    \begin{bmatrix}
        \omega^{0\cdot k}\\
        \omega^{1\cdot k}\\
        \vdots\\
        \omega^{(l-1)\cdot k}
    \end{bmatrix}

    \odot

    \begin{bmatrix}
        \sum_{i=k}^{n-1} e_{i} \omega^{0 \cdot (i-k)}\\
        \sum_{i=k}^{n-1} e_{i} \omega^{1 \cdot (i-k)}\\
        \vdots\\
        \sum_{i=k}^{n-1} e_{i} \omega^{(l-1) \cdot (i-k)}
    \end{bmatrix}

\\ &=
    \begin{bmatrix}
        s_m(\omega^0)\\
        s_m(\omega^1)\\
        \vdots\\
        s_m(\omega^{l-1})
    \end{bmatrix}
    +
    \begin{bmatrix}
        \sum_{i=k}^{n-1} e_{i} \omega^{0 \cdot i}\\
        \sum_{i=k}^{n-1} e_{i} \omega^{1 \cdot i}\\
        \vdots\\
        \sum_{i=k}^{n-1} e_{i} \omega^{(l-1) \cdot i}
    \end{bmatrix}

\\ &=
    \begin{bmatrix}
        s_m(\omega^0)\\
        s_m(\omega^1)\\
        \vdots\\
        s_m(\omega^{l-1})
    \end{bmatrix}
    +
    \begin{bmatrix}
        s_c(\omega^0)\\
        s_c(\omega^1)\\
        \vdots\\
        s_c(\omega^{l-1})
    \end{bmatrix}

\\ &= S_m + S_c
\\ &= S
\end{align*}
$$


---


$$
\begin{align*}
\mathbf V \cdot E_{c0} &=
\begin{bmatrix}
    \omega^{0\cdot 0}&\omega ^{0\cdot 1}&\cdots &\omega ^{0\cdot (n-1)}\\
    \omega^{1\cdot 0}&\omega^{1\cdot 1}&\cdots &\omega ^{1\cdot (n-1)}\\
    \vdots &\vdots &\ddots &\vdots\\
    \omega^{(n-1)\cdot 0}&\omega ^{(n-1)\cdot 1}&\cdots &\omega^{(n-1)\cdot (n-1)}
\end{bmatrix}
\begin{bmatrix}
    a_0\\
    a_0\\
    a_1\\
    a_1\\
    \vdots\\
    a_{l-1}\\
    a_{l-1}\\
    0\\
    0\\
    \vdots\\
\end{bmatrix}
\\&=
\begin{bmatrix}
    \sum_{i=0}^{l-1} a_{i} (\omega^{0 \cdot 2i} + \omega^{0 \cdot (2i + 1)})\\
    \sum_{i=0}^{l-1} a_{i} (\omega^{1 \cdot 2i} + \omega^{1 \cdot (2i + 1)})\\
    \vdots\\
    \sum_{i=0}^{l-1} a_{i} (\omega^{(n-1) \cdot 2i} + \omega^{(n-1) \cdot (2i + 1)})
\end{bmatrix}
\\&=
\begin{bmatrix}
    \sum_{i=0}^{l-1} a_{i} \omega^{0 \cdot 2i} (1 + \omega^0)\\
    \sum_{i=0}^{l-1} a_{i} \omega^{1 \cdot 2i} (1 + \omega^1)\\
    \vdots\\
    \sum_{i=0}^{l-1} a_{i} \omega^{(n-1) \cdot 2i} (1 + \omega^{n-1})
\end{bmatrix}
\end{align*}
$$