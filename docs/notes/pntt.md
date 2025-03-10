
Given the Vandermonde matrix:
$$
\mathbf V =
\begin{bmatrix}
    \alpha^{0\cdot 0}&\alpha ^{0\cdot 1}&\cdots &\alpha ^{0\cdot (n-1)}\\
    \alpha^{1\cdot 0}&\alpha^{1\cdot 1}&\cdots &\alpha ^{1\cdot (n-1)}\\
    \vdots &\vdots &\ddots &\vdots\\
    \alpha^{(n-1)\cdot 0}&\alpha ^{(n-1)\cdot 1}&\cdots &\alpha^{(n-1)\cdot (n-1)}
\end{bmatrix}
$$


If $n$ is a primary root of unity then:
$$
\mathbf V^{-1} =
\begin{bmatrix}
    \alpha^{-0\cdot 0}&\alpha ^{-0\cdot 1}&\cdots &\alpha ^{-0\cdot (n-1)}\\
    \alpha^{-1\cdot 0}&\alpha^{-1\cdot 1}&\cdots &\alpha ^{-1\cdot (n-1)}\\
    \vdots &\vdots &\ddots &\vdots\\
    \alpha^{-(n-1)\cdot 0}&\alpha ^{-(n-1)\cdot 1}&\cdots &\alpha^{-(n-1)\cdot (n-1)}
\end{bmatrix}
$$


$$
C =
\begin{bmatrix}
    c_0\\
    c_1\\
    \vdots\\
    c_k\\
    \vdots\\
    c_n\\
\end{bmatrix}

\qquad

c_k\dots c_n = 0
$$

$$
\begin{align*}
\mathbf V \cdot C &=
    \begin{bmatrix}
        \alpha^{0\cdot 0}&\alpha ^{0\cdot 1}&\cdots &\alpha ^{0\cdot (n-1)}\\
        \alpha^{1\cdot 0}&\alpha^{1\cdot 1}&\cdots &\alpha ^{1\cdot (n-1)}\\
        \vdots &\vdots &\ddots &\vdots\\
        \alpha^{(k-1)\cdot 0}&\alpha ^{(k-1)\cdot 1}&\dots &\alpha^{(k-1)\cdot (k-1)}\\
        \vdots &\vdots &\ddots &\vdots\\
        \alpha^{(n-1)\cdot 0}&\alpha ^{(n-1)\cdot 1}&\cdots &\alpha^{(n-1)\cdot (n-1)}\\
    \end{bmatrix}
    \begin{bmatrix}
        c_0\\
        c_1\\
        \vdots\\
        c_k\\
        \vdots\\
        c_n\\
    \end{bmatrix}

=
    \begin{bmatrix}
        \sum_{i=0}^{k-1} c_i\alpha^{0\cdot i}\\
        \sum_{i=0}^{k-1} c_i\alpha^{1\cdot i}\\
        \vdots\\
        \sum_{i=0}^{k-1} c_i\alpha^{(n-1)\cdot i}\\
    \end{bmatrix}
\end{align*}
$$