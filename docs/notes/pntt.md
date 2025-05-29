
Given the Vandermonde matrix:
$$
\mathbf V =
\begin{bmatrix}
    \omega^{0\cdot 0}&\omega ^{0\cdot 1}&\cdots &\omega ^{0\cdot (n-1)}\\
    \omega^{1\cdot 0}&\omega^{1\cdot 1}&\cdots &\omega ^{1\cdot (n-1)}\\
    \vdots &\vdots &\ddots &\vdots\\
    \omega^{(n-1)\cdot 0}&\omega ^{(n-1)\cdot 1}&\cdots &\omega^{(n-1)\cdot (n-1)}\\
\end{bmatrix}
$$


If $\omega$ is a primary $n$-th root of unity then:
$$
\mathbf V^{-1} =
\begin{bmatrix}
    \omega^{-0\cdot 0}&\omega ^{-0\cdot 1}&\cdots &\omega ^{-0\cdot (n-1)}\\
    \omega^{-1\cdot 0}&\omega^{-1\cdot 1}&\cdots &\omega ^{-1\cdot (n-1)}\\
    \vdots &\vdots &\ddots &\vdots\\
    \omega^{-(n-1)\cdot 0}&\omega ^{-(n-1)\cdot 1}&\cdots &\omega^{-(n-1)\cdot (n-1)}\\
\end{bmatrix}
$$


---


$$
\mathbf a =
\begin{bmatrix}
    a_0\\
    a_1\\
    \vdots\\
    a_n\\
\end{bmatrix}
$$

---


Circulant matrix
================

$$

\mathbf{C} =
\begin{bmatrix}
c_0      & c_{n-1} & \cdots  & c_2     & c_1     \\
c_1      & c_0     & c_{n-1} &         & c_2     \\
\vdots   & c_1     & c_0     & \ddots  & \vdots  \\
c_{n-2}  &         & \ddots  & \ddots  & c_{n-1} \\
c_{n-1}  & c_{n-2} & \cdots  & c_1     & c_0     \\
\end{bmatrix}

$$


$$
\begin{bmatrix}
    \omega^0 & \omega^0 & \omega^{4} & \omega^{3} & \omega^{2} & \omega^{1} \\
    \omega^{1} & \omega^0 & \omega^0 & \omega^{4} & \omega^{3} & \omega^{2} \\
    \omega^{2} & \omega^{1} & \omega^0 & \omega^0 & \omega^{4} & \omega^{3} \\
    \omega^{3} & \omega^{2} & \omega^{1} & \omega^0 & \omega^0 & \omega^{4} \\
    \omega^{4} & \omega^{3} & \omega^{2} & \omega^{1} & \omega^0 & \omega^0 \\
    \omega^0 & \omega^{4} & \omega^{3} & \omega^{2} & \omega^{1} & \omega^0 \\
\end{bmatrix}
\odot
\begin{bmatrix}
    \omega^0 \\ \omega^{-1} \\ \omega^{-2} \\ \omega^{-3} \\ \omega^{-4} \\ \omega^{0} \\
\end{bmatrix}
=
\begin{bmatrix}
    \omega^0 & \omega^0 & \omega^{4} & \omega^{3} & \omega^{2} & \omega^{1} \\
    \omega^{0} & \omega^4 & \omega^4 & \omega^{3} & \omega^{2} & \omega^{1} \\
    \omega^{0} & \omega^{4} & \omega^3 & \omega^3 & \omega^{2} & \omega^{1} \\
    \omega^{0} & \omega^{4} & \omega^{3} & \omega^2 & \omega^2 & \omega^{1} \\
    \omega^{0} & \omega^{4} & \omega^{3} & \omega^{2} & \omega^1 & \omega^1 \\
    \omega^0 & \omega^{4} & \omega^{3} & \omega^{2} & \omega^{1} & \omega^0 \\
\end{bmatrix}
$$



$$

\begin{bmatrix}
    \omega^0 & \omega^3 & \omega^2 & \omega^1 \\
    \omega^1 & \omega^0 & \omega^3 & \omega^2 \\
    \omega^2 & \omega^1 & \omega^0 & \omega^3 \\
    \omega^3 & \omega^2 & \omega^1 & \omega^0 \\
\end{bmatrix}
\odot
\begin{bmatrix}
    \omega^0 \\ \omega^{-1} \\ \omega^{-2} \\ \omega^{-3} \\
\end{bmatrix}
=
\begin{bmatrix}
    \omega^0 & \omega^3 & \omega^2 & \omega^1 \\
    \omega^0 & \omega^4 & \omega^2 & \omega^1 \\
    \omega^0 & \omega^4 & \omega^3 & \omega^1 \\
    \omega^0 & \omega^4 & \omega^3 & \omega^2 \\
\end{bmatrix}
$$

---


$$
\mathbf{C} =
\begin{bmatrix}
    \omega^0 & \omega^1 & \omega^2 & \omega^3 & \omega^4 \\
    \omega^1 & \omega^2 & \omega^3 & \omega^4 & \omega^0 \\
    \omega^2 & \omega^3 & \omega^4 & \omega^0 & \omega^1 \\
    \omega^3 & \omega^4 & \omega^0 & \omega^1 & \omega^2 \\
    \omega^4 & \omega^0 & \omega^1 & \omega^2 & \omega^3 \\
\end{bmatrix}
\begin{bmatrix}
    a_0\\a_1\\a_2\\a_3\\a_4\\
\end{bmatrix}
=
\begin{bmatrix}
    \omega^0 a_0 &+& \omega^1 a_1 &+& \omega^2 a_2 &+& \omega^3 a_3 &+& \omega^4 a_4\\
    \omega^1 a_0 &+& \omega^2 a_1 &+& \omega^3 a_2 &+& \omega^4 a_3 &+& \omega^0 a_4\\
    \omega^2 a_0 &+& \omega^3 a_1 &+& \omega^4 a_2 &+& \omega^0 a_3 &+& \omega^1 a_4\\
    \omega^3 a_0 &+& \omega^4 a_1 &+& \omega^0 a_2 &+& \omega^1 a_3 &+& \omega^2 a_4\\
    \omega^4 a_0 &+& \omega^0 a_1 &+& \omega^1 a_2 &+& \omega^2 a_3 &+& \omega^3 a_4\\
\end{bmatrix}
$$



---

$$
NTT^\psi(\mathbf a) = \hat a_j = \sum_{i=0}^{n-1} \psi^{2ij+i} a_i
$$

$$
\mathbf{\hat a} =
\begin{bmatrix}
    1 & \psi^1 & \psi^2 & \psi^3 & \psi^4\\
    1 & \psi^3 & \psi^6 & \psi^9 & \psi^{12}\\
    1 & \psi^5 & \psi^{10} & \psi^{15} & \psi^{20}\\
    1 & \psi^7 & \psi^{14} & \psi^{21} & \psi^{28}\\
    1 & \psi^9 & \psi^{18} & \psi^{27} & \psi^{36}\\
\end{bmatrix}
\begin{bmatrix}
    a_0\\a_1\\a_2\\a_3\\a_4\\
\end{bmatrix}
=
\begin{bmatrix}
    1 & \psi^1 & \psi^2 & \psi^3 & \psi^4\\
    1 & \psi^3 & \psi^1 & \psi^4 & \psi^2\\
    1 & \psi^0 & \psi^0 & \psi^0 & \psi^0\\
    1 & \psi^2 & \psi^4 & \psi^1 & \psi^3\\
    1 & \psi^4 & \psi^3 & \psi^2 & \psi^1\\
\end{bmatrix}
\begin{bmatrix}
    a_0\\a_1\\a_2\\a_3\\a_4\\
\end{bmatrix}
$$
---

$$
NTT(\mathbf a) = \hat{\mathbf a} = \mathbf V \cdot \mathbf a
$$

$$
\hat a_j = \sum_{i=0}^{n-1} \omega^{ij} a_i
$$


$$
\mathbf V =
\begin{bmatrix}
    1 & 1 & 1 & 1 & 1\\
    1 & \omega^1 & \omega^2 & \omega^3 & \omega^4\\
    1 & \omega^2 & \omega^4 & \omega^1 & \omega^3\\
    1 & \omega^3 & \omega^1 & \omega^4 & \omega^2\\
    1 & \omega^4 & \omega^3 & \omega^2 & \omega^1\\
\end{bmatrix}

\quad\Longrightarrow\quad
\begin{bmatrix}
    a_0 &+& a_1 &+& a_2 &+& a_3 &+& a_4\\
    a_0 &+& \omega^1 a_1 &+& \omega^2 a_2 &+& \omega^3 a_3 &+& \omega^4 a_4\\
    a_0 &+& \omega^2 a_1 &+& \omega^4 a_2 &+& \omega^1 a_3 &+& \omega^3 a_4\\
    a_0 &+& \omega^3 a_1 &+& \omega^1 a_2 &+& \omega^4 a_3 &+& \omega^2 a_4\\
    a_0 &+& \omega^4 a_1 &+& \omega^3 a_2 &+& \omega^2 a_3 &+& \omega^1 a_4\\
\end{bmatrix}
$$

$$
\quad\Longrightarrow\quad
\begin{bmatrix}
    a_0 &+& a_1 &+& a_2 &+& a_3 &+& a_4\\
    a_0 &+& \omega^1 a_1 &+& \omega^2 a_2 &+& \omega^{-2} a_3 &+& \omega^{-1} a_4\\
    a_0 &+& \omega^2 a_1 &+& \omega^{-1} a_2 &+& \omega^1 a_3 &+& \omega^{-2} a_4\\
    a_0 &+& \omega^{-2} a_1 &+& \omega^1 a_2 &+& \omega^{-1} a_3 &+& \omega^2 a_4\\
    a_0 &+& \omega^{-1} a_1 &+& \omega^{-2} a_2 &+& \omega^2 a_3 &+& \omega^1 a_4\\
\end{bmatrix}
$$



$$
\begin{align*}
\hat{\mathbf a} &=
    \begin{bmatrix}
        \omega^{0\cdot 0}&\omega ^{0\cdot 1}&\cdots &\omega ^{0\cdot (n-1)}\\
        \omega^{1\cdot 0}&\omega^{1\cdot 1}&\cdots &\omega ^{1\cdot (n-1)}\\
        \vdots &\vdots &\ddots &\vdots\\
        \omega^{(n-1)\cdot 0}&\omega ^{(n-1)\cdot 1}&\cdots &\omega^{(n-1)\cdot (n-1)}\\
    \end{bmatrix}
    \begin{bmatrix}
        a_0\\
        a_1\\
        \vdots\\
        a_n\\
    \end{bmatrix}

=
    \begin{bmatrix}
    \sum_{i=0}^{n-1}\omega^{i\cdot 0} a_i\\
    \sum_{i=0}^{n-1}\omega^{i\cdot 1} a_i\\
    \vdots\\
    \sum_{i=0}^{n-1}\omega^{i\cdot (n-1)} a_i\\
    \end{bmatrix}
\end{align*}
$$


---

$$
\omega^{k+i \cdot n} = \omega^k , \quad i \in \mathbb{N}
$$

For $GF(2^k)$, with $k \in \mathbb{N^*}$, $n$ will always be odd.


---



$$
\begin{align*}
\hat{\mathbf a} &=

    \begin{bmatrix}
    \sum_{i=0}^{n-1}\omega^{i\cdot 0} a_i\\
    \sum_{i=0}^{n-1}\omega^{i\cdot 1} a_i\\
    \vdots\\
    \sum_{i=0}^{n-1}\omega^{i\cdot (n-1)} a_i\\
    \end{bmatrix}

\\&=
    \begin{bmatrix}
    \sum_{i=0}^{n-1}\omega^{i\cdot 0} a_i\\
    \sum_{i=0}^{n-1}\omega^{i\cdot 1} a_i\\
    \vdots\\
    \sum_{i=0}^{n-1}\omega^{i\cdot (n-1)} a_i\\
    \end{bmatrix}
\end{align*}
$$







<br><br><br><br>
<br><br><br><br>



Cooley-Tukey
============

$$
\begin{align*}
\hat a_j &= \sum_{i=0}^{n-1} \omega^{ij} a_i
\\&= \sum_{\substack{i = 0 \\ i \text{ even}}}^{n-1} \omega^{ij} a_i + \sum_{\substack{i = 1 \\ i \text{ odd}}}^{n-1} \omega^{ij} a_i
\\&= \sum_{i = 0}^{(n-1) / 2} \omega^{2ij} a_{2i} + \sum_{i = 0}^{(n-1)/2 - 1} \omega^{(2i+1)j} a_{2i+1}
\\&= \sum_{i = 0}^{(n-1) / 2 - 1} \omega^{2ij} a_{2i} + \omega^{(n-1) j} + \omega^j\sum_{i = 0}^{(n-1)/2 - 1} \omega^{2ij} a_{2i+1}
\end{align*}
$$

---

$$
\begin{align*}
\hat a_{j + (n-1)/2} &= \sum_{i = 0}^{(n-1) / 2 - 1} \omega^{2i(j + (n-1)/2)} a_{2i} + \omega^{(n-1) (j + (n-1)/2)} + \omega^{j + (n-1)/2}\sum_{i = 0}^{(n-1)/2 - 1} \omega^{2i(j + (n-1)/2)} a_{2i+1}
\\ &= \sum_{i = 0}^{(n-1) / 2 - 1} \omega^{2ij + in-i} a_{2i} + \omega^{(n-1) (j + (n-1)/2)} + \omega^{j + (n-1)/2}\sum_{i = 0}^{(n-1)/2 - 1} \omega^{2ij + in-i} a_{2i+1}
\\ &= \sum_{i = 0}^{(n-1) / 2 - 1} \omega^{2ij -i} a_{2i} + \omega^{(n-1) (j + (n-1)/2)} + \omega^{j + (n-1)/2}\sum_{i = 0}^{(n-1)/2 - 1} \omega^{2ij -i} a_{2i+1}
\end{align*}
$$

Deduction for Cooley-Tukey fails because $n$ is odd, resulting in $\omega^{2ij -i}$ inside the summations for $\hat a_{j + (n-1)/2}$

---

<!-- 
---

$$
C =
\begin{bmatrix}
    a_0\\
    a_1\\
    \vdots\\
    a_k\\
    \vdots\\
    a_n\\
\end{bmatrix}

\qquad

a_k\dots a_n = 0
$$

$$
\begin{align*}
\mathbf V \cdot C &=
    \begin{bmatrix}
        \omega^{0\cdot 0}&\omega ^{0\cdot 1}&\cdots &\omega ^{0\cdot (n-1)}\\
        \omega^{1\cdot 0}&\omega^{1\cdot 1}&\cdots &\omega ^{1\cdot (n-1)}\\
        \vdots &\vdots &\ddots &\vdots\\
        \omega^{(k-1)\cdot 0}&\omega ^{(k-1)\cdot 1}&\dots &\omega^{(k-1)\cdot (k-1)}\\
        \vdots &\vdots &\ddots &\vdots\\
        \omega^{(n-1)\cdot 0}&\omega ^{(n-1)\cdot 1}&\cdots &\omega^{(n-1)\cdot (n-1)}\\
    \end{bmatrix}
    \begin{bmatrix}
        a_0\\
        a_1\\
        \vdots\\
        a_k\\
        \vdots\\
        a_n\\
    \end{bmatrix}

=
    \begin{bmatrix}
        \sum_{i=0}^{k-1}\omega^{0\cdot i} a_i\\
        \sum_{i=0}^{k-1}\omega^{1\cdot i}a_i\\
        \vdots\\
        \sum_{i=0}^{k-1}\omega^{(n-1)\cdot i} a_i\\
    \end{bmatrix}
\end{align*}
$$
 -->
