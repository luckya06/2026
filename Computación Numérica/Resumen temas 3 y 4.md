# CN — Resumen Temas 3 y 4

---

# TEMA 3 — APROXIMACIÓN DE FUNCIONES

> **Interpolación** → la función pasa EXACTAMENTE por los puntos (datos exactos)  
> **Ajuste / Regresión** → minimiza el error (datos con ruido)  
> Con $n+1$ nodos distintos existe un único polinomio de grado $\leq n$ que interpola.

---

## 3.1 Interpolación de Lagrange

### Polinomios fundamentales

$$L_i(x) = \prod_{j=0,\, j\neq i}^{n} \frac{x - x_j}{x_i - x_j}$$

Propiedad: $L_i(x_j) = \begin{cases} 1 & \text{si } j = i \\ 0 & \text{si } j \neq i \end{cases}$

### Polinomio interpolante

$$P_n(x) = \sum_{i=0}^{n} f(x_i)\, L_i(x)$$

### Error de interpolación

$$E(x) = f(x) - P_n(x) = f^{(n+1)}(c)\,\frac{(x-x_0)(x-x_1)\cdots(x-x_n)}{(n+1)!}$$

**Cota práctica** ($c$ desconocido, se acota $|f^{(n+1)}|$ con $M$):

$$|E(x)| \leq \frac{M}{(n+1)!}\,|(x-x_0)(x-x_1)\cdots(x-x_n)|$$

> ⚠️ Si hay $n+1$ nodos → derivada de orden $n+1$. Acotar $|\sin| \leq 1$, $|\cos| \leq 1$, etc.

---

## 3.1b Forma de Newton — Diferencias divididas

**Tabla** (los $c_k$ son la primera fila de cada columna):

| $x$ | $f[x]$ | orden 1 | orden 2 | orden 3 |
|-----|--------|---------|---------|---------|
| $x_0$ | $c_0$ | | | |
| $x_1$ | $f[x_1]$ | $c_1 = \dfrac{f[x_1]-f[x_0]}{x_1-x_0}$ | | |
| $x_2$ | $f[x_2]$ | $\dfrac{f[x_2]-f[x_1]}{x_2-x_1}$ | $c_2 = \dfrac{[\cdot]-[\cdot]}{x_2-x_0}$ | |
| $x_3$ | $f[x_3]$ | $\cdots$ | $\cdots$ | $c_3$ |

**Polinomio:**

$$P_n(x) = c_0 + c_1(x-x_0) + c_2(x-x_0)(x-x_1) + c_3(x-x_0)(x-x_1)(x-x_2) + \cdots$$

> Ventaja: al añadir un nodo nuevo solo se añade un término.

---

## 3.2 Interpolación a trozos — Spline cúbico

En cada subintervalo $[x_i, x_{i+1}]$: $S_i(x) = a_i + b_i x + c_i x^2 + d_i x^3$

Para $n+1$ nodos → $n$ subintervalos → $4n$ incógnitas:

| Condición | Ecuaciones |
|-----------|-----------|
| Interpolación en cada nodo | $2n$ |
| Continuidad de $S'$ en interiores | $n-1$ |
| Continuidad de $S''$ en interiores | $n-1$ |
| **Faltan 2** → condición de contorno | $2$ |

**Tipos:**
- **Natural:** $S''(x_0) = 0$ y $S''(x_n) = 0$
- **Sujeto (clamped):** $S'(x_0) = f'(x_0)$ y $S'(x_n) = f'(x_n)$

---

## 3.3 Mínimos cuadrados — Recta de regresión

Minimiza $E = \sum_{k=1}^{m}(a_0 + a_1 x_k - y_k)^2$. Sistema normal:

$$\begin{pmatrix} m & \sum x_k \\ \sum x_k & \sum x_k^2 \end{pmatrix} \begin{pmatrix} a_0 \\ a_1 \end{pmatrix} = \begin{pmatrix} \sum y_k \\ \sum x_k y_k \end{pmatrix}$$

### Linealizaciones frecuentes

| Función | Transformación | Forma lineal |
|---------|---------------|--------------|
| $y = b\,e^{ax}$ | $Y = \ln y,\quad B = \ln b$ | $Y = B + ax$ |
| $y = b\,x^a$ | $Y = \ln y,\quad X = \ln x,\quad B = \ln b$ | $Y = B + aX$ |
| $y = a\,e^{bx}$ | $Y = \ln y,\quad A = \ln a$ | $Y = A + bx$ |
| $y = \frac{a}{x+b}$ | $Y = 1/y$ | sistema $2\times2$ directo |

> ⚠️ Transforma TODOS los $y_i$ (ej: $Y_i = \ln(y_i)$, no solo la fórmula abstracta).

---

## 3.4–3.5 Bases polinómicas y polinomios ortogonales

Sistema normal general: $A^T A\,\mathbf{a} = A^T \mathbf{y}$, donde $A_{ij} = \varphi_j(x_i)$.

Si la base es **ortogonal** respecto al producto escalar $\langle f,g \rangle$:

$$a_j = \frac{\langle f,\, \varphi_j \rangle}{\langle \varphi_j,\, \varphi_j \rangle}$$

- Discreto: $\langle f, g \rangle = \sum_i f(x_i)\,g(x_i)$
- Continuo: $\langle f, g \rangle = \int f(x)\,g(x)\,dx$

### Fourier (base trigonométrica, periodo $T$)

$$a_0 = \frac{1}{T}\int f(x)\,dx \qquad a_k = \frac{2}{T}\int f(x)\cos(k\omega x)\,dx \qquad b_k = \frac{2}{T}\int f(x)\sin(k\omega x)\,dx$$

---

# TEMA 4 — DERIVACIÓN E INTEGRACIÓN NUMÉRICA

---

## 4.1 Derivación numérica — una variable

| Fórmula | Expresión | Error | Orden |
|---------|-----------|-------|-------|
| Progresiva | $f'(x_0) \approx \dfrac{f(x_0+h)-f(x_0)}{h}$ | $-\dfrac{h}{2}f''(c)$ | $O(h)$ |
| Regresiva | $f'(x_0) \approx \dfrac{f(x_0)-f(x_0-h)}{h}$ | $+\dfrac{h}{2}f''(c)$ | $O(h)$ |
| **Centrada** | $f'(x_0) \approx \dfrac{f(x_0+h)-f(x_0-h)}{2h}$ | $-\dfrac{h^2}{6}f'''(c)$ | $\mathbf{O(h^2)}$ |

Segunda derivada centrada:

$$f''(x_0) \approx \frac{f(x_0+h) - 2f(x_0) + f(x_0-h)}{h^2} \qquad \text{error } O(h^2)$$

Orden de la fórmula: si $E_h = Kh^n$ → $n$ es el orden. Con $h=0.1$: orden 1 → error $\sim 0.1$, orden 2 → error $\sim 0.01$.

---

## 4.2 Derivadas parciales

$$\frac{\partial f}{\partial x}(x_0,y_0) \approx \frac{f(x_0+h,y_0)-f(x_0-h,y_0)}{2h}$$

$$\frac{\partial f}{\partial y}(x_0,y_0) \approx \frac{f(x_0,y_0+k)-f(x_0,y_0-k)}{2k}$$

$$\frac{\partial^2 f}{\partial x\,\partial y} \approx \frac{f(x+h,y+k)-f(x+h,y-k)-f(x-h,y+k)+f(x-h,y-k)}{4hk}$$

---

## 4.3 Integración — Newton-Cotes

### Trapecio simple $(h = b-a)$

$$\int_a^b f(x)\,dx \approx \frac{h}{2}\bigl(f(a)+f(b)\bigr) \qquad \text{error: } -\frac{h^3}{12}f''(c)$$

Exacta para grado $\leq 1$.

### Simpson simple $(h = \frac{b-a}{2},\; m = \frac{a+b}{2})$

$$\int_a^b f(x)\,dx \approx \frac{b-a}{6}\left(f(a)+4f(m)+f(b)\right) \qquad \text{error: } -\frac{h^5}{90}f^{(4)}(c)$$

Exacta para grado $\leq 3$ (aunque usa polinomio de grado 2).

### Trapecio compuesta $(m$ subintervalos, $h=\frac{b-a}{m})$

$$\int_a^b f(x)\,dx \approx \frac{h}{2}\bigl[f(x_0)+2f(x_1)+2f(x_2)+\cdots+2f(x_{m-1})+f(x_m)\bigr]$$

$$\text{Error total: } -\frac{(b-a)h^2}{12}f''(c) \qquad \Rightarrow \quad m > \sqrt{\frac{(b-a)^3 M_2}{12\,\varepsilon}}$$

### Simpson compuesta $(m$ par$)$

$$\int_a^b f(x)\,dx \approx \frac{h}{3}\bigl[f(x_0)+4f(x_1)+2f(x_2)+4f(x_3)+\cdots+4f(x_{m-1})+f(x_m)\bigr]$$

$$\text{Error total: } -\frac{(b-a)h^4}{180}f^{(4)}(c) \qquad \Rightarrow \quad m > \sqrt[4]{\frac{(b-a)^5 M_4}{180\,\varepsilon}}$$

> Patrón trapecio: $1,2,2,\ldots,2,1$ — Patrón Simpson: $1,4,2,4,2,\ldots,4,1$

---

## 4.4 Cuadratura de Gauss-Legendre

Con $n$ puntos, exacta para grado $\leq 2n-1$.

| $n$ | Nodos $x_i$ en $[-1,1]$ | Pesos $w_i$ | Exacta grado $\leq$ |
|-----|------------------------|------------|---------------------|
| 1 | $0$ | $2$ | $1$ |
| 2 | $\pm\frac{1}{\sqrt{3}} \approx \pm 0.5774$ | $1,\; 1$ | $3$ |
| 3 | $0,\; \pm\sqrt{\frac{3}{5}} \approx \pm 0.7746$ | $\frac{8}{9},\; \frac{5}{9},\;\frac{5}{9}$ | $5$ |
| 4 | $\pm 0.3399,\; \pm 0.8611$ | $0.6521,\; 0.3479$ | $7$ |

**Fórmula de 2 nodos:**
$$\int_{-1}^{1} f(x)\,dx \approx f\!\left(-\tfrac{1}{\sqrt{3}}\right) + f\!\left(\tfrac{1}{\sqrt{3}}\right)$$

**Fórmula de 3 nodos:**
$$\int_{-1}^{1} f(x)\,dx \approx \frac{5}{9}f\!\left(-\sqrt{\tfrac{3}{5}}\right) + \frac{8}{9}f(0) + \frac{5}{9}f\!\left(\sqrt{\tfrac{3}{5}}\right)$$

**Cambio de variable para $[a,b]$:**

$$y_i = \frac{b-a}{2}\,x_i + \frac{a+b}{2} \qquad \Rightarrow \qquad \int_a^b f(x)\,dx \approx \frac{b-a}{2}\sum_{i} w_i\, f(y_i)$$

> ⚠️ No olvides el factor $\dfrac{b-a}{2}$ delante.

---

## Resumen de errores

| Método | Error | Orden |
|--------|-------|-------|
| Deriv. progresiva/regresiva | $\pm\dfrac{h}{2}f''(c)$ | $O(h)$ |
| Deriv. centrada 1ª | $-\dfrac{h^2}{6}f'''(c)$ | $O(h^2)$ |
| Deriv. 2ª centrada | $\dfrac{h^2}{12}f^{(4)}(c)$ | $O(h^2)$ |
| Interp. Lagrange ($n+1$ nodos) | $\dfrac{M}{(n+1)!}\|\omega(x)\|$ | depende de $M$ |
| Trapecio simple | $-\dfrac{h^3}{12}f''(c)$ | $O(h^3)$ |
| Simpson simple | $-\dfrac{h^5}{90}f^{(4)}(c)$ | $O(h^5)$ |
| Trapecio compuesta | $-\dfrac{(b-a)h^2}{12}f''(c)$ | $O(h^2)$ |
| Simpson compuesta | $-\dfrac{(b-a)h^4}{180}f^{(4)}(c)$ | $O(h^4)$ |

---

## ⚠️ Errores típicos de examen

- **Lagrange:** el denominador de $L_i$ es un número fijo (sin $x$), solo depende de los nodos.
- **Newton:** los $c_k$ son la **primera fila** de cada columna, no la diagonal.
- **Error:** con $n+1$ nodos → derivada de orden $n+1$.
- **Simpson compuesta:** $m$ debe ser **PAR**.
- **Gauss:** multiplicar siempre por $\dfrac{b-a}{2}$ al cambiar de $[-1,1]$ a $[a,b]$.
- **Linealización:** transforma TODOS los $y_i$, no solo la fórmula.
- **Spline natural:** condición $S''=0$ en extremos, **no** $S'=0$.
- **Centrada:** denominador es $2h$, no $h$.

---

## Orden óptimo de problemas

| # | Prob. | Concepto | Prioridad |
|---|-------|----------|-----------|
| 1 | 3.1 | Lagrange básico — 3 nodos | 🔴 ESENCIAL |
| 2 | 3.4 | Interpolación lineal en datos reales | 🔴 ESENCIAL |
| 3 | 3.5 | Interpolación cuadrática | 🔴 ESENCIAL |
| 4 | 3.7 | Newton — diferencias divididas, tabla log | 🔴 ESENCIAL |
| 5 | 3.11 | Newton — datos en tabla | 🔴 ESENCIAL |
| 6 | 3.8 | Cota de error de interpolación | 🔴 ESENCIAL |
| 7 | 4.1 | Las 3 fórmulas de derivación + orden | 🔴 ESENCIAL |
| 8 | 4.8 | Trapecio simple | 🔴 ESENCIAL |
| 9 | 4.12 | Simpson  | 🔴 ESENCIAL | 
| 10 | 3.12 | Spline cúbico natural completo | 🟡 IMPORTANTE |
| 11 | 3.14 | Recta de regresión (mínimos cuadrados) | 🟡 IMPORTANTE |
| 12 | 3.15 | Linealización + mínimos cuadrados | 🟡 IMPORTANTE |
| 13 | 4.9 | Simpson y trapecio compuestas | 🟡 IMPORTANTE |
| 14 | 4.11 | Error compuestas — calcular nº subintervalos | 🟡 IMPORTANTE |
| 15 | 3.16 | Ajuste curva no lineal → linealización | 🟠 ÚTIL |
| 16 | 4.12 | Gauss-Legendre 2 puntos | 🟠 ÚTIL |
| 17 | 3.18 | Mínimos cuadrados con base polinómica | 🟠 ÚTIL |
| 18 | 3.20 | Aprox. continua con base $\{1,x^2\}$ | 🟠 ÚTIL |
| 19 | 4.3 | Derivadas parciales numéricas | 🟢 SI SOBRA |
| 20 | 3.24 | Desarrollo de Fourier discreto | 🟢 SI SOBRA |
