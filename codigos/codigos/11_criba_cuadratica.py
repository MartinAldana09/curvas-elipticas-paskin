"""
11_criba_cuadratica.py
======================

Criba cuadratica basica (Pomerance, 1981) para factorizar enteros.

Idea: congruencia de cuadrados
------------------------------
Si encontramos X, Y con  X^2 = Y^2 (mod N)  y  X != +-Y (mod N), entonces
N divide a (X - Y)(X + Y) sin dividir a ninguno de los dos factores, y
mcd(X - Y, N) es un factor propio de N (es la misma idea de la
proposicion sobre raices de 1 modulo N = pq).

Como encontrar X e Y
--------------------
1. Se toman x cerca de sqrt(N) y se calcula Q(x) = x^2 - N, que es "pequeno".
2. Se buscan los Q(x) que se factorizan por completo con primos pequenos
   (la "base de factores"). Solo sirven primos p con N cuadrado modulo p,
   porque p | Q(x) significa x^2 = N (mod p).
3. Cada Q(x) suave da un vector de exponentes modulo 2. Con mas relaciones
   que primos, el algebra lineal sobre F_2 encuentra un subconjunto cuyo
   producto tiene todos los exponentes pares: un cuadrado perfecto Y^2.
   Entonces X = prod(x) cumple X^2 = Y^2 (mod N).
4. La "criba": en lugar de factorizar cada Q(x), para cada primo p se
   marcan de un golpe todos los x con p | Q(x) (son x = +-r mod p, donde
   r^2 = N mod p), sumando log p. Los x con suma grande son candidatos.

Ejemplo pequeno (resultado esperado)
------------------------------------
N = 1649:  41^2 - N = 32 = 2^5,   43^2 - N = 200 = 2^3 * 5^2.
32 * 200 = 6400 = 80^2,  X = 41 * 43 = 1763 = 114 (mod N),
mcd(114 - 80, 1649) = 17,  y 1649 = 17 * 97.

AVISO: version didactica, sin las optimizaciones de las implementaciones
reales (polinomios multiples, variacion de primos grandes, etc.).
"""
import math


# ---------------------------------------------------------------------------
# Herramientas
# ---------------------------------------------------------------------------
def criba_eratostenes(n):
    es = bytearray([1]) * (n + 1)
    es[0:2] = b"\x00\x00"
    for i in range(2, math.isqrt(n) + 1):
        if es[i]:
            es[i * i::i] = bytearray(len(es[i * i::i]))
    return [i for i in range(n + 1) if es[i]]


def legendre(a, p):
    """Simbolo de Legendre (a/p) por el criterio de Euler."""
    r = pow(a % p, (p - 1) // 2, p)
    return -1 if r == p - 1 else r


def raiz_modular(n, p):
    """Tonelli-Shanks: r con r^2 = n (mod p), p primo impar, (n/p) = 1."""
    n %= p
    if p == 2:
        return n
    if p % 4 == 3:
        return pow(n, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0:
        q, s = q // 2, s + 1
    z = 2
    while legendre(z, p) != -1:
        z += 1
    m, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2, i = t2 * t2 % p, i + 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c, t, r = i, b * b % p, t * b * b % p, r * b % p
    return r


# ---------------------------------------------------------------------------
# Paso 1: base de factores
# ---------------------------------------------------------------------------
def base_de_factores(N, B):
    """-1 y los primos p <= B con N cuadrado modulo p."""
    base = [-1, 2]
    for p in criba_eratostenes(B)[1:]:
        if legendre(N, p) == 1:
            base.append(p)
    return base


def factorizar_en_base(q, base):
    """Vector de exponentes de q sobre la base, o None si q no es suave."""
    exps = [0] * len(base)
    if q < 0:
        exps[0], q = 1, -q
    for i, p in enumerate(base[1:], start=1):
        while q % p == 0:
            q //= p
            exps[i] += 1
    return exps if q == 1 else None


# ---------------------------------------------------------------------------
# Paso 2: criba
# ---------------------------------------------------------------------------
def cribar(N, base, M):
    """Busca relaciones (x, Q(x), exponentes) con x en [sqrt(N) - M, sqrt(N) + M]."""
    raiz = math.isqrt(N)
    inicio = raiz - M
    logs = [0.0] * (2 * M + 1)
    for p in base[1:]:
        lp = math.log(p)
        if p == 2:
            raices = [x for x in (0, 1) if (x * x - N) % 2 == 0]
        else:
            r = raiz_modular(N, p)
            raices = [r, p - r]
        for r in raices:
            primero = (r - inicio) % p
            for j in range(primero, 2 * M + 1, p):
                logs[j] += lp
    relaciones = []
    for j, valor in enumerate(logs):
        x = inicio + j
        q = x * x - N
        if q == 0:
            continue
        if valor > math.log(abs(q)) - 2.5 * math.log(base[-1]):   # umbral tolerante
            exps = factorizar_en_base(q, base)
            if exps is not None:
                relaciones.append((x, q, exps))
    return relaciones


# ---------------------------------------------------------------------------
# Paso 3: algebra lineal sobre F_2
# ---------------------------------------------------------------------------
def dependencias_F2(vectores):
    """Eliminacion gaussiana en F_2. Cada vector se guarda como un entero
    (bit i = paridad del exponente del primo i). Se lleva registro de que
    relaciones se combinaron; cuando una fila se anula, esa combinacion es
    un cuadrado perfecto."""
    filas = []
    for idx, v in enumerate(vectores):
        bits = sum(1 << i for i, e in enumerate(v) if e % 2)
        filas.append([bits, 1 << idx])
    pivotes = {}
    dependencias = []
    for fila in filas:
        bits, combo = fila
        while bits:
            alto = bits.bit_length() - 1
            if alto in pivotes:
                pb, pc = pivotes[alto]
                bits, combo = bits ^ pb, combo ^ pc
            else:
                pivotes[alto] = (bits, combo)
                break
        if bits == 0:
            dependencias.append(combo)
    return dependencias


# ---------------------------------------------------------------------------
# Paso 4: congruencia de cuadrados
# ---------------------------------------------------------------------------
def factor_desde_dependencia(N, relaciones, base, combo):
    X, exps_total = 1, [0] * len(base)
    for i, (x, q, exps) in enumerate(relaciones):
        if combo >> i & 1:
            X = X * x % N
            exps_total = [a + b for a, b in zip(exps_total, exps)]
    Y = 1
    for p, e in zip(base[1:], exps_total[1:]):
        Y = Y * pow(p, e // 2, N) % N
    g = math.gcd(X - Y, N)
    return g if 1 < g < N else None


def criba_cuadratica(N, B=None, M=None, mostrar=False):
    if N % 2 == 0:
        return 2
    r = math.isqrt(N)
    if r * r == N:
        return r
    if B is None:                                     # parametros heuristicos
        L = math.exp(math.sqrt(math.log(N) * math.log(math.log(N))))
        B = max(30, int(L ** 0.6))
    base = base_de_factores(N, B)
    if M is None:
        M = 5 * B
    # Se agranda el intervalo hasta tener mas relaciones que elementos en la base
    while True:
        relaciones = cribar(N, base, M)
        if len(relaciones) > len(base) + 5 or M > 4 * 10 ** 6:
            break
        M *= 2
    if mostrar:
        print(f"  base de factores: {len(base)} elementos (primos hasta {B});"
              f" intervalo de {2 * M + 1} valores; relaciones: {len(relaciones)}")
    for combo in dependencias_F2([e for _, _, e in relaciones]):
        f = factor_desde_dependencia(N, relaciones, base, combo)
        if f:
            return f
    return None


# ---------------------------------------------------------------------------
# Pruebas
# ---------------------------------------------------------------------------
def pruebas():
    assert raiz_modular(10, 13) ** 2 % 13 == 10
    assert raiz_modular(2, 41) ** 2 % 41 == 2          # 41 = 1 mod 8: caso dificil
    # Ejemplo pequeno: 1649
    assert (41 * 41 - 1649, 43 * 43 - 1649) == (32, 200)
    assert math.gcd(41 * 43 % 1649 - 80, 1649) == 17
    assert criba_cuadratica(1649) in (17, 97)
    # Numeros mas grandes (semiprimos)
    for p, q in [(1000003, 1000033), (2147483647, 4294967291), (999999999989, 1000000000039)]:
        f = criba_cuadratica(p * q)
        assert f in (p, q), (p, q, f)
    print("11_criba_cuadratica: todas las pruebas pasaron.")


if __name__ == "__main__":
    N = 1649
    print("Ejemplo pequeno, N = 1649:")
    for x in (41, 42, 43):
        print(f"  {x}^2 - N = {x * x - N}")
    X, Y = 41 * 43 % N, 80
    print(f"  32 * 200 = 6400 = 80^2;  X = 41*43 mod N = {X};  mcd(X - Y, N) = {math.gcd(X - Y, N)}")
    print()
    import time
    for p, q in [(2147483647, 4294967291), (999999999989, 1000000000039)]:
        N = p * q
        t = time.perf_counter()
        print(f"N = {N} ({len(str(N))} cifras)")
        f = criba_cuadratica(N, mostrar=True)
        print(f"  factor: {f}, cofactor: {N // f}  ({time.perf_counter() - t:.2f} s)")
    print()
    pruebas()
