"""
10_coppersmith_lll.py
=====================

Reduccion de reticulos (LLL) y el teorema de Coppersmith para encontrar
raices pequenas de polinomios modulo N. Aplicaciones a RSA con e = 3:

  A. Mensaje "estereotipado": se conoce todo el mensaje salvo un pedazo x.
  B. Ataque de Hastad con relleno lineal fijo (Boneh, Teorema 6).

Fundamento matematico
---------------------
Teorema (Coppersmith). Sea f en Z[x] monico de grado delta y X = N^(1/delta - eps).
Se pueden hallar eficientemente todos los x0 con |x0| < X y f(x0) = 0 (mod N).

Idea (version de Howgrave-Graham):
  1. Se construyen polinomios que se anulan en x0 modulo N^m:
        g_{i,j}(x) = x^j * N^(m-i) * f(x)^i     (0 <= i < m, 0 <= j < delta)
        h_j(x)     = x^j * f(x)^m               (0 <= j < t)
  2. Con sus coeficientes evaluados en xX se arma un reticulo (matriz
     triangular) y se reduce con LLL.
  3. El primer vector reducido da un polinomio h con coeficientes pequenos.
     Si ||h(xX)|| < N^m / sqrt(w), entonces |h(x0)| < N^m y, como
     h(x0) = 0 (mod N^m), en realidad h(x0) = 0 EN LOS ENTEROS.
  4. Se buscan las raices enteras de h (aqui: mcd de dos polinomios cortos,
     o busqueda por biseccion).

LLL: se usa la version entera de de Weger (Cohen, "A Course in Computational
Algebraic Number Theory", algoritmo 2.6.7), que solo opera con enteros.

Resultados esperados
--------------------
* LLL en el reticulo {(1, 0), (1000, 1)}... -> vectores cortos (ver pruebas).
* A: RSA de 512 bits, e = 3: con un reticulo de dimension 10 se recuperan
     120 bits desconocidos del mensaje (con 140 bits ya no alcanza: acercarse
     al limite teorico N^(1/3) ~ 170 bits exige reticulos mas grandes).
* B: tres destinatarios con relleno fijo M_i = M + i * 2^450 y M de 300 bits
     (modulos de 512 bits): se recupera M con un reticulo de dimension 7.

Equivalente en SageMath (para comparar):
    P.<x> = PolynomialRing(Zmod(N))
    f = (prefijo + x)^3 - c
    f.small_roots(X=2^120, beta=1)

La ejecucion completa tarda menos de un minuto (LLL en Python puro).

AVISO: codigo educativo. Los ataques solo funcionan contra RSA mal usado
(exponente pequeno sin relleno aleatorio).

Prueba con tus propios datos
----------------------------
    python 10_coppersmith_lll.py "Gato2026"
o ejecuta el archivo sin nada mas y escribe tus datos cuando te los pida.
Con  python 10_coppersmith_lll.py --ejemplo  se ve el ejemplo de clase y las pruebas.
"""
import math
import random
from fractions import Fraction


# ---------------------------------------------------------------------------
# Herramientas de RSA (para generar ejemplos)
# ---------------------------------------------------------------------------
def es_primo_probable(n, rondas=30):
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29):
        if n % p == 0:
            return n == p
    s, r = 0, n - 1
    while r % 2 == 0:
        s, r = s + 1, r // 2
    for _ in range(rondas):
        x = pow(random.randrange(2, n - 1), r, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True


def primo_aleatorio(bits):
    while True:
        c = random.getrandbits(bits) | (3 << (bits - 2)) | 1
        if es_primo_probable(c):
            return c


def modulo_rsa(bits, e=3):
    """N = p*q con mcd(e, phi(N)) = 1."""
    while True:
        p, q = primo_aleatorio(bits // 2), primo_aleatorio(bits // 2)
        if p != q and math.gcd(e, (p - 1) * (q - 1)) == 1:
            return p * q


# ---------------------------------------------------------------------------
# 1. LLL entero (de Weger / Cohen 2.6.7)
# ---------------------------------------------------------------------------
def lll(base, delta=(3, 4)):
    """Reduce una base de un reticulo (vectores enteros linealmente independientes).

    Mantiene, en enteros, los datos de Gram-Schmidt:
        d[i]       = det de la matriz de Gram de los primeros i vectores,
        lam[k][j]  = d[j+1] * mu_{k,j}.
    Condicion de Lovasz (con delta = 3/4): si
        ||b*_k||^2 < (delta - mu_{k,k-1}^2) ||b*_{k-1}||^2
    se intercambian b_k y b_{k-1}.
    """
    b = [list(v) for v in base]
    n = len(b)
    punto = lambda u, v: sum(x * y for x, y in zip(u, v))
    d = [0] * (n + 1)
    lam = [[0] * n for _ in range(n)]
    num, den = delta
    d[0], d[1] = 1, punto(b[0], b[0])
    k, kmax = 1, 0

    def reducir(k, l):
        if 2 * abs(lam[k][l]) > d[l + 1]:
            q = (2 * lam[k][l] + d[l + 1]) // (2 * d[l + 1])
            b[k] = [x - q * y for x, y in zip(b[k], b[l])]
            lam[k][l] -= q * d[l + 1]
            for i in range(l):
                lam[k][i] -= q * lam[l][i]

    def intercambiar(k):
        b[k], b[k - 1] = b[k - 1], b[k]
        for j in range(k - 1):
            lam[k][j], lam[k - 1][j] = lam[k - 1][j], lam[k][j]
        l = lam[k][k - 1]
        B = (d[k - 1] * d[k + 1] + l * l) // d[k]
        for i in range(k + 1, kmax + 1):
            t = lam[i][k]
            lam[i][k] = (d[k + 1] * lam[i][k - 1] - l * t) // d[k]
            lam[i][k - 1] = (B * t + l * lam[i][k]) // d[k + 1]
        d[k] = B

    while k < n:
        if k > kmax:                       # Gram-Schmidt incremental
            kmax = k
            for j in range(k + 1):
                u = punto(b[k], b[j])
                for i in range(j):
                    u = (d[i + 1] * u - lam[k][i] * lam[j][i]) // d[i]
                if j < k:
                    lam[k][j] = u
                else:
                    if u == 0:
                        raise ValueError("los vectores son linealmente dependientes")
                    d[k + 1] = u
        reducir(k, k - 1)
        if den * d[k + 1] * d[k - 1] < num * d[k] * d[k] - den * lam[k][k - 1] ** 2:
            intercambiar(k)
            k = max(1, k - 1)
        else:
            for l in range(k - 2, -1, -1):
                reducir(k, l)
            k += 1
    return b


# ---------------------------------------------------------------------------
# 2. Polinomios (listas de coeficientes, c[i] acompana a x^i)
# ---------------------------------------------------------------------------
def p_mul(a, b):
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            r[i + j] += x * y
    return r


def p_pot(a, k):
    r = [1]
    for _ in range(k):
        r = p_mul(r, a)
    return r


def p_eval(c, x):
    r = 0
    for coef in reversed(c):
        r = r * x + coef
    return r


def p_limpiar(c):
    c = list(c)
    while len(c) > 1 and c[-1] == 0:
        c.pop()
    return c


def p_es_cero(c):
    return len(c) == 1 and c[0] == 0


def p_resto(a, b):
    """Resto de dividir a entre b en Q[x]."""
    a = list(a)
    while len(a) >= len(b) and not p_es_cero(a):
        factor = a[-1] / b[-1]
        desplazamiento = len(a) - len(b)
        for i, coef in enumerate(b):
            a[i + desplazamiento] -= factor * coef
        a = p_limpiar(a[:-1]) if len(a) > 1 else [Fraction(0)]
    return a


def p_mcd_racional(a, b):
    """mcd monico de dos polinomios en Q[x] (algoritmo de Euclides, como en Z)."""
    a = p_limpiar([Fraction(x) for x in a])
    b = p_limpiar([Fraction(x) for x in b])
    while not p_es_cero(b):
        a, b = b, p_resto(a, b)
    return [x / a[-1] for x in a]


def raices_enteras(c, lo, hi):
    """Raices enteras de c en [lo, hi] por biseccion entre puntos criticos."""
    c = p_limpiar(c)
    if len(c) <= 1:
        return []
    if len(c) == 2:
        return [-c[0] // c[1]] if c[0] % c[1] == 0 and lo <= -c[0] // c[1] <= hi else []
    derivada = [i * c[i] for i in range(1, len(c))]
    criticos = _intervalos(derivada, lo, hi)
    puntos = sorted({lo, hi} | {t for t in criticos if lo <= t <= hi} |
                    {t + 1 for t in criticos if lo <= t + 1 <= hi})
    raices = set()
    for a, b in zip(puntos, puntos[1:]):
        fa, fb = p_eval(c, a), p_eval(c, b)
        for x, fx in ((a, fa), (b, fb)):
            if fx == 0:
                raices.add(x)
        if fa * fb < 0:
            while b - a > 1:
                m = (a + b) // 2
                fm = p_eval(c, m)
                if fm == 0:
                    raices.add(m)
                    break
                if (fm < 0) == (fa < 0):
                    a, fa = m, fm
                else:
                    b = m
    return sorted(raices)


def _intervalos(c, lo, hi):
    """Enteros t tales que hay una raiz real de c en [t, t+1]."""
    c = p_limpiar(c)
    if len(c) <= 1:
        return []
    derivada = [i * c[i] for i in range(1, len(c))]
    criticos = _intervalos(derivada, lo, hi) if len(c) > 2 else []
    puntos = sorted({lo, hi} | {t for t in criticos if lo <= t <= hi} |
                    {t + 1 for t in criticos if lo <= t + 1 <= hi})
    salida = set()
    for a, b in zip(puntos, puntos[1:]):
        fa, fb = p_eval(c, a), p_eval(c, b)
        if fa == 0:
            salida.add(a)
        if fb == 0:
            salida.add(b)
        if fa * fb < 0:
            while b - a > 1:
                m = (a + b) // 2
                fm = p_eval(c, m)
                if (fm < 0) == (fa < 0):
                    a, fa = m, fm
                else:
                    b = m
            salida.add(a)
    return sorted(salida)


# ---------------------------------------------------------------------------
# 3. Coppersmith / Howgrave-Graham
# ---------------------------------------------------------------------------
def coppersmith(f, N, X, m=2, t=1):
    """Raices x0 de f (monico, coeficientes mod N) con |x0| < X.

    f: lista de coeficientes [f0, f1, ..., f_delta] con f_delta = 1.
    Devuelve (raices, dimension del reticulo).
    """
    f = [x % N for x in f]
    delta = len(f) - 1
    assert f[-1] == 1, "f debe ser monico"
    polinomios = []
    for i in range(m):
        fi = p_pot(f, i)
        for j in range(delta):
            polinomios.append([0] * j + [coef * N ** (m - i) for coef in fi])
    fm = p_pot(f, m)
    for j in range(t):
        polinomios.append([0] * j + fm)
    w = len(polinomios)                              # dimension = delta*m + t
    base = []
    for g in polinomios:
        fila = [0] * w
        for k, coef in enumerate(g):
            fila[k] = coef * X ** k                  # coeficientes de g(xX)
        base.append(fila)
    reducida = lll(base)
    # Volvemos a polinomios en x dividiendo por X^k (division exacta)
    cortos = [[v[k] // X ** k for k in range(w)] for v in reducida[:3]]
    candidatos = set()
    g = p_mcd_racional(cortos[0], cortos[1])
    if len(g) == 2:                                   # mcd lineal: x + g0
        r = -g[0]
        if r.denominator == 1:
            candidatos.add(int(r))
    if not candidatos:
        candidatos.update(raices_enteras(cortos[0], -X, X))
    raices = sorted(r for r in candidatos if abs(r) < X and p_eval(f, r) % N == 0)
    return raices, w


# ---------------------------------------------------------------------------
# 4. Aplicaciones
# ---------------------------------------------------------------------------
def ataque_mensaje_estereotipado(N, c, prefijo, bits_desconocidos, m=3, t=1):
    """M = prefijo * 2^k + x con x < 2^k desconocido y c = M^3 mod N."""
    P = prefijo << bits_desconocidos
    # f(x) = (P + x)^3 - c = x^3 + 3P x^2 + 3P^2 x + P^3 - c
    f = [(P ** 3 - c) % N, (3 * P * P) % N, (3 * P) % N, 1]
    raices, w = coppersmith(f, N, 2 ** bits_desconocidos, m=m, t=t)
    return (P + raices[0] if raices else None), w


def ataque_hastad_relleno(modulos, cifrados, rellenos, X, m=2, t=1):
    """Cada destinatario recibe (a_i M + b_i)^3 mod N_i.

    Se hace monico cada g_i(x) = (a_i x + b_i)^3 - c_i y se combinan con los
    coeficientes del TCR T_i (T_i = 1 mod N_i, 0 mod N_j) en un solo polinomio
    G(x) = sum T_i g_i(x)  modulo N = N_1 N_2 N_3, que tiene la raiz pequena M.
    """
    N = math.prod(modulos)
    G = [0, 0, 0, 0]
    for Ni, ci, (ai, bi) in zip(modulos, cifrados, rellenos):
        gi = p_mul(p_mul([bi, ai], [bi, ai]), [bi, ai])
        gi[0] -= ci
        inv = pow(ai ** 3, -1, Ni)                    # hacer monico modulo N_i
        gi = [(coef * inv) % Ni for coef in gi]
        Mi = N // Ni
        Ti = Mi * pow(Mi, -1, Ni)                     # 1 mod N_i, 0 mod N_j
        G = [(Gk + Ti * gk) % N for Gk, gk in zip(G, gi)]
    raices, w = coppersmith(G, N, X, m=m, t=t)
    return (raices[0] if raices else None), w


# ---------------------------------------------------------------------------
# Pruebas
# ---------------------------------------------------------------------------
def pruebas():
    # LLL en un reticulo de dimension 2 (resultado clasico: base casi ortogonal)
    red = lll([[1, 0], [1000, 1]])
    assert sorted(abs(x) for v in red for x in v) == [0, 0, 1, 1]
    # LLL encuentra una relacion entera: 3*a - 2*b = 0 con a = 2, b = 3
    red = lll([[1, 0, 2 * 10 ** 6], [0, 1, 3 * 10 ** 6]])
    assert [abs(x) for x in red[0]] == [3, 2, 0]
    # Raices enteras de un polinomio
    assert raices_enteras(p_mul(p_mul([-5, 1], [7, 1]), [-123456789, 1]), -10 ** 12, 10 ** 12) == [-7, 5, 123456789]
    random.seed(10)
    # A. Mensaje estereotipado: RSA de 512 bits, e = 3, 120 bits desconocidos
    N = modulo_rsa(512)
    prefijo = int.from_bytes(b"La clave de hoy es: ", "big")
    secreto = random.getrandbits(120)
    M = (prefijo << 120) + secreto
    c = pow(M, 3, N)
    recuperado, _ = ataque_mensaje_estereotipado(N, c, prefijo, 120)
    assert recuperado == M
    # B. Hastad con relleno lineal: tres modulos de 512 bits, M de 300 bits
    modulos = [modulo_rsa(512) for _ in range(3)]
    M = random.getrandbits(300)
    rellenos = [(1, i << 450) for i in (1, 2, 3)]       # M_i = M + i * 2^450
    cifrados = [pow(a * M + b, 3, Ni) for Ni, (a, b) in zip(modulos, rellenos)]
    assert ataque_hastad_relleno(modulos, cifrados, rellenos, 2 ** 300)[0] == M
    print("10_coppersmith_lll: todas las pruebas pasaron.")


def ejemplo_de_clase():
    """Lo que se muestra sin datos: el ejemplo de clase del informe."""
    print("LLL sobre la base {(1, 0, 2000000), (0, 1, 3000000)}:")
    print("  ", lll([[1, 0, 2 * 10 ** 6], [0, 1, 3 * 10 ** 6]]),
          " (el primer vector expresa la relacion 3*2000000 - 2*3000000 = 0)")
    random.seed(2026)
    N = modulo_rsa(512)
    prefijo = int.from_bytes(b"La clave de hoy es: ", "big")
    for k in (120, 140):
        secreto = random.getrandbits(k)
        M = (prefijo << k) + secreto
        c = pow(M, 3, N)
        r, w = ataque_mensaje_estereotipado(N, c, prefijo, k)
        print(f"A. RSA-512, e = 3, {k} bits desconocidos (N^(1/3) ~ 2^{N.bit_length() // 3}):"
              f" dimension {w}, {'RECUPERADO' if r == M else 'no recuperado'}")
    modulos = [modulo_rsa(512) for _ in range(3)]
    M = random.getrandbits(300)
    rellenos = [(1, i << 450) for i in (1, 2, 3)]
    cifrados = [pow(a * M + b, 3, Ni) for Ni, (a, b) in zip(modulos, rellenos)]
    r, w = ataque_hastad_relleno(modulos, cifrados, rellenos, 2 ** 300)
    print(f"B. Hastad con relleno fijo M_i = M + i*2^450, tres modulos de 512 bits:"
          f" dimension {w}, {'RECUPERADO' if r == M else 'no recuperado'}")


# ---------------------------------------------------------------------------
# Prueba con tus propios datos
# ---------------------------------------------------------------------------
def leer_datos(pregunta):
    """Devuelve los datos que escribe el usuario (texto) o "" si no hay.

    * python ARCHIVO.py tus datos   -> usa lo que va despues del nombre del archivo.
    * python ARCHIVO.py             -> si hay teclado, hace la pregunta; con solo
                                       Enter se muestra el ejemplo de clase.
    * python ARCHIVO.py --ejemplo   -> ejemplo de clase y pruebas, sin preguntar.
    """
    import builtins
    import os
    import sys
    argumentos = sys.argv[1:]
    # Si el codigo se pego en una celda de Jupyter o Colab, sys.argv trae los
    # datos del "kernel" (-f archivo.json), que no son del usuario.
    if "launcher" in os.path.basename(sys.argv[0] if sys.argv else "") or argumentos[:1] == ["-f"]:
        argumentos = []
    if argumentos:
        return "" if argumentos == ["--ejemplo"] else " ".join(argumentos)
    hay_teclado = getattr(builtins, "__IPYTHON__", False) or (
        sys.stdin is not None and sys.stdin.isatty())
    if not hay_teclado:
        return ""
    try:
        return input(pregunta).strip()
    except EOFError:
        return ""


def corto(n, cifras=12):
    """Muestra un entero muy largo como 'primeras...ultimas cifras (k cifras)'."""
    s = str(n)
    if len(s) <= 2 * cifras + 3:
        return s
    return f"{s[:cifras]}...{s[-cifras:]} ({len(s)} cifras)"


def leer_enteros(texto):
    """Convierte '240 46 13' o '240, 46, 13' en la lista [240, 46, 13]."""
    return [int(t) for t in texto.replace(",", " ").split()]


PREGUNTA = ("Escribe una clave secreta corta (hasta 15 caracteres) que Eve intentara\n"
            "descubrir con el metodo de Coppersmith. Solo Enter = ejemplo de clase: ")
USO = 'python 10_coppersmith_lll.py "Gato2026"'


def con_tus_datos(secreto):
    import time
    s = secreto.encode("utf-8")
    if len(s) > 15:
        print("La clave tiene mas de 15 bytes (120 bits). Con el reticulo de dimension 10 de")
        print("este programa el metodo alcanza unos 120 bits desconocidos (el limite teorico")
        print("para RSA-512 con e = 3 es cerca de 170 bits). Usa una clave mas corta, por ejemplo:")
        print("  " + USO)
        return
    texto_fijo = b"La clave de hoy es: "
    k = 8 * len(s)
    prefijo = int.from_bytes(texto_fijo, "big")
    x0 = int.from_bytes(s, "big")
    M = (prefijo << k) + x0
    print(f'Mensaje que Alice cifra: "{(texto_fijo + s).decode("utf-8", errors="replace")}"')
    print(f"Eve conoce la parte fija \"{texto_fijo.decode()}\" pero no tu clave ({k} bits desconocidos).")

    N = modulo_rsa(512)
    c = pow(M, 3, N)
    print("\n1. Lo que viaja por la red (RSA-512 con e = 3, sin relleno):")
    print(f"   N = {corto(N)}")
    print(f"   c = M^3 mod N = {corto(c)}")

    print("\n2. Eve escribe M = P + x, con P = (parte fija)*2^k conocido y x < 2^k desconocido:")
    print("   f(x) = (P + x)^3 - c = 0 (mod N), un polinomio de grado 3 con una raiz PEQUENA.")
    print(f"   Coppersmith: arma un reticulo con multiplos de f y potencias de N, lo reduce con LLL")
    print("   y obtiene un polinomio con la misma raiz pero SIN modulo, que se resuelve normalmente.")
    print("   Se prueba primero con un reticulo pequeno y se agranda si no alcanza:")
    recuperado = None
    for m_ in (1, 2, 3):
        t = time.perf_counter()
        recuperado, w = ataque_mensaje_estereotipado(N, c, prefijo, k, m=m_)
        t = time.perf_counter() - t
        print(f"   dimension {w:2d}: {'raiz encontrada' if recuperado is not None else 'no alcanza'}"
              f"  ({t:.1f} s)", flush=True)
        if recuperado is not None:
            break
    if recuperado is None:
        print("   No se encontro la raiz (prueba con una clave mas corta).")
        return
    x = recuperado - (prefijo << k)
    descubierta = x.to_bytes(len(s), "big").decode("utf-8", errors="replace")
    print(f"\n3. Raiz encontrada: x = {x}  ->  \"{descubierta}\"")
    print(f"   Es tu clave? {recuperado == M}")
    print(f"   Probar todas las claves posibles habria exigido hasta 2^{k} = {corto(2 ** k)} intentos.")
    print("\nDefensa: usar relleno aleatorio (OAEP); asi no hay partes conocidas del mensaje.")


if __name__ == "__main__":
    datos = leer_datos(PREGUNTA)
    if datos:
        con_tus_datos(datos)
        print("\n(Para ver el ejemplo de clase y las pruebas:  python 10_coppersmith_lll.py --ejemplo)")
    else:
        ejemplo_de_clase()
        print()
        pruebas()
