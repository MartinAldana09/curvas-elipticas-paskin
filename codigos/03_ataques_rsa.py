"""
03_ataques_rsa.py
=================

Ataques clasicos a RSA estudiados en Boneh, "Twenty Years of Attacks on the
RSA Cryptosystem" (Notices of the AMS, 1999) y en Trappe-Washington, cap. 6.

Contenido
---------
A. Conocer d permite factorizar N (Hecho 1 de Boneh).
B. Modulo comun.
C. Firma "ciega" (blinding).
D. Ataque de Wiener (exponente privado pequeno, fracciones continuas).
E. Ataque de difusion de Hastad (e = 3, Teorema Chino del Resto).
F. Metodo p - 1 de Pollard (factorizacion).
G. Ataque por fallo en el descifrado con TCR (Boneh-DeMillo-Lipton).

Ninguno de estos ataques "rompe" RSA bien usado: todos aprovechan una mala
eleccion de parametros o una implementacion descuidada.

AVISO: codigo educativo, para entender por que existen las reglas de uso
seguro de RSA. No debe usarse contra sistemas ajenos.
"""
import math
import random

# ---------------------------------------------------------------------------
# Herramientas
# ---------------------------------------------------------------------------
def euclides_extendido(a, b):
    r0, r1, x0, x1, y0, y1 = a, b, 1, 0, 0, 1
    while r1:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        x0, x1 = x1, x0 - q * x1
        y0, y1 = y1, y0 - q * y1
    return r0, x0, y0


def inverso_modular(a, n):
    g, x, _ = euclides_extendido(a % n, n)
    if g != 1:
        raise ValueError("no invertible")
    return x % n


def tcr(residuos, modulos):
    M = math.prod(modulos)
    x = 0
    for a_i, m_i in zip(residuos, modulos):
        M_i = M // m_i
        x += a_i * M_i * inverso_modular(M_i, m_i)
    return x % M


def raiz_entera(n, k):
    """Mayor entero r con r^k <= n (metodo de Newton con enteros)."""
    if n < 2:
        return n
    r = 1 << ((n.bit_length() + k - 1) // k)      # cota superior inicial
    while True:
        s = ((k - 1) * r + n // r ** (k - 1)) // k
        if s >= r:
            break
        r = s
    while r ** k > n:
        r -= 1
    while (r + 1) ** k <= n:
        r += 1
    return r


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
        c = random.getrandbits(bits) | (1 << (bits - 1)) | 1
        if es_primo_probable(c):
            return c


def claves_rsa(bits=512, e=65537):
    while True:
        p, q = primo_aleatorio(bits // 2), primo_aleatorio(bits // 2)
        phi = (p - 1) * (q - 1)
        if p != q and math.gcd(e, phi) == 1:
            return p * q, e, inverso_modular(e, phi), p, q


# ---------------------------------------------------------------------------
# A. Conocer d  ==>  factorizar N
# ---------------------------------------------------------------------------
def factorizar_con_d(N, e, d, intentos=100):
    """Factoriza N conociendo (e, d).

    k = ed - 1 es multiplo de phi(N). Se escribe k = 2^t * r con r impar.
    Para g al azar, la sucesion g^r, g^{2r}, ..., g^{2^t r} termina en 1.
    El ultimo termino distinto de 1, llamemoslo x, cumple x^2 = 1 (mod N).
    Si ademas x != -1, es una raiz cuadrada NO trivial de 1 y
    mcd(x - 1, N) es un factor propio de N (por el TCR hay 4 raices de 1).
    """
    k = e * d - 1
    t, r = 0, k
    while r % 2 == 0:
        t, r = t + 1, r // 2
    for _ in range(intentos):
        g = random.randrange(2, N - 1)
        x = pow(g, r, N)
        if x in (1, N - 1):
            continue
        for _ in range(t):
            y = pow(x, 2, N)
            if y == 1:
                if x != N - 1:
                    p = math.gcd(x - 1, N)
                    return p, N // p
                break
            if y == N - 1:
                break
            x = y
    return None


# ---------------------------------------------------------------------------
# B. Modulo comun
# ---------------------------------------------------------------------------
def ataque_modulo_comun(N, e1, e2, c1, c2):
    """Recupera m si el MISMO mensaje se cifro con (N, e1) y (N, e2), mcd(e1,e2)=1.

    Bezout: a*e1 + b*e2 = 1  ==>  c1^a * c2^b = m^{a e1 + b e2} = m (mod N).
    Un exponente negativo se maneja con el inverso modular.
    """
    g, a, b = euclides_extendido(e1, e2)
    assert g == 1
    parte1 = pow(c1, a, N) if a >= 0 else pow(inverso_modular(c1, N), -a, N)
    parte2 = pow(c2, b, N) if b >= 0 else pow(inverso_modular(c2, N), -b, N)
    return parte1 * parte2 % N


# ---------------------------------------------------------------------------
# C. Firma ciega (blinding)
# ---------------------------------------------------------------------------
def firma_ciega(N, e, d, M):
    """Muestra que firmar M' = r^e M permite obtener la firma de M.

    S' = (r^e M)^d = r^{ed} M^d = r * S (mod N),  luego S = S' * r^{-1}.
    Como ataque: alguien logra que se firme M' (que parece inocente) y
    obtiene la firma de M. Como defensa: el firmante ciega internamente
    sus calculos para que tiempo y consumo no dependan de M.
    """
    while True:
        r = random.randrange(2, N - 1)
        if math.gcd(r, N) == 1:
            break
    M_ciego = pow(r, e, N) * M % N
    S_ciega = pow(M_ciego, d, N)            # lo unico que calcula el firmante
    S = S_ciega * inverso_modular(r, N) % N
    return S


# ---------------------------------------------------------------------------
# D. Ataque de Wiener
# ---------------------------------------------------------------------------
def fraccion_continua(num, den):
    coef = []
    while den:
        q = num // den
        coef.append(q)
        num, den = den, num - q * den
    return coef


def convergentes(coef):
    h_ant, h, k_ant, k = 0, 1, 1, 0
    for a in coef:
        h_ant, h = h, a * h + h_ant
        k_ant, k = k, a * k + k_ant
        yield h, k


def wiener(N, e, mostrar=False):
    """Recupera (p, q, d) si d es pequeno (Wiener: d < N^{1/4}/3 garantiza exito).

    De ed - k phi(N) = 1 se obtiene |e/N - k/d| < 1/(2 d^2), y por la teoria
    de fracciones continuas k/d es un convergente de e/N. Para cada
    convergente k/d se prueba phi = (ed - 1)/k y se resuelve
        x^2 - (N - phi + 1) x + N = 0,
    cuyas raices serian p y q.
    """
    for i, (k, d) in enumerate(convergentes(fraccion_continua(e, N))):
        if k == 0 or (e * d - 1) % k:
            continue
        phi = (e * d - 1) // k
        s = N - phi + 1                       # s = p + q
        disc = s * s - 4 * N
        if disc < 0:
            continue
        raiz = math.isqrt(disc)
        if raiz * raiz == disc:
            p, q = (s + raiz) // 2, (s - raiz) // 2
            if p * q == N:
                if mostrar:
                    print(f"  exito en el convergente numero {i}: k = {k}, d = {d}")
                return p, q, d
    return None


# ---------------------------------------------------------------------------
# E. Ataque de Hastad (difusion con e pequeno)
# ---------------------------------------------------------------------------
def hastad(cifrados, modulos, e=3):
    """Si el mismo M se cifra con exponente e para e destinatarios,
    el TCR da C = M^e (mod N1...Ne). Como M^e < N1...Ne, C = M^e en Z
    y basta sacar la raiz e-esima entera. Hacen falta al menos e cifrados.
    """
    if len(cifrados) < e:
        raise ValueError(f"se necesitan al menos {e} cifrados")
    C = tcr(cifrados[:e], modulos[:e])
    M = raiz_entera(C, e)
    return M if M ** e == C else None


# ---------------------------------------------------------------------------
# F. Metodo p - 1 de Pollard
# ---------------------------------------------------------------------------
def pollard_p_menos_1(n, B=100, a=2):
    """Encuentra un factor p de n si p - 1 solo tiene factores primos <= B.

    Si (p - 1) | B!, por Fermat a^{B!} = 1 (mod p), asi que p | mcd(a^{B!} - 1, n).
    """
    b = a
    for j in range(2, B + 1):
        b = pow(b, j, n)
    g = math.gcd(b - 1, n)
    return g if 1 < g < n else None


# ---------------------------------------------------------------------------
# G. Fallo en el descifrado con TCR
# ---------------------------------------------------------------------------
def ataque_por_fallo(N, e, d, p, q, M):
    """Simula una firma S = M^d calculada con TCR donde falla la parte modulo q.

    La firma defectuosa S' es correcta modulo p pero no modulo q, asi que
    p | (S'^e - M) y q no lo divide:  mcd(S'^e - M, N) = p.
    """
    S_p = pow(M, d % (p - 1), p)
    S_q = (pow(M, d % (q - 1), q) + 1) % q      # <- error de "hardware"
    h = inverso_modular(q, p) * (S_p - S_q) % p
    S_defectuosa = S_q + h * q
    return math.gcd(pow(S_defectuosa, e, N) - M, N)


# ---------------------------------------------------------------------------
# Pruebas con datos fijos (resultados verificados)
# ---------------------------------------------------------------------------
EJEMPLO_WIENER = dict(N=1966981193543797, e=323815174542919)
EJEMPLO_HASTAD = dict(
    modulos=[2469247531693, 11111502225583, 44444222221411],
    cifrados=[359335245251, 10436363975495, 5135984059593],
)
EJEMPLO_POLLARD = 376875575426394855599989992897873239


def pruebas():
    random.seed(7)
    # A
    N, e, d, p, q = claves_rsa(512)
    assert set(factorizar_con_d(N, e, d)) == {p, q}
    # B
    e1, e2 = 65537, 17
    m = 123456789
    assert ataque_modulo_comun(N, e1, e2, pow(m, e1, N), pow(m, e2, N)) == m
    # C
    M = 424242
    assert firma_ciega(N, e, d, M) == pow(M, d, N)
    # D: ejemplo fijo; la clave verdadera es d = 3511
    assert wiener(**EJEMPLO_WIENER) == (52783789, 37264873, 3511)
    #    y una clave generada con d pequeno a proposito
    while True:
        p, q = primo_aleatorio(256), primo_aleatorio(256)
        if q < p < 2 * q:
            break
    N, phi = p * q, (p - 1) * (q - 1)
    while True:
        d = random.randrange(3, raiz_entera(N, 4) // 3)
        if math.gcd(d, phi) == 1:
            break
    e = inverso_modular(d, phi)
    assert wiener(N, e)[2] == d
    # E: los modulos son productos de dos primos; el mensaje es 805121215
    assert hastad(**EJEMPLO_HASTAD) == 805121215
    #    con claves generadas
    msg = 2026
    mods = [claves_rsa(512, e=3)[0] for _ in range(3)]
    assert hastad([pow(msg, 3, n) for n in mods], mods, 3) == msg
    # F: el factor encontrado es 430553161739796481
    f = pollard_p_menos_1(EJEMPLO_POLLARD, B=100)
    assert f == 430553161739796481 and EJEMPLO_POLLARD % f == 0
    # G
    N, e, d, p, q = claves_rsa(512)
    assert ataque_por_fallo(N, e, d, p, q, 99999) == p
    print("03_ataques_rsa: todas las pruebas pasaron.")


if __name__ == "__main__":
    print("D. Wiener sobre N =", EJEMPLO_WIENER["N"])
    p, q, d = wiener(**EJEMPLO_WIENER, mostrar=True)
    print(f"  p = {p}, q = {q}, d = {d}")
    print(f"  (notese que d = {d} > N^(1/4)/3 = {EJEMPLO_WIENER['N'] ** 0.25 / 3:.1f}: "
          "la cota de Wiener es suficiente, no necesaria)")
    print("E. Hastad con e = 3:")
    print("  mensaje recuperado:", hastad(**EJEMPLO_HASTAD))
    print("F. Pollard p - 1 con B = 100:")
    f = pollard_p_menos_1(EJEMPLO_POLLARD, B=100)
    print(f"  factor = {f},  cofactor = {EJEMPLO_POLLARD // f}")
    print()
    pruebas()
