"""
05_multiplicacion_escalar.py
============================

La operacion "facil" de la criptografia de curvas elipticas es calcular kP;
la operacion "dificil" es, dados P y Q = kP, encontrar k (ECDLP).

Contenido
---------
1. Duplicar y sumar (analogo a los cuadrados sucesivos de 01).
2. Escalera de Montgomery: misma cantidad de operaciones para cualquier k
   (defensa contra ataques de canal lateral por tiempo).
3. Conteo de operaciones: unas 2 log2(k) en lugar de k - 1.
4. ECDLP por fuerza bruta y por "paso de bebe - paso de gigante"
   (O(sqrt(n)) operaciones), midiendo como crece el trabajo con n.

Resultado esperado de la escalera para k = 13 = 1101_2 sobre la curva de
trabajo (pares (R0, R1)):  (O, G) -> (G, 2G) -> (3G, 4G) -> (6G, 7G) -> (13G, 14G).

AVISO: codigo educativo; no es de tiempo constante a nivel de maquina.
"""
import math
import random
import time

INFINITO = None


class Curva:
    def __init__(self, a, b, p, contar=False):
        self.a, self.b, self.p = a, b, p
        self.operaciones = 0

    def suma(self, P, Q):
        self.operaciones += 1
        if P is INFINITO:
            return Q
        if Q is INFINITO:
            return P
        (x1, y1), (x2, y2), p = P, Q, self.p
        if x1 == x2 and (y1 + y2) % p == 0:
            return INFINITO
        if P == Q:
            m = (3 * x1 * x1 + self.a) * pow(2 * y1, -1, p) % p
        else:
            m = (y2 - y1) * pow(x2 - x1, -1, p) % p
        x3 = (m * m - x1 - x2) % p
        return (x3, (m * (x1 - x3) - y1) % p)

    def negativo(self, P):
        return INFINITO if P is INFINITO else (P[0], (-P[1]) % self.p)


# ---------------------------------------------------------------------------
# 1. Duplicar y sumar
# ---------------------------------------------------------------------------
def duplicar_y_sumar(E, k, P, mostrar=False):
    """Recorre los bits de k de izquierda a derecha: D = duplicar, S = sumar.
    El primer bit (siempre 1) solo copia P; luego, por cada bit: D, y S si el bit es 1."""
    if k == 0:
        return INFINITO
    R = P
    traza = []
    for bit in bin(k)[3:]:
        R = E.suma(R, R)
        traza.append("D")
        if bit == "1":
            R = E.suma(R, P)
            traza.append("S")
    if mostrar:
        print(f"  k = {k} = {bin(k)[2:]}_2  ->  operaciones: {' '.join(traza)}")
    return R


# ---------------------------------------------------------------------------
# 2. Escalera de Montgomery
# ---------------------------------------------------------------------------
def escalera_montgomery(E, k, P, mostrar=False):
    """Mantiene R1 - R0 = P. En cada bit se hace UNA suma y UNA duplicacion,
    sin importar si el bit es 0 o 1 (el patron de operaciones no revela k)."""
    R0, R1 = INFINITO, P
    pasos = [("inicio", R0, R1)]
    for bit in bin(k)[2:]:
        if bit == "0":
            R1 = E.suma(R0, R1)
            R0 = E.suma(R0, R0)
        else:
            R0 = E.suma(R0, R1)
            R1 = E.suma(R1, R1)
        pasos.append((bit, R0, R1))
    if mostrar:
        fmt = lambda P: "O" if P is INFINITO else f"({P[0]},{P[1]})"
        for bit, a, b in pasos:
            print(f"  {bit:>6}: R0 = {fmt(a)},  R1 = {fmt(b)}")
    return R0, pasos


# ---------------------------------------------------------------------------
# 4. ECDLP
# ---------------------------------------------------------------------------
def log_discreto_fuerza_bruta(E, G, Q, n):
    R = INFINITO
    for k in range(n):
        if R == Q:
            return k
        R = E.suma(R, G)
    return None


def paso_bebe_paso_gigante(E, G, Q, n):
    """Escribe k = i*m + j con m = ceil(sqrt(n)). Guarda jG (pasos de bebe)
    y recorre Q - i(mG) (pasos de gigante) hasta encontrar una coincidencia."""
    m = math.isqrt(n) + 1
    tabla = {}
    R = INFINITO
    for j in range(m):
        tabla.setdefault(R, j)
        R = E.suma(R, G)
    mG_neg = E.negativo(duplicar_y_sumar(E, m, G))
    gamma = Q
    for i in range(m):
        if gamma in tabla:
            return (i * m + tabla[gamma]) % n
        gamma = E.suma(gamma, mG_neg)
    return None


def orden_por_fuerza_bruta(E, P):
    n, Q = 1, P
    while Q is not INFINITO:
        Q = E.suma(Q, P)
        n += 1
    return n


def es_primo(n):
    return n > 1 and all(n % d for d in range(2, math.isqrt(n) + 1))


def curva_con_punto_de_orden_primo(bits, semilla=0):
    """Busca (por azar) una curva sobre un primo de 'bits' bits con un punto de
    orden primo grande. Solo para experimentos pequenos."""
    rnd = random.Random(semilla)
    while True:
        p = rnd.randrange(2 ** (bits - 1), 2 ** bits)
        if not es_primo(p) or p % 4 != 3:
            continue
        a = rnd.randrange(p)
        x = rnd.randrange(p)
        y = rnd.randrange(1, p)
        b = (y * y - x ** 3 - a * x) % p          # obligamos a que (x, y) este en la curva
        if (4 * a ** 3 + 27 * b * b) % p == 0:
            continue
        E = Curva(a, b, p)
        n = orden_por_fuerza_bruta(E, (x, y))
        if es_primo(n) and n > p // 2:
            return E, (x, y), n


# ---------------------------------------------------------------------------
# Pruebas
# ---------------------------------------------------------------------------
def pruebas():
    E = Curva(24, 13, 29)
    G = (1, 3)
    # Escalera para k = 13
    R, pasos = escalera_montgomery(E, 13, G)
    multiplo = {INFINITO: 0}
    Q = INFINITO
    for k in range(1, 20):
        Q = E.suma(Q, G)
        multiplo[Q] = k % 19
    assert [(multiplo[a], multiplo[b]) for _, a, b in pasos] == [(0, 1), (1, 2), (3, 4), (6, 7), (13, 14)]
    # Ambos metodos coinciden con la suma repetida
    for k in range(60):
        A = duplicar_y_sumar(E, k, G)
        B = escalera_montgomery(E, k, G)[0]
        assert A == B and multiplo[A] == k % 19
    # ECDLP en la curva de trabajo
    for k in range(19):
        Q = duplicar_y_sumar(E, k, G)
        assert log_discreto_fuerza_bruta(E, G, Q, 19) == k
        assert paso_bebe_paso_gigante(E, G, Q, 19) == k
    print("05_multiplicacion_escalar: todas las pruebas pasaron.")


if __name__ == "__main__":
    E = Curva(24, 13, 29)
    G = (1, 3)
    print("Duplicar y sumar con k = 13:")
    duplicar_y_sumar(E, 13, G, mostrar=True)
    print("Escalera de Montgomery con k = 13:")
    escalera_montgomery(E, 13, G, mostrar=True)

    print("\nCosto de calcular kP para k de 256 bits (conteo de sumas/duplicaciones):")
    k = random.getrandbits(256) | (1 << 255)
    E.operaciones = 0
    duplicar_y_sumar(E, k, G)
    print(f"  duplicar y sumar: {E.operaciones} operaciones (sumar G una por una: ~2^256)")

    print("\nECDLP: el trabajo crece como sqrt(n) incluso con el mejor metodo generico")
    print("  bits de p | orden n   | fuerza bruta (s) | paso bebe-gigante (s)")
    for bits in (10, 14, 18):
        E2, G2, n = curva_con_punto_de_orden_primo(bits, semilla=bits)
        k = n - 1                       # peor caso para la fuerza bruta
        Q = duplicar_y_sumar(E2, k, G2)
        t0 = time.perf_counter(); k1 = log_discreto_fuerza_bruta(E2, G2, Q, n); t1 = time.perf_counter()
        k2 = paso_bebe_paso_gigante(E2, G2, Q, n); t2 = time.perf_counter()
        assert k1 == k2 == k
        print(f"  {bits:9d} | {n:9d} | {t1 - t0:16.4f} | {t2 - t1:20.4f}")
    print()
    pruebas()
