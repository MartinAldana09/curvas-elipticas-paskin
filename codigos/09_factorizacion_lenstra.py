"""
09_factorizacion_lenstra.py
===========================

Metodo de curvas elipticas de Lenstra (ECM) para factorizar enteros.

Idea
----
Se hacen las formulas de la suma "modulo n" como si n fuera primo. Para sumar
hay que invertir un denominador modulo n. Si en algun momento ese
denominador NO es invertible, entonces mcd(denominador, n) es un factor de n.
Por el TCR, la curva modulo n = p*q se comporta como dos curvas: una modulo p
y otra modulo q. Si los multiplos de P llegan a O antes modulo p que modulo q,
el denominador se anula modulo p pero no modulo q: el mcd separa a p.

Es la version eliptica del metodo p - 1 de Pollard (ver 03_ataques_rsa.py):
alli se necesita que p - 1 sea "suave"; aqui basta que #E(F_p) lo sea para
ALGUNA curva, y si una curva no sirve se cambia por otra.

Resultados esperados (Trappe-Washington, seccion 16.3):
    n = 2773,  E: y^2 = x^3 + 4x + 4,  P = (1, 3):  2P = (1771, 705),
    al calcular 3P hay que invertir 1770 y mcd(1770, 2773) = 59;  2773 = 59 * 47.

AVISO: version didactica; hay implementaciones mucho mas eficientes (GMP-ECM).
"""
import math
import random


class FactorEncontrado(Exception):
    def __init__(self, factor):
        self.factor = factor


INFINITO = None


def inverso_o_factor(a, n):
    g = math.gcd(a % n, n)
    if g == 1:
        return pow(a, -1, n)
    raise FactorEncontrado(g)


def suma_mod_n(P, Q, a, n):
    if P is INFINITO:
        return Q
    if Q is INFINITO:
        return P
    (x1, y1), (x2, y2) = P, Q
    if x1 == x2 and (y1 + y2) % n == 0:
        return INFINITO
    if P == Q:
        m = (3 * x1 * x1 + a) * inverso_o_factor(2 * y1, n) % n
    else:
        m = (y2 - y1) * inverso_o_factor(x2 - x1, n) % n
    x3 = (m * m - x1 - x2) % n
    return (x3, (m * (x1 - x3) - y1) % n)


def multiplicar_mod_n(k, P, a, n):
    R, A = INFINITO, P
    while k:
        if k & 1:
            R = suma_mod_n(R, A, a, n)
        A = suma_mod_n(A, A, a, n)
        k >>= 1
    return R


def ecm_una_curva(n, a, P, B):
    """Calcula B! P modulo n; devuelve un factor si aparece."""
    try:
        Q = P
        for j in range(2, B + 1):
            Q = multiplicar_mod_n(j, Q, a, n)
        return None
    except FactorEncontrado as f:
        return f.factor if f.factor != n else None


def ecm(n, B=200, curvas=200, semilla=None):
    """Prueba curvas aleatorias: se elige P = (x0, y0) y a, y luego b
    queda determinado para que P este en la curva."""
    rnd = random.Random(semilla)
    for _ in range(curvas):
        x0, y0, a = rnd.randrange(n), rnd.randrange(n), rnd.randrange(n)
        f = ecm_una_curva(n, a, (x0, y0), B)
        if f:
            return f
    return None


def pruebas():
    n, a, P = 2773, 4, (1, 3)
    P2 = suma_mod_n(P, P, a, n)
    assert P2 == (1771, 705)
    try:
        suma_mod_n(P2, P, a, n)
        raise AssertionError("se esperaba encontrar un factor")
    except FactorEncontrado as f:
        assert f.factor == 59 and n == 59 * 47
    # Ejemplo de Trappe-Washington (apendice): 193279 = 347 * 557, B = 12
    assert ecm_una_curva(193279, -10, (2, 4), 12) in (347, 557)
    # Un numero de 30 cifras = primo de 10 cifras * primo de 20 cifras
    p, q = 1000000007, 10000000000000000051
    f = ecm(p * q, B=300, curvas=300, semilla=1)
    assert f in (p, q)
    print("09_factorizacion_lenstra: todas las pruebas pasaron.")


if __name__ == "__main__":
    n, a, P = 2773, 4, (1, 3)
    print("n = 2773, E: y^2 = x^3 + 4x + 4, P = (1,3)")
    P2 = suma_mod_n(P, P, a, n)
    print("  2P =", P2)
    try:
        suma_mod_n(P2, P, a, n)
    except FactorEncontrado as f:
        print(f"  calculando 3P: mcd(1770, 2773) = {f.factor}  ->  2773 = {f.factor} * {n // f.factor}")
    m = 1000000007 * 10000000000000000051
    print("\nECM sobre", m)
    f = ecm(m, B=300, curvas=300, semilla=1)
    print("  factor encontrado:", f, " cofactor:", m // f)
    print()
    pruebas()
