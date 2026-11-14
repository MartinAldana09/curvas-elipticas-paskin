"""
04_curvas_elipticas.py
======================

Curvas elipticas sobre un cuerpo finito F_p (p primo, p > 3) en forma corta
de Weierstrass:

        E :  y^2 = x^3 + a x + b,      con  4 a^3 + 27 b^2 != 0 (mod p).

Contenido
---------
1. Discriminante y no singularidad.
2. Generacion y validacion de puntos (lista completa de E(F_p)).
3. Teorema de Hasse:  | #E(F_p) - (p + 1) | <= 2 sqrt(p).
4. Elemento neutro (punto en el infinito), inverso y suma de puntos.
5. Orden de un punto y tabla de multiplos kG.
6. Verificacion EXHAUSTIVA de los axiomas de grupo en una curva pequena.

Formulas de la suma (P = (x1, y1), Q = (x2, y2), P != -Q)
--------------------------------------------------------
    m = (y2 - y1) / (x2 - x1)          si P != Q   (recta secante)
    m = (3 x1^2 + a) / (2 y1)          si P == Q   (recta tangente)
    x3 = m^2 - x1 - x2,     y3 = m (x1 - x3) - y1.
Dividir modulo p = multiplicar por el inverso modular.

Curva de trabajo del proyecto
-----------------------------
    y^2 = x^3 + 24 x + 13  sobre F_29,  G = (1, 3),  #E = 38,  orden(G) = 19.

AVISO: parametros pequenos, solo para aprender. No son seguros.
"""
import math

INFINITO = None        # representamos el punto en el infinito O con None


class CurvaEliptica:
    def __init__(self, a, b, p):
        self.a, self.b, self.p = a % p, b % p, p
        if self.discriminante() == 0:
            raise ValueError("la curva es singular: 4a^3 + 27b^2 = 0 (mod p)")

    def __repr__(self):
        return f"y^2 = x^3 + {self.a}x + {self.b} sobre F_{self.p}"

    # 1. Discriminante -----------------------------------------------------
    def discriminante(self):
        """4a^3 + 27b^2 mod p. Es 0 exactamente cuando x^3+ax+b tiene raiz repetida."""
        return (4 * self.a ** 3 + 27 * self.b ** 2) % self.p

    # 2. Puntos ------------------------------------------------------------
    def contiene(self, P):
        if P is INFINITO:
            return True
        x, y = P
        return (y * y - (x ** 3 + self.a * x + self.b)) % self.p == 0

    def puntos(self):
        """Lista de todos los puntos afines de E(F_p) (busqueda directa, p pequeno)."""
        p = self.p
        raices = {}
        for y in range(p):
            raices.setdefault(y * y % p, []).append(y)
        lista = []
        for x in range(p):
            for y in raices.get((x ** 3 + self.a * x + self.b) % p, []):
                lista.append((x, y))
        return lista

    def numero_de_puntos(self):
        """#E(F_p), contando el punto en el infinito."""
        return len(self.puntos()) + 1

    def cumple_hasse(self):
        N = self.numero_de_puntos()
        return (N - (self.p + 1)) ** 2 <= 4 * self.p   # |N - p - 1| <= 2 sqrt(p)

    # 4. Operaciones de grupo ---------------------------------------------
    def negativo(self, P):
        if P is INFINITO:
            return INFINITO
        x, y = P
        return (x, (-y) % self.p)

    def suma(self, P, Q):
        if P is INFINITO:
            return Q
        if Q is INFINITO:
            return P
        x1, y1 = P
        x2, y2 = Q
        p = self.p
        if x1 == x2 and (y1 + y2) % p == 0:
            return INFINITO                               # P + (-P) = O
        if P == Q:
            m = (3 * x1 * x1 + self.a) * pow(2 * y1, -1, p) % p
        else:
            m = (y2 - y1) * pow(x2 - x1, -1, p) % p
        x3 = (m * m - x1 - x2) % p
        y3 = (m * (x1 - x3) - y1) % p
        return (x3, y3)

    def resta(self, P, Q):
        return self.suma(P, self.negativo(Q))

    def multiplicar(self, k, P):
        """kP por duplicar y sumar (ver 05_multiplicacion_escalar.py)."""
        if k < 0:
            return self.multiplicar(-k, self.negativo(P))
        R, A = INFINITO, P
        while k:
            if k & 1:
                R = self.suma(R, A)
            A = self.suma(A, A)
            k >>= 1
        return R

    # 5. Orden ------------------------------------------------------------
    def orden(self, P):
        """Menor n >= 1 con nP = O."""
        n, Q = 1, P
        while Q is not INFINITO:
            Q = self.suma(Q, P)
            n += 1
        return n

    def tabla_de_multiplos(self, G):
        n = self.orden(G)
        return [(k, self.multiplicar(k, G)) for k in range(1, n + 1)]


def formato(P):
    return "O" if P is INFINITO else f"({P[0]},{P[1]})"


# ---------------------------------------------------------------------------
# Pruebas
# ---------------------------------------------------------------------------
def pruebas():
    # Ejemplo de Trappe-Washington: y^2 = x^3 + 2x - 1 mod 5 tiene 7 puntos
    E5 = CurvaEliptica(2, -1, 5)
    assert sorted(E5.puntos()) == [(0, 2), (0, 3), (2, 1), (2, 4), (4, 1), (4, 4)]
    # Curva singular: y^2 = x^3 (a = b = 0)
    try:
        CurvaEliptica(0, 0, 29)
        raise AssertionError("deberia ser singular")
    except ValueError:
        pass
    # Curva de trabajo
    E = CurvaEliptica(24, 13, 29)
    G = (1, 3)
    assert E.contiene(G) and not E.contiene((1, 4))
    assert E.discriminante() == 3
    assert E.numero_de_puntos() == 38 and E.cumple_hasse()
    assert E.orden(G) == 19
    esperado = [(1, 3), (11, 10), (23, 28), (0, 10), (19, 7), (18, 19), (15, 6),
                (20, 24), (4, 12), (4, 17), (20, 5), (15, 23), (18, 10), (19, 22),
                (0, 19), (23, 1), (11, 19), (1, 26), INFINITO]
    assert [Q for _, Q in E.tabla_de_multiplos(G)] == esperado
    # Axiomas de grupo verificados con TODOS los puntos (38 puntos -> 54 872 ternas)
    todos = E.puntos() + [INFINITO]
    for P in todos:
        assert E.suma(P, INFINITO) == P                         # neutro
        assert E.suma(P, E.negativo(P)) is INFINITO             # inverso
        assert E.orden(P) and 38 % E.orden(P) == 0              # Lagrange
        for Q in todos:
            R = E.suma(P, Q)
            assert E.contiene(R)                                # clausura
            assert R == E.suma(Q, P)                            # conmutatividad
    for P in todos:
        for Q in todos:
            PQ = E.suma(P, Q)
            for R in todos:
                assert E.suma(PQ, R) == E.suma(P, E.suma(Q, R))  # asociatividad
    # Hasse en muchas curvas
    for p in [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61]:
        for a in range(p):
            for b in range(p):
                if (4 * a ** 3 + 27 * b ** 2) % p:
                    assert CurvaEliptica(a, b, p).cumple_hasse()
    print("04_curvas_elipticas: todas las pruebas pasaron.")


if __name__ == "__main__":
    E = CurvaEliptica(24, 13, 29)
    G = (1, 3)
    print("Curva:", E)
    print("Discriminante 4a^3 + 27b^2 mod 29 =", E.discriminante(), "(distinto de 0: no singular)")
    N = E.numero_de_puntos()
    print(f"#E(F_29) = {N};  Hasse: |{N} - 30| = {abs(N - 30)} <= 2*sqrt(29) = {2 * math.sqrt(29):.2f}")
    print("Puntos afines:", " ".join(formato(P) for P in E.puntos()))
    print("Multiplos de G = (1,3):")
    for k, Q in E.tabla_de_multiplos(G):
        print(f"  {k:2d}G = {formato(Q)}")
    print("Ejemplo: (1,3) + (11,10) =", formato(E.suma((1, 3), (11, 10))))
    print()
    pruebas()
