"""
07_elgamal_eliptico.py
======================

Cifrado de ElGamal sobre curvas elipticas y representacion de mensajes como
puntos (metodo de Koblitz). Basado en Trappe-Washington, secciones 16.2.3 y 16.5.

1. Mensajes como puntos (Koblitz)
---------------------------------
Para guardar un entero m en la coordenada x se prueban x = mK + j,
j = 0, 1, ..., K - 1, hasta que x^3 + ax + b sea un cuadrado modulo p.
Cada intento falla con probabilidad cercana a 1/2, asi que todo el proceso
falla con probabilidad cercana a 1/2^K. Se recupera m = floor(x / K).

2. ElGamal eliptico
-------------------
Bob tiene clave privada d y publica Q = dG. Para cifrar el punto M,
Alice elige k al azar y envia (C1, C2) = (kG, M + kQ).
Bob descifra:  C2 - d C1 = M + k(dG) - d(kG) = M.

Resultados esperados (curva de clase, d = 7, Q = (15, 6), M = (4, 12), k = 3):
    C1 = (23, 28),  kQ = (11, 10),  C2 = (20, 5),  descifrado (4, 12).

AVISO: esquema didactico; no usar para proteger informacion real.
"""
import random

INFINITO = None


class Curva:
    def __init__(self, a, b, p):
        self.a, self.b, self.p = a % p, b % p, p

    def suma(self, P, Q):
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

    def multiplicar(self, k, P):
        R, A = INFINITO, P
        while k:
            if k & 1:
                R = self.suma(R, A)
            A = self.suma(A, A)
            k >>= 1
        return R

    def raiz_cuadrada(self, v):
        """Raiz cuadrada modulo p por busqueda (p pequeno) o formula si p = 3 mod 4."""
        v %= self.p
        if self.p % 4 == 3:
            y = pow(v, (self.p + 1) // 4, self.p)
            return y if y * y % self.p == v else None
        for y in range(self.p):
            if y * y % self.p == v:
                return y
        return None


def koblitz_codificar(E, m, K):
    for j in range(K):
        x = m * K + j
        y = E.raiz_cuadrada(x ** 3 + E.a * x + E.b)
        if y is not None:
            return (x, y)
    raise ValueError("no se encontro punto; aumente K")


def koblitz_decodificar(P, K):
    return P[0] // K


def elgamal_cifrar(E, G, Q, M, k=None):
    if k is None:
        k = random.randrange(1, E.p)
    return E.multiplicar(k, G), E.suma(M, E.multiplicar(k, Q))


def elgamal_descifrar(E, d, C1, C2):
    return E.suma(C2, E.negativo(E.multiplicar(d, C1)))


def pruebas():
    E = Curva(24, 13, 29)
    G, d = (1, 3), 7
    Q = E.multiplicar(d, G)
    assert Q == (15, 6)
    # Koblitz con K = 4 en la curva de clase
    assert koblitz_codificar(E, 1, 4) == (4, 12)        # x = 4 funciona al primer intento
    assert koblitz_codificar(E, 2, 4)[0] == 9        # x = 8 falla, x = 9 funciona
    # Ejemplo de Trappe-Washington: p = 179, y^2 = x^3 + 2x + 7, K = 10, m = 5 -> (51, +-11)
    E179 = Curva(2, 7, 179)
    P5 = koblitz_codificar(E179, 5, 10)
    assert P5[0] == 51 and P5[1] in (11, 168) and koblitz_decodificar(P5, 10) == 5
    # ElGamal de clase
    M = (4, 12)
    C1, C2 = elgamal_cifrar(E, G, Q, M, k=3)
    assert C1 == (23, 28) and E.multiplicar(3, Q) == (11, 10) and C2 == (20, 5)
    assert elgamal_descifrar(E, d, C1, C2) == M
    # Ejemplo de Trappe-Washington: y^2 = x^3 + 3x + 45 mod 8831, G = (4, 11), 3G = (413, 1808)
    E8831 = Curva(3, 45, 8831)
    assert E8831.multiplicar(3, (4, 11)) == (413, 1808)
    # Ida y vuelta aleatoria
    random.seed(3)
    for _ in range(50):
        k = random.randrange(1, 19)
        assert elgamal_descifrar(E, d, *elgamal_cifrar(E, G, Q, M, k)) == M
    print("07_elgamal_eliptico: todas las pruebas pasaron.")


if __name__ == "__main__":
    E = Curva(24, 13, 29)
    G, d = (1, 3), 7
    Q = E.multiplicar(d, G)
    print("Clave de Bob: d = 7, Q = 7G =", Q)
    print("Koblitz con K = 4:")
    for m in (1, 2):
        P = koblitz_codificar(E, m, 4)
        print(f"  m = {m}  ->  punto {P}  ->  m = floor({P[0]}/4) = {koblitz_decodificar(P, 4)}")
    M = (4, 12)
    C1, C2 = elgamal_cifrar(E, G, Q, M, k=3)
    print("Cifrado de M = (4,12) con k = 3:  C1 =", C1, " C2 =", C2)
    print("Descifrado: C2 - d*C1 =", elgamal_descifrar(E, d, C1, C2))
    print()
    pruebas()
