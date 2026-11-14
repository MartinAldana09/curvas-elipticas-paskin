"""
08_firmas_ECDSA.py
==================

Firma digital ECDSA (FIPS 186) y el ataque por nonce repetido.

Firmar (clave privada d, publica Q = dG, G de orden primo n, e = H(m)):
    1. elegir k al azar en [1, n-1]          (k = "nonce", secreto y unico)
    2. R = kG,   r = x(R) mod n              (si r = 0, otro k)
    3. s = k^{-1} (e + d r) mod n            (si s = 0, otro k)
    Firma: (r, s).

Verificar:
    w = s^{-1},  u1 = e w,  u2 = r w  (mod n),   X = u1 G + u2 Q.
    Aceptar si x(X) mod n == r.
Por que funciona:  X = w(e + r d) G = w (k s) G = kG = R.

Nonce repetido: si dos firmas (r, s1), (r, s2) de e1 != e2 usan el mismo k,
    k = (e1 - e2) / (s1 - s2)   y   d = (s1 k - e1) / r    (mod n).

Resultados esperados en la curva de clase (n = 19, d = 7, e = 5, k = 3):
    firma (4, 11);  con e' = 7 la verificacion falla;  firma de e2 = 8: (4, 12);
    el ataque recupera k = 3 y d = 7.

AVISO: educativo. En la practica k se genera de forma determinista (RFC 6979)
o con un generador aleatorio criptografico, y se usan bibliotecas auditadas.
"""
import hashlib
import secrets

INFINITO = None


class Curva:
    def __init__(self, a, b, p, n, G):
        self.a, self.b, self.p, self.n, self.G = a % p, b % p, p, n, G

    def contiene(self, P):
        x, y = P
        return (y * y - x ** 3 - self.a * x - self.b) % self.p == 0

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

    def multiplicar(self, k, P):
        R, A = INFINITO, P
        while k:
            if k & 1:
                R = self.suma(R, A)
            A = self.suma(A, A)
            k >>= 1
        return R


CURVA_CLASE = Curva(24, 13, 29, n=19, G=(1, 3))

P256 = Curva(
    a=-3,
    b=0x5AC635D8AA3A93E7B3EBBD55769886BC651D06B0CC53B0F63BCE3C3E27D2604B,
    p=2**256 - 2**224 + 2**192 + 2**96 - 1,
    n=0xFFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC632551,
    G=(0x6B17D1F2E12C4247F8BCE6E563A440F277037D812DEB33A0F4A13945D898C296,
       0x4FE342E2FE1A7F9B8EE7EB4A7C0F9E162BCE33576B315ECECBB6406837BF51F5),
)

# secp256k1 (la curva de Bitcoin, SEC 2): y^2 = x^3 + 7
SECP256K1 = Curva(
    a=0, b=7,
    p=2**256 - 2**32 - 977,
    n=0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141,
    G=(0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798,
       0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8),
)


def hash_a_entero(mensaje, n):
    """e = SHA-256(mensaje) reducido a la longitud en bits de n."""
    h = int.from_bytes(hashlib.sha256(mensaje.encode("utf-8")).digest(), "big")
    exceso = 256 - n.bit_length()
    return h >> exceso if exceso > 0 else h


def firmar(E, d, e, k=None):
    while True:
        kk = k if k is not None else secrets.randbelow(E.n - 1) + 1
        R = E.multiplicar(kk, E.G)
        r = R[0] % E.n
        s = pow(kk, -1, E.n) * (e + d * r) % E.n
        if r and s:
            return r, s
        if k is not None:
            raise ValueError("ese k produce r = 0 o s = 0")


def verificar(E, Q, e, firma):
    r, s = firma
    if not (1 <= r < E.n and 1 <= s < E.n):
        return False
    w = pow(s, -1, E.n)
    u1, u2 = e * w % E.n, r * w % E.n
    X = E.suma(E.multiplicar(u1, E.G), E.multiplicar(u2, Q))
    return X is not INFINITO and X[0] % E.n == r


def ataque_nonce_repetido(n, e1, firma1, e2, firma2):
    (r1, s1), (r2, s2) = firma1, firma2
    assert r1 == r2, "las firmas no comparten nonce"
    k = (e1 - e2) * pow(s1 - s2, -1, n) % n
    d = (s1 * k - e1) * pow(r1, -1, n) % n
    return k, d


def pruebas():
    E = CURVA_CLASE
    d = 7
    Q = E.multiplicar(d, E.G)
    assert Q == (15, 6)
    f1 = firmar(E, d, 5, k=3)
    assert f1 == (4, 11) and verificar(E, Q, 5, f1)
    assert not verificar(E, Q, 7, f1)                # mensaje alterado
    f2 = firmar(E, d, 8, k=3)
    assert f2 == (4, 12)
    assert ataque_nonce_repetido(E.n, 5, f1, 8, f2) == (3, 7)
    # Parametros estandar correctos
    for C in (P256, SECP256K1):
        assert C.contiene(C.G) and C.multiplicar(C.n, C.G) is INFINITO
    # Firma real con P-256 y SHA-256
    d = secrets.randbelow(P256.n - 1) + 1
    Q = P256.multiplicar(d, P256.G)
    e = hash_a_entero("Pasquin Matematico", P256.n)
    f = firmar(P256, d, e)
    assert verificar(P256, Q, e, f)
    assert not verificar(P256, Q, hash_a_entero("Pasquin matematico", P256.n), f)
    # Nonce repetido a tamano real: se recupera d de 256 bits
    k = secrets.randbelow(P256.n - 1) + 1
    e1, e2 = hash_a_entero("mensaje 1", P256.n), hash_a_entero("mensaje 2", P256.n)
    k_rec, d_rec = ataque_nonce_repetido(P256.n, e1, firmar(P256, d, e1, k), e2, firmar(P256, d, e2, k))
    assert (k_rec, d_rec) == (k, d)
    print("08_firmas_ECDSA: todas las pruebas pasaron.")


if __name__ == "__main__":
    E = CURVA_CLASE
    d = 7
    Q = E.multiplicar(d, E.G)
    print("Curva de clase: n = 19, d = 7, Q =", Q)
    f1 = firmar(E, d, 5, k=3)
    print("  firma de e = 5 con k = 3:", f1, " valida:", verificar(E, Q, 5, f1))
    print("  la misma firma para e' = 7:", verificar(E, Q, 7, f1))
    f2 = firmar(E, d, 8, k=3)
    print("  firma de e = 8 con el MISMO k:", f2)
    print("  ataque -> (k, d) =", ataque_nonce_repetido(E.n, 5, f1, 8, f2))

    print("\nP-256 con SHA-256:")
    d = secrets.randbelow(P256.n - 1) + 1
    Q = P256.multiplicar(d, P256.G)
    e = hash_a_entero("Hola", P256.n)
    f = firmar(P256, d, e)
    print("  firma valida:", verificar(P256, Q, e, f))
    k = secrets.randbelow(P256.n - 1) + 1
    e1, e2 = hash_a_entero("transaccion A", P256.n), hash_a_entero("transaccion B", P256.n)
    _, d_rec = ataque_nonce_repetido(P256.n, e1, firmar(P256, d, e1, k), e2, firmar(P256, d, e2, k))
    print("  nonce repetido -> clave privada recuperada:", d_rec == d)
    print()
    pruebas()
