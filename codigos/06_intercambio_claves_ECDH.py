"""
06_intercambio_claves_ECDH.py
=============================

Intercambio de claves de Diffie-Hellman sobre curvas elipticas (ECDH).

Idea
----
Datos publicos: la curva E, un punto base G de orden primo n.
    Alice elige a secreto y publica A = aG.
    Bob   elige b secreto y publica B = bG.
    Alice calcula aB = a(bG);  Bob calcula bA = b(aG).
Ambos obtienen el mismo punto S = (ab)G porque la multiplicacion escalar
cumple a(bG) = (ab)G = b(aG). Eve ve A y B, pero para obtener S tendria que
resolver (en general) el ECDLP. De la coordenada x de S se deriva la clave
simetrica (por ejemplo, para AES) con una funcion de derivacion (KDF).

Contenido
---------
1. Ejemplo de clase en y^2 = x^3 + 24x + 13 sobre F_29:  a = 3, b = 5 -> S = 15G = (0, 19).
2. ECDH con la curva estandar P-256 (NIST) y comprobacion de sus parametros.
3. Hombre en el medio: por que las claves publicas deben ir firmadas.

AVISO: implementacion educativa (sin tiempo constante ni validaciones completas).
Para usos reales: bibliotecas como 'cryptography' (X25519 / ECDH con P-256).
"""
import hashlib
import secrets

INFINITO = None


class Curva:
    def __init__(self, a, b, p, n=None, G=None):
        self.a, self.b, self.p, self.n, self.G = a % p, b % p, p, n, G

    def contiene(self, P):
        if P is INFINITO:
            return True
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


# Curva de trabajo (didactica)
CURVA_CLASE = Curva(24, 13, 29, n=19, G=(1, 3))

# P-256 = secp256r1 (NIST FIPS 186-4 / SP 800-186)
P256 = Curva(
    a=-3,
    b=0x5AC635D8AA3A93E7B3EBBD55769886BC651D06B0CC53B0F63BCE3C3E27D2604B,
    p=2**256 - 2**224 + 2**192 + 2**96 - 1,
    n=0xFFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC632551,
    G=(0x6B17D1F2E12C4247F8BCE6E563A440F277037D812DEB33A0F4A13945D898C296,
       0x4FE342E2FE1A7F9B8EE7EB4A7C0F9E162BCE33576B315ECECBB6406837BF51F5),
)


def generar_par(E):
    d = secrets.randbelow(E.n - 1) + 1
    return d, E.multiplicar(d, E.G)


def secreto_compartido(E, mi_privada, su_publica):
    if not E.contiene(su_publica):                 # validacion de la clave recibida
        raise ValueError("el punto recibido no esta en la curva")
    return E.multiplicar(mi_privada, su_publica)


def derivar_clave(S, longitud=16):
    """KDF muy simple: SHA-256 de la coordenada x (solo ilustrativo)."""
    x = S[0]
    return hashlib.sha256(x.to_bytes((x.bit_length() + 7) // 8 or 1, "big")).digest()[:longitud]


def pruebas():
    E = CURVA_CLASE
    A = E.multiplicar(3, E.G)
    B = E.multiplicar(5, E.G)
    assert A == (23, 28) and B == (19, 7)
    assert secreto_compartido(E, 3, B) == secreto_compartido(E, 5, A) == (0, 19)
    assert E.multiplicar(15, E.G) == (0, 19)
    # Parametros de P-256: G esta en la curva y nG = O
    assert P256.contiene(P256.G)
    assert P256.multiplicar(P256.n, P256.G) is INFINITO
    # ECDH real con P-256
    a, A = generar_par(P256)
    b, B = generar_par(P256)
    assert derivar_clave(secreto_compartido(P256, a, B)) == derivar_clave(secreto_compartido(P256, b, A))
    # Validacion de puntos
    try:
        secreto_compartido(P256, a, (1, 1))
        raise AssertionError("debio rechazar un punto fuera de la curva")
    except ValueError:
        pass
    print("06_intercambio_claves_ECDH: todas las pruebas pasaron.")


if __name__ == "__main__":
    E = CURVA_CLASE
    print("Ejemplo de clase:", "G =", E.G, " n =", E.n)
    A, B = E.multiplicar(3, E.G), E.multiplicar(5, E.G)
    print("  Alice: a = 3, envia A = 3G =", A)
    print("  Bob:   b = 5, envia B = 5G =", B)
    print("  Alice calcula 3B =", secreto_compartido(E, 3, B), "; Bob calcula 5A =", secreto_compartido(E, 5, A))
    print("  secreto: 15G =", E.multiplicar(15, E.G))

    print("\nECDH con P-256:")
    a, A = generar_par(P256)
    b, B = generar_par(P256)
    k_alice = derivar_clave(secreto_compartido(P256, a, B))
    k_bob = derivar_clave(secreto_compartido(P256, b, A))
    print("  clave AES-128 de Alice:", k_alice.hex())
    print("  clave AES-128 de Bob:  ", k_bob.hex())

    print("\nHombre en el medio (Eve intercepta y reemplaza A y B):")
    e, Eve_pub = generar_par(P256)
    k_alice_eve = derivar_clave(secreto_compartido(P256, a, Eve_pub))
    k_eve_alice = derivar_clave(secreto_compartido(P256, e, A))
    print("  Alice cree hablar con Bob, pero comparte clave con Eve:", k_alice_eve == k_eve_alice)
    print("  Defensa: firmar las claves publicas (ECDSA, ver 08_firmas_ECDSA.py).")
    print()
    pruebas()
