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

Prueba con tus propios datos
----------------------------
    python 07_elgamal_eliptico.py "Martin Aldana"
o ejecuta el archivo sin nada mas y escribe tus datos cuando te los pida.
Con  python 07_elgamal_eliptico.py --ejemplo  se ve el ejemplo de clase y las pruebas.
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


def ejemplo_de_clase():
    """Lo que se muestra sin datos: el ejemplo de clase del informe."""
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


PREGUNTA = "Escribe tu mensaje (por ejemplo, tu nombre). Solo Enter = ejemplo de clase: "
USO = 'python 07_elgamal_eliptico.py "Martin Aldana"'

# P-256 (NIST): y^2 = x^3 - 3x + b sobre F_p, con p = 3 (mod 4)
P256_P = 2**256 - 2**224 + 2**192 + 2**96 - 1
P256_B = 0x5AC635D8AA3A93E7B3EBBD55769886BC651D06B0CC53B0F63BCE3C3E27D2604B
P256_N = 0xFFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC632551
P256_G = (0x6B17D1F2E12C4247F8BCE6E563A440F277037D812DEB33A0F4A13945D898C296,
          0x4FE342E2FE1A7F9B8EE7EB4A7C0F9E162BCE33576B315ECECBB6406837BF51F5)


def punto_corto(P):
    return "O" if P is INFINITO else f"(x = {corto(P[0], 8)}, y = {corto(P[1], 8)})"


def con_tus_datos(mensaje):
    azar = random.SystemRandom()
    E = Curva(-3, P256_B, P256_P)
    G, n, K = P256_G, P256_N, 100
    print(f'Tu mensaje: "{mensaje}"')
    print("Curva P-256 (la de NIST), Koblitz con K = 100.")
    d = azar.randrange(1, n)
    Q = E.multiplicar(d, G)
    print("\n1. Claves de Bob:")
    print(f"   privada d = {corto(d)}")
    print(f"   publica Q = dG = {punto_corto(Q)}")

    datos = mensaje.encode("utf-8")
    bloques = [datos[i:i + 30] for i in range(0, len(datos), 30)]
    print(f"\n2. El texto se parte en {len(bloques)} bloque(s) de hasta 30 bytes; cada bloque es un numero m.")
    recuperado = b""
    for i, bloque in enumerate(bloques, 1):
        m = int.from_bytes(bloque, "big")
        print(f"\n   Bloque {i}: {bloque.decode('utf-8', errors='replace')!r}  ->  m = {corto(m)}")
        M = koblitz_codificar(E, m, K)
        j = M[0] - m * K
        print(f"   Koblitz: se prueba x = 100m + j con j = 0, 1, ...; sirvio j = {j}"
              f" ({j + 1} intento(s)): M = {punto_corto(M)}")
        k = azar.randrange(1, n)
        C1, C2 = elgamal_cifrar(E, G, Q, M, k)
        print(f"   Alice elige k al azar y envia C1 = kG = {punto_corto(C1)}")
        print(f"                              C2 = M + kQ = {punto_corto(C2)}")
        M2 = elgamal_descifrar(E, d, C1, C2)
        m2 = koblitz_decodificar(M2, K)
        print(f"   Bob calcula C2 - d*C1 = M (igual: {M2 == M})  ->  m = floor(x/100) = {corto(m2)}")
        recuperado += m2.to_bytes(len(bloque), "big")
    print(f'\n3. Bob junta los bloques: "{recuperado.decode("utf-8", errors="replace")}"')

    M = koblitz_codificar(E, int.from_bytes(bloques[0], "big"), K)
    C1a, _ = elgamal_cifrar(E, G, Q, M, azar.randrange(1, n))
    C1b, _ = elgamal_cifrar(E, G, Q, M, azar.randrange(1, n))
    print("\n4. Si se cifra dos veces el mismo bloque, el cifrado cambia (porque k cambia):")
    print(f"   C1 la primera vez: {punto_corto(C1a)}")
    print(f"   C1 la segunda vez: {punto_corto(C1b)}")
    print("   Asi un espia no puede saber si dos mensajes cifrados son iguales.")


if __name__ == "__main__":
    datos = leer_datos(PREGUNTA)
    if datos:
        con_tus_datos(datos)
        print("\n(Para ver el ejemplo de clase y las pruebas:  python 07_elgamal_eliptico.py --ejemplo)")
    else:
        ejemplo_de_clase()
        print()
        pruebas()
