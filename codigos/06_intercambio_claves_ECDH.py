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

Prueba con tus propios datos
----------------------------
    python 06_intercambio_claves_ECDH.py 3 5      (a = 3 para Alice, b = 5 para Bob)
o ejecuta el archivo sin nada mas y escribe tus datos cuando te los pida.
Con  python 06_intercambio_claves_ECDH.py --ejemplo  se ve el ejemplo de clase y las pruebas.
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


def ejemplo_de_clase():
    """Lo que se muestra sin datos: el ejemplo de clase del informe."""
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


PREGUNTA = ("Escribe el numero secreto de Alice y el de Bob (dos enteros >= 1, ejemplo: 3 5).\n"
            "Solo Enter = ejemplo de clase: ")
USO = "python 06_intercambio_claves_ECDH.py 3 5      (a = 3 para Alice, b = 5 para Bob)"


def con_tus_datos(texto):
    try:
        a, b = leer_enteros(texto)[:2]
    except ValueError:
        print("Escribe dos numeros enteros. Ejemplo:\n  " + USO)
        return
    if a < 1 or b < 1:
        print("Los secretos deben ser >= 1. Ejemplo:\n  " + USO)
        return
    E = CURVA_CLASE
    fmt = lambda P: "O" if P is INFINITO else f"({P[0]},{P[1]})"
    print(f"Secretos: Alice a = {a}, Bob b = {b}.  Curva de clase: G = (1,3), orden n = 19.")
    A = E.multiplicar(a, E.G)
    B = E.multiplicar(b, E.G)
    print("\n1. Cada uno calcula su clave publica y la envia por un canal publico:")
    print(f"   Alice: A = aG = {a}G = {fmt(A)}")
    print(f"   Bob:   B = bG = {b}G = {fmt(B)}")
    if A is INFINITO or B is INFINITO:
        print("   Ojo: un secreto multiplo de 19 da el punto O, que no esconde nada. Usa otro numero.")
        return
    S_alice = E.multiplicar(a, B)
    S_bob = E.multiplicar(b, A)
    print("\n2. Cada uno combina SU secreto con la clave publica del otro:")
    print(f"   Alice calcula aB = {a}*{fmt(B)} = {fmt(S_alice)}")
    print(f"   Bob   calcula bA = {b}*{fmt(A)} = {fmt(S_bob)}")
    print(f"   Iguales: {S_alice == S_bob}, porque aB = a(bG) = (ab)G = b(aG) = bA.")
    if a * b < 19:
        print(f"   Secreto compartido: (ab)G = {a * b}G = {fmt(S_alice)}")
    else:
        print(f"   Secreto compartido: (ab)G = {a * b}G = {(a * b) % 19}G = {fmt(S_alice)}"
              f"   (porque 19G = O, se reduce {a * b} mod 19)")
    if S_alice is not INFINITO:
        print(f"   Clave de 128 bits derivada (SHA-256 de la x): {derivar_clave(S_alice).hex()}")

    print("\n3. Eve ve G, A y B. En esta curva tan pequena puede probar todos los multiplos de G:")
    R, k = INFINITO, 0
    while R != A:
        R, k = E.suma(R, E.G), k + 1
    print(f"   encuentra {k}G = A, asi que conoce a mod 19 = {k} y calcula el secreto: {fmt(E.multiplicar(k, B))}")
    print("   Por eso en la practica el orden de G tiene unos 256 bits.")

    print("\n4. Lo mismo con P-256 (claves aleatorias de 256 bits, como en la practica):")
    a2, A2 = generar_par(P256)
    b2, B2 = generar_par(P256)
    print(f"   a = {corto(a2)}   A = (x = {corto(A2[0])}, y = ...)")
    print(f"   b = {corto(b2)}   B = (x = {corto(B2[0])}, y = ...)")
    k_alice = derivar_clave(secreto_compartido(P256, a2, B2))
    k_bob = derivar_clave(secreto_compartido(P256, b2, A2))
    print(f"   clave de Alice: {k_alice.hex()}")
    print(f"   clave de Bob:   {k_bob.hex()}   (iguales: {k_alice == k_bob})")
    print("   Encontrar a desde A exigiria del orden de 2^128 operaciones.")


if __name__ == "__main__":
    datos = leer_datos(PREGUNTA)
    if datos:
        con_tus_datos(datos)
        print("\n(Para ver el ejemplo de clase y las pruebas:  python 06_intercambio_claves_ECDH.py --ejemplo)")
    else:
        ejemplo_de_clase()
        print()
        pruebas()
