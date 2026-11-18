"""
02_rsa.py
=========

Implementacion didactica del criptosistema RSA (Rivest, Shamir, Adleman, 1978).

Fundamento matematico
---------------------
1. Se eligen dos primos p, q y se calcula N = p*q.
2. phi(N) = (p - 1)(q - 1).
3. Se elige e con mcd(e, phi(N)) = 1 (es comun e = 65537).
4. d = e^{-1} mod phi(N)  (Euclides extendido).
5. Clave publica (N, e); clave privada (p, q, d).
6. Cifrar:    c = m^e mod N.     Descifrar: m = c^d mod N.

Correctitud: ed = 1 + k*phi(N), y por el teorema de Euler
    c^d = m^{ed} = m * (m^{phi(N)})^k = m (mod N)   si mcd(m, N) = 1;
el caso mcd(m, N) != 1 se resuelve trabajando modulo p y modulo q (TCR).

Datos de entrada / resultados esperados
---------------------------------------
Ejemplo de clase: p = 5, q = 11, e = 3  ->  N = 55, phi = 40, d = 27,
m = 7  ->  c = 13  ->  descifrado 7.

AVISO: este RSA "de libro" (sin relleno aleatorio OAEP) es INSEGURO.
Solo sirve para entender las matematicas.

Prueba con tus propios datos
----------------------------
    python 02_rsa.py "Martin Aldana"
o ejecuta el archivo sin nada mas y escribe tus datos cuando te los pida.
Con  python 02_rsa.py --ejemplo  se ve el ejemplo de clase y las pruebas.
"""
import random

# ---------------------------------------------------------------------------
# Herramientas (copiadas de 01 para que el archivo funcione por si solo)
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
        raise ValueError(f"{a} no es invertible modulo {n}")
    return x % n


# ---------------------------------------------------------------------------
# Primalidad: prueba de Miller-Rabin
# ---------------------------------------------------------------------------
def es_primo_probable(n, rondas=40):
    """Prueba probabilistica de Miller-Rabin.

    Escribe n - 1 = 2^s * r con r impar. Si n es primo, para todo a
    la sucesion a^r, a^{2r}, ..., a^{2^s r} empieza en 1 o pasa por -1.
    Un a que no cumpla esto es un "testigo" de que n es compuesto.
    """
    if n < 2:
        return False
    pequenos = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
    for p in pequenos:
        if n % p == 0:
            return n == p
    s, r = 0, n - 1
    while r % 2 == 0:
        s += 1
        r //= 2
    for _ in range(rondas):
        a = random.randrange(2, n - 1)
        x = pow(a, r, n)
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
    """Primo aleatorio de 'bits' bits con los dos bits mas altos en 1
    (asi el producto de dos de ellos tiene exactamente 2*bits bits)."""
    while True:
        candidato = random.getrandbits(bits) | (3 << (bits - 2)) | 1
        if es_primo_probable(candidato):
            return candidato


# ---------------------------------------------------------------------------
# RSA
# ---------------------------------------------------------------------------
def generar_claves(bits=1024, e=65537):
    """Genera (clave_publica, clave_privada) con N de 'bits' bits."""
    while True:
        p = primo_aleatorio(bits // 2)
        q = primo_aleatorio(bits // 2)
        if p == q:
            continue
        phi = (p - 1) * (q - 1)
        if euclides_extendido(e, phi)[0] == 1:
            break
    N = p * q
    d = inverso_modular(e, phi)
    return (N, e), (p, q, d)


def claves_desde_primos(p, q, e):
    """Construye las claves a partir de primos dados (ejemplos de clase)."""
    N, phi = p * q, (p - 1) * (q - 1)
    d = inverso_modular(e, phi)
    return (N, e), (p, q, d)


def cifrar(m, clave_publica):
    N, e = clave_publica
    if not 0 <= m < N:
        raise ValueError("el mensaje debe ser un entero entre 0 y N - 1")
    return pow(m, e, N)


def descifrar(c, clave_privada):
    p, q, d = clave_privada
    return pow(c, d, p * q)


def descifrar_tcr(c, clave_privada):
    """Descifrado acelerado con el TCR: dos potencias pequenas en vez de una grande.

    m_p = c^{d mod (p-1)} mod p,  m_q = c^{d mod (q-1)} mod q,
    y luego se combinan con el TCR (formula de Garner).
    """
    p, q, d = clave_privada
    m_p = pow(c, d % (p - 1), p)
    m_q = pow(c, d % (q - 1), q)
    h = (inverso_modular(q, p) * (m_p - m_q)) % p
    return m_q + h * q


def texto_a_entero(texto):
    return int.from_bytes(texto.encode("utf-8"), "big")


def entero_a_texto(m):
    return m.to_bytes((m.bit_length() + 7) // 8, "big").decode("utf-8")


# ---------------------------------------------------------------------------
# Pruebas
# ---------------------------------------------------------------------------
def pruebas():
    # Ejemplo de clase
    pub, priv = claves_desde_primos(5, 11, 3)
    assert pub == (55, 3) and priv[2] == 27
    c = cifrar(7, pub)
    assert c == 13 and descifrar(c, priv) == 7
    # Todos los mensajes posibles modulo 55 (incluye mcd(m, N) != 1)
    for m in range(55):
        assert descifrar(cifrar(m, pub), priv) == m
    # Claves de tamano realista
    random.seed(2026)
    pub, priv = generar_claves(1024)
    N, e = pub
    assert N.bit_length() == 1024
    mensaje = "Hola, Pasquin Matematico"
    m = texto_a_entero(mensaje)
    c = cifrar(m, pub)
    assert descifrar(c, priv) == m == descifrar_tcr(c, priv)
    assert entero_a_texto(descifrar(c, priv)) == mensaje
    # Primalidad
    assert es_primo_probable(2**127 - 1) and not es_primo_probable(2**128 + 1)
    assert not es_primo_probable(561)              # numero de Carmichael
    print("02_rsa: todas las pruebas pasaron.")


def ejemplo_de_clase():
    """Lo que se muestra sin datos: el ejemplo de clase del informe."""
    pub, priv = claves_desde_primos(5, 11, 3)
    print("Ejemplo de clase: clave publica", pub, " d =", priv[2])
    c = cifrar(7, pub)
    print("  m = 7  ->  c =", c, " ->  descifrado:", descifrar(c, priv))
    random.seed(1)
    pub, priv = generar_claves(512)
    m = texto_a_entero("ECC")
    c = cifrar(m, pub)
    print("RSA de 512 bits: N tiene", pub[0].bit_length(), "bits")
    print("  'ECC' -> c =", hex(c)[:34], "...  -> descifrado:", entero_a_texto(descifrar_tcr(c, priv)))


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
USO = 'python 02_rsa.py "Martin Aldana"'


def con_tus_datos(mensaje):
    print(f'Tu mensaje: "{mensaje}"')

    # ---------------- Parte A: version de juguete -------------------------
    p, q, e = 61, 53, 17
    (N, _), (_, _, d) = claves_desde_primos(p, q, e)
    print("\n=== Parte A. RSA de juguete, letra por letra (p = 61, q = 53) ===")
    print(f"N = 61*53 = {N},  phi = 60*52 = {(p - 1) * (q - 1)},  e = {e},  d = e^(-1) mod phi = {d}")
    print("Cada caracter se cambia por su codigo (byte UTF-8, un numero de 0 a 255) y se cifra:")
    print("   caracter | codigo m | cifrado c = m^17 mod 3233 | descifrado c^2753 mod 3233")
    vistos = {}
    for byte in mensaje.encode("utf-8")[:24]:
        c = pow(byte, e, N)
        letra = repr(chr(byte)) if byte < 128 else "(byte)"
        print(f"   {letra:>8} | {byte:8d} | {c:25d} | {pow(c, d, N):10d}")
        vistos.setdefault(byte, []).append(c)
    repetidos = [chr(b) for b, lista in vistos.items() if len(lista) > 1 and b < 128]
    if len(mensaje.encode("utf-8")) > 24:
        print("   (solo se muestran los primeros 24 caracteres)")
    if len(repetidos) == 1:
        print(f"Fijate: {repetidos[0]!r} se repite y SIEMPRE da el mismo cifrado.")
    elif repetidos:
        print(f"Fijate: {', '.join(repr(r) for r in repetidos)} se repiten y cada una SIEMPRE da el mismo cifrado.")
    print("Cifrar letra por letra es inseguro: letras iguales dan cifrados iguales")
    print("(analisis de frecuencias, como en el cifrado de Cesar) y N = 3233 se factoriza al instante.")

    # ---------------- Parte B: tamano real ---------------------------------
    m = texto_a_entero(mensaje)
    bits = max(512, (m.bit_length() // 64 + 2) * 64)      # N debe ser mayor que m
    print(f"\n=== Parte B. RSA de {bits} bits con el mensaje completo ===")
    print("(los primos se eligen al azar: si ejecutas otra vez, todos los numeros cambian,")
    print(" pero el descifrado siempre devuelve tu mensaje)")
    (N, e), (p, q, d) = generar_claves(bits)
    phi = (p - 1) * (q - 1)
    print(f"\nPaso 1. Bob genera dos primos de {bits // 2} bits (probados con Miller-Rabin):")
    print(f"   p = {p}")
    print(f"   q = {q}")
    print(f"Paso 2. N = p*q ({N.bit_length()} bits, {len(str(N))} cifras):")
    print(f"   N = {N}")
    print(f"Paso 3. phi(N) = (p - 1)(q - 1) = {phi}")
    print(f"Paso 4. Exponente publico e = {e}  (mcd(e, phi) = 1)")
    print(f"Paso 5. Exponente privado d = e^(-1) mod phi(N) = {d}")
    print(f"   comprobacion: e*d mod phi(N) = {e * d % phi}")
    print("   Clave publica (la conoce todo el mundo): (N, e).  Clave privada (solo Bob): d, p, q.")
    datos = mensaje.encode("utf-8")
    print("Paso 6. Alice convierte el texto en un numero (bytes UTF-8 en base 256):")
    print(f"   bytes: {' '.join(f'{b:02X}' for b in datos)}   (hexadecimal)")
    print(f"   m = {m}   (m < N: si)")
    c = cifrar(m, (N, e))
    print("Paso 7. Alice cifra con la clave publica: c = m^e mod N =")
    print(f"   c = {c}")
    m2 = descifrar(c, (p, q, d))
    print("Paso 8. Bob descifra con su clave privada: c^d mod N =")
    print(f"   {m2}   -> igual a m: {m2 == m}")
    print(f"   (con el TCR, que es mas rapido, da lo mismo: {descifrar_tcr(c, (p, q, d)) == m})")
    print(f'Paso 9. Bob convierte el numero en texto: "{entero_a_texto(m2)}"')
    print("\nAVISO: este RSA 'de libro' no usa relleno aleatorio (OAEP); es solo para aprender.")


if __name__ == "__main__":
    datos = leer_datos(PREGUNTA)
    if datos:
        con_tus_datos(datos)
        print("\n(Para ver el ejemplo de clase y las pruebas:  python 02_rsa.py --ejemplo)")
    else:
        ejemplo_de_clase()
        print()
        pruebas()
