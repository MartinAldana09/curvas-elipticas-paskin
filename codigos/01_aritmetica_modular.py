"""
01_aritmetica_modular.py
========================

Herramientas de aritmetica modular que se usan en RSA y en las curvas elipticas.

Contenido
---------
1. Algoritmo de Euclides (maximo comun divisor).
2. Algoritmo de Euclides extendido (identidad de Bezout).
3. Inverso modular.
4. Exponenciacion modular rapida (cuadrados sucesivos).
5. Teorema Chino del Resto (TCR).
6. Funcion phi de Euler.
7. Fracciones continuas y convergentes (se usan en el ataque de Wiener).

Fundamento matematico
---------------------
* a = b (mod n)  significa que n divide a (a - b).
* Bezout: existen enteros x, y con a*x + b*y = mcd(a, b).
* a tiene inverso modulo n  <=>  mcd(a, n) = 1. El inverso sale de Bezout.
* Euler: si mcd(a, n) = 1 entonces a^phi(n) = 1 (mod n).
* TCR: si los modulos son coprimos dos a dos, el sistema x = a_i (mod m_i)
  tiene solucion unica modulo el producto de los m_i.

Como ejecutar
-------------
    python 01_aritmetica_modular.py

AVISO: codigo con fines educativos. Para proteger informacion real se deben
usar bibliotecas criptograficas reconocidas y parametros estandarizados.

Prueba con tus propios datos
----------------------------
    python 01_aritmetica_modular.py 240 46 13      (a = 240, n = 46, e = 13; e es opcional)
o ejecuta el archivo sin nada mas y escribe tus datos cuando te los pida.
Con  python 01_aritmetica_modular.py --ejemplo  se ve el ejemplo de clase y las pruebas.
"""


# ---------------------------------------------------------------------------
# 1. Euclides
# ---------------------------------------------------------------------------
def mcd(a, b, mostrar=False):
    """Maximo comun divisor por divisiones sucesivas.

    Ejemplo: mcd(240, 46) = 2.
    """
    a, b = abs(a), abs(b)
    while b != 0:
        q, r = divmod(a, b)
        if mostrar:
            print(f"  {a} = {q}*{b} + {r}")
        a, b = b, r
    return a


# ---------------------------------------------------------------------------
# 2. Euclides extendido
# ---------------------------------------------------------------------------
def euclides_extendido(a, b):
    """Devuelve (g, x, y) con a*x + b*y = g = mcd(a, b).

    Se mantienen dos combinaciones lineales de a y b que se van
    actualizando igual que los restos del algoritmo de Euclides.
    Ejemplo: euclides_extendido(240, 46) = (2, -9, 47).
    """
    r0, r1 = a, b
    x0, x1 = 1, 0
    y0, y1 = 0, 1
    while r1 != 0:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        x0, x1 = x1, x0 - q * x1
        y0, y1 = y1, y0 - q * y1
    return r0, x0, y0


# ---------------------------------------------------------------------------
# 3. Inverso modular
# ---------------------------------------------------------------------------
def inverso_modular(a, n):
    """Inverso de a modulo n (existe si y solo si mcd(a, n) = 1).

    Ejemplo: inverso_modular(3, 40) = 27, porque 3*27 = 81 = 2*40 + 1.
    """
    g, x, _ = euclides_extendido(a % n, n)
    if g != 1:
        raise ValueError(f"{a} no tiene inverso modulo {n} (mcd = {g})")
    return x % n


# ---------------------------------------------------------------------------
# 4. Exponenciacion modular rapida
# ---------------------------------------------------------------------------
def potencia_modular(base, exponente, modulo, mostrar=False):
    """Calcula base^exponente mod modulo con cuadrados sucesivos.

    Se recorre el exponente en binario: en cada paso se eleva al cuadrado
    y, si el bit vale 1, se multiplica por la base. Hace O(log e)
    multiplicaciones en lugar de e - 1.
    Ejemplo: potencia_modular(3, 13, 100) = 23.
    """
    resultado = 1
    cuadrado = base % modulo
    e = exponente
    paso = 0
    while e > 0:
        if e & 1:
            resultado = (resultado * cuadrado) % modulo
        if mostrar:
            print(f"  bit {paso}: {e & 1}  base^(2^{paso}) = {cuadrado}  acumulado = {resultado}")
        cuadrado = (cuadrado * cuadrado) % modulo
        e >>= 1
        paso += 1
    return resultado


# ---------------------------------------------------------------------------
# 5. Teorema Chino del Resto
# ---------------------------------------------------------------------------
def tcr(residuos, modulos):
    """Resuelve x = residuos[i] (mod modulos[i]) con modulos coprimos dos a dos.

    Formula constructiva: con M = prod(m_i) y M_i = M / m_i,
        x = sum( a_i * M_i * (M_i^{-1} mod m_i) )  (mod M).
    Ejemplo: tcr([2, 3], [3, 5]) = 8.
    """
    M = 1
    for m in modulos:
        M *= m
    x = 0
    for a_i, m_i in zip(residuos, modulos):
        M_i = M // m_i
        x += a_i * M_i * inverso_modular(M_i, m_i)
    return x % M


# ---------------------------------------------------------------------------
# 6. Funcion phi de Euler
# ---------------------------------------------------------------------------
def factorizar_por_tanteo(n):
    """Factorizacion por division de prueba (solo para n pequenos)."""
    factores = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            factores[d] = factores.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        factores[n] = factores.get(n, 0) + 1
    return factores


def phi_euler(n):
    """phi(n) = n * prod(1 - 1/p) sobre los primos p que dividen a n.

    Ejemplos: phi(12) = 4, phi(p) = p - 1, phi(p*q) = (p-1)(q-1).
    """
    resultado = n
    for p in factorizar_por_tanteo(n):
        resultado = resultado // p * (p - 1)
    return resultado


# ---------------------------------------------------------------------------
# 7. Fracciones continuas
# ---------------------------------------------------------------------------
def fraccion_continua(num, den):
    """Coeficientes [a0; a1, a2, ...] de la fraccion continua de num/den."""
    coef = []
    while den != 0:
        q = num // den
        coef.append(q)
        num, den = den, num - q * den
    return coef


def convergentes(coef):
    """Convergentes h_i/k_i a partir de los coeficientes.

    Recurrencias: h_i = a_i h_{i-1} + h_{i-2},  k_i = a_i k_{i-1} + k_{i-2},
    con h_{-1} = 1, h_{-2} = 0, k_{-1} = 0, k_{-2} = 1.
    """
    h_ant, h = 0, 1
    k_ant, k = 1, 0
    lista = []
    for a in coef:
        h_ant, h = h, a * h + h_ant
        k_ant, k = k, a * k + k_ant
        lista.append((h, k))
    return lista


# ---------------------------------------------------------------------------
# Pruebas (resultados esperados)
# ---------------------------------------------------------------------------
def pruebas():
    # Congruencias
    assert (38 - 3) % 7 == 0                       # 38 = 3 (mod 7)
    # Euclides y Bezout
    assert mcd(240, 46) == 2
    g, x, y = euclides_extendido(240, 46)
    assert (g, x, y) == (2, -9, 47) and 240 * x + 46 * y == 2
    # Inverso modular
    assert inverso_modular(3, 40) == 27
    try:
        inverso_modular(4, 40)
        raise AssertionError("4 no deberia tener inverso modulo 40")
    except ValueError:
        pass
    # Exponenciacion rapida
    assert potencia_modular(3, 13, 100) == 23 == pow(3, 13, 100)
    for b, e, m in [(7, 3, 55), (13, 27, 55), (2, 10**6, 10**9 + 7)]:
        assert potencia_modular(b, e, m) == pow(b, e, m)
    # TCR
    assert tcr([2, 3], [3, 5]) == 8
    assert tcr([1, 2, 3], [5, 7, 9]) % 5 == 1
    # Euler
    assert phi_euler(12) == 4 and phi_euler(55) == 40 and phi_euler(97) == 96
    # Teorema de Euler: a^phi(n) = 1 (mod n) si mcd(a, n) = 1
    n = 55
    for a in range(1, n):
        if mcd(a, n) == 1:
            assert pow(a, phi_euler(n), n) == 1
    # Fracciones continuas: 415/93 = [4; 2, 6, 7]
    assert fraccion_continua(415, 93) == [4, 2, 6, 7]
    assert convergentes([4, 2, 6, 7])[-1] == (415, 93)
    print("01_aritmetica_modular: todas las pruebas pasaron.")


def ejemplo_de_clase():
    """Lo que se muestra sin datos: el ejemplo de clase del informe."""
    print("Euclides para mcd(240, 46):")
    print("  mcd =", mcd(240, 46, mostrar=True))
    print("Bezout:", euclides_extendido(240, 46), "  (2 = 240*(-9) + 46*47)")
    print("Inverso de 3 modulo 40:", inverso_modular(3, 40))
    print("3^13 mod 100 paso a paso:")
    print("  resultado =", potencia_modular(3, 13, 100, mostrar=True))
    print("TCR: x = 2 (mod 3), x = 3 (mod 5)  ->  x =", tcr([2, 3], [3, 5]))
    print("phi(12) =", phi_euler(12))


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


PREGUNTA = ("Escribe dos enteros positivos a y n, y si quieres un exponente e\n"
            "(ejemplo: 240 46 13). Solo Enter = ejemplo de clase: ")
USO = 'python 01_aritmetica_modular.py 240 46 13      (a = 240, n = 46, e = 13; e es opcional)'


def con_tus_datos(texto):
    try:
        numeros = leer_enteros(texto)
    except ValueError:
        print("Escribe solo numeros enteros. Ejemplo:\n  " + USO)
        return
    if len(numeros) < 2 or numeros[0] < 0 or numeros[1] < 2:
        print("Hacen falta dos enteros a >= 0 y n >= 2. Ejemplo:\n  " + USO)
        return
    a, n = numeros[0], numeros[1]
    e = numeros[2] if len(numeros) > 2 else None
    print(f"Tus datos: a = {a}, n = {n}" + (f", e = {e}" if e is not None else ""))

    print(f"\n1. Congruencia: {a} = {a // n}*{n} + {a % n},  asi que  {a} = {a % n} (mod {n}).")

    print(f"\n2. Algoritmo de Euclides para mcd({a}, {n}) (cada fila: dividendo = cociente*divisor + resto):")
    g = mcd(a, n, mostrar=True)
    print(f"   mcd({a}, {n}) = {g}   (el ultimo resto distinto de 0)")

    g, x, y = euclides_extendido(a, n)
    print(f"\n3. Identidad de Bezout (Euclides extendido):  {a}*({x}) + {n}*({y}) = {a * x + n * y}")

    print(f"\n4. Inverso de {a} modulo {n}:")
    if g == 1:
        inv = inverso_modular(a, n)
        print(f"   {a}^(-1) = {inv} (mod {n})  (sale del coeficiente {x} de Bezout: {x} mod {n} = {inv})")
        print(f"   comprobacion: {a}*{inv} = {a * inv} = {a * inv // n}*{n} + {a * inv % n}")
    else:
        print(f"   NO existe, porque mcd({a}, {n}) = {g} es distinto de 1.")

    print(f"\n5. Funcion phi de Euler:")
    if n < 10 ** 12:
        factores = factorizar_por_tanteo(n)
        escrito = " * ".join(f"{p}^{k}" if k > 1 else str(p) for p, k in factores.items())
        ph = phi_euler(n)
        print(f"   {n} = {escrito},  phi({n}) = {ph}"
              "   (cuantos numeros entre 1 y n no comparten factores con n)")
        if g == 1:
            print(f"   Teorema de Euler: {a}^phi({n}) = {a}^{ph} mod {n} = {pow(a, ph, n)}  (siempre da 1)")
    else:
        print("   n es demasiado grande para factorizarlo por tanteo; se omite.")

    if e is not None:
        if e < 0:
            print("\n6. El exponente debe ser >= 0; se omite la potencia.")
        else:
            print(f"\n6. {a}^{e} mod {n} con cuadrados sucesivos; {e} en binario = {bin(e)[2:]}.")
            print("   En cada paso se eleva al cuadrado; si el bit es 1, se multiplica al acumulado:")
            r = potencia_modular(a, e, n, mostrar=e.bit_length() <= 64)
            print(f"   resultado: {a}^{e} mod {n} = {r}   (comprobacion con pow de Python: {pow(a, e, n)})")

    coef = fraccion_continua(a, n)
    print(f"\n7. Fraccion continua de {a}/{n} = [{coef[0]}; {', '.join(map(str, coef[1:]))}]")
    conv = convergentes(coef)
    print("   convergentes (fracciones que se acercan a a/n):",
          ", ".join(f"{h}/{k}" for h, k in conv[:10]) + (" ..." if len(conv) > 10 else ""))


if __name__ == "__main__":
    datos = leer_datos(PREGUNTA)
    if datos:
        con_tus_datos(datos)
        print("\n(Para ver el ejemplo de clase y las pruebas:  python 01_aritmetica_modular.py --ejemplo)")
    else:
        ejemplo_de_clase()
        print()
        pruebas()
