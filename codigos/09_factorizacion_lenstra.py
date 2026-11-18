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

Prueba con tus propios datos
----------------------------
    python 09_factorizacion_lenstra.py 1000000016000000063
o ejecuta el archivo sin nada mas y escribe tus datos cuando te los pida.
Con  python 09_factorizacion_lenstra.py --ejemplo  se ve el ejemplo de clase y las pruebas.
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


def ejemplo_de_clase():
    """Lo que se muestra sin datos: el ejemplo de clase del informe."""
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


PREGUNTA = ("Escribe un entero n >= 2 para factorizarlo con curvas elipticas\n"
            "(ejemplo: 2773 o 1000000016000000063). Solo Enter = ejemplo de clase: ")
USO = "python 09_factorizacion_lenstra.py 1000000016000000063"


def es_primo_probable(n):
    """Miller-Rabin con bases fijas (exacto para n < 3.3 * 10^24) y algunas al azar."""
    if n < 2:
        return False
    bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41]
    for p in bases:
        if n % p == 0:
            return n == p
    s, r = 0, n - 1
    while r % 2 == 0:
        s, r = s + 1, r // 2
    for a in bases + [random.randrange(2, n - 1) for _ in range(10)]:
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


def ecm_con_detalle(n, max_curvas=600):
    """Prueba curvas al azar con cotas B crecientes; devuelve (factor, curva, B, a, P)."""
    rnd = random.Random()
    B = 100
    for curva in range(1, max_curvas + 1):
        if curva % 40 == 0 and B < 20000:
            B *= 2
        x0, y0, a = rnd.randrange(n), rnd.randrange(n), rnd.randrange(n)
        f = ecm_una_curva(n, a, (x0, y0), B)
        if f:
            return f, curva, B, a, (x0, y0)
    return None


def con_tus_datos(texto):
    import time
    try:
        n = leer_enteros(texto)[0]
    except (ValueError, IndexError):
        print("Escribe un numero entero. Ejemplo:\n  " + USO)
        return
    if n < 2:
        print("n debe ser >= 2. Ejemplo:\n  " + USO)
        return
    print(f"Tu numero: n = {n}  ({len(str(n))} cifras)")
    if len(str(n)) > 40:
        print("Aviso: con mas de 40 cifras esta version didactica puede tardar mucho.")
    inicio = time.perf_counter()
    factores = {}

    resto = n
    pequenos = []
    for p in (2, 3, 5, 7):
        while resto % p == 0:
            factores[p] = factores.get(p, 0) + 1
            pequenos.append(p)
            resto //= p
    print(f"\n1. Division por 2, 3, 5 y 7 (las curvas modulo n necesitan un n sin esos factores): "
          + (f"se encontraron {', '.join(map(str, pequenos))}" if pequenos else "ninguno divide a n"))

    pendientes = [resto] if resto > 1 else []
    paso = 2
    while pendientes:
        m = pendientes.pop()
        if es_primo_probable(m):
            factores[m] = factores.get(m, 0) + 1
            continue
        r = math.isqrt(m)
        if r * r == m:
            print(f"\n{paso}. {m} es un cuadrado perfecto: {m} = {r}^2")
            pendientes += [r, r]
            paso += 1
            continue
        print(f"\n{paso}. ECM sobre {m}:")
        resultado = ecm_con_detalle(m)
        if resultado is None:
            print("   No se encontro factor con estas curvas; se deja como esta.")
            factores[m] = factores.get(m, 0) + 1
            continue
        f, curva, B, a, (x0, y0) = resultado
        b = (y0 * y0 - x0 ** 3 - a * x0) % m
        print(f"   curvas probadas: {curva}.  La que funciono: y^2 = x^3 + a x + b (mod {m}) con")
        print(f"   a = {corto(a)}, b = {corto(b)}, P = ({corto(x0)}, {corto(y0)}), B = {B}.")
        print(f"   Al calcular B! P (multiplicar P por 2, luego por 3, ..., hasta B) aparecio")
        print(f"   un denominador NO invertible modulo {m}:")
        print(f"   su mcd con {m} es {f}  ->  {m} = {f} * {m // f}")
        pendientes += [f, m // f]
        paso += 1

    escrito = " * ".join(f"{p}^{k}" if k > 1 else str(p) for p, k in sorted(factores.items()))
    producto = math.prod(p ** k for p, k in factores.items())
    print(f"\nFactorizacion: {n} = {escrito}")
    print(f"comprobacion (multiplicar): {producto == n};  tiempo: {time.perf_counter() - inicio:.2f} s")
    if len(factores) == 1 and list(factores.values()) == [1]:
        print("(n es primo: no tiene factores propios)")


if __name__ == "__main__":
    datos = leer_datos(PREGUNTA)
    if datos:
        con_tus_datos(datos)
        print("\n(Para ver el ejemplo de clase y las pruebas:  python 09_factorizacion_lenstra.py --ejemplo)")
    else:
        ejemplo_de_clase()
        print()
        pruebas()
