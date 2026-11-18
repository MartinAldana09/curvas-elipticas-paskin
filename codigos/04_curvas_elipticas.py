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

Prueba con tus propios datos
----------------------------
    python 04_curvas_elipticas.py 2 3 97      (a = 2, b = 3, p = 97; puedes agregar un punto: 2 3 97 3 6)
o ejecuta el archivo sin nada mas y escribe tus datos cuando te los pida.
Con  python 04_curvas_elipticas.py --ejemplo  se ve el ejemplo de clase y las pruebas.
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


def ejemplo_de_clase():
    """Lo que se muestra sin datos: el ejemplo de clase del informe."""
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


PREGUNTA = ("Escribe a, b y p (p primo mayor que 3) para la curva y^2 = x^3 + ax + b,\n"
            "y si quieres un punto x y (ejemplo: 2 3 97  o  24 13 29 1 3).\n"
            "Solo Enter = ejemplo de clase: ")
USO = "python 04_curvas_elipticas.py 2 3 97      (a = 2, b = 3, p = 97; puedes agregar un punto: 2 3 97 3 6)"


def es_primo_pequeno(n):
    return n > 1 and all(n % d for d in range(2, math.isqrt(n) + 1))


def explicar_suma(E, P, Q):
    """Imprime la cuenta de P + Q con las formulas de la secante o la tangente."""
    p = E.p
    if P is INFINITO or Q is INFINITO or (P[0] == Q[0] and (P[1] + Q[1]) % p == 0):
        print(f"   {formato(P)} + {formato(Q)} = O")
        return
    (x1, y1), (x2, y2) = P, Q
    if P == Q:
        num, den = (3 * x1 * x1 + E.a) % p, (2 * y1) % p
        print(f"   Tangente: m = (3x1^2 + a)/(2y1) = (3*{x1}^2 + {E.a})/(2*{y1}) = {num}/{den} (mod {p})")
    else:
        num, den = (y2 - y1) % p, (x2 - x1) % p
        print(f"   Secante: m = (y2 - y1)/(x2 - x1) = ({y2} - {y1})/({x2} - {x1}) = {num}/{den} (mod {p})")
    inv = pow(den, -1, p)
    m = num * inv % p
    print(f"   dividir es multiplicar por el inverso: {den}^(-1) = {inv} (mod {p}),  m = {num}*{inv} mod {p} = {m}")
    x3 = (m * m - x1 - x2) % p
    y3 = (m * (x1 - x3) - y1) % p
    print(f"   x3 = m^2 - x1 - x2 = {m}^2 - {x1} - {x2} = {x3} (mod {p})")
    print(f"   y3 = m(x1 - x3) - y1 = {m}*({x1} - {x3}) - {y1} = {y3} (mod {p})")
    print(f"   resultado: {formato(P)} + {formato(Q)} = ({x3},{y3})")


def con_tus_datos(texto):
    try:
        numeros = leer_enteros(texto)
    except ValueError:
        print("Escribe solo numeros enteros. Ejemplo:\n  " + USO)
        return
    if len(numeros) < 3:
        print("Hacen falta tres numeros: a, b y p. Ejemplo:\n  " + USO)
        return
    a, b, p = numeros[:3]
    if p <= 3 or not es_primo_pequeno(p):
        print(f"p = {p} no sirve: debe ser un primo mayor que 3 (por ejemplo 29, 97, 101, 1009).")
        return
    if p > 200000:
        print("p es demasiado grande para contar los puntos uno por uno (usa p < 200000).")
        return
    print(f"Tu curva: y^2 = x^3 + {a % p}x + {b % p} sobre F_{p}")
    D = (4 * a ** 3 + 27 * b ** 2) % p
    print(f"\n1. Discriminante: 4a^3 + 27b^2 mod {p} = {D}")
    if D == 0:
        print("   Es 0: la curva es SINGULAR (tiene un pico o un nudo) y sus puntos no forman")
        print("   un grupo como el de las curvas elipticas. Cambia a o b.")
        return
    print("   Es distinto de 0: la curva es no singular (es una curva eliptica).")
    E = CurvaEliptica(a, b, p)
    puntos = E.puntos()
    N = len(puntos) + 1
    print(f"\n2. Puntos: se prueba cada x de 0 a {p - 1} y se busca y con y^2 = x^3 + ax + b (mod {p}).")
    if len(puntos) <= 120:
        print("   " + " ".join(formato(P) for P in puntos))
    else:
        print("   primeros: " + " ".join(formato(P) for P in puntos[:20]) + " ...")
    print(f"   #E(F_{p}) = {len(puntos)} puntos (x, y) + el punto en el infinito O = {N}")
    cota = 2 * math.sqrt(p)
    print(f"\n3. Teorema de Hasse: |#E - (p + 1)| = |{N} - {p + 1}| = {abs(N - p - 1)}"
          f" <= 2*sqrt({p}) = {cota:.2f}   {'(se cumple)' if abs(N - p - 1) <= cota else '(?)'}")
    if es_primo_pequeno(N):
        print(f"   #E = {N} es primo: todo punto distinto de O genera la curva entera.")
    else:
        divisores = [k for k in range(2, N) if N % k == 0]
        print(f"   #E = {N} no es primo; por Lagrange, el orden de cada punto divide a {N}"
              f" (divisores: 1, {', '.join(map(str, divisores))}, {N}).")

    G = None
    if len(numeros) >= 5:
        G = (numeros[3] % p, numeros[4] % p)
        if not E.contiene(G):
            print(f"\n   Ojo: ({numeros[3]},{numeros[4]}) no esta en la curva:"
                  f" {G[1]}^2 mod {p} = {G[1] ** 2 % p}, pero x^3 + ax + b = {(G[0] ** 3 + E.a * G[0] + E.b) % p}.")
            G = None
    if G is None:
        G = next((P for P in puntos if P[1] != 0), None)
        if G is None:
            print("\nLa curva no tiene puntos con y distinto de 0; no hay tabla de multiplos.")
            return
    n = E.orden(G)
    print(f"\n4. Punto G = {formato(G)}: su orden es {n} (el menor n con nG = O).")
    print("   Tabla de multiplos (cada fila es la anterior + G):")
    tabla = E.tabla_de_multiplos(G)
    filas = tabla if n <= 40 else tabla[:12] + [None] + tabla[-3:]
    for fila in filas:
        if fila is None:
            print("    ...")
        else:
            k, Q = fila
            print(f"   {k:4d}G = {formato(Q)}")

    print("\n5. Como se calcula 2G (recta tangente):")
    G2 = E.suma(G, G)
    explicar_suma(E, G, G)
    if G2 is not INFINITO and n > 3:
        print("   Y 3G = 2G + G (recta secante):")
        explicar_suma(E, G2, G)


if __name__ == "__main__":
    datos = leer_datos(PREGUNTA)
    if datos:
        con_tus_datos(datos)
        print("\n(Para ver el ejemplo de clase y las pruebas:  python 04_curvas_elipticas.py --ejemplo)")
    else:
        ejemplo_de_clase()
        print()
        pruebas()
