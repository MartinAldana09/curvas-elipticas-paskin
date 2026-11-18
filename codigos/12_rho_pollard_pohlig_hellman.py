"""
12_rho_pollard_pohlig_hellman.py
================================

Dos ataques al problema del logaritmo discreto eliptico (ECDLP):
dado Q = kG, encontrar k.

1. Rho de Pollard (1978)
------------------------
Se hace una "caminata pseudoaleatoria" R_{i+1} = R_i + M_j, donde el salto
M_j = c_j G + d_j Q depende de la coordenada x de R_i. Cada punto se lleva
escrito como R = a G + b Q. Como el grupo es finito, la caminata termina
en un ciclo (forma de letra rho). Con el metodo de la tortuga y la liebre
(Floyd) se encuentra una coincidencia
        a G + b Q = a' G + b' Q   =>   k = (a - a') / (b' - b)  (mod n).
Por la paradoja del cumpleanos se esperan unos sqrt(pi n / 2) pasos y casi
nada de memoria.

2. Pohlig-Hellman
-----------------
Si el orden n de G se factoriza como prod q^e con primos q pequenos, se
calcula k modulo cada q^e trabajando en subgrupos pequenos y se junta todo
con el Teorema Chino del Resto. El costo depende del MAYOR primo de n, no
de n. Por eso en criptografia se exige que G tenga orden primo grande.

Resultados esperados
--------------------
* Curva de trabajo y^2 = x^3 + 24x + 13 sobre F_29 (#E = 38 = 2 * 19):
  ambos metodos recuperan k.
* Curvas de ~32 bits: rho necesita del orden de 10^5 pasos si el orden es
  primo, mientras que Pohlig-Hellman resuelve al instante si el orden es suave.

AVISO: educativo. Con n de 256 bits, rho necesitaria unas 2^128 operaciones.

Prueba con tus propios datos
----------------------------
    python 12_rho_pollard_pohlig_hellman.py 123456789
o ejecuta el archivo sin nada mas y escribe tus datos cuando te los pida.
Con  python 12_rho_pollard_pohlig_hellman.py --ejemplo  se ve el ejemplo de clase y las pruebas.
"""
import math
import random
import time

INFINITO = None


class Curva:
    def __init__(self, a, b, p):
        self.a, self.b, self.p = a % p, b % p, p
        self.sumas = 0

    def contiene(self, P):
        if P is INFINITO:
            return True
        x, y = P
        return (y * y - x ** 3 - self.a * x - self.b) % self.p == 0

    def suma(self, P, Q):
        self.sumas += 1
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
        if k < 0:
            return self.multiplicar(-k, self.negativo(P))
        R, A = INFINITO, P
        while k:
            if k & 1:
                R = self.suma(R, A)
            A = self.suma(A, A)
            k >>= 1
        return R


# ---------------------------------------------------------------------------
# Herramientas de teoria de numeros
# ---------------------------------------------------------------------------
def factorizar(n):
    f, d = {}, 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def es_primo(n):
    return n > 1 and factorizar(n) == {n: 1}


def tcr(residuos, modulos):
    M = math.prod(modulos)
    return sum(a * (M // m) * pow(M // m, -1, m) for a, m in zip(residuos, modulos)) % M


# ---------------------------------------------------------------------------
# 1. Rho de Pollard
# ---------------------------------------------------------------------------
def rho_pollard(E, G, Q, n, particiones=16, semilla=0):
    """Resuelve Q = kG con G de orden primo n."""
    rnd = random.Random(semilla)
    while True:
        saltos = []
        for _ in range(particiones):
            c, d = rnd.randrange(n), rnd.randrange(n)
            saltos.append((E.suma(E.multiplicar(c, G), E.multiplicar(d, Q)), c, d))

        def paso(R, a, b):
            j = 0 if R is INFINITO else R[0] % particiones
            M, c, d = saltos[j]
            return E.suma(R, M), (a + c) % n, (b + d) % n

        a0, b0 = rnd.randrange(n), rnd.randrange(n)
        R0 = E.suma(E.multiplicar(a0, G), E.multiplicar(b0, Q))
        tortuga = liebre = (R0, a0, b0)
        pasos = 0
        while True:
            tortuga = paso(*tortuga)
            liebre = paso(*paso(*liebre))
            pasos += 1
            if tortuga[0] == liebre[0]:
                break
        (_, a1, b1), (_, a2, b2) = tortuga, liebre
        if (b2 - b1) % n:
            return (a1 - a2) * pow(b2 - b1, -1, n) % n, pasos
        # coincidencia inutil (b1 = b2): se reinicia con otros saltos


# ---------------------------------------------------------------------------
# 2. Pohlig-Hellman
# ---------------------------------------------------------------------------
def log_en_subgrupo_pequeno(E, G, Q, q):
    """k con Q = kG, G de orden primo q pequeno (busqueda directa)."""
    R = INFINITO
    for k in range(q):
        if R == Q:
            return k
        R = E.suma(R, G)
    raise ValueError("Q no esta en el subgrupo generado por G")


def pohlig_hellman(E, G, Q, n):
    """Resuelve Q = kG con ord(G) = n, usando la factorizacion de n."""
    residuos, modulos = [], []
    for q, e in factorizar(n).items():
        G0 = E.multiplicar(n // q, G)                 # orden q
        k_qe = 0
        for i in range(e):                            # digitos de k en base q
            resto = E.suma(Q, E.negativo(E.multiplicar(k_qe, G)))
            Qi = E.multiplicar(n // q ** (i + 1), resto)
            k_qe += log_en_subgrupo_pequeno(E, G0, Qi, q) * q ** i
        residuos.append(k_qe)
        modulos.append(q ** e)
    return tcr(residuos, modulos)


# ---------------------------------------------------------------------------
# Busqueda de curvas para los experimentos
# ---------------------------------------------------------------------------
def orden_de_la_curva(E, P):
    """#E(F_p) por paso de bebe - paso de gigante en el intervalo de Hasse.
    Se buscan los m en [p+1-2sqrt(p), p+1+2sqrt(p)] con mP = O; si el orden de P
    es mayor que 4 sqrt(p), hay uno solo y es #E(F_p)."""
    p = E.p
    lo = p + 1 - 2 * math.isqrt(p) - 2
    ancho = 4 * math.isqrt(p) + 5
    s = math.isqrt(ancho) + 1
    bebes = {}
    R = INFINITO
    for j in range(s):
        bebes.setdefault(R, []).append(j)
        R = E.suma(R, P)
    paso = E.negativo(E.multiplicar(s, P))
    gamma = E.negativo(E.multiplicar(lo, P))         # -(lo) P
    candidatos = []
    for i in range(s + 1):
        if gamma in bebes:                           # (lo + i s + j) P = O
            candidatos += [lo + i * s + j for j in bebes[gamma]]
        gamma = E.suma(gamma, paso)
    candidatos = [m for m in candidatos if abs(m - (p + 1)) <= 2 * math.isqrt(p) + 1]
    return candidatos[0] if len(set(candidatos)) == 1 else None


def buscar_curva(bits, condicion, semilla):
    """Curva y^2 = x^3 + ax + b sobre F_p (p primo de 'bits' bits, p = 3 mod 4)
    y un punto G cuyo orden N = #E cumple condicion(N)."""
    rnd = random.Random(semilla)
    while True:
        p = rnd.randrange(2 ** (bits - 1), 2 ** bits) | 3
        if not es_primo(p):
            continue
        a, b = rnd.randrange(p), rnd.randrange(p)
        if (4 * a ** 3 + 27 * b * b) % p == 0:
            continue
        E = Curva(a, b, p)
        x = rnd.randrange(p)
        v = (x ** 3 + a * x + b) % p
        y = pow(v, (p + 1) // 4, p)
        if y * y % p != v or y == 0:
            continue
        G = (x, y)
        N = orden_de_la_curva(E, G)
        if N is None or not condicion(N):
            continue
        if all(E.multiplicar(N // q, G) is not INFINITO for q in factorizar(N)):
            return E, G, N                           # G genera todo E(F_p)


# ---------------------------------------------------------------------------
# Pruebas
# ---------------------------------------------------------------------------
def pruebas():
    E = Curva(24, 13, 29)
    # Curva de trabajo: G = (1,3) tiene orden 19; P = G + (24,0) tiene orden 38
    G = (1, 3)
    for k in range(19):
        assert rho_pollard(E, G, E.multiplicar(k, G), 19, particiones=4, semilla=k)[0] == k
    P = E.suma(G, (24, 0))
    assert E.multiplicar(38, P) is INFINITO and E.multiplicar(19, P) is not INFINITO
    for k in range(38):
        assert pohlig_hellman(E, P, E.multiplicar(k, P), 38) == k
    # Curvas de ~24 bits
    E1, G1, n1 = buscar_curva(24, es_primo, semilla=1)
    k = 1234567 % n1
    assert rho_pollard(E1, G1, E1.multiplicar(k, G1), n1)[0] == k
    suave = lambda N: max(factorizar(N)) < 200
    E2, G2, n2 = buscar_curva(24, suave, semilla=2)
    k = 7654321 % n2
    assert pohlig_hellman(E2, G2, E2.multiplicar(k, G2), n2) == k
    print("12_rho_pollard_pohlig_hellman: todas las pruebas pasaron.")


def ejemplo_de_clase():
    """Lo que se muestra sin datos: el ejemplo de clase del informe."""
    E = Curva(24, 13, 29)
    G = (1, 3)
    k, pasos = rho_pollard(E, G, (15, 6), 19, particiones=4)
    print(f"Curva de trabajo: Q = (15,6) = kG  ->  rho de Pollard da k = {k} ({pasos} iteraciones)")
    P = E.suma(G, (24, 0))
    print(f"P = G + (24,0) = {P} tiene orden 38 = 2*19; Pohlig-Hellman con Q = 25P:",
          pohlig_hellman(E, P, E.multiplicar(25, P), 38))

    print("\nCurvas de 32 bits:")
    E1, G1, n1 = buscar_curva(32, es_primo, semilla=7)
    print(f"  (a) orden PRIMO   n = {n1} = {factorizar(n1)}")
    k = random.randrange(1, n1)
    Q = E1.multiplicar(k, G1)
    E1.sumas = 0
    t = time.perf_counter()
    k1, pasos = rho_pollard(E1, G1, Q, n1)
    print(f"      rho de Pollard: k correcto = {k1 == k}, {pasos} iteraciones"
          f" (sqrt(pi n/2) ~ {int(math.sqrt(math.pi * n1 / 2))}), {E1.sumas} sumas, "
          f"{time.perf_counter() - t:.2f} s")
    suave = lambda N: max(factorizar(N)) < 1000
    E2, G2, n2 = buscar_curva(32, suave, semilla=8)
    print(f"  (b) orden SUAVE   n = {n2} = {factorizar(n2)}")
    k = random.randrange(1, n2)
    Q = E2.multiplicar(k, G2)
    E2.sumas = 0
    t = time.perf_counter()
    k2 = pohlig_hellman(E2, G2, Q, n2)
    print(f"      Pohlig-Hellman: k correcto = {k2 == k}, {E2.sumas} sumas, "
          f"{time.perf_counter() - t:.3f} s")


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


PREGUNTA = ("Escribe tu numero secreto k (un entero >= 1). Eve solo vera el punto Q = kG\n"
            "e intentara recuperar k. Solo Enter = ejemplo de clase: ")
USO = "python 12_rho_pollard_pohlig_hellman.py 123456789"


def pohlig_hellman_con_detalle(E, G, Q, n):
    """Igual que pohlig_hellman, pero imprime k modulo cada potencia de primo."""
    residuos, modulos = [], []
    for q, e in factorizar(n).items():
        G0 = E.multiplicar(n // q, G)
        k_qe = 0
        for i in range(e):
            resto = E.suma(Q, E.negativo(E.multiplicar(k_qe, G)))
            Qi = E.multiplicar(n // q ** (i + 1), resto)
            k_qe += log_en_subgrupo_pequeno(E, G0, Qi, q) * q ** i
        print(f"      subgrupo de orden {q ** e:>4}:  k = {k_qe} (mod {q ** e})")
        residuos.append(k_qe)
        modulos.append(q ** e)
    return tcr(residuos, modulos)


def con_tus_datos(texto):
    try:
        k = leer_enteros(texto)[0]
    except (ValueError, IndexError):
        print("Escribe un numero entero. Ejemplo:\n  " + USO)
        return
    if k < 1:
        print("k debe ser >= 1. Ejemplo:\n  " + USO)
        return
    print(f"Tu secreto: k = {k}")

    E = Curva(24, 13, 29)
    G = (1, 3)
    print("\n1. Curva de clase (G = (1,3) de orden 19):")
    if k % 19 == 0:
        print("   k es multiplo de 19, asi que Q = kG = O: no esconde nada. Prueba otro k.")
    else:
        Q = E.multiplicar(k, G)
        k1, pasos = rho_pollard(E, G, Q, 19, particiones=4)
        print(f"   Q = kG = {Q}.  Rho de Pollard encuentra k = {k1} en {pasos} pasos.")
        print(f"   Como 19G = O, Eve solo puede (y solo necesita) saber k mod 19 = {k % 19}: {k1 == k % 19}")

    print("\n2. Curva de 32 bits con orden PRIMO (asi deben ser las curvas en la practica):")
    E1, G1, n1 = buscar_curva(32, es_primo, semilla=7)
    k_1 = k % n1
    if k_1 == 0:
        print("   k es multiplo del orden; prueba otro k.")
    else:
        Q1 = E1.multiplicar(k_1, G1)
        print(f"   orden n = {n1} (primo).  Q = kG = {Q1}")
        E1.sumas = 0
        t = time.perf_counter()
        k1, pasos = rho_pollard(E1, G1, Q1, n1)
        print(f"   Rho de Pollard: k = {k1}  (correcto: {k1 == k_1});  {pasos} pasos,"
              f" ~sqrt(pi n/2) = {int(math.sqrt(math.pi * n1 / 2))};  {time.perf_counter() - t:.2f} s")

    print("\n3. Curva de 32 bits con orden SUAVE (producto de primos pequenos: mala eleccion):")
    suave = lambda N: max(factorizar(N)) < 1000
    E2, G2, n2 = buscar_curva(32, suave, semilla=8)
    k_2 = k % n2
    Q2 = E2.multiplicar(k_2, G2)
    escrito = " * ".join(f"{p}^{e}" if e > 1 else str(p) for p, e in factorizar(n2).items())
    print(f"   orden n = {n2} = {escrito}.  Q = kG = {Q2}")
    print("   Pohlig-Hellman resuelve un problema pequeno en cada subgrupo:")
    t = time.perf_counter()
    k2 = pohlig_hellman_con_detalle(E2, G2, Q2, n2)
    print(f"   y el TCR los junta: k = {k2}  (correcto: {k2 == k_2});  {time.perf_counter() - t:.3f} s")
    if k >= n2:
        print(f"   (Eve obtiene k mod n = {k_2}, que es todo lo que importa porque nG = O)")

    print("\nCon un orden primo de 256 bits, rho necesitaria unas 2^128 operaciones: imposible hoy.")


if __name__ == "__main__":
    datos = leer_datos(PREGUNTA)
    if datos:
        con_tus_datos(datos)
        print("\n(Para ver el ejemplo de clase y las pruebas:  python 12_rho_pollard_pohlig_hellman.py --ejemplo)")
    else:
        ejemplo_de_clase()
        print()
        pruebas()
