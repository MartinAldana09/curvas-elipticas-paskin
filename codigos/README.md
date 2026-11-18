# Códigos del proyecto: de RSA a las curvas elípticas

Material computacional del informe de auxiliar de investigación del proyecto
**Grupo de Mordell y Cota de Chabauty-Coleman (2026)** — Grupo PROMENTE,
Fundación Universitaria Konrad Lorenz.

Autor: Martin Alexander Aldana Usma · Investigador principal: Wilmar Reinero Bolaños Chávez.

> **Aviso.** Todos los programas son **educativos**. Los ejemplos con parámetros
> pequeños (por ejemplo, la curva $y^2=x^3+24x+13$ sobre $\mathbb{F}_{29}$) **no son
> seguros**. Para proteger información real se deben usar bibliotecas criptográficas
> auditadas (por ejemplo `cryptography` u OpenSSL) y parámetros estandarizados
> (NIST SP 800-186, RFC 7748).

## Archivos (de menor a mayor complejidad)

| Archivo | Tema | Resultado que se verifica |
|---|---|---|
| `01_aritmetica_modular.py` | Euclides, Bézout, inverso, potencias rápidas, TCR, φ, fracciones continuas | `mcd(240,46)=2`, `3^13 mod 100 = 23` |
| `02_rsa.py` | RSA: claves, cifrado, descifrado, descifrado con TCR | `p=5, q=11, e=3 → d=27`, `7 → 13 → 7` |
| `03_ataques_rsa.py` | d ⇒ factorizar, módulo común, firma ciega, Wiener, Håstad, Pollard p−1, fallos | Wiener: `d=3511`; Håstad: `M=805121215` |
| `04_curvas_elipticas.py` | Discriminante, puntos de E(F_p), Hasse, suma, inverso, orden | `#E(F_29)=38`, `ord(G)=19`, axiomas de grupo |
| `05_multiplicacion_escalar.py` | Duplicar y sumar, escalera de Montgomery, ECDLP | escalera para `k=13` |
| `06_intercambio_claves_ECDH.py` | ECDH (curva de clase y P-256), hombre en el medio | secreto `15G=(0,19)` |
| `07_elgamal_eliptico.py` | Mensajes como puntos (Koblitz) y ElGamal elíptico | `C1=(23,28)`, `C2=(20,5)` |
| `08_firmas_ECDSA.py` | ECDSA, verificación, ataque por nonce repetido (también en P-256) | firma `(4,11)`, se recupera `d=7` |
| `09_factorizacion_lenstra.py` | Factorización con curvas elípticas (ECM) | `2773 = 59·47` |
| `10_coppersmith_lll.py` | Reducción LLL y teorema de Coppersmith: mensaje estereotipado y Håstad con relleno fijo | RSA-512, `e=3`: se recuperan 120 bits desconocidos |
| `11_criba_cuadratica.py` | Criba cuadrática: base de factores, criba, álgebra lineal sobre F₂ | `1649 = 17·97`; números de 19 y 25 cifras |
| `12_rho_pollard_pohlig_hellman.py` | Ataques al ECDLP: rho de Pollard y Pohlig–Hellman | orden primo vs. orden suave |

Cada archivo funciona **solo** (no depende de los demás) y trae comentarios con la
matemática. Tiene dos modos:

* **Con tus datos**: hace todo el proceso paso a paso con lo que tú escribas
  (tu nombre, tus números, tu curva...).
* **Ejemplo de clase**: muestra el ejemplo fijo del informe y ejecuta la función
  `pruebas()`, que comprueba los resultados esperados. Por eso este modo imprime
  siempre lo mismo: sus datos están escritos dentro del programa.

## Tres formas de ejecutar un programa

```bash
python 02_rsa.py "Martin Aldana"   # 1. tus datos después del nombre del archivo
python 02_rsa.py                   # 2. el programa te pide los datos (Enter = ejemplo de clase)
python 02_rsa.py --ejemplo         # 3. ejemplo de clase y pruebas, sin preguntar
python ejecutar_todo.py            #    los doce programas en modo ejemplo, con sus pruebas
```

Las comillas son necesarias cuando el texto tiene espacios. En Spyder, VS Code o
Jupyter (`%run 02_rsa.py`) funciona la forma 2: el programa pregunta y tú escribes.
Si pegas el código en una celda de Colab, también te pregunta.

## Qué escribir y qué muestra cada programa

| Programa | Qué escribes (ejemplo) | Qué te muestra con tus datos |
|---|---|---|
| `01_aritmetica_modular.py` | dos enteros `a n` y, si quieres, un exponente `e`: `240 46 13` | congruencia, Euclides fila por fila, Bézout, inverso de `a` módulo `n` (o por qué no existe), φ(n) y Euler, `a^e mod n` bit por bit, fracción continua de `a/n` |
| `02_rsa.py` | un texto: `"Martin Aldana"` | (A) RSA de juguete letra por letra con `N = 3233`, y por qué es inseguro; (B) RSA de 512 bits: `p`, `q`, `N`, φ, `e`, `d`, el texto convertido en número `m`, el cifrado `c`, el descifrado y el texto recuperado. Los primos son aleatorios: **cada ejecución da números distintos, pero siempre se recupera tu mensaje** |
| `03_ataques_rsa.py` | un texto de hasta 50 caracteres: `"Martin"` | cómo una atacante recupera tu mensaje o la clave con Håstad, módulo común, `d` filtrado, Wiener (`d` pequeño) y un fallo de hardware |
| `04_curvas_elipticas.py` | `a b p` y, si quieres, un punto `x y`: `2 3 97` o `24 13 29 1 3` | discriminante, todos los puntos, `#E`, Hasse, orden del punto, tabla `kG` y la cuenta completa de `2G` (tangente) y `3G` (secante) |
| `05_multiplicacion_escalar.py` | un entero `k`: `13` | `k` en binario, duplicar y sumar paso a paso, escalera de Montgomery y número de operaciones |
| `06_intercambio_claves_ECDH.py` | el secreto de Alice y el de Bob: `3 5` | `A = aG`, `B = bG`, el secreto común `abG`, cómo lo rompe Eve en la curva pequeña y una versión con P-256 |
| `07_elgamal_eliptico.py` | un texto: `"Martin Aldana"` | claves en P-256, el texto convertido en puntos (Koblitz), `C1`, `C2`, descifrado y por qué cifrar dos veces da resultados distintos |
| `08_firmas_ECDSA.py` | el texto a firmar | huella SHA-256, firma `(r, s)`, verificación, qué pasa si se cambia un carácter y cómo un nonce repetido entrega la clave privada |
| `09_factorizacion_lenstra.py` | un entero: `1000000016000000063` | qué curva encontró el factor, con qué cota `B`, y la factorización completa |
| `10_coppersmith_lll.py` | una clave de hasta 15 caracteres: `"Gato2026"` | Alice cifra "La clave de hoy es: tu clave" con `e = 3`; Coppersmith + LLL la descubre |
| `11_criba_cuadratica.py` | un entero: `1649` o `1000000016000000063` | la idea (valores `x² − N` suaves y un cuadrado), la criba y la factorización completa |
| `12_rho_pollard_pohlig_hellman.py` | un entero secreto `k`: `123456789` | rho de Pollard y Pohlig–Hellman recuperan tu `k` a partir de `Q = kG` |

El modo ejemplo del archivo 10 tarda cerca de un minuto (la reducción LLL se hace
en Python puro); con tus datos tarda segundos. Los demás tardan pocos segundos.

## Cómo ejecutar

**En el computador** (Python 3.8 o superior, sin paquetes adicionales): abra una
terminal en la carpeta `codigos` y use los comandos de arriba.

**En Google Colab** (sin instalar nada): abra el cuaderno con el enlace

```
https://colab.research.google.com/github/MartinAldana09/curvas-elipticas-paskin/blob/main/codigos/Experimentos_ECC.ipynb
```

y ejecute las celdas en orden. La primera celda descarga los programas
(en ella, `REPO` debe tener la dirección del repositorio). En cada celda hay una
variable para cambiar, por ejemplo:

```python
MENSAJE = "Martin Aldana"        # <- cambia esto
!python 02_rsa.py "{MENSAJE}"
```

También puede usar una celda nueva:

```python
!git clone https://github.com/MartinAldana09/curvas-elipticas-paskin repo
%cd repo/codigos
!python 02_rsa.py "Martin Aldana"
```

## Referencias principales

* W. Trappe y L. C. Washington, *Introduction to Cryptography with Coding Theory*, 2.ª ed., Pearson, 2006 (caps. 3, 6 y 16).
* D. Boneh, "Twenty Years of Attacks on the RSA Cryptosystem", *Notices of the AMS* 46(2), 1999.
* J. W. Bos *et al.*, "Elliptic Curve Cryptography in Practice", *Financial Cryptography* 2014.
* D. Coppersmith, "Small solutions to polynomial equations, and low exponent RSA vulnerabilities", *J. Cryptology* 10(4), 1997.
* C. Pomerance, "A Tale of Two Sieves", *Notices of the AMS* 43(12), 1996.
* L. C. Washington, *Elliptic Curves: Number Theory and Cryptography*, 2.ª ed., CRC, 2008.
