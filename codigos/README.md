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

Cada archivo funciona **solo** (no depende de los demás), trae comentarios con la
matemática, imprime un ejemplo y termina con una función `pruebas()` que comprueba
los resultados esperados. El archivo 10 tarda menos de un minuto
(la reducción LLL se hace en Python puro); los demás, pocos segundos.

## Cómo ejecutar

**En el computador** (Python 3.8 o superior, sin paquetes adicionales):

```bash
python 04_curvas_elipticas.py      # un archivo
python ejecutar_todo.py            # los doce, con sus pruebas
```

**En Google Colab** (sin instalar nada): abra el cuaderno con el enlace

```
https://colab.research.google.com/github/MartinAldana09/curvas-elipticas-paskin/blob/main/codigos/Experimentos_ECC.ipynb
```

y ejecute las celdas en orden. La primera celda descarga los programas
(en ella, `REPO` debe tener la dirección del repositorio). También puede usar
una celda nueva:

```python
!git clone https://github.com/MartinAldana09/curvas-elipticas-paskin repo
%cd repo/codigos
!python 04_curvas_elipticas.py
```

## Referencias principales

* W. Trappe y L. C. Washington, *Introduction to Cryptography with Coding Theory*, 2.ª ed., Pearson, 2006 (caps. 3, 6 y 16).
* D. Boneh, "Twenty Years of Attacks on the RSA Cryptosystem", *Notices of the AMS* 46(2), 1999.
* J. W. Bos *et al.*, "Elliptic Curve Cryptography in Practice", *Financial Cryptography* 2014.
* D. Coppersmith, "Small solutions to polynomial equations, and low exponent RSA vulnerabilities", *J. Cryptology* 10(4), 1997.
* C. Pomerance, "A Tale of Two Sieves", *Notices of the AMS* 43(12), 1996.
* L. C. Washington, *Elliptic Curves: Number Theory and Cryptography*, 2.ª ed., CRC, 2008.
