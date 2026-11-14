"""
ejecutar_todo.py
================
Ejecuta los doce programas en orden y reporta si sus pruebas pasaron.

    python ejecutar_todo.py
"""
import pathlib
import subprocess
import sys
import time

carpeta = pathlib.Path(__file__).resolve().parent
archivos = sorted(carpeta.glob("[0-9][0-9]_*.py"))
fallos = 0
for archivo in archivos:
    inicio = time.perf_counter()
    r = subprocess.run([sys.executable, str(archivo)], capture_output=True, text=True)
    duracion = time.perf_counter() - inicio
    ok = r.returncode == 0 and "todas las pruebas pasaron" in r.stdout
    fallos += not ok
    print(f"{'OK   ' if ok else 'FALLO'} {archivo.name:38s} {duracion:6.1f} s")
    if not ok:
        print(r.stdout[-1500:], r.stderr[-1500:])
print(f"\n{len(archivos) - fallos} de {len(archivos)} programas correctos.")
sys.exit(1 if fallos else 0)
