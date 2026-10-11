"""Comparar resultados de MATLAB/Octave con las ecuaciones Python verificadas.

Ejecutar después de matlab/paso_a_paso.m:
    uv run --locked python pruebas/verificar_matlab.py
"""
from pathlib import Path
import numpy as np
from feal.modelo import aluminio, hierro, energia, diferencia, coeficientes_gamma


def main():
    # La ruta se obtiene desde este archivo, no desde una carpeta personal.
    raiz = Path(__file__).resolve().parents[1]
    salida = raiz / "output" / "matlab"
    # skiprows=1 omite los nombres de columna y lee solo los números.
    datos = np.loadtxt(salida / "malla.csv", delimiter=",", skiprows=1)
    # Cada fila es una composición y cada columna una magnitud declarada.
    assert datos.shape == (10001, 21), "Dimensiones incorrectas en malla.csv"
    x = np.linspace(0, 0.3872, 10001)
    np.testing.assert_allclose(datos[:, 0], x, rtol=0, atol=1e-16)
    # Comparamos las cuatro filas matemáticas guardadas ahora como columnas.
    for orden in range(4):
        np.testing.assert_allclose(datos[:, 1+orden], aluminio(x, orden), rtol=1e-13, atol=2e-14)
        np.testing.assert_allclose(datos[:, 5+orden], hierro(x, orden), rtol=1e-11, atol=2e-14)
        np.testing.assert_allclose(datos[:, 9+orden], energia(x, orden), rtol=1e-11, atol=2e-9)
        np.testing.assert_allclose(datos[:, 13+orden], np.exp(aluminio(x, orden)), rtol=1e-12, atol=1e-14)
    b = coeficientes_gamma(12)
    coeficientes = np.loadtxt(salida / "coeficientes_gamma.csv", delimiter=",", skiprows=1)
    np.testing.assert_array_equal(coeficientes[:, 0], np.arange(13))
    np.testing.assert_allclose(coeficientes[:, 1], b, rtol=1e-12, atol=1e-13)
    errores_infinitos = []
    for posicion, grado in enumerate((3, 5, 8, 12)):
        esperado = np.polynomial.polynomial.polyval(x, b[:grado+1])
        np.testing.assert_allclose(datos[:, 17+posicion], esperado, rtol=1e-12, atol=1e-14)
        errores_infinitos.append([grado, np.max(100*abs(esperado/np.exp(aluminio(x))-1))])
    metricas_infinitas = np.loadtxt(salida / "taylor_gamma.csv", delimiter=",", skiprows=1)
    np.testing.assert_allclose(metricas_infinitas, errores_infinitos, rtol=1e-10, atol=1e-11)
    esperadas = []
    for orden in range(4):
        delta = diferencia(x, orden)
        esperadas.append([orden, np.max(abs(delta)), np.sqrt(np.mean(delta**2)),
                          np.max(100*abs(np.expm1(delta)))])
    metricas = np.loadtxt(salida / "metricas.csv", delimiter=",", skiprows=1)
    np.testing.assert_allclose(metricas, esperadas, rtol=1e-12, atol=1e-12)
    print("MATLAB/Octave coincide con Python: 10001 puntos, 20 columnas calculadas,")
    print("13 coeficientes de la serie infinita y las dos tablas de métricas.")


if __name__ == "__main__":
    main()
