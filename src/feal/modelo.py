"""Funciones puras; x es fracción molar y las energías se expresan en J/mol."""
import numpy as np

A, B = -4.6095, 0.741
R, T = 8.314462618, 1873.0
X_MAX = 0.3872
C = np.array([-3.8685, 4.773, 2.06, -2.964])
C_RK = np.array([A+B, -2*A-6*B, A+9*B, -4*B])


def validar(x, orden):
    x = np.asarray(x, dtype=float)
    if not isinstance(orden, (int, np.integer)) or not 0 <= orden <= 3:
        raise ValueError("El orden debe ser un entero entre 0 y 3.")
    if np.any(~np.isfinite(x)) or np.any((x < 0) | (x >= 1)):
        raise ValueError("La fracción molar debe ser finita y cumplir 0 <= x < 1.")
    return x


def aluminio(x, orden=3, coeficientes=C):
    """P_N = ln(gamma_Al); el orden N retiene N+1 términos."""
    x = validar(x, orden)
    return np.polynomial.polynomial.polyval(x, coeficientes[:orden+1])


def resto_logaritmico(x, k):
    """ln(1-x) + sum(x**j/j), evitando cancelación cerca de cero."""
    directo = np.log1p(-x) + sum(x**j/j for j in range(1, k+1))
    # Para x < 0.1, la cola converge rápidamente; el error omitido es < 1e-32.
    pequeno = np.where(x < 0.1, x, 0.0)
    cola = -sum(pequeno**j/j for j in range(k+1, 33))
    return np.where(x < 0.1, cola, directo)


def hierro(x, orden=3, coeficientes=C):
    """H_N obtenido por Gibbs-Duhem con H_N(0)=0."""
    x = validar(x, orden)
    return sum((k*coeficientes[k]*resto_logaritmico(x, k)
                for k in range(1, orden+1)), start=np.zeros_like(x))


def energia(x, orden=3):
    x = validar(x, orden)
    return R*T*(x*aluminio(x, orden)+(1-x)*hierro(x, orden))


def diferencia(x, orden):
    """P_N-F_3 calculado con los términos omitidos, sin restar valores próximos."""
    x = validar(x, orden)
    return -sum((C[k]*x**k for k in range(orden+1, 4)), start=np.zeros_like(x))


def coeficientes_gamma(grado):
    """Recurrencia del Taylor directo de exp(F_3), no exp(P_N)."""
    if not isinstance(grado, int) or grado < 0:
        raise ValueError("El grado debe ser un entero no negativo.")
    b = np.zeros(grado+1)
    b[0] = np.exp(C[0])
    for n in range(1, grado+1):
        b[n] = sum(k*C[k]*b[n-k] for k in range(1, min(3,n)+1))/n
    return b


def limite_uno_por_ciento(orden):
    """Primera frontera del intervalo contiguo desde cero con error <= 1 %."""
    lo, hi = 0.0, X_MAX
    for _ in range(70):
        medio = (lo+hi)/2
        if 100*abs(np.expm1(diferencia(medio, orden))) <= 1:
            lo = medio
        else:
            hi = medio
    return lo
