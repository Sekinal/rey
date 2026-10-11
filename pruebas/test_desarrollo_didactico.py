"""Comprobar que la explicación ejecutable conserva los resultados científicos."""
import csv
import importlib.util
from pathlib import Path
import numpy as np
import pytest
from feal.modelo import aluminio, hierro, energia, diferencia, coeficientes_gamma
from feal.desglose import generar_desglose


@pytest.fixture(scope="module")
def ejemplo():
    ruta = Path(__file__).resolve().parents[1] / "ejemplos/paso_a_paso.py"
    especificacion = importlib.util.spec_from_file_location("ejemplo", ruta)
    modulo = importlib.util.module_from_spec(especificacion)
    especificacion.loader.exec_module(modulo)
    return modulo.calcular()


@pytest.mark.parametrize("orden", range(4))
def test_programa_comentado_contra_modelo(ejemplo, orden):
    x = ejemplo["x"]
    # Métodos de evaluación diferentes pueden redondear en distinto orden.
    np.testing.assert_allclose(ejemplo["p"][orden], aluminio(x, orden), atol=2e-15)
    np.testing.assert_allclose(ejemplo["h"][orden], hierro(x, orden), atol=2e-15)
    np.testing.assert_allclose(ejemplo["g"][orden], energia(x, orden), atol=2e-11)
    np.testing.assert_allclose(ejemplo["delta"][orden], diferencia(x, orden), atol=1e-15)
    np.testing.assert_allclose(ejemplo["gamma"][orden], np.exp(aluminio(x, orden)), atol=1e-15)


def test_extension_comentada(ejemplo):
    np.testing.assert_allclose(ejemplo["b"], coeficientes_gamma(12), atol=1e-14)
    for grado, valores in ejemplo["q"].items():
        esperado = np.polynomial.polynomial.polyval(ejemplo["x"], ejemplo["b"][:grado+1])
        np.testing.assert_allclose(valores, esperado, atol=1e-15)


@pytest.fixture(scope="module")
def desgloses(tmp_path_factory):
    salida = tmp_path_factory.mktemp("desgloses")
    (salida / "datos").mkdir()
    (salida / "tablas").mkdir()
    generar_desglose(salida)
    return salida


def test_operaciones_en_todas_las_filas(desgloses):
    with (desgloses / "datos/desglose_malla.csv").open() as archivo:
        filas = list(csv.DictReader(archivo))
    assert len(filas) == 40004
    # Cada fila debe ser reconstruible con sus cuatro productos explícitos.
    for fila in filas:
        orden = int(fila["orden"])
        p = sum(float(fila[f"t{k}"]) for k in range(orden+1))
        assert p == pytest.approx(float(fila["P_N"]), abs=2e-15)
        delta = float(fila["delta"])
        assert delta**2 == pytest.approx(float(fila["error_cuadrado"]), abs=1e-15)
    for orden in range(4):
        seleccion = [f for f in filas if int(f["orden"]) == orden]
        cuadrados = np.array([float(f["error_cuadrado"]) for f in seleccion])
        x = np.array([float(f["x"]) for f in seleccion])
        np.testing.assert_allclose(cuadrados, diferencia(x,orden)**2, atol=1e-15)


def test_recurrencia_desglosada(desgloses):
    with (desgloses / "datos/coeficientes_gamma.csv").open() as archivo:
        filas = list(csv.DictReader(archivo))
    assert len(filas) == 13
    for fila in filas[1:]:
        suma = sum(float(fila[f"producto_{k}"]) for k in (1,2,3))
        assert suma == pytest.approx(float(fila["suma_o_b0"]),abs=1e-13)
        assert suma/int(fila["n"]) == pytest.approx(float(fila["b_n"]),abs=1e-13)


def test_biseccion_registrada(desgloses):
    with (desgloses / "datos/biseccion.csv").open() as archivo:
        filas = list(csv.DictReader(archivo))
    assert len(filas) == 140
    for orden in (1,2):
        seleccion = [f for f in filas if int(f["orden"]) == orden]
        for previo, siguiente in zip(seleccion, seleccion[1:]):
            if float(previo["error_porcentaje"]) <= 1:
                assert siguiente["inferior"] == previo["medio"]
                assert siguiente["superior"] == previo["superior"]
            else:
                assert siguiente["superior"] == previo["medio"]
                assert siguiente["inferior"] == previo["inferior"]
