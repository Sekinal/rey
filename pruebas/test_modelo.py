"""Verificaciones independientes de las ecuaciones y resultados del solucionario."""
import numpy as np
import pytest
import sympy as sp
from feal.modelo import (C, C_RK, X_MAX, aluminio, hierro, energia, diferencia,
                          coeficientes_gamma, limite_uno_por_ciento)


def test_derivacion_parcial_molar_y_redondeo():
    x = sp.Symbol("x")
    A, B = sp.Rational("-4.6095"), sp.Rational("0.741")
    f = x*(1-x)*(A+B*(1-2*x))
    al = sp.expand(f+(1-x)*sp.diff(f,x))
    fe = sp.expand(f-x*sp.diff(f,x))
    assert al == sum(sp.Rational(v)*x**k for k,v in enumerate(
        ["-3.8685","4.773","2.0595","-2.964"]))
    assert fe == -sp.Rational("2.3865")*x**2-sp.Rational("2.964")*x**3
    assert al.subs(x,1) == 0
    malla = np.linspace(0,X_MAX,100)
    np.testing.assert_allclose(aluminio(malla)-aluminio(malla,coeficientes=C_RK),0.0005*malla**2,atol=1e-15)


@pytest.mark.parametrize("orden",range(4))
def test_gibbs_duhem_simbolico(orden):
    x=sp.Symbol("x")
    c=list(map(sp.Rational,["-3.8685","4.773","2.06","-2.964"]))
    p=sum(c[k]*x**k for k in range(orden+1))
    h=sum(k*c[k]*(sp.log(1-x)+sum(x**j/j for j in range(1,k+1))) for k in range(1,orden+1))
    assert sp.simplify((1-x)*sp.diff(h,x)+x*sp.diff(p,x)) == 0
    assert hierro(0,orden) == 0


@pytest.mark.parametrize("orden",[1,2,3])
def test_hierro_con_cuadratura_independiente(orden):
    nodos,pesos=np.polynomial.legendre.leggauss(64)
    for x in [1e-9,1e-4,0.1,X_MAX]:
        t=(nodos+1)*x/2
        dp=sum(k*C[k]*t**(k-1) for k in range(1,orden+1))
        integral=-x/2*np.sum(pesos*t/(1-t)*dp)
        np.testing.assert_allclose(hierro(x,orden),integral,rtol=2e-12,atol=1e-30)


def test_metricas_publicadas():
    x=np.linspace(0,X_MAX,10001)
    esperados=[(1.984887,1.141205,86.2604),(0.136781,0.074787,12.7839),
               (0.172062,0.065041,18.7751),(0,0,0)]
    for n,(emax,rmse,egamma) in enumerate(esperados):
        d=diferencia(x,n)
        assert np.max(abs(d)) == pytest.approx(emax,abs=5e-7)
        assert np.sqrt(np.mean(d*d)) == pytest.approx(rmse,abs=5e-7)
        assert np.max(100*abs(np.expm1(d))) == pytest.approx(egamma,abs=5e-5)
    assert energia(X_MAX)/1000 == pytest.approx(-16.41466,abs=5e-6)


def test_recurrencia_exponencial_con_serie_simbolica():
    x=sp.Symbol("x")
    # Se separa la constante para simplificar el desarrollo simbólico.
    serie=sp.exp(sum(sp.Rational(str(C[k]))*x**k for k in range(1,4))).series(x,0,9).removeO()
    esperado=[float(serie.coeff(x,k))*np.exp(C[0]) for k in range(9)]
    np.testing.assert_allclose(coeficientes_gamma(8),esperado,rtol=1e-12,atol=1e-14)


@pytest.mark.parametrize("orden",[1,2])
def test_frontera_tolerancia(orden):
    limite=limite_uno_por_ciento(orden)
    assert 100*abs(np.expm1(diferencia(limite,orden))) == pytest.approx(1,abs=1e-12)
    assert 100*abs(np.expm1(diferencia(limite+1e-6,orden))) > 1


@pytest.mark.parametrize("x",[-0.1,1,np.nan,np.inf])
def test_dominio(x):
    with pytest.raises(ValueError):
        hierro(x)
