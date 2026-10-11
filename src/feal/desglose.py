"""Cuentas intermedias para leer el desarrollo y auditar cada resultado.

No cambia los parámetros: reutiliza el modelo verificado. Los CSV conservan
los valores numéricos; las cadenas LaTeX se redondean solo para mostrarlas.
"""
import csv
from pathlib import Path
import numpy as np
from .modelo import C, R, T, X_MAX, aluminio, hierro, energia, diferencia, coeficientes_gamma


def generar_desglose(salida):
    """Crear tablas explicativas y registros completos de operaciones repetidas."""
    salida = Path(salida)
    datos = salida / "datos"
    tablas = salida / "tablas"

    def tex(nombre, texto):
        # UTF-8 conserva los acentos de las explicaciones.
        (tablas / nombre).write_text(texto, encoding="utf-8")

    def tabla(nombre, titulos, filas, formato):
        # Cada & separa columnas; dos barras inversas terminan una fila LaTeX.
        contenido = "\\begin{tabular}{" + formato + "}\\toprule\n"
        contenido += " & ".join(titulos) + r" \\ \midrule" + "\n"
        contenido += "\n".join(" & ".join(map(str, f)) + r" \\" for f in filas)
        contenido += "\n\\bottomrule\\end{tabular}\n"
        tex(nombre, contenido)

    def csv_tabla(nombre, titulos, filas):
        # newline evita líneas vacías extra al abrir el archivo en otros sistemas.
        with (datos / nombre).open("w", encoding="utf-8", newline="") as archivo:
            escritor = csv.writer(archivo)
            escritor.writerow(titulos)
            escritor.writerows(filas)

    puntos = [0.01, 0.1, 0.2, 0.3, X_MAX]
    filas = [[f"{x:.4f}"] + [f"{C[k]*x**k:.9f}" for k in range(4)] for x in puntos]
    tabla("productos.tex", ["$x$", "$c_0$", "$c_1x$", "$c_2x^2$", "$c_3x^3$"], filas, "rrrrr")

    x = np.linspace(0, X_MAX, 10001)
    # El generador entrega una fila a la vez para no guardar 40004 listas en memoria.
    def filas_malla():
        for valor in x:
            terminos = C * valor**np.arange(4)
            for orden in range(4):
                p = aluminio(valor, orden)
                d = diferencia(valor, orden)
                yield [valor, orden, valor**2, valor**3, *terminos, p, aluminio(valor),
                       d, abs(d), d*d, np.exp(p), 100*abs(np.expm1(d))]

    csv_tabla("desglose_malla.csv", ["x", "orden", "x2", "x3", "t0", "t1", "t2", "t3",
              "P_N", "F3", "delta", "error_absoluto", "error_cuadrado", "gamma_N",
              "error_gamma_porcentaje"], filas_malla())
    sumas = []
    for orden in range(4):
        delta = diferencia(x, orden)
        suma = float(np.sum(delta**2))
        promedio = suma / len(x)
        sumas.append([orden, f"{suma:.10f}", f"{promedio:.10f}", f"{np.sqrt(promedio):.10f}"])
    tabla("sumas_error.tex", ["Orden", "$S_N$", "$S_N/10001$", "$\\sqrt{S_N/10001}$"], sumas, "rrrr")

    lineas = []
    for orden in range(4):
        delta = float(diferencia(0.1, orden))
        exponencial = np.exp(delta)
        error = 100*abs(np.expm1(delta))
        lineas.append(r"\noindent Orden " + str(orden) + ":\n" + r"\["
            + rf"e_{{\gamma,{orden}}}=100|e^{{{delta:.6f}}}-1|"
            + rf"\simeq100|{exponencial:.10f}-1|\simeq {error:.8f}\,\%." + r"\]" + "\n")
    tex("errores_ejemplo.tex", "".join(lineas))

    # Desglose de la integración en una composición que se puede revisar a mano.
    valor = 0.1
    logaritmo = np.log1p(-valor)
    s1 = logaritmo + valor
    s2 = s1 + valor**2/2
    s3 = s2 + valor**3/3
    contribuciones = [C[1]*s1, 2*C[2]*s2, 3*C[3]*s3]
    h = float(hierro(valor))
    p = float(aluminio(valor))
    u = valor*p
    v = (1-valor)*h
    texto = rf"""
Primero, $1-x=1-0.1=0.9$ y $\ln(1-x)=\ln(0.9)\simeq{logaritmo:.12f}$.
Los tres corchetes, calculados uno después de otro, son
\begin{{align*}}
s_1&=\ln(0.9)+0.1\simeq{s1:.12f},\\
s_2&=s_1+0.1^2/2=s_1+0.005\simeq{s2:.12f},\\
s_3&=s_2+0.1^3/3=s_2+0.000333333333\simeq{s3:.12f}.
\end{{align*}}
Multiplicamos por los factores exteriores:
\[
4.773s_1\simeq{contribuciones[0]:.12f},\quad
4.12s_2\simeq{contribuciones[1]:.12f},\quad
-8.892s_3\simeq{contribuciones[2]:.12f}.
\]
Sumamos esas tres contribuciones: $H_3\simeq{h:.12f}$. Ahora ponderamos:
\[
xP_3=(0.1)({p:.6f})={u:.10f}.
\]\[
(1-x)H_3=(0.9)({h:.12f})\simeq{v:.12f}.
\]
La suma es $xP_3+(1-x)H_3\simeq{u+v:.12f}$. Finalmente,
\[
g_3^E=(15572.986483514)({u+v:.12f})\simeq{energia(valor):.8f}\ \mathrm{{J/mol}}.
\]
Dividimos entre 1000: $g_3^E\simeq{energia(valor)/1000:.8f}$ kJ/mol.
"""
    tex("hierro_ejemplo.tex", texto)

    valor = X_MAX
    terminos = C*valor**np.arange(4)
    parciales = np.cumsum(terminos)
    texto = rf"""
\paso{{1. Potencias y productos en $x={valor:.4f}$.}}
\[
x^2=({valor:.4f})({valor:.4f})={valor**2:.8f},\quad
x^3=({valor**2:.8f})({valor:.4f})={valor**3:.12f}.
\]
\begin{{align*}}
t_0&=-3.8685,\\
t_1&=4.773({valor:.4f})={terminos[1]:.7f},\\
t_2&=2.06({valor**2:.8f})={terminos[2]:.10f},\\
t_3&=-2.964({valor**3:.12f})={terminos[3]:.15f}.
\end{{align*}}
\paso{{2. Sumas acumuladas.}}
\begin{{align*}}
P_0&={parciales[0]:.7f},\\
P_1&={parciales[0]:.7f}+{terminos[1]:.7f}={parciales[1]:.7f},\\
P_2&={parciales[1]:.7f}+{terminos[2]:.10f}={parciales[2]:.10f},\\
P_3&={parciales[2]:.10f}+({terminos[3]:.15f})\simeq{parciales[3]:.12f}.
\end{{align*}}
\paso{{3. Exponencial y actividad.}}
\[
\gamma_3=e^{{{parciales[3]:.12f}}}\simeq{np.exp(parciales[3]):.10f},\quad
a_{{\Al,3}}=({valor:.4f})({np.exp(parciales[3]):.10f})\simeq{valor*np.exp(parciales[3]):.10f}.
\]
\paso{{4. Hierro y energía.}}
$1-x=0.6128$ y $\ln(0.6128)\simeq{np.log1p(-valor):.12f}$. Usando la forma agrupada,
\begin{{align*}}
H_3&=0.001[\ln(0.6128)+0.3872]-2.386({valor**2:.8f})-2.964({valor**3:.12f})\\
&\simeq{0.001*(np.log1p(-valor)+valor):.12f}-{2.386*valor**2:.12f}-{2.964*valor**3:.12f}\\
&\simeq{hierro(valor):.12f}.
\end{{align*}}
Entonces $xP_3\simeq{valor*parciales[3]:.12f}$ y $(1-x)H_3\simeq{(1-valor)*hierro(valor):.12f}$.
Su suma es ${valor*parciales[3]+(1-valor)*hierro(valor):.12f}$. Multiplicando por $RT$:
\[
g_3^E\simeq{energia(valor):.8f}\ \mathrm{{J/mol}}={energia(valor)/1000:.8f}\ \mathrm{{kJ/mol}}.
\]
"""
    tex("extremo.tex", texto)

    iteraciones = []
    for orden in (1, 2):
        inferior, superior = 0.0, X_MAX
        for paso in range(1, 71):
            medio = (inferior+superior)/2
            error = 100*abs(np.expm1(diferencia(medio, orden)))
            aceptado = error <= 1
            iteraciones.append([orden, paso, inferior, superior, medio, float(error),
                                "subir inferior" if aceptado else "bajar superior"])
            if aceptado:
                inferior = medio
            else:
                superior = medio
    csv_tabla("biseccion.csv", ["orden", "paso", "inferior", "superior", "medio", "error_porcentaje", "decision"], iteraciones)
    tabla("biseccion.tex", ["Paso", "$l$", "$u$", "$m$", "Error (\\%)", "Decisión"],
          [[paso, f"{l:.8f}", f"{u:.8f}", f"{m:.8f}", f"{e:.5f}", d]
           for orden,paso,l,u,m,e,d in iteraciones if orden == 1 and paso <= 6], "rrrrrl")

    b = coeficientes_gamma(12)
    desarrollo = [rf"$b_0=e^{{-3.8685}}\simeq{b[0]:.12f}$. Este es el punto de partida." + "\n"]
    filas_coeficientes = [[0, 0., 0., 0., b[0], b[0]]]
    for n in range(1, 13):
        # Para n=1 o n=2 se omiten los productos con índices negativos.
        productos = [k*C[k]*b[n-k] for k in range(1, min(3,n)+1)]
        desarrollo.append(r"\paso{Coeficiente $b_{"+str(n)+r"}$:}" + "\n")
        desarrollo.append(r"\begin{align*}" + "\n")
        for k, producto in enumerate(productos, start=1):
            desarrollo.append(rf"u_{k}&=({k})({C[k]:g})b_{{{n-k}}}"
                + rf"\simeq({k})({C[k]:g})({b[n-k]:.12f})\simeq {producto:.12f},\\" + "\n")
        desarrollo.append(rf"s&="+"+".join(f"({p:.12f})" for p in productos)
                          +rf"\simeq {sum(productos):.12f},\\"+"\n")
        desarrollo.append(rf"b_{{{n}}}&=s/{n}\simeq({sum(productos):.12f})/{n}\simeq {b[n]:.12f}."+"\n")
        desarrollo.append(r"\end{align*}"+"\n")
        if n in (3, 6, 9):
            desarrollo.append(r"\clearpage"+"\n")
        rellenos = productos + [0.]*(3-len(productos))
        filas_coeficientes.append([n, *rellenos, sum(productos), b[n]])
    tex("recurrencia_desglosada.tex", "".join(desarrollo))
    csv_tabla("coeficientes_gamma.csv", ["n", "producto_1", "producto_2", "producto_3", "suma_o_b0", "b_n"], filas_coeficientes)
    terminos = b*X_MAX**np.arange(13)
    sumas = np.cumsum(terminos)
    tabla("terminos_gamma.tex", ["$n$", "$b_n$", "$x^n$", "$b_nx^n$", "$Q_n$"],
          [[n, f"{b[n]:.10f}", f"{X_MAX**n:.10f}", f"{terminos[n]:.10f}", f"{sumas[n]:.10f}"]
           for n in range(13)], "rrrrr")
    csv_tabla("terminos_gamma_extremo.csv", ["n", "b_n", "x_potencia_n", "termino", "suma_acumulada"],
              [[n,b[n],X_MAX**n,terminos[n],sumas[n]] for n in range(13)])
