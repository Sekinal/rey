"""Programa completo y comentado: Taylor de Fe-Al, tablas y cuatro figuras.

Ejecutar desde la raíz: uv run --locked python ejemplos/paso_a_paso.py
No requiere importar feal ni instalar LaTeX para generar sus resultados.
"""

# Path permite construir rutas sin depender de las barras de cada sistema.
from pathlib import Path
# csv escribe tablas que se pueden abrir en Excel, MATLAB o un editor.
import csv
# NumPy trabaja con muchos números a la vez; lo llamaremos np.
import numpy as np
# Matplotlib dibuja las gráficas; primero elegimos un modo sin ventana.
import matplotlib
# Agg permite guardar imágenes incluso en una computadora sin escritorio.
matplotlib.use("Agg")
# pyplot contiene las instrucciones para crear y guardar figuras.
import matplotlib.pyplot as plt

# Cada posición guarda un coeficiente: C[0], C[1], C[2], C[3].
C = np.array([-3.8685, 4.773, 2.06, -2.964])
# La constante R convierte la energía adimensional a energía por mol.
R = 8.314462618
# Todos los cálculos corresponden a esta temperatura en kelvin.
T = 1873.0
# El extremo procede de la tabla 2; no significa 38.72 % en masa.
X_MAX = 0.3872


# def define una función; no ejecuta todavía las líneas de su interior.
def calcular_logaritmo(x, orden):
    """Sumar, sin ocultar las potencias, P_N = suma de c_k*x**k."""
    # Creamos un arreglo de ceros con una entrada por composición.
    resultado = np.zeros_like(x, dtype=float)
    # range termina antes de su límite: orden+1 incluye el propio orden.
    for k in range(orden + 1):
        # ** es potencia; el cálculo se hace para todas las entradas de x.
        potencia = x**k
        # C[k] selecciona el coeficiente asociado a esa potencia.
        termino = C[k] * potencia
        # Acumulamos el nuevo término sin borrar los anteriores.
        resultado = resultado + termino
    # return entrega el arreglo calculado a quien llamó a la función.
    return resultado


def calcular_diferencia(x, orden):
    """Calcular P_N-F_3 como el negativo de los términos omitidos."""
    # Al principio todavía no hemos sumado ninguna diferencia.
    delta = np.zeros_like(x, dtype=float)
    # Los omitidos comienzan en orden+1 y terminan en k=3.
    for k in range(orden + 1, 4):
        # Restamos cada término ausente: delta = menos el resto de Taylor.
        delta = delta - C[k] * x**k
    # Si orden=3, el ciclo no se ejecuta y devuelve ceros.
    return delta


def calcular_hierro(x, orden):
    """Integrar Gibbs-Duhem con la condición H_N(0)=0."""
    # Habrá un resultado H por cada composición.
    h = np.zeros_like(x, dtype=float)
    # El término constante no contribuye porque su derivada es cero.
    for k in range(1, orden + 1):
        # log1p(-x) significa logaritmo natural de 1-x.
        corchete = np.log1p(-x)
        # Añadimos x, x**2/2, ... hasta x**k/k.
        for j in range(1, k + 1):
            # / divide; esta línea añade un solo término del corchete.
            corchete = corchete + x**j / j
        # Marcamos las posiciones en las que conviene evitar cancelación.
        pequeno = x < 0.1
        # La cola estable solo se evaluará en las composiciones marcadas.
        cola = np.zeros_like(x[pequeno])
        # range(k+1,33) incluye las potencias desde k+1 hasta 32.
        for j in range(k + 1, 33):
            # La cola tiene signo negativo: -suma de x**j/j.
            cola = cola - x[pequeno]**j / j
        # Sustituimos únicamente las posiciones pequeñas por la cola.
        corchete[pequeno] = cola
        # k*C[k] procede de derivar c_k*x**k antes de integrar.
        h = h + k * C[k] * corchete
    # Para orden 0, el ciclo queda vacío y h permanece igual a cero.
    return h


def calcular_coeficientes_gamma(grado):
    """Obtener b_0 hasta b_grado mediante la recurrencia deducida."""
    # Necesitamos grado+1 espacios porque empezamos a contar desde cero.
    b = np.zeros(grado + 1)
    # np.exp calcula e elevado al argumento; aquí da gamma en cero.
    b[0] = np.exp(C[0])
    # Comenzamos en 1 porque b[0] ya está calculado.
    for n in range(1, grado + 1):
        # Reiniciamos la suma para cada coeficiente nuevo.
        suma = 0.0
        # min impide usar índices negativos cuando n vale 1 o 2.
        for k in range(1, min(3, n) + 1):
            # b[n-k] ya existe porque su índice es menor que n.
            producto = k * C[k] * b[n-k]
            # Incorporamos uno de los productos de la recurrencia.
            suma = suma + producto
        # La división entre n viene del término n*b_n de la derivada.
        b[n] = suma / n
    # Entregamos todos los coeficientes, no solo el último.
    return b


def calcular():
    """Calcular todos los números antes de guardarlos o dibujarlos."""
    # La malla incluye ambos extremos y contiene 10001 números.
    x = np.linspace(0.0, X_MAX, 10001)
    # Una matriz de 4 filas almacena los órdenes 0, 1, 2 y 3.
    p = np.zeros((4, len(x)))
    # h tendrá el logaritmo del hierro, con la misma organización.
    h = np.zeros_like(p)
    # g guardará la energía de exceso en J/mol.
    g = np.zeros_like(p)
    # gamma guarda el coeficiente del Al, distinto de su logaritmo.
    gamma = np.zeros_like(p)
    # delta guarda el error con signo frente a F_3.
    delta = np.zeros_like(p)
    # error_gamma guarda porcentajes, no fracciones.
    error_gamma = np.zeros_like(p)
    # Las métricas tienen columnas: orden, máximo, RMSE, máximo porcentual.
    metricas = []
    # Ejecutamos exactamente cuatro iteraciones: orden=0,1,2,3.
    for orden in range(4):
        # Una fila completa recibe P_N para todos los x.
        p[orden] = calcular_logaritmo(x, orden)
        # Esta fila recibe H_N compatible con el mismo orden.
        h[orden] = calcular_hierro(x, orden)
        # exp se aplica a cada elemento y recupera el coeficiente.
        gamma[orden] = np.exp(p[orden])
        # Multiplicamos cada especie por su fracción y después por RT.
        g[orden] = R * T * (x*p[orden] + (1-x)*h[orden])
        # Calculamos la diferencia mediante los términos omitidos.
        delta[orden] = calcular_diferencia(x, orden)
        # expm1(z) evalúa exp(z)-1 con precisión para z pequeño.
        error_gamma[orden] = 100 * np.abs(np.expm1(delta[orden]))
        # abs elimina el signo y max elige el mayor valor de la fila.
        maximo = np.max(np.abs(delta[orden]))
        # Primero elevamos al cuadrado cada error.
        cuadrados = delta[orden]**2
        # Sumamos todos los cuadrados.
        suma = np.sum(cuadrados)
        # len(x) es el número de puntos, 10001, no el número de intervalos.
        promedio = suma / len(x)
        # sqrt calcula la raíz cuadrada del promedio.
        rmse = np.sqrt(promedio)
        # Este máximo corresponde a porcentajes del coeficiente.
        maximo_gamma = np.max(error_gamma[orden])
        # append añade una fila al final de la lista de métricas.
        metricas.append([orden, maximo, rmse, maximo_gamma])
    # Para la extensión infinita calculamos 13 coeficientes, del 0 al 12.
    b = calcular_coeficientes_gamma(12)
    # Un diccionario relacionará cada grado con su arreglo de valores.
    q = {}
    # Estos son los cuatro grados que compararemos.
    for grado in (3, 5, 8, 12):
        # Cada polinomio comienza con una suma de ceros.
        aproximacion = np.zeros_like(x)
        # Se incluyen todos los términos desde n=0 hasta n=grado.
        for n in range(grado + 1):
            # Es el Taylor directo de gamma, no la exponencial de P_N.
            aproximacion = aproximacion + b[n] * x**n
        # Guardamos el arreglo bajo su grado: q[3], q[5], q[8] o q[12].
        q[grado] = aproximacion
    # Entregamos un diccionario para poder verificar cada resultado.
    return {"x": x, "p": p, "h": h, "g": g, "gamma": gamma,
            "delta": delta, "error_gamma": error_gamma,
            "metricas": metricas, "b": b, "q": q}


def guardar_tabla(ruta, cabecera, filas):
    """Escribir una tabla CSV con nombres de columnas en español."""
    # with cierra el archivo incluso si ocurre un error durante la escritura.
    with ruta.open("w", encoding="utf-8", newline="") as archivo:
        # writer prepara un escritor que separa las columnas mediante comas.
        escritor = csv.writer(archivo)
        # La primera fila contiene los nombres de las columnas.
        escritor.writerow(cabecera)
        # Después escribimos una fila de datos por cada elemento de filas.
        escritor.writerows(filas)


def graficar_pares(x, izquierda, derecha, nombre, etiquetas, salida):
    """Dibujar dos paneles; cada entrada de la lista será una curva."""
    # subplots crea una figura con una fila y dos columnas de ejes.
    figura, ejes = plt.subplots(1, 2, figsize=(10, 4), layout="constrained")
    # enumerate entrega el número de orden y su arreglo de valores.
    for orden, valores in enumerate(izquierda):
        # ejes[0] es el panel izquierdo; plot une pares de puntos.
        ejes[0].plot(x, valores, label=f"Orden {orden}")
    # Repetimos para el panel derecho, ejes[1].
    for orden, valores in enumerate(derecha):
        # La etiqueta se usará en la leyenda de colores.
        ejes[1].plot(x, valores, label=f"Orden {orden}")
    # Asignamos un texto diferente a cada eje vertical.
    for indice in range(2):
        # Ambos paneles comparten el significado de x.
        ejes[indice].set_xlabel("Fracción molar de Al")
        # etiquetas[0] y etiquetas[1] describen magnitudes diferentes.
        ejes[indice].set_ylabel(etiquetas[indice])
        # La cuadrícula ayuda a leer valores sin modificar los datos.
        ejes[indice].grid(alpha=0.2)
        # Mostramos la relación entre colores y órdenes.
        ejes[indice].legend()
    # / combina una carpeta con un nombre de archivo usando Path.
    figura.savefig(salida / nombre, dpi=180)
    # Cerramos esta figura para liberar memoria antes de crear otra.
    plt.close(figura)


def main():
    """Ejecutar cálculos, guardar tablas y producir todas las gráficas."""
    # __file__ es este archivo; parents[1] es la raíz del repositorio.
    raiz = Path(__file__).resolve().parents[1]
    # Los resultados didácticos no sustituyen los de la edición del PDF.
    salida = raiz / "output" / "didactico"
    # parents crea carpetas intermedias; exist_ok permite repetir el comando.
    salida.mkdir(parents=True, exist_ok=True)
    # Ahora sí ejecutamos las operaciones definidas en calcular().
    resultados = calcular()
    # Extraemos la malla del diccionario usando su clave de texto.
    x = resultados["x"]
    # Guardamos las métricas de todos los órdenes en un archivo.
    guardar_tabla(salida / "metricas.csv",
                  ["orden", "error_maximo", "rmse", "error_gamma_porcentaje"],
                  resultados["metricas"])
    # zip forma parejas (índice,coeficiente): una por cada b_n.
    guardar_tabla(salida / "coeficientes_gamma.csv", ["n", "b_n"],
                  zip(range(13), resultados["b"]))
    # El primer panel muestra logaritmos y el segundo sus exponenciales.
    graficar_pares(x, resultados["p"], resultados["gamma"], "actividad.png",
                   ["ln gamma del Al", "gamma del Al"], salida)
    # El error absoluto se obtiene tomando abs de la matriz de diferencias.
    errores = np.abs(resultados["delta"])
    # Aquí se ven el error del logaritmo y el porcentaje del coeficiente.
    graficar_pares(x, errores, resultados["error_gamma"], "errores.png",
                   ["Error absoluto de ln gamma", "Error de gamma (%)"], salida)
    # Dividir entre 1000 convierte las energías de J/mol a kJ/mol.
    graficar_pares(x, resultados["h"], resultados["g"] / 1000,
                   "termodinamica.png", ["ln gamma del Fe", "gE (kJ/mol)"], salida)
    # Para la serie infinita usamos una figura adicional con dos paneles.
    figura, ejes = plt.subplots(1, 2, figsize=(10, 4), layout="constrained")
    # La fila 3 de gamma es exp(F_3), nuestra referencia.
    referencia = resultados["gamma"][3]
    # Dibujamos primero esa referencia en el panel izquierdo.
    ejes[0].plot(x, referencia, color="black", label="Referencia exp(F3)")
    # items entrega cada pareja grado,arreglo del diccionario q.
    for grado, valores in resultados["q"].items():
        # Mostramos el polinomio directo Q_m en el primer panel.
        ejes[0].plot(x, valores, label=f"Grado {grado}")
        # El error directo usa cociente; valores ya aproxima gamma.
        error = 100 * np.abs(valores / referencia - 1)
        # El segundo panel muestra su error relativo porcentual.
        ejes[1].plot(x, error, label=f"Grado {grado}")
    # Añadimos etiquetas horizontales, cuadrícula y leyenda a ambos paneles.
    for eje in ejes:
        # eje representa uno de los paneles en cada iteración.
        eje.set_xlabel("Fracción molar de Al")
        # Una cuadrícula tenue no oculta las curvas.
        eje.grid(alpha=0.2)
        # legend utiliza los textos label de las llamadas anteriores.
        eje.legend(fontsize=8)
    # El primer eje vertical es el valor del coeficiente.
    ejes[0].set_ylabel("Taylor directo de gamma")
    # El segundo eje vertical es su error porcentual.
    ejes[1].set_ylabel("Error relativo (%)")
    # Guardamos la cuarta imagen; no hace falta abrir una ventana.
    figura.savefig(salida / "serie_infinita.png", dpi=180)
    # Liberamos también esta figura.
    plt.close(figura)
    # Imprimimos una tabla pequeña en la terminal para comprobar la ejecución.
    for fila in resultados["metricas"]:
        # print no modifica los cálculos: solo muestra una fila.
        print(fila)
    # Indicamos dónde quedaron todos los archivos generados.
    print(f"Tablas y cuatro gráficas guardadas en: {salida}")


# Esta condición es cierta al ejecutar el archivo y falsa al importarlo.
if __name__ == "__main__":
    # Llamamos al procedimiento que realiza el trabajo completo.
    main()
