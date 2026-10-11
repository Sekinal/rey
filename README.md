# Taylor y termodinámica de Fe-Al: solucionario reproducible

Desarrollo en español de los incisos a–g del primer parcial de Matemáticas Aplicadas. Reconstruye el contenido matemático del solucionario proporcionado, con una nueva maquetación editable y una explicación desde cero, paso a paso. Calcula las aproximaciones de Taylor, sus errores, la respuesta del hierro mediante Gibbs–Duhem y la energía de Gibbs de exceso. Incluye una extensión al Taylor directo de la exponencial.

**Documento listo para leer:** [solucionario generado, 60 páginas](output/pdf/solucionario.pdf). Incluye procedencia de cada dato, reglas básicas, todas las derivadas de Taylor, sustituciones, factoriales, cuentas de error, integración de Gibbs–Duhem y los programas completos en Python y MATLAB comentados instrucción por instrucción.

**Para empezar con el código:** [Python didáctico completo](ejemplos/paso_a_paso.py) · [MATLAB completo](matlab/paso_a_paso.m).

## Reproducir

Se necesita `uv` (la ejecución original usó 0.11.15) y `pdflatex`. `uv` descarga Python 3.12.12 y usa las versiones fijadas en `uv.lock`. La primera ejecución requiere conexión para descargar dependencias.

En Debian o Ubuntu, instalar LaTeX y las herramientas de revisión:

```bash
sudo apt-get update
sudo apt-get install -y texlive-latex-base texlive-latex-recommended texlive-fonts-recommended lmodern poppler-utils
```

Desde la raíz de este repositorio:

```bash
uv sync --locked
uv run --locked pytest
uv run --locked feal
```

El último comando regenera todo, incluido `output/pdf/solucionario.pdf`. Para calcular y producir figuras sin instalar LaTeX:

```bash
uv run --locked feal --solo-datos
```

No se requiere recuperar archivos de `/home/ieqr/Downloads`: las tres entradas se incluyen en `referencias/`. Los comandos pueden ejecutarse de nuevo sin intervención manual. Se deben ejecutar desde un clon o copia del repositorio; el programa localiza la raíz a partir de su archivo fuente.

## Ejecutar los programas didácticos

Python, desde la raíz del repositorio:

```bash
uv sync --locked
uv run --locked python ejemplos/paso_a_paso.py
```

Genera tablas y cuatro gráficas en `output/didactico/`, sin necesitar LaTeX. Todas las funciones están incluidas; el script no importa el paquete interno `feal`.

MATLAB, con la raíz del repositorio como carpeta actual:

```matlab
addpath('matlab');
resultados = paso_a_paso;
resultados.metricas
```

Genera tablas, la malla completa y cuatro gráficas en `output/matlab/`. El archivo incluye todas las funciones auxiliares y no requiere toolboxes adicionales. La función devuelve también una estructura para explorar los valores en la consola. Los índices de MATLAB empiezan en 1: `C(k+1)` representa el coeficiente matemático `c_k`.

La compatibilidad se validó con **GNU Octave 7.3.0** en un contenedor Debian; no se ejecutó el motor propietario MATLAB, que no está instalado en este entorno. Para repetir esa comprobación:

```bash
docker build -f Dockerfile.octave -t rey-octave .
docker run --rm --user "$(id -u):$(id -g)" -v "$PWD:/proyecto" rey-octave
uv run --locked python pruebas/verificar_matlab.py
```

El verificador contrasta los 10 001 puntos y las 20 columnas calculadas, los 13 coeficientes de la serie infinita y las dos tablas de métricas. Las tolerancias están declaradas en el archivo del verificador. En sistemas sin `id`, puede omitirse `--user`; las instrucciones de MATLAB no requieren Docker.

## Archivos

| Ruta | Contenido |
|---|---|
| `referencias/` | Los tres PDF originales, sin modificaciones, y sus huellas SHA-256 |
| `src/feal/modelo.py` | Funciones numéricas y parámetros con unidades |
| `src/feal/reproducir.py` | Generación de tablas, gráficas y PDF |
| `src/feal/desglose.py` | Operaciones intermedias, recurrencia y registro de bisección |
| `ejemplos/paso_a_paso.py` | Programa Python completo, comentado por instrucción |
| `matlab/paso_a_paso.m` | Programa MATLAB completo, comentado por instrucción |
| `documentos/solucionario.tex` y `documentos/capitulos/` | Desarrollo completo en español, editable |
| `pruebas/test_modelo.py` | Comprobaciones algebraicas y numéricas independientes |
| `output/datos/` | Malla completa, métricas, composiciones, extensión y resumen JSON |
| `output/figuras/` | Figuras vectoriales PDF y vistas PNG |
| `output/tablas/` | Tablas y valores generados para LaTeX |
| `output/manifiesto.json` | Huellas SHA-256 de los resultados |

Los resultados intermedios se regeneran y se excluyen de Git. El PDF final se conserva en Git para facilitar su lectura.

## Convenciones científicas

- Temperatura: 1873 K. Composición: fracción molar de Al, `0 <= x <= 0.3872`. Malla de métricas: 10 001 puntos, con ambos extremos.
- Parámetros Redlich–Kister: `A = -4.6095`, `B = 0.741`, tabla 1 del artículo. La derivación exacta da `c2 = 2.0595`.
- Referencia de comparación: ecuación (26), `F3 = -3.8685 + 4.773*x + 2.06*x**2 - 2.964*x**3`. El redondeo publicado cambia la referencia en `0.0005*x**2`.
- Orden N significa N+1 términos. Los truncamientos conservan los mismos coeficientes: no se reajusta cada orden.
- Los errores se miden frente a un modelo, no frente a observaciones experimentales. El modelo binario no resuelve el equilibrio ternario ni su cinética.
- El 0.07 % de Al mencionado en el artículo es porcentaje en masa, no fracción molar. La validez física no se extiende automáticamente a toda la malla matemática.
- Cerca de cero se evita la cancelación en Gibbs–Duhem usando la cola de la serie del logaritmo. Los límites del 1 % se calculan mediante bisección.

## Validación

Las pruebas contrastan la derivación parcial molar con álgebra simbólica exacta, verifican Gibbs–Duhem para todos los órdenes y comparan la fórmula del hierro con cuadratura de Gauss–Legendre independiente. También comprueban las métricas del solucionario original, la recurrencia de la exponencial y las fronteras de tolerancia.

Validación de esta edición: 24 pruebas aprobadas, ejecución de los programas didácticos Python y MATLAB compatible, y comparación de todas las columnas de MATLAB/Octave con Python. El PDF se revisa visualmente después de compilar. La reproducción editorial se verifica también desde una copia limpia y mediante una segunda ejecución. La receta `Dockerfile.octave` se ejecutó; el `Dockerfile` de generación editorial se conserva como alternativa documentada.

Resultados de control:

| Orden | Error máximo de ln γAl | RMSE | Error relativo máximo de γAl |
|---|---:|---:|---:|
| 0 | 1.984887 | 1.141205 | 86.2604 % |
| 1 | 0.136781 | 0.074787 | 12.7839 % |
| 2 | 0.172062 | 0.065041 | 18.7751 % |
| 3 | 0 | 0 | 0 % |

Para inspeccionar visualmente el documento regenerado:

```bash
mkdir -p tmp/pdfs
pdftoppm -scale-to 1400 -png output/pdf/solucionario.pdf tmp/pdfs/pagina
```

## Auditar cada cuenta

Después de `uv run --locked feal`:

- `output/datos/desglose_malla.csv`: 40 004 filas, una por composición y orden, con potencias, cuatro productos, suma, referencia, diferencia, cuadrado del error y porcentaje.
- `output/datos/coeficientes_gamma.csv`: productos, suma y resultado de la recurrencia para los coeficientes del 0 al 12.
- `output/datos/terminos_gamma_extremo.csv`: cada potencia, término y suma acumulada del Taylor directo en `x = 0.3872`.
- `output/datos/biseccion.csv`: 70 iteraciones por orden, con extremos, punto medio, error y decisión.
- Las tablas insertadas en el PDF se generan desde estos mismos cálculos. Los números se redondean solo para mostrarlos.

La malla se elige para evaluar ecuaciones: sus puntos **no son mediciones experimentales**. Los coeficientes de la serie infinita se deducen del polinomio publicado mediante la recurrencia; no se toman de una tabla adicional del artículo.

## Reproducibilidad y automatización

Python y sus dependencias están fijados. Se fija `SOURCE_DATE_EPOCH` y se eliminan fechas variables de las figuras PDF; dos ejecuciones con las mismas herramientas producen los mismos resultados. `output/manifiesto.json` permite contrastar sus huellas. Excluye las carpetas `output/didactico/` y `output/matlab/`, que contienen ejecuciones complementarias; se contrastan numéricamente, no por bytes.

La versión y las fuentes de TeX Live dependen del sistema. No se garantiza igualdad byte por byte entre sistemas operativos o distribuciones diferentes de LaTeX; sí se fijan el modelo, los datos de entrada, la malla y las dependencias numéricas. El `Dockerfile` proporciona una receta Linux para instalar las herramientas; la imagen base y los paquetes del sistema no están fijados por digest, por lo que tampoco constituye una imagen archivada bit a bit.

```bash
docker build -t solucionario-feal .
docker run --rm -v "$PWD/output:/proyecto/output" solucionario-feal
```

La integración continua en `.github/workflows/verificar.yml` ejecuta pruebas, regenera el documento, compara dos ejecuciones y guarda los resultados como artefacto. Se activará al subir el repositorio a GitHub; no se configura ni publica un remoto automáticamente.

## Procedencia

1. Ishfaq, M.; Pande, M. M. (2025). *Thermodynamic analysis of Fe-Al-O system using Redlich-Kister (R-K) polynomial and Wagner’s Taylor series expansion*. The Journal of Chemical Thermodynamics, 202, 107415. DOI: [10.1016/j.jct.2024.107415](https://doi.org/10.1016/j.jct.2024.107415). Parámetros: tabla 1, p. 5; ecuaciones 16–26, p. 6; tabla 2, p. 7.
2. Universidad Veracruzana. *Primer Examen Parcial: Matemáticas Aplicadas*, 6 de octubre de 2026.
3. *Solucionario razonado: Series de Taylor y termodinámica del sistema Fe-Al-O*, PDF proporcionado, fechado el 7 de octubre de 2026.

Los PDF de referencia conservan los derechos de sus respectivos titulares. Su inclusión documenta las entradas de este repositorio local; no se les asigna una licencia nueva.
