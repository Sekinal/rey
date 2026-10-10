# Taylor y termodinámica de Fe-Al: solucionario reproducible

Desarrollo en español de los incisos a–g del primer parcial de Matemáticas Aplicadas. Reconstruye el contenido matemático del solucionario proporcionado, con una nueva maquetación editable. Calcula las aproximaciones de Taylor, sus errores, la respuesta del hierro mediante Gibbs–Duhem y la energía de Gibbs de exceso. Incluye una extensión al Taylor directo de la exponencial.

**Documento listo para leer:** [solucionario generado](output/pdf/solucionario.pdf).

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

## Archivos

| Ruta | Contenido |
|---|---|
| `referencias/` | Los tres PDF originales, sin modificaciones, y sus huellas SHA-256 |
| `src/feal/modelo.py` | Funciones numéricas y parámetros con unidades |
| `src/feal/reproducir.py` | Generación de tablas, gráficas y PDF |
| `documentos/solucionario.tex` | Desarrollo completo en español, editable |
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

Validación realizada: 16 pruebas aprobadas; PDF de siete páginas revisado visualmente; reconstrucción desde una copia limpia con un entorno virtual nuevo. Todas las huellas de los resultados coincidieron, incluido el PDF. La receta Docker y el flujo remoto de GitHub Actions se incluyen, pero no se ejecutaron durante esta validación local.

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

## Reproducibilidad y automatización

Python y sus dependencias están fijados. Se fija `SOURCE_DATE_EPOCH` y se eliminan fechas variables de las figuras PDF; dos ejecuciones con las mismas herramientas producen los mismos resultados. `output/manifiesto.json` permite contrastar sus huellas.

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
