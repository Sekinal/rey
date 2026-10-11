function resultados = paso_a_paso(carpeta_salida)
% SOLUCIONARIO FE-AL: programa completo, comentado y sin toolboxes adicionales.
% Desde la raiz del repositorio: addpath('matlab'); paso_a_paso
% Devuelve resultados y guarda CSV y cuatro figuras PNG.
% Los comentarios explican cada instruccion; el PDF desarrolla las ecuaciones.

% nargin cuenta los argumentos recibidos; cero significa usar la ruta habitual.
if nargin == 0
    % mfilename localiza este archivo; fileparts extrae la carpeta que lo contiene.
    carpeta_codigo = fileparts(mfilename('fullpath'));
    % La carpeta superior a matlab es la raiz del repositorio.
    raiz = fileparts(carpeta_codigo);
    % fullfile une nombres con el separador apropiado para cada sistema.
    carpeta_salida = fullfile(raiz, 'output', 'matlab');
% end cierra el bloque if.
end
% exist(...,'dir') comprueba si ya existe esa carpeta.
if ~exist(carpeta_salida, 'dir')
    % mkdir crea la carpeta; ~ significa negacion logica.
    mkdir(carpeta_salida);
end

% MATLAB cuenta desde 1: C(1)=c_0, C(2)=c_1, C(3)=c_2, C(4)=c_3.
C = [-3.8685, 4.773, 2.06, -2.964];
% Constante de los gases en J/(mol K).
R = 8.314462618;
% Temperatura fija en kelvin.
T = 1873;
% El limite de composicion procede de la tabla 2 del articulo.
x_max = 0.3872;
% linspace incluye ambos extremos; x es un vector fila de 10001 posiciones.
x = linspace(0, x_max, 10001);
% numel cuenta todas las posiciones del vector.
M = numel(x);
% Cada fila de P representa un orden y cada columna una composicion.
P = zeros(4, M);
% H almacenara los logaritmos del coeficiente del hierro.
H = zeros(4, M);
% G almacenara energia de exceso en J/mol.
G = zeros(4, M);
% Gamma almacenara el coeficiente del aluminio, no su logaritmo.
Gamma = zeros(4, M);
% Delta guardara P_N-F_3 con su signo.
Delta = zeros(4, M);
% Error_gamma guardara el error relativo porcentual del coeficiente.
Error_gamma = zeros(4, M);
% Columnas: orden, error maximo, RMSE y maximo error porcentual.
metricas = zeros(4, 4);

% 0:3 produce 0,1,2,3; aqui el orden matematico sigue empezando en cero.
for orden = 0:3
    % La fila fisica del arreglo es orden+1 porque MATLAB empieza en uno.
    fila = orden + 1;
    % Dos puntos selecciona todas las columnas de esa fila.
    P(fila, :) = logaritmo_aluminio(x, orden, C);
    % Integramos Gibbs-Duhem para el mismo orden, no siempre para el cubico.
    H(fila, :) = logaritmo_hierro(x, orden, C);
    % exp se aplica a cada componente del vector.
    Gamma(fila, :) = exp(P(fila, :));
    % .* multiplica elemento por elemento; RT convierte a J/mol.
    G(fila, :) = R*T*(x.*P(fila, :) + (1-x).*H(fila, :));
    % Calculamos los terminos omitidos para evitar restar numeros casi iguales.
    Delta(fila, :) = diferencia_logaritmo(x, orden, C);
    % expm1(z)=exp(z)-1; abs elimina el signo de cada resultado.
    Error_gamma(fila, :) = 100*abs(expm1(Delta(fila, :)));
    % max escoge el mayor error absoluto entre todas las composiciones.
    error_maximo = max(abs(Delta(fila, :)));
    % .^2 eleva cada error al cuadrado; ^2 seria una potencia matricial.
    cuadrados = Delta(fila, :).^2;
    % sum suma los 10001 cuadrados del vector.
    suma = sum(cuadrados);
    % Dividimos entre el numero de puntos, no entre el de intervalos.
    promedio = suma/M;
    % sqrt obtiene la raiz cuadrada del promedio.
    rmse = sqrt(promedio);
    % Este maximo corresponde al error de gamma en porcentaje.
    error_porcentual = max(Error_gamma(fila, :));
    % Los corchetes forman una fila de cuatro numeros.
    metricas(fila, :) = [orden, error_maximo, rmse, error_porcentual];
end

% Calculamos b_0,...,b_12; el resultado tendra 13 posiciones.
b = coeficientes_exponencial(12, C);
% Estos grados incluyen, respectivamente, 4,6,9 y 13 terminos.
grados = [3, 5, 8, 12];
% Cada fila de Q representa uno de esos cuatro grados.
Q = zeros(4, M);
% Guardaremos el grado y su error relativo maximo.
metricas_infinita = zeros(4, 2);
% La fila 4 corresponde al orden 3; es exp(F_3), la referencia.
referencia = Gamma(4, :);
% indice indica la posicion dentro de la lista grados.
for indice = 1:numel(grados)
    % Recuperamos el grado concreto: 3,5,8 o 12.
    grado = grados(indice);
    % Recorremos todos los terminos, incluido el constante n=0.
    for n = 0:grado
        % b(n+1) representa b_n; .^ aplica la potencia a cada composicion.
        Q(indice, :) = Q(indice, :) + b(n+1)*x.^n;
    end
    % ./ divide entrada por entrada; Q ya aproxima gamma directamente.
    error = 100*abs(Q(indice, :)./referencia - 1);
    % Una fila guarda el grado y el mayor porcentaje obtenido.
    metricas_infinita(indice, :) = [grado, max(error)];
end

% Construimos la ruta del primer CSV y escribimos sus columnas.
escribir_csv(fullfile(carpeta_salida, 'metricas.csv'), ...
    'orden,error_maximo,rmse,error_gamma_porcentaje', metricas);
% (0:12).' convierte un vector fila a columna sin conjugacion compleja.
escribir_csv(fullfile(carpeta_salida, 'coeficientes_gamma.csv'), ...
    'n,b_n', [(0:12).', b.']);
% Guardamos tambien los errores de la expansion directa de gamma.
escribir_csv(fullfile(carpeta_salida, 'taylor_gamma.csv'), ...
    'grado,error_max_porcentaje', metricas_infinita);
% Transponemos cada matriz para que una fila del CSV sea una composicion.
datos = [x.', P.', H.', G.', Gamma.', Q.'];
% La cabecera enumera las columnas en el mismo orden que datos.
cabecera = ['x,P0,P1,P2,P3,H0,H1,H2,H3,G0,G1,G2,G3,' ...
    'gamma0,gamma1,gamma2,gamma3,Q3,Q5,Q8,Q12'];
% Escribimos las 10001 filas calculadas para compararlas con Python.
escribir_csv(fullfile(carpeta_salida, 'malla.csv'), cabecera, datos);

% Las llaves crean una lista de textos, llamada arreglo de celdas.
leyendas = {'Orden 0', 'Orden 1', 'Orden 2', 'Orden 3'};
% Primera figura: logaritmos a la izquierda, coeficientes a la derecha.
graficar_par(x, P, Gamma, 'ln gamma del Al', 'gamma del Al', ...
    leyendas, fullfile(carpeta_salida, 'actividad.png'));
% Segunda figura: error absoluto y error relativo porcentual.
graficar_par(x, abs(Delta), Error_gamma, 'Error absoluto de ln gamma', ...
    'Error de gamma (%)', leyendas, fullfile(carpeta_salida, 'errores.png'));
% Tercera figura: dividir G entre 1000 convierte J/mol a kJ/mol.
graficar_par(x, H, G/1000, 'ln gamma del Fe', 'gE (kJ/mol)', ...
    leyendas, fullfile(carpeta_salida, 'termodinamica.png'));
% La cuarta figura necesita leyendas por grado, no por orden del logaritmo.
leyendas_q = {'Grado 3', 'Grado 5', 'Grado 8', 'Grado 12'};
% Creamos una matriz para los porcentajes de cada Q.
errores_q = zeros(size(Q));
% Recorremos las cuatro filas ya calculadas.
for indice = 1:4
    % El cociente se calcula por componentes, por eso se usa ./.
    errores_q(indice, :) = 100*abs(Q(indice, :)./referencia - 1);
end
% Dibujamos los Q y sus errores; la referencia se usa en el calculo del error.
graficar_par(x, Q, errores_q, 'Taylor directo de gamma', ...
    'Error relativo (%)', leyendas_q, ...
    fullfile(carpeta_salida, 'serie_infinita.png'));

% struct crea un contenedor con campos identificados por nombre.
resultados = struct();
% El punto accede a un campo del contenedor, no multiplica.
resultados.x = x;
% Guardamos cada magnitud para poder consultarla en la consola.
resultados.P = P;
resultados.H = H;
resultados.G = G;
resultados.Gamma = Gamma;
resultados.Q = Q;
resultados.b = b;
resultados.metricas = metricas;
% disp muestra valores sin alterar el calculo.
disp('Columnas: orden, error maximo, RMSE, error de gamma (%)');
disp(metricas);
% fprintf inserta la ruta en %s y termina con un salto de linea \n.
fprintf('Tablas y cuatro figuras guardadas en: %s\n', carpeta_salida);
% Fin de la funcion principal; las funciones auxiliares siguen en el mismo archivo.
end

function p = logaritmo_aluminio(x, orden, C)
% Implementa P_N=sumatoria de c_k*x^k.
% size(x) conserva el numero de posiciones y la orientacion de x.
p = zeros(size(x));
% Incluimos el termino k=0, que es la constante.
for k = 0:orden
    % .^ significa potencia por elemento.
    potencia = x.^k;
    % La posicion k+1 guarda el coeficiente matematico c_k.
    termino = C(k+1)*potencia;
    % Acumulamos los terminos sin sustituir la suma anterior.
    p = p + termino;
end
end

function delta = diferencia_logaritmo(x, orden, C)
% Implementa el negativo de la suma de terminos omitidos.
% La diferencia comienza en cero para cada composicion.
delta = zeros(size(x));
% (orden+1):3 es vacio si orden=3; no se ejecuta ninguna iteracion.
for k = (orden+1):3
    % Restar el termino cubico negativo produce una contribucion positiva.
    delta = delta - C(k+1)*x.^k;
end
end

function h = logaritmo_hierro(x, orden, C)
% Implementa la integral analitica de Gibbs-Duhem, con estabilizacion diluida.
% h empieza en cero, compatible con H_0=0.
h = zeros(size(x));
% No se incluye k=0 porque la derivada de una constante es cero.
for k = 1:orden
    % log1p(-x) calcula el logaritmo natural de 1-x.
    corchete = log1p(-x);
    % Añadimos los primeros k terminos positivos.
    for j = 1:k
        % Cada potencia se divide por su indice j.
        corchete = corchete + x.^j/j;
    end
    % La mascara logica marca donde x es menor que 0.1.
    pequeno = x < 0.1;
    % Seleccionamos unicamente esas posiciones para la cola estable.
    xs = x(pequeno);
    % Reservamos un cero para cada posicion seleccionada.
    cola = zeros(size(xs));
    % Sumamos hasta j=32; el PDF acota el resto omitido.
    for j = (k+1):32
        % Todos los terminos de esta cola son negativos o cero.
        cola = cola - xs.^j/j;
    end
    % Reemplazamos solo los puntos pequeños, conservando el resto.
    corchete(pequeno) = cola;
    % El factor k procede de derivar c_k*x^k.
    h = h + k*C(k+1)*corchete;
end
end

function b = coeficientes_exponencial(grado, C)
% Implementa la recurrencia para el Taylor directo de exp(F_3).
% El vector tiene grado+1 posiciones porque incluye b_0.
b = zeros(1, grado+1);
% La posicion 1 corresponde a b_0, no a b_1.
b(1) = exp(C(1));
% El contador n conserva el indice matematico del coeficiente buscado.
for n = 1:grado
    % Reiniciamos la suma para este nuevo coeficiente.
    suma = 0;
    % min evita pedir coeficientes de indice matematico negativo.
    for k = 1:min(3, n)
        % C(k+1)=c_k y b(n-k+1)=b_(n-k).
        producto = k*C(k+1)*b(n-k+1);
        % Incorporamos un producto a la suma.
        suma = suma + producto;
    end
    % La posicion n+1 corresponde a b_n; dividimos entre el n matematico.
    b(n+1) = suma/n;
end
end

function escribir_csv(ruta, cabecera, datos)
% Guarda numeros con 17 cifras significativas, utiles para comparaciones.
% fopen abre el archivo para escritura; -1 indica que no pudo abrirlo.
archivo = fopen(ruta, 'w');
% Detectamos un error de permisos o de ruta antes de escribir.
if archivo == -1
    % error interrumpe la ejecucion con un mensaje legible.
    error('No se pudo abrir el archivo: %s', ruta);
end
% onCleanup garantiza el cierre al salir, incluso si ocurre un error.
cerrar = onCleanup(@() fclose(archivo)); %#ok<NASGU>
% Escribimos los nombres de columnas y un salto de linea.
fprintf(archivo, '%s\n', cabecera);
% size(datos,1) cuenta filas; recorremos una por composicion o metrica.
for fila = 1:size(datos, 1)
    % El primer numero se escribe sin coma delante.
    fprintf(archivo, '%.17g', datos(fila, 1));
    % Las columnas siguientes si llevan una coma antes del numero.
    for columna = 2:size(datos, 2)
        % %.17g conserva 17 cifras significativas de cada numero.
        fprintf(archivo, ',%.17g', datos(fila, columna));
    end
    % Terminamos la fila despues de escribir todas sus columnas.
    fprintf(archivo, '\n');
end
end

function graficar_par(x, izquierda, derecha, titulo_izq, titulo_der, etiquetas, ruta)
% Dibuja y guarda dos paneles; no requiere interaccion con ventanas.
% Para curvas decrecientes, la leyenda abajo evita ocultar el inicio.
if strcmp(titulo_izq, 'ln gamma del Fe')
    % southwest significa esquina inferior izquierda.
    posicion_leyenda = 'southwest';
else
    % Para las curvas crecientes dejamos libre la zona inferior.
    posicion_leyenda = 'northwest';
end
% figure crea una figura invisible de tamaño suficiente para leer las etiquetas.
figura = figure('Visible', 'off', 'Position', [100, 100, 1100, 440]);
% Primer panel: una fila, dos columnas, posicion 1.
subplot(1, 2, 1);
% .' transpone sin conjugar; cada columna pasa a representar una curva.
plot(x, izquierda.', 'LineWidth', 1.2);
% xlabel indica el significado del eje horizontal.
xlabel('Fraccion molar de Al');
% ylabel indica la magnitud calculada del panel izquierdo.
ylabel(titulo_izq);
% legend asocia un texto con cada curva en el mismo orden.
legend(etiquetas, 'Location', posicion_leyenda);
% grid on agrega lineas de referencia visual.
grid on;
% Activamos el segundo panel de la misma figura.
subplot(1, 2, 2);
% Dibujamos los datos correspondientes al panel derecho.
plot(x, derecha.', 'LineWidth', 1.2);
% Ambos paneles tienen la misma variable horizontal.
xlabel('Fraccion molar de Al');
% El eje vertical cambia segun la propiedad.
ylabel(titulo_der);
% Esta leyenda usa los mismos ordenes o grados que la primera.
legend(etiquetas, 'Location', posicion_leyenda);
% La cuadricula no modifica las curvas.
grid on;
% print exporta PNG; -r150 establece la resolucion en puntos por pulgada.
print(figura, ruta, '-dpng', '-r150');
% Cerramos la figura para liberar memoria.
close(figura);
end
