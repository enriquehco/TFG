classdef dataControler < handle
    %Esta clase es la que actúa como controlador, actualizando el modelo y
    %la vista.
    
    properties
        red;
        cargadoDatosRed;
        imagenBool;
        textArea;
        directorioActual;
    end
    
    methods
        %Constructor
        function obj = dataControler()
            obj.red = BrainNetwork(" ");
            obj.textArea = OutputBox("> Es necesario cargar unos datos...");
            obj.directorioActual = obj.red.directorioActual;
            
            if(~isempty(obj.red.filas))
                obj.cargadoDatosRed = 1;
            else
                obj.cargadoDatosRed = 0;
            end
        end
        
        %Funcion de carga de nuevos datos
        function cargados = cargaNuevosDatos(obj)
            redAuxiliar = obj.red;
            obj.red = BrainNetwork(obj.directorioActual);
            
            if(~isempty(obj.red.filas))
                obj.cargadoDatosRed = 1;
                obj.directorioActual = obj.red.directorioActual;
            else
                obj.cargadoDatosRed = 0;
                obj.red = redAuxiliar;
            end
            
            cargados = obj.cargadoDatosRed;
                
        end
        
        %Función para aplicar el umbral seleccionado a la red
        function  source = aplicaUmbral(obj,umbral)
            if(umbral >= 0 && umbral < 1)
                id = obj.centinelaImagen();
                source = obj.red.aplicaUmbral(umbral,id);
            else
                obj.textArea.addText("> El valor del umbral no es correcto");
            end
        end
        
        %Función para pintar la matriz
        function source = pintaMatriz(obj,id)
            if(~isempty(obj.red.matriz))
                source = obj.red.pintaMatriz(id);
            end
        end
        
        %Función para calcular las medidas de una red
        function calculaMedidas(obj)
            obj.red.calculaMedidas();
        end
        
        %Función para añadir texto al textArea
        function text = addText(obj,texto)
            if(~isempty(texto))
                text = obj.textArea.addText(texto);
            end
        end
        
        %Función para cargar el archivo de coordenadas 3D
        function cargaCoordenadas3D(obj)
             obj.red.cargaCoordenadas3D();
        end
        
        %Función para iniciar los campos de nodos y enlaces
        function campos = iniciaCamposNodosEnlaces(obj)
            campos = obj.red.iniciaCamposNodosEnlaces();
        end
        
        %Función para comprobar si se han cargado unas coordenadas
        function cargadas = coordenadasCargadas(obj)
            cargadas = ~isempty(obj.red.coordenadas3D);
        end
        
        %Función para visualizar la red en 2D
        function visualiza2D(obj)
            obj.red.visualiza2D();
        end
        
        %Función para visualizar la red en 3D
        function visualiza3D(obj)
            obj.red.visualiza3D();
        end
        
        %Función para obtener el nombre de los nodos
        function names = getNombreNodos(obj)
            names = obj.red.getNombreNodos();
        end
        
        %Calculo del grado medio de la red
        function mediumDegree = calculaGradoFeatureGlobal(obj)
            mediumDegree = obj.red.calculaGradoMedio();
        end
        
        %Calculo del coeficiente de clustering
        function coefClustering = calculaCoefClustering(obj)
            coefClustering = obj.red.calculaCoefClustering();
        end
        
        %Calculo del indice de modularidad
        function modularidad = calculaModularidad(obj)
            modularidad = obj.red.calculaModularidad();
        end
        
        %Calculo de la dimension fractal
        function dimFractal = calculaDimensionFractal(obj)
            dimFractal = obj.red.calculaDimensionFractal();
        end
        
        %Calculo del grado de un nodo
        function degree = calculaGradoNodo(obj,nodo)
            degree = obj.red.calculaGradoNodo(nodo);
        end
        
        %Calculo de la cercanía de un nodo
        function cercania = calculaCercaniaNodo(obj,nodo)
            cercania = obj.red.calculaCercaniaNodo(nodo);
        end
        
        %Calculo de la intermediación de un nodo
        function intermediacion = calculaIntermediacionNodo(obj,nodo)
            intermediacion = obj.red.calculaIntermediacionNodo(nodo);
        end
        
        %Calculo de la centralidad de vector propio de un nodo
        function vectorPropio = calculaVectorPropio(obj,nodo)
            vectorPropio = obj.red.calculaVectorPropio(nodo);
        end
        
        %Calculo de la escala de dimension fractal
        function calculaCompactBoxBurning(obj)
             tStart = cputime;
             obj.red.calculaCompactBoxBurning();
             tEnd = cputime - tStart;
        end
        
        %Calculo de la escala de dimension fractal con Greedy Coloring
        function calculaGreedyColoring(obj)
            obj.red.calculaGreedyColoring();
        end
        
        %Calculo de la escala de dimension fractal con Random Sequential
        function calculaRandomSequential(obj)
            obj.red.calculaRandomSequential();
        end
        
        %Calculo de la escala de dimension fractal con Merge algorithm
        function calculaMergeAlgorithm(obj)
            obj.red.calculaMergeAlgorithm();
        end
        
        %Calculo de la escala de dimension fractal con OBC algorithm
        function calculaOBCA(obj)
            obj.red.calculaOBCA();
        end
        
        %Calculo de la escala de dimension fractal con MEMB algorithm
        function calculaMEMB(obj)
            obj.red.calculaMEMB();
        end
        
        %Calculo de la escala de dimension fractal con REMCC algorithm
        function calculaREMCC(obj)
            obj.red.calculaREMCC();
        end
        
        %Calculo de la escala de dimension fractal con PSO algorithm
        function calculaPSO(obj)
            obj.red.calculaPSO();
        end
        
        %Calculo de la escala de dimension fractal con DE algorithm
        function calculaDE(obj)
            obj.red.calculaDE();
        end
        
        %Calculo de todos los algoritmos
        function calculaTodos(obj)
            obj.red.calculaTodos();
        end
        
        %Calculo de la dimensión fractal
        function dimFractal = calculaPendienteEscalaCajas(obj,lb_min,lb_max)
            dimFractal = obj.red.calculaPendienteEscalaCajas(lb_min,lb_max);
        end
        
        %Función para guardar las medidas globales
        function exportGlobalFeaturesToCSV(obj)
            obj.red.exportGlobalFeaturesToCSV();
        end
        
        %Función para guardar las medidas locales
        function exportLocalFeaturesToCSV(obj)
            obj.red.exportLocalFeaturesToCSV();
        end
        
        %Función para devolver el mensaje de informacion adicional
        function mensaje = informacionAdicional(obj)
            mensaje = obj.red.informacionAdicional();
        end
        
        %Función para devolver el mensaje de contacto
        function mensaje = contactaConmigo(obj)
            mensaje = obj.red.contactaConmigo();
        end
        
        %Función para restablecer los datos
        function image = restableceMatriz(obj)
             id = obj.centinelaImagen();
             image = obj.red.restableceMatriz(id);
        end
        
        %Función para controlar el display de las imagenes
        function id = centinelaImagen(obj)
            if(obj.imagenBool == 0)
               obj.imagenBool = 1; 
            elseif(obj.imagenBool == 1)
               obj.imagenBool = 2;
            else
                obj.imagenBool = 0; 
            end
            
            id = obj.imagenBool;
        end
        
        %Funcion para borrar las imagenes al cerrar la aplcicación
        function borraImagenes(obj)
            currentFolderArray = pwd;
            currentFolder = string(currentFolderArray);
            
            pos = strfind(currentFolder,"\");
            limit = size(pos);
            max = pos(limit(2));
            currentFolder = extractBetween(currentFolder,1,max);
            
            image0 = currentFolder + '\matrizAdyacencia0.png';
            image1 = currentFolder + '\matrizAdyacencia1.png';
            image2 = currentFolder + '\matrizAdyacencia2.png';
            
            if(isfile(image0))
                delete(image0);
            end
            
            if(isfile(image1))
                delete(image1);
            end
            
            if(isfile(image2))
                delete(image2);
            end
            
        end
        
        function restableceDatosCargados(obj)
            obj.cargadoDatosRed = 1;
        end
    end
end

