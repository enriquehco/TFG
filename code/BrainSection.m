classdef BrainSection < handle
    %Esta clase representa a cada uno de los nodos de la matriz de
    %adyacencia que a su vez representan a las distintas partes del cerebro
    
    properties
        id;
        features;
        coorX;
        coorY;
        coorZ;
        nombre;
    end
    
    methods
        %Constructor
        function obj = BrainSection(identificador,datos)
            obj.id = identificador;
            obj.features = NodeFeatures(identificador);
            
        end
        
        %Función para obtener el id
        function id = get.id(obj)
            id = obj.id;
        end
        
        %Función para cargar las coordenadas X, Y, Z además del nombre
        function cargaCoordenada3D(obj, datos)
            indiceDatosFila = obj.id + 1;
            obj.coorX = datos(indiceDatosFila,3);
            obj.coorY = datos(indiceDatosFila,4);
            obj.coorZ = datos(indiceDatosFila,5);
            obj.nombre = datos(indiceDatosFila,2);
        end
        
        %Función para obtener el grado medio del nodo
        function degree = calculaGradoNodo(obj,datos)
            degree = obj.features.calculaGrado(datos);
        end
        
        %Función para obtener la cercania de un nodo
        function cercania = calculaCercaniaNodo(obj,datos)
            cercania = obj.features.calculaCercaniaNodo(datos);
        end
        
        %Función para obtener la intermediación de un nodo
        function intermediacion = calculaIntermediacionNodo(obj,datos)
            intermediacion = obj.features.calculaIntermediacionNodo(datos);
        end
        
        %Función para obtener la centraliad de vector propio de un nodo
        function vectorPropio = calculaVectorPropio(obj,datos)
            vectorPropio = obj.features.calculaVectorPropio(datos);
        end
    end
end

