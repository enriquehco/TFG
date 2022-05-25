classdef NodeFeatures < Features
    %Esta clase representa las medidas que se pueden tomar de cada uno de
    %los nodos que componen la red.
    
    properties
        grado;
        cercania;
        intermediacion;
        centrVectorPropio;
        id;
    end
    
    methods
        
        %Constructor, aquí usamos la propiedad de la clase padre abstracta
        %id para tener identificado cada nodo.
        function obj = NodeFeatures(identificador)
            obj.id = identificador;
        end
        
        %GRADO
        function grado = calculaGrado(obj,datos)
            obj.grado = 0;
            numEnlaces = 0;
            filas = size(datos);
            
            for i=obj.id:obj.id
                for j=1:filas(2)
                    if(datos(i,j) > 0)
                        numEnlaces = numEnlaces + 1;
                    end
                end
            end
            obj.grado = numEnlaces;
            grado = obj.grado;
            
        end
        
        %CERCANÍA
        function cercania = calculaCercaniaNodo(obj,datos)
            grafo = graph(datos);
            cercaniaNodos = centrality(grafo,'closeness');
            cercania = cercaniaNodos(obj.id);
            obj.cercania = cercania;
        end
        
        %INTERMEDIACIÓN
        function intermediacion = calculaIntermediacionNodo(obj,datos)
            grafo = graph(datos);
            intermediacionNodos = centrality(grafo,'betweenness');
            intermediacion = intermediacionNodos(obj.id);
            obj.intermediacion = intermediacion;
        end
        
        %CENTRALIDAD DE VECTOR PROPIO
        function vectorPropioCentralidad = calculaVectorPropio(obj,datos)
            grafo = graph(datos);
            centrVectorPropio = centrality(grafo,'eigenvector');
            obj.centrVectorPropio = centrVectorPropio(obj.id);
            vectorPropioCentralidad = centrVectorPropio(obj.id);
        end
        
    end
end

