classdef BrainNetwork < handle
    % Esta clase representa la matriz de adyacencia construida a partir
    % de los datos EEG.
    
    properties
        filas
        columnas
        matriz;
        matrizOriginal;
        features;
        nombreRed;
        brainSections;
        coordenadas3D;
        umbralImagen;
        directorioActual;
    end
    
    methods
        %Constructor de red aleatoria con valores entre 0 y 1
        function obj = BrainNetwork(directorio)
            if(directorio == " ")
                [filename,pathname] = uigetfile('*.csv', 'Selecciona el archivo de datos');
            else
                [filename,pathname] = uigetfile('*.csv', 'Selecciona el archivo de datos',directorio);
            end
            
            obj.directorioActual = pathname;
            rutaArchivo = strcat(pathname,filename);
            ruta = convertCharsToStrings(rutaArchivo);
            obj.umbralImagen = 0;
            
            if(filename ~= 0)
                obj.matriz = readtable(ruta);
                
                [obj.filas,obj.columnas] = size(obj.matriz);
            
                obj.matriz = obj.matriz{:,:};
                
                for i=1:obj.filas
                    brainSec = BrainSection(i,obj.matriz);
                    obj.brainSections = [obj.brainSections,brainSec]; 
                end
            
                obj.nombreRed = filename;
                obj.features = NetworkFeatures(filename);
                obj.matrizOriginal = obj.matriz;
            end
            
            
        end
        
        %Método para obtener matriz
        function data = get.matriz(obj)
           data = obj.matriz; 
        end
        
        %Método para obtener número de nodos
        function numNodos = get.filas(obj)
           numNodos = obj.filas; 
        end
        
        %Método para aplicar un umbral a la matriz 
        function source = aplicaUmbral(obj,valorUmbral,id)
            
            for i=1:obj.filas
                for j=1:obj.columnas
                    if(obj.matriz(i,j) < valorUmbral)
                        obj.matriz(i,j) = 0;
                    end
                end
            end
            
            source = obj.pintaMatriz(id);
           
        end
        
        %Método para pintar la matriz de adyacencia
        function source = pintaMatriz(obj,id)
            set(0,'DefaultFigureVisible','off');
            
            f1 = figure;
            imagesc(obj.matriz);
            axis square;    
            colorbar;
            
            currentFolder = pwd;
            pos = strfind(currentFolder,"\");
            limit = size(pos);
            currentFolder = string(currentFolder);
            max = pos(limit(2));
            currentFolder = extractBetween(currentFolder,1,max);
            
            if(id == 0)
                saveas(f1,fullfile(currentFolder, 'matrizAdyacencia0'),'png');
                source = strcat(currentFolder,"\matrizAdyacencia0.png");
                obj.umbralImagen = 1;
                pause(2);
            elseif(id==1)
                saveas(f1,fullfile(currentFolder, 'matrizAdyacencia1'),'png');
                source = strcat(currentFolder,"\matrizAdyacencia1.png");
                obj.umbralImagen = 0;
                pause(2);
            else
                saveas(f1,fullfile(currentFolder, 'matrizAdyacencia2'),'png');
                source = strcat(currentFolder,"\matrizAdyacencia2.png");
                obj.umbralImagen = 1;
                pause(2);
            end
            
            
            set(0,'DefaultFigureVisible','on');
        end
        
        %Función para calcular el coeficiente de clustering
        function medida = calculaCoefClustering(obj)
            medida = obj.features.calculaCoefClustering(obj.matriz);
        end
        
        %Función para calcular las medidas de la red
        function medida = calculaGradoMedio(obj)
            medida = obj.features.calculaGradoMedio(obj.matriz);
        end
        
        %Función para calcular el indice de modularidad
        function modularidad = calculaModularidad(obj)
            modularidad = obj.features.calculaModularidad(obj.matriz);
        end
        
        %Función para calcular la dimension fractal
        function calculaCompactBoxBurning(obj)
            nodos = obj.getNombreNodos();
            limit = size(nodos);
            for i=1:limit(2)
                nombreNodo = nodos(1,i);
                nodos(2,i) = obj.correspondenciaNodoId(nombreNodo);
            end
            
            obj.features.calculaCompactBoxBurning(obj.matriz,nodos);
        end
        
        %Funcion para calcular la dimension fractal (con Coloreado greedy)
        function calculaGreedyColoring(obj)
            nodos = obj.getNombreNodos();
            limit = size(nodos);
            nodosfin = [1:limit(2)];
            nodosfin = nodosfin.';
%             for i=1:limit(2)
%                 nombreNodo = nodos(1,i);
%                 nodos(2,i) = obj.correspondenciaNodoId(nombreNodo);
%             end
            obj.features.calculaGreedyColoring(obj.matriz,nodosfin);
        end
        
        %Funcion para calcular la dimension fractal (con Random Sequential)
        function calculaRandomSequential(obj)
            nodos = obj.getNombreNodos();
            limit = size(nodos);
            for i=1:limit(2)
                nombreNodo = nodos(1,i);
                nodos(2,i) = obj.correspondenciaNodoId(nombreNodo);
            end
            
            obj.features.calculaRandomSequential(obj.matriz,nodos);
        end
        
        %Funcion para calcular la dimension fractal (con Merge algorithm)
        function calculaMergeAlgorithm(obj)
            nodos = obj.getNombreNodos();
            limit = size(nodos);
            othernodes = [1:1:limit(2)];
            for i=1:limit(2)
                nombreNodo = nodos(1,i);
                nodos(2,i) = i;
            end
            obj.features.calculaMergeAlgorithm(obj.matriz,othernodes);
        end
        
        %Funcion para calcular la dimension fractal (con OBCA)
        function calculaOBCA(obj)
            nodos = obj.getNombreNodos();
            limit = size(nodos);
            othernodes = [1:1:limit(2)];
            obj.features.calculaOBCA(obj.matriz,othernodes);
        end
        
        %Funcion para calcular la dimension fractal (con OBCA)
        function calculaMEMB(obj)
            nodos = obj.getNombreNodos();
            limit = size(nodos);
            othernodes = [1:1:limit(2)];
            obj.features.calculaMEMB(obj.matriz,othernodes);
        end
        
        %Función para cargar el archivo de coordenadas 3D
        function cargaCoordenadas3D(obj)
            if(isempty(obj.coordenadas3D))
                [filename,pathname] = uigetfile('*.pos', 'Abre un archivo de coordenadas 3D',obj.directorioActual);
                rutaArchivo = strcat(pathname,filename);
                
                if(rutaArchivo ~= "")
                    ruta = convertCharsToStrings(rutaArchivo);
                    fileID = fopen(ruta,'r');
                    coordenadasArchivo = textscan(fileID,'%d8 %s %f32 %f32 %f32');
                    fclose(fileID);
                    limit = size(coordenadasArchivo);
                    limit = limit(2);
                
                    for i=1:limit
                        limit2 = size(coordenadasArchivo{:,i});
                        limitCoordenadas = limit2(1);
                    
                        for j=1:limitCoordenadas
                            aux = coordenadasArchivo{:,i}; 
                            obj.coordenadas3D{j,i} = aux(j,1);
                        end
                    end
                
                    limit3 = size(obj.brainSections);
                    limitBrainSections = limit3(2);
                    for i=1:limitBrainSections
                        obj.brainSections(1,i).cargaCoordenada3D(obj.coordenadas3D);
                    end
                end
                
            else
                errordlg('Ya hay unas coordenadas 3D cargadas','Fatal Error');
            end
            
        end
        
        %Función para contar el número de enlaces
        function enlaces = cuentaEnlaces(obj)
            enlaces = 0;
            limit = 0;
            for i=1:obj.filas
                for j=1:limit
                    if(obj.matriz(i,j) > 0)
                        enlaces = enlaces + 1;
                    end
                end
                limit = limit + 1;
            end
            
        end
        
        %Función para contar el número de nodos
        function nodos = cuentaNodos(obj)
            nodos = 0;
            for i=1:obj.filas
                for j=1:obj.columnas
                    if(obj.matriz(i,j) ~= 0)
                        nodos = nodos + 1;
                        break;
                    end
                end
            end
            
        end
        
        %Función para iniciar los campos de nodos y enlaces
        function campos = iniciaCamposNodosEnlaces(obj)
            campos(1,1) = obj.cuentaNodos;
            campos(1,2) = obj.cuentaEnlaces;
        end
        
        function rgb = double2rgb(obj,img, map, bounds, varargin)

            if nargin < 2
                error('Need to specify a colormap');
            end

            % ensure map is a numeric array
            if ischar(map)
                map = feval(map, 256);
            end

            % extract background value
            bgColor = [1 1 1];
            if ~isempty(varargin)
                bgColor = parseColor(varargin{1});
            end

            % get valid pixels (finite value)
            valid = isfinite(img);
            if ~exist('bounds', 'var') || isempty(bounds)
               bounds = [min(img(valid)) max(img(valid))];
            end

            % convert finite values to indices between 1 and map length
            n = size(map, 1);
            inds = (img(valid) - bounds(1)) / (bounds(end) - bounds(1)) * (n-1);
            inds = floor(min(max(inds, 0), n-1))+1;

            % compute the 3 bands
            dim = size(img);
            r = ones(dim) * bgColor(1); r(valid) = map(inds, 1);
            g = ones(dim) * bgColor(2); g(valid) = map(inds, 2);
            b = ones(dim) * bgColor(3); b(valid) = map(inds, 3);

            % concatenate the 3 bands to form an rgb image
            if length(dim) == 2
                % case of 2D image
                rgb = cat(3, r, g, b);
    
            else
                % case of 3D image: need to play with channels
                dim2 = [dim(1:2) 3 dim(3:end)];
                rgb = zeros(dim2, class(map));
                rgb(:,:,1,:) = r;
                rgb(:,:,2,:) = g;
                rgb(:,:,3,:) = b;
            end
        end

        function color = parseColor(color)

            if ischar(color)
               switch(color)
                  case 'k'
                     color = [0 0 0];
                 case 'w'
                     color = [1 1 1];
                 case 'r'
                     color = [1 0 0];
                  case 'g'
                      color = [0 1 0];
                 case 'b'
                     color = [0 0 1];
                 case 'c'
                     color = [0 1 1];
                 case 'm'
                       color = [1 0 1];
                 case 'y'
                       color = [1 1 0];
               otherwise 
                 error('Unknown color string');
               end
            end
        end
        
        %Función para visualizar la red en 2D
        function visualiza2D(obj)
           
            grafo = graph(obj.matriz);
            links = table2array(grafo.Edges);
            wij = links(:,3);
            a = 1;
            b = 3;
            wij = ( (wij - min(wij)).*(b-a) )/(max(wij) - min(wij));
            wij = wij+0.5;
            rgb = squeeze(obj.double2rgb(wij, colormap(jet))); close;
            
            a = 1;            
            b = 3;    
            wij = ( (wij - min(wij)).*(b-a) )/(max(wij) - min(wij));
            wij = wij+0.5;
            
            figure
            hold on;
            
            for ll = 1 : length(links)  
                node1X = obj.brainSections(links(ll,1)).coorX{1,1};
                node1Y = obj.brainSections(links(ll,1)).coorY{1,1};
                node2X = obj.brainSections(links(ll,2)).coorX{1,1};
                node2Y = obj.brainSections(links(ll,2)).coorY{1,1};
                line([node1X,node2X], [node1Y,node2Y], 'LineWidth', wij(ll), 'Color', rgb(ll,:));
            end
            
            colormap(jet)
            colorbar
            
            for nn = 1 : obj.filas
                x = obj.brainSections(nn).coorX{1,1};
                y = obj.brainSections(nn).coorY{1,1};
                name = string(obj.brainSections(nn).nombre);
                grado = obj.brainSections(nn).calculaGradoNodo(obj.matriz);
                
                if(grado ~= 0)
                   plot(x, y, 'ro', 'MarkerSize', grado, 'MarkerFaceColor', 'y')
                   text('Position', [x+0.5, y], 'String', name, 'FontSize', 11); 
                end
                
            end
            
            view([270 90]);            
            axis equal; axis off;
            set(gcf, 'units','normalized','outerposition',[0 0 1 1]) 
        end
        
        %Función para visualizar la red en 3D
        function visualiza3D(obj)
           
            grafo = graph(obj.matriz);
            links = table2array(grafo.Edges);
            wij = links(:,3);
            a = 1;
            b = 3;
            wij = ( (wij - min(wij)).*(b-a) )/(max(wij) - min(wij));
            wij = wij+0.5;
            rgb = squeeze(obj.double2rgb(wij, colormap(jet))); close;
            
            a = 1;            
            b = 3;    
            wij = ( (wij - min(wij)).*(b-a) )/(max(wij) - min(wij));
            wij = wij+0.5;
            
            figure
            hold on;
            
            for ll = 1 : length(links)  
                node1X = obj.brainSections(links(ll,1)).coorX{1,1};
                node1Y = obj.brainSections(links(ll,1)).coorY{1,1};
                node1Z = obj.brainSections(links(ll,1)).coorZ{1,1};
                node2X = obj.brainSections(links(ll,2)).coorX{1,1};
                node2Y = obj.brainSections(links(ll,2)).coorY{1,1};
                node2Z = obj.brainSections(links(ll,2)).coorZ{1,1};
                line([node1X,node2X], [node1Y,node2Y], [node1Z,node2Z], 'LineWidth', wij(ll), 'Color', rgb(ll,:));
            end
            
            colormap(jet)
            colorbar
            
            for nn = 1 : obj.filas
                x = obj.brainSections(nn).coorX{1,1};
                y = obj.brainSections(nn).coorY{1,1};
                z = obj.brainSections(nn).coorZ{1,1};
                name = string(obj.brainSections(nn).nombre);
                grado = obj.brainSections(nn).calculaGradoNodo(obj.matriz);
                
                if(grado ~= 0)
                   plot3(x, y, z, 'ro', 'MarkerSize', grado, 'MarkerFaceColor', 'y')
                   text('Position', [x+0.5, y, z], 'String', name, 'FontSize', 11);
                end
                
            end
            
            view([270 90]);            
            axis equal; axis off;
            set(gcf, 'units','normalized','outerposition',[0 0 1 1]) 
        end
        
        %Función para obtener el nombre de los nodos
        function names = getNombreNodos(obj)
            aux = size(obj.brainSections);
            limit = aux(2);
            for nn = 1 : limit
                names(nn) = string(obj.brainSections(nn).nombre);
            end
        end
           
        
        %Función para obtener el id de un nodo según su nombre
        function idNodo = correspondenciaNodoId(obj,nodo)
            nodo = string(nodo);
            
            for i=1:obj.filas
                nameCell = obj.brainSections(i).nombre;
                nameNodo = string(nameCell);
                if(nameNodo == nodo)
                    idNodo = obj.brainSections(i).id;
                end
            end
        end
        
        %Calculo del grado de un nodo
        function degree = calculaGradoNodo(obj,nodo)
            idNodo = obj.correspondenciaNodoId(nodo);
            degree = obj.brainSections(idNodo).calculaGradoNodo(obj.matriz);
        end
        
        %Calculo de la cercanía de un nodo
        function cercania = calculaCercaniaNodo(obj,nodo)
            idNodo = obj.correspondenciaNodoId(nodo);
            cercania = obj.brainSections(idNodo).calculaCercaniaNodo(obj.matriz);
        end
        
        %Calculo de la intermediación de un nodo
        function intermediacion = calculaIntermediacionNodo(obj,nodo)
            idNodo = obj.correspondenciaNodoId(nodo);
            intermediacion = obj.brainSections(idNodo).calculaIntermediacionNodo(obj.matriz);
        end
        
        %Calculo de la centralidad de vector propio de un nodo
        function vectorPropio = calculaVectorPropio(obj,nodo)
            idNodo = obj.correspondenciaNodoId(nodo);
            vectorPropio = obj.brainSections(idNodo).calculaVectorPropio(obj.matriz);
        end
        
        %Función para guardar las medidas globales
        function exportGlobalFeaturesToCSV(obj)
            grado = obj.calculaGradoMedio();
            clustering = obj.calculaCoefClustering();
            modularidad = obj.calculaModularidad();
            GradoMedio = [grado];
            CoefClustering = [clustering];
            Modularidad = [modularidad];
            T = table(GradoMedio,CoefClustering,Modularidad);
            dname = uigetdir('C:\');
            fileName = strcat(dname,'\globalFeatures.csv');
            
            if(fileName(1) ~= '\globalFeatures.csv')
                writetable(T,fileName,'Delimiter',',');
            end
            
        end
        
        %Función para guardar las medidas locales
        function exportLocalFeaturesToCSV(obj)
            for i=1:obj.filas
                name = obj.brainSections(i).nombre;
                Nombres(i,1) = string(name);
                Grado(i,1) = obj.calculaGradoNodo(name);
                Cercania(i,1) = obj.calculaCercaniaNodo(name);
                Intermediacion(i,1) = obj.calculaIntermediacionNodo(name);
                Ctr_VectorProprio(i,1) = obj.calculaVectorPropio(name);
            end
            
            T = table(Nombres,Grado,Cercania,Intermediacion,Ctr_VectorProprio,'RowNames',Nombres);
            dname = uigetdir('C:\');
            fileName = strcat(dname,'\localFeatures.csv');
            writetable(T,fileName,'Delimiter',',');
        end
        
        %Función para devolver el mensaje de informacion adicional
        function mensaje = informacionAdicional(obj)
            mensaje1 = "Esta aplicación ha sido el proyecto de trabajo de fin de grado de Daniel Cruz Andrades";
            mensaje2 = ", alumno de la Escuela Técnica Superior de Ingeniería Informática y Telecomunicaciones de la UGR";
            mensaje = strcat(mensaje1,mensaje2);
        end
        
        %Función para devolver el mensaje de contacto
        function mensaje = contactaConmigo(obj)
            mensaje = "e-mail de contacto : danielcruz16@correo.ugr.es";
        end
        
        %Función para restablecer los datos
        function image = restableceMatriz(obj,id)
            obj.matriz = obj.matrizOriginal;
            image = obj.pintaMatriz(id);
        end
        
        %Calculo de la dimensión fractal
        function dimFractal = calculaPendienteEscalaCajas(obj,lb_min,lb_max)
            dimFractal = obj.features.calculaPendienteEscalaCajas(lb_min,lb_max,obj.matriz);
        end
        
    end
end

