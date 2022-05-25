classdef NetworkFeatures < Features
    %Esta clase representa las distintas medidas que se pueden tomar de la
    %red.
    
    properties
        gradoMedio;
        coefClustering;
        modularidad;
        boxCovering; 
        id;
    end
    
    methods
        
        %Constructor, aquí usamos la propiedad de la clase padre abstracta
        %id para tener identificado el nombre de la red.
        function obj = NetworkFeatures(identificador)
            obj.id = identificador;
        end
        
        % GRADO
        function grado = calculaGradoMedio(obj,datos)
            limit = 0;
            numEnlaces = 0;
            numNodos = size(datos);
            
            for i=1:numNodos(1)
                for j=1:limit
                    if(datos(i,j) > 0)
                        numEnlaces = numEnlaces + 1;
                    end
                end
                limit = limit + 1;
            end
            
            obj.gradoMedio = (2*numEnlaces)/numNodos(1);
            grado = obj.gradoMedio;
        end
        
        %COEFICIENTE DE CLUSTERING
        function coefClusteringMedio = calculaCoefClustering(obj,datos)
            numNodos = size(datos);
            enlaceVecino = 0;
            
            for i=1:numNodos(1)
                
                %primero vemos el numero de enlaces que tiene el nodo
                cont = 0;
                for j=1:numNodos(1)
                   if(datos(i,j) > 0)
                       cont = cont + 1;
                       nodosVecinos(i,cont) = j;
                   end
                end
                
                %contabilizamos el numero de enlaces de sus vecinos.
                for j=1:cont-1
                    for k = j:cont
                        if(datos(nodosVecinos(i,j),nodosVecinos(i,k)) > 0)
                            enlaceVecino = enlaceVecino + 1;
                        end
                    end
                end
                
                if(enlaceVecino ~= 0)
                    coeficienteNodos(i) = (2*enlaceVecino)/(cont*(cont-1));    
                else
                    coeficienteNodos(i) = 0;
                end
                
                enlaceVecino = 0;
                
            end
            
            %Calculamos el coeficiente de clustering medio
            C = 0;
            limit = size(coeficienteNodos);
            for i=1:limit(2)
                C = C + coeficienteNodos(i);
            end
            
            coefClusteringMedio = C/numNodos(1);
            obj.coefClustering = coefClusteringMedio;
        end
        
        %ÍNDICE DE MODULARIDAD
        function modularidad = calculaModularidad(obj,datos)
            modularidad = community_louvain(datos,1);
            obj.modularidad = modularidad;
        end
        
        %BOX-COVERING
        function  calculaCompactBoxBurning(obj,matriz,nodos)
            limit = size(matriz);
            valores = [];
            
            for i=1:limit(1)
                iteraciones(i) = i;
                for j=1:limit(1)
                    if(matriz(i,j) ~= 0)
                       matriz(i,j) = 1; 
                    end
                end
            end
            
            G = graph(matriz);
            distancias = distances(G);
            nodos 
            
            %%ALGORITMO COMPACT BOX BURNING%%
            
                for k=1:limit(1)
                    U = nodos(1,:);
                    tamBox = k;
                    boxes = {};
                    cont = 1;
                    
                    while(~isempty(U))
                        box = [];
                        C = U;
                        
                        while(~isempty(C))
                            p = randsample(C,1);
                            id_p = obj.nameToIdNodo(p,nodos);
                            box = [box,p];
                            limit_3 = size(U);
                            
                            for j=1:limit_3(2)
                                p_aux = U(j);
                                id_p_aux = obj.nameToIdNodo(p_aux,nodos);
                                if(obj.distanciaMayor(id_p,id_p_aux,tamBox,distancias) && p_aux ~= p)
                                    C = C(C~=p_aux);
                                end
                            end
                            
                            C = C(C~=p);
                            
                        end
                    
                        boxes{1,cont} = box;
                        cont = cont + 1;
                
                        aux = size(box);
                        limit_2 = aux(2);
                
                        for l=1:limit_2
                            U = U(U~=box(l));
                        end
                    end
                    
                    num_boxes = size(boxes);
                    valores(k) = num_boxes(2);
                    
                    
                    
                end
                
                obj.boxCovering = valores;
                
           
           x = iteraciones;
           y = obj.boxCovering(iteraciones);
           
           Lb = x;
           Nb = y;
           
           N = log10(Nb)';
           R = log10(Lb)';
           figure;
           set(gcf,'color',[1 1 1]);
           set(gca, 'FontSize',15);
           p = scatter(R, N,'filled');
           title("Escala de cajas");
           hold on;

           dtt = p.DataTipTemplate;
           dtt.DataTipRows(1).Label = "Lb";
           dtt.DataTipRows(1).Value = Lb;
           dtt.DataTipRows(2).Label = "Nb";
           dtt.DataTipRows(2).Value = Nb;
           
        end
        
        %Cálculo de dimension fractal con el algoritmo Greedy Coloring
        function calculaGreedyColoring(obj,matriz,nodos)
            limit = size(matriz);
            valores = [];
%             iteraciones = [];
            
            for i=1:limit(1)
                iteraciones(i) = i;
            end

            for k=1:limit(1)
                boxes = {};
                cont = 1;
                tamBox = k;

                G_b = obj.calculaGrafoAuxiliar(matriz,nodos,tamBox);

                nodesSorted = nodos(:,1);
                tam1 = size(nodesSorted);
                color(1:tam1(1)) = -1;
        
                edges = G_b.Edges;
                edges = edges{:,:};
        
                %datasample coge una muestra aleatoria de tamaño tam1 de
                %nodesSorted sin repetición
                nodesSorted = datasample(nodesSorted,tam1(1),'Replace',false);

                %maximo de colores usados para este tamaño de caja
                color_max=1;

                %para cada nodo del grafo
                for l=1:tam1(1)
                    %lista de colores adyacentes a este color
                    adj_col = [];


                    %lista de nodos adyacentes
                    adj_nod = edges(find(edges(:,1)==l),2);

                    %lista de nodos adyacentes (en la otra direccion del edge)
                    adj_nod2 = edges(find(edges(:,2)==l),1);
                    %este array indica los colores de los nodos adjacentes
                    adj_col = color(adj_nod);
                    adj_col2 = color(adj_nod2);

                    %ciclar los posibles colores hasta encontrar el menor posible a
                    %asignar
                    possible_col=1;
                    cond_met = false;
                    while (~cond_met && possible_col <= tam1(1)) 
                        if(~ismember(possible_col,adj_col) &&  ~ismember(possible_col,adj_col2))
                            color(l) = possible_col;
                            if(possible_col > color_max)
                                color_max = possible_col;
                            end
                            cond_met = true;
                        end
                        possible_col = possible_col+1;
                    end     
                end
%                 color_max = color_max+1;

                for i=1:color_max
                    boxes{cont,1} = nodos(color(:)==i);
                    cont = cont+1;
                end

                %este condicional es para mostrar que se crean las cajas
                %correctamente, se demuestra con k==2 cuando se realiza un coloreo
                %de 3 colores.
%                 if(k==2)
%                     plot(G_b);
%                     fprintf("tamaño de caja %d",k);
%                     for j=1:color_max
%                         boxes{j,1}
%                     end
%                 end

                num_boxes = size(boxes);
                valores(k) = num_boxes(1);

            end
             
            obj.boxCovering = valores;

            x = iteraciones
            y = obj.boxCovering(iteraciones)

            Lb = x;
            Nb = y;

            N = log10(Nb)';
            R = log10(Lb)';
            figure;
            set(gcf,'color',[1 1 1]);
            set(gca, 'FontSize',15);
            p = scatter(R, N,'filled');
            title("Escala de cajas");
            hold on;

            dtt = p.DataTipTemplate;
            dtt.DataTipRows(1).Label = "Lb";
            dtt.DataTipRows(1).Value = Lb;
            dtt.DataTipRows(2).Label = "Nb";
            dtt.DataTipRows(2).Value = Nb;
        end
        
        function calculaMEMB(obj,matriz,nodos)
            limit = size(matriz);
            valores = [];
            iteraciones = [];
            
            for i=1:limit(1)
                iteraciones(i)=i;
            end
            
            G = graph(matriz);
            distancias = distances(G);
            
            for k=1:limit(1)
               U = nodos(1,:);
               NCN = nodos(1,:);
               CN = {};
               max_exc_mass = 0;
               max_exc_node;
               %%%%
               tamBox = k;
               boxes = {};
               cont = 1;
               %%%%%
               while(~isempty(U) && ~isempty(NCN))
                  for l=1:limit(1)  
                    p = randsample(C,1);
                    id_p = obj.nameToIdNodo(p,nodos);
                    if (ex_mass(id_p) > max_exc_mass)
                        max_exc_mass = ex_mass(id_p);
                        max_exc_node = id_p;
                    end
                  end
                  CN = [CN max_exc_node];
                  NCN = NCN(NCN~= max_exc_node);
                  cov_list = calc_chem_distance(nodes,max_exc_node);
                  U = U(U~= cov_list);
               end
               %once all nodes have been processed
               central_dists = {};
               box_ids = {};
               for cn=1:size(CN)
                   box_ids(cn) = cn;
               end
               for n=1:size(nodos(1,:))
                   central_dists = [central_dists calc_central_dist(nameToIdNodo(n,nodos))];
               end
               %sort NCN according to central_dists
            end
        end
        %%%%%%%%%
        %Función para calcular el box-covering con el método Merge
        %Algorithm
        function calculaMergeAlgorithm(obj,matriz,nodos)
           limit = size(matriz);
           valores = [];
           iteraciones = 1:1:size(1);
           
           for k=1:limit(1)
                %todo nodo es una caja
                boxes = {};
                limit2 = size(nodos);
                for n=1:limit2(2)
                    boxes{end+1,1} = nodos(1,n);
                end
                l = 1;
                while l < k
                    boxes_vis = {};
                    while (~isempty(boxes))
                        bnum = randi(limit2(2));
                        b = boxes{bnum,1};
                        adj = {};
                        cell_adj = {};
                        for i=1:limit(1)
                            compound_box = {};
                            if i~=bnum
                                subbox1 = b;
                                subbox2 = boxes{i};
                                compound_box = [subbox1 subbox2];   
                                cell_box = {subbox1 subbox2};
                                adj{end+1} = compound_box; 
                                cell_adj{end+1} = cell_box;
                            end
                        end
                        if ~isempty(adj)
                            n_adj = size(adj);
                            choices = {};
                            cell_choices = {};
                            for c=1:n_adj(2)
                                b_size = size(adj{c});
                                if b_size(2) <= l
                                    choices{end+1} = adj{c};
                                    cell_choices{end+1} = cell_adj{c};
                                end
                            end
                            %función aislada, para comparación entre celdas
                            limit3 = size(choices);
                            if limit3(2) > 0
                                cnum = randi(limit3(2));
                                choice = choices{cnum,1};
                                cell_choice = cell_choices{cnum,1};
                                boxes_vis{end+1} = choice{1};
                                f1 = @(a) a == cell_choice{1};
                                f2 = @(a) a == cell_choice{2};
                                members1 = cellfun(f1,boxes,'UniformOutput',false);
                                members2 = cellfun(f2,boxes,'UniformOutput',false);
                                m_size = size(members1);
                                for q=1:m_size(1)
                                    if ~ismember(0,members1{q})
                                        boxes(q) = [];
                                        members1(q) = [];
                                        members2(q) = [];
                                        m_size(1) = m_size(1) - 1;
                                    end
                                    if ~ismember(0,members2{q})
                                        boxes(q) = [];
                                        members1(q) = [];
                                        members2(q) = [];
                                        m_size(1) = m_size(1) - 1;
                                    end
                                end
                            end
                        else
                            boxes_vis{end+1} = b;
                            cand_size = size(boxes);
                            for u=1:cand_size(1)
                                if isequal(b,boxes{u})
                                    boxes(u) = [];
                                    cand_size = cand_size -1;
                                end
                            end
                        end
                    end
                    boxes = boxes_vis;
                    l = l + 1;
                end
                num_boxes = size(boxes);
                valores(k) = num_boxes(2);
           end
           obj.boxCovering = valores;
    
            x = iteraciones;
            y = obj.boxCovering(iteraciones);
    
            Lb = x;
            Nb = y;
    
            N = log10(Nb)';
            R = log10(Lb)';
            figure;
            set(gcf,'color',[1 1 1]);
            set(gca, 'FontSize',15);
            p = scatter(R,N,'filled');
            title("Escala de cajas");
            hold on;
    
            dtt = p.DataTipTemplate;
            dtt.DataTipRows(1).Label = "Lb";
            dtt.DataTipRows(1).Value = Lb;
            dtt.DataTipRows(2).Label = "Nb";
            dtt.DataTipRows(2).Value = Nb;
        end
        %%%%%%%%%
        
        
        
        %%%%%%%%%
        %Funcion para calcular el box-covering con el método random sequential
        function calculaRandomSequential(obj,matriz,nodos)
        %RANDOMSEQUENTIAL Summary of this function goes here
        %   Detailed explanation goes here
            limit = size(matriz);
            valores = [];
            iteraciones = [];
            
            for i=1:limit(1)
                iteraciones(i)=i;
            end
    
            G = obj.normalizaGrafo(matriz);
            G = graph(G);
            distancias = distances(G);
    
            edges = G.Edges(:,1);
            edges = edges{:,:};
    
            for k=1:limit(1)
                boxes = {};
                cont = 1;
                tamBox = k;
       
                U = nodos(1,:);
       
                while (~isempty(U))
                    limit2 = size(U);
%                   c = U(randi(limit2,1),:);
                    c = randsample(U,1);
                    id_c = obj.nameToIdNodo(c,nodos);
                    box = [];
                    box = [box,c];
                    for i=1:limit2(2)
                        d = U(i);
                        id_d = obj.nameToIdNodo(d,nodos);
                        if(~obj.distanciaMayor(id_c,id_d,tamBox,distancias) && c~=d)
                            box = [box,d];
                        end
                    end
                    
                    aux = size(box);
                    limit3 = aux(2);
                    
                    for i=1:limit3
                        U = U(U~=box(i));
                    end
                    
                    boxes{1,cont} = box;
                    cont = cont+1;
                    fprintf("la caja %d contiene \n", k);
                    box
                end
                
                num_boxes = size(boxes);
                valores(k) = num_boxes(2);
                fprintf("numcajas en iteracion %d es %d \n",k,num_boxes(2));
            end
            
            obj.boxCovering = valores;
    
            x = iteraciones;
            y = obj.boxCovering(iteraciones);
    
            Lb = x;
            Nb = y;
    
            N = log10(Nb)';
            R = log10(Lb)';
            figure;
            set(gcf,'color',[1 1 1]);
            set(gca, 'FontSize',15);
            p = scatter(R,N,'filled');
            title("Escala de cajas");
            hold on;
    
            dtt = p.DataTipTemplate;
            dtt.DataTipRows(1).Label = "Lb";
            dtt.DataTipRows(1).Value = Lb;
            dtt.DataTipRows(2).Label = "Nb";
            dtt.DataTipRows(2).Value = Nb;
        end
        %%%%%%%%%
        
        %Cálculo de la matriz (grafo) auxiliar
        function newgraph = calculaGrafoAuxiliar(obj,matriz,nodos,lb)
            matriz = obj.normalizaGrafo(matriz);
            edges = [1 1];
            G = graph(matriz);
            distancias = distances(G);
            limit2 = size(nodos);
    
            for v1=1:limit2
                %id1 = nameToIdNodo(v1,nodos);
                for v2=1:limit2
                    %id2 = nameToIdNodo(v2,nodos);
                    if(obj.distanciaMayor(v1,v2,lb,distancias) && v1~=v2)
                        if(~ismember([v2 v1],edges,'rows'))
                            edges = [edges; [v1 v2]];
                        end
                    end
                end
            end

            edges = edges(2:end,:);
            G2 = graph;
            G2 = addnode(G2,limit2(1));
            s = edges(:,1);
            t = edges(:,2);

            G2 = addedge(G2,s,t);

            newgraph = G2;
        end
        %%%%%%%%%
        
        %Transforma la matriz de adyacencia a valores 0 o 1 (si el valor
        %original es 0, entonces el valor resultante es 0, cualquier otro
        %valor equivale a 1)
        function newmat = normalizaGrafo(obj,matriz)
            limit = size(matriz);
            for i=1:limit(1)
                for j=1:limit(1)
                    if(matriz(i,j) ~= 0)
                        newmat(i,j) = 1;
                    end
                end
            end
        end
        %%%%%%%%%
        
        %Funcion para comprobar si un nodo tiene mayor distancia con
        %respecto a otros que la indicada por argumento
        function result = distanciaMayor(obj,p,p_aux,tamBox,distancias)
            result = 0;
            if(isstring(p) && isstring(p_aux))
                p = str2num(p);
                p_aux = str2num(p_aux);
            end
            if(distancias(p,p_aux) >= tamBox)
                result = 1;
            end
        end
        
        %Funcion de correspondencia nodo to id
        function id = nameToIdNodo(obj,p,nodos)
            aux = size(nodos);
            limit = aux(2);
            for i=1:limit
                if(nodos(1,i) == p)
                    id = nodos(2,i);
                end
            end
        end
        
        %Calculo de la dimensión fractal
        function fd_titulo = calculaPendienteEscalaCajas(obj,lb_min,lb_max,matriz)
           limit = size(matriz);
            
           for i=1:limit(1)
               iteraciones(i) = i;
           end
           
           x = iteraciones;
           y = obj.boxCovering(iteraciones);
           
           Lb = x;
           Nb = y;
           
           N = log10(Nb)';
           R = log10(Lb)';
           
           % linear regression computation
           Rr = R(lb_min : lb_max); % R limited to the selected range of boxes
           Nr = N(lb_min : lb_max); % N limited to the selected range of boxes
           x = [ones(length(Rr), 1) Rr]; % adds a column of ones to Rr
           b = x \ Nr; % b(1) is the y-intercept and b(2) is the slope of the line 
           y = x * b; % linear regression
           corr = 1 - sum((Nr - y).^2) / sum((Nr - mean(Nr)).^2); % correlation

           figure;
           hold on;
           % plots the regression line and the fractal dimension results
           plot(Rr, y, 'r', 'LineWidth', 2);
           xlabel('log Lb');
           ylabel('log Nb');
           fd_titulo = num2str(abs(b(2)),'%.2f');% 3dfd proveniente de b(2) pero redondeada a 2 decimales
           title(sprintf('Correlación = %f', corr));
           p = scatter(R, N);
           
           dtt = p.DataTipTemplate;
           dtt.DataTipRows(1).Label = "Lb";
           dtt.DataTipRows(1).Value = Lb;
           dtt.DataTipRows(2).Label = "Nb";
           dtt.DataTipRows(2).Value = Nb;
           
           
        end
        
    end
end

