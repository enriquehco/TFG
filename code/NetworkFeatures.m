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
        
        %%%%%%%%%
        %Función para calcular el box-covering con el método Merge
        %Algorithm
        function calculaMergeAlgorithm(obj,matriz,nodos)
           limit = size(matriz);
           valores = [];
           iteraciones = 1:1:limit(1);
           n_nodos = size(nodos);
           G = graph(matriz);
           distancias = distances(G)
           for k=1:limit(1)
                boxes = {};
                for i=1:n_nodos(2)
                    boxes{end+1} = nodos(i);
                end
                %todo nodo es una caja
                l = 1;
                while l <= k
                    boxes_vis = {};
                    n_boxes = size(boxes);
                    while ~isempty(boxes)
                        bnum = randi(n_boxes(2));
                        b = boxes(bnum);
                        choice = obj.find_adjacents(b,boxes,l,distancias);
                        if ~isempty(choice)
                            %boxes_vis{end+1} = [choice{1} choice{2}];
                            boxes_vis{end+1} = [choice{:,1} choice{:,2}];
                            f1 = @(a) isempty(setdiff(choice{1},a)) & isempty(setdiff(a,choice{1}));
                            f2 = @(a) isempty(setdiff(choice{2},a)) & isempty(setdiff(a,choice{2}));
                            matches1 = cellfun(f1,boxes,'UniformOutput',false);
                            matches2 = cellfun(f2,boxes,'UniformOutput',false);
                            u = n_boxes(2);
                            while u ~= 0
                                if isequal(matches1{u},1)
                                    boxes(u) = [];
                                    n_boxes(2) = n_boxes(2) - 1;
                                elseif isequal(matches2{u},1)
                                    boxes(u) = [];
                                    n_boxes(2) = n_boxes(2) - 1;
                                end
                                u = u - 1;
                            end
                        else
                            boxes_vis{end+1} = b{:};
                            boxes(bnum) = [];
                            n_boxes(2) = n_boxes(2) - 1;
                        end
                    end
                    boxes = boxes_vis;
                    l = l + 1;
                end
                boxes
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
        
        %%%%%%%%%%
        %Funcion para calcular el box-covering con el método overlapping
        %box-covering algorithm
        
        function calculaOBCA(obj,matriz,nodos)
            limit = size(matriz);
            valores = [];
            iteraciones = 1:1:limit(1);
            n_nodos = size(nodos);
            G = graph(matriz);
            distancias = distances(G);
            for k=1:limit(1)
                boxes = {};
                n_nodes = size(nodos);
                f = zeros(1,n_nodes(2));
                nd = obj.get_degree(matriz);
                aux_m = [nodos;f;nd];
                [temp,idx] = sort(aux_m(3,:));
                full_m = aux_m(:,idx);
                
                for n=1:n_nodes(2)
                    if full_m(2,n) > 0
                        continue;
                    end
                    %node has relative frequency = 0, so it's a center node for a
                    %candidate box, build said box around the node with distance k
                    box = obj.build_c_box(n,k,full_m,distancias);
                    [temp2,idx2] = sort(box(2,:));
                    C_box = box(:,idx2); %sort by f (2nd row)
                    tam_cb = size(C_box);
                    for i=1:tam_cb(2)
                        if ~ismember(C_box(1,i),box(1))
                            continue;
                        end
                        for j=i+1:tam_cb(2)
                            if distancias(C_box(1,i),C_box(1,j)) > k
                                if C_box(1,j) == full_m(1,n)
                                    box = box(:,box(1,:)~=C_box(1,i));
                                    break
                                else
                                    box = box(:,box(1,:)~=C_box(1,j));
                                end
                            end
                        end
                    end
                    for p=1:tam_cb(2)
                        fidx = find(full_m(1,:) == box(1,p));
                        full_m(2,fidx) = full_m(2,fidx) + 1;
                    end
                    boxes{end+1} = double(box(1,:));
                end
                numbox = size(boxes);
                for l=1:numbox(2)
                    uniquebox = false;
                    curr_b = boxes{l};
                    scb = size(curr_b);
                    for j=1:scb(2)
                        bidx = find(full_m(1,:) == curr_b(j));
                        if full_m(2,bidx) < 2
                            uniquebox = true;
                            break;
                        end
                    end
                    if ~uniquebox
                        boxes(l) = [];
                        for j=1:scb(2)
                            bidx2 = find(full_m(1,:) == curr_b(j));
                            full_m(2,bidx2) = full_m(2,bidx2) - 1;
                        end
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
            p = scatter(R,N,'filled');
            title("Escala de cajas");
            hold on;
    
            dtt = p.DataTipTemplate;
            dtt.DataTipRows(1).Label = "Lb";
            dtt.DataTipRows(1).Value = Lb;
            dtt.DataTipRows(2).Label = "Nb";
            dtt.DataTipRows(2).Value = Nb;
        end
        %%%%%%%%%%
        
        %%%%%%%%%%
        %Función para calcular el box-covering con el método maximal
        %excluded mass burning.
        function calculaMEMB(obj,matriz,nodos)
            matriz(matriz~=0) = 1;
            G = graph(matriz);
            distancias = distances(G);
            tammat = size(matriz);
            valores = [];
            iteraciones = 1:1:tammat(1);
            
            for k=1:tammat(2)
                U = nodos; %list of uncovered nodes
                C = []; %list of center nodes
                boxes = {};
                NCN = nodos;
                while ~isempty(U)
                    [p, p_box] = obj.maxexcmass(U,NCN,distancias,k); 
                    C(end+1) = p;
                    boxes{end+1} = p;
                    NCN = NCN(NCN~=p);
                    U = setdiff(U,p_box);
                end
                cdist = obj.calc_dists(nodos,C,distancias,k); 
                C_bis = cdist(:,cdist(3,:)~=0);
                cbtam = size(C_bis);
                ctam = size(C);
                for i=1:cbtam(2)
                    smallers = cdist(:,cdist(3,:)<C_bis(3,i));
                    stam = size(smallers);
                    smaller = smallers(:,randi(stam(2)));
                    for j=1:ctam(2)
                        if ismember(smaller(2,:),boxes{j})
                            boxes{j}(end+1) = C_bis(1,i);
                        end
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
            p = scatter(R,N,'filled');
            title("Escala de cajas");
            hold on;
    
            dtt = p.DataTipTemplate;
            dtt.DataTipRows(1).Label = "Lb";
            dtt.DataTipRows(1).Value = Lb;
            dtt.DataTipRows(2).Label = "Nb";
            dtt.DataTipRows(2).Value = Nb;
        end
        %%%%%%%%%%
        
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
        
        %%%%%%%%%
        function [node,mbox] = maxexcmass(obj,U,NCN,distancias,k)
            tamncn = size(NCN);
            tamu = size(U);
            max_mass = 0;
            node = 0;
            for i=1:tamncn(2)
                ex_mass = 0;
                curr_box = [];
                curr_box(end+1) = NCN(i);
                for j=1:tamu(2)
                    if distancias(NCN(i),U(j)) <= k && U(j) ~= node
                        ex_mass = ex_mass + 1;
                        curr_box(end+1) = U(j);
                    end
                end
                if ex_mass > max_mass
                    node = NCN(i);
                    max_mass = ex_mass;
                    mbox = curr_box;
                end
            end
        end
        %%%%%%%%%
        
        %%%%%%%%%
        function mindists = calc_dists(obj,nodos,C,distancias,k)
            tn = size(nodos);
            tc = size(C);
            mindists = [nodos;zeros(1,tn(2));zeros(1,tn(2))+k];
            for i=1:tn(2)
                for j=1:tc(2)
                    if distancias(nodos(i),C(j)) <= mindists(3,i)
                        mindists(3,i) = distancias(nodos(i),C(j));
                        mindists(2,i) = C(j);
                    end
                end
            end
            mindists = sortrows(mindists',3)';
        end
        %%%%%%%%%
        
        %%%%%%%%%
        %Devuelve el grado de los nodos
        function degrees = get_degree(obj,matriz)
            tamm = size(matriz);
            degrees = zeros(1,tamm(1));
            
            for i=1:tamm(1)
                degrees(1,i) = nnz(matriz(i,:));
            end
        end
        %%%%%%%%%
        
        %%%%%%%%%
        %Construye una caja según distancia a un nodo central
        function box = build_c_box(obj,n,k,full_m,distancias)
            box = int16.empty(3,0);
            n_el = full_m(1,n);
            box(:,end+1) = full_m(:,n);
            tami = size(full_m);
            for i=1:tami(2)
                ni = full_m(1,i);
                if n_el~=ni && distancias(n_el,ni) <= k
                    box(:,end+1) = full_m(:,i);
                end
            end
        end
        
        %%%%%%%%%
        %Cálculo de las posibles cajas adyacentes para algoritmo merge
        function choice = find_adjacents(obj,b,boxes,l,distances)
            num_boxes = size(boxes);
            tamb = size(b{:});
            valid_boxes = {};
            bwascell = false;
            if isa(b,'cell') || tamb(2) > 1
                bwascell = true;
                b = b{:};
            end
            for i=1:num_boxes(2)
                tamc = size(boxes{i});
                mergeable = true;
                cwascell = false;
                compbox = boxes(i);
                
                if isa(compbox,'cell') || tamc(2) > 1
                    compbox = compbox{:};
                    cwascell = true;
                end
                
                for j=1:tamb(2)
                    for k=1:tamc(2)
                        if bwascell
                            bvalue = b(j);
                        elseif ~bwascell
                            bvalue = b{j};
                        end
                        
                        if cwascell
                            cvalue = compbox(k);
                        elseif ~cwascell
                            cvalue = compbox{k};
                        end
                        if cvalue == bvalue || distances(cvalue,bvalue) > l
                            mergeable = false;
                        end
                    end
                end
                
                if ~cwascell && isa(compbox,'cell')
                    compbox = compbox{:};
                end
                if ~bwascell && isa(b,'cell')
                    b = b{:};
                end
                if mergeable
                    valid_boxes{end+1} = {compbox b};
                end
            end
            tam_vb = size(valid_boxes);
            if tam_vb(2) > 0
                choicenum = randi(tam_vb(2));
                choice = valid_boxes{choicenum};
            else
                choice = {};
            end
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

