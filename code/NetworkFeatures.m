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
        function  calculaCompactBoxBurning(obj,matriz,nodos,combo)
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
            if combo == 0
                obj.draw_box_plot(iteraciones,valores);
            else
                obj.draw_combo_plot(iteraciones,valores);
            end
        end
        
        %Cálculo de dimension fractal con el algoritmo Greedy Coloring
        function devuelve = calculaGreedyColoring(obj, matriz, nodos,combo, tb, aux_g, neigh)
            limit = size(matriz);
            threshold = 0.5;
            valores = [];
            iteraciones = [1:1:limit(1)];
            matriz(matriz>=threshold)=1;matriz(matriz<threshold)=0; 
            nargin
            if nargin == 4
                onetime=false;
            else
                if ~exist('neigh','var')
                    onetime = true;
                    [G_b,neighbors] = obj.auxiliary_graph(matriz,tb);
                else
                    onetime=true;
                    G_b = aux_g;
                    neighbors = neigh;
                end
            end
            for k=1:limit(1)
                if ~onetime
                    [G_b,neighbors] = obj.auxiliary_graph(matriz,k);
                else
                    k = tb;
                end
                boxes = {};
                cont = 1;
                colors = zeros(1,limit(1));
                for i=1:limit(1)
                    curr_neigh = neighbors{i};
                    for j=1:limit(1)
                        if ~ismember(j,colors(curr_neigh))
                            colors(i) = j;
                            break
                        end
                    end
                end
                max_col = max(colors);
                valores(k) = max_col;
                if onetime
                    devuelve = colors;
                    break
                end
            end
            if ~onetime
                if combo == 0
                    obj.draw_box_plot(iteraciones,valores);
                elseif  combo == 1
                    obj.draw_combo_plot(iteraciones,valores);
                end
            end
        end
        
        %%%%%%%%%
        %Función para calcular el box-covering con el método Merge
        %Algorithm
        function calculaMergeAlgorithm(obj,matriz,nodos,combo)
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
           if combo == 0
               obj.draw_box_plot(iteraciones,valores);
           else
               obj.draw_combo_plot(iteraciones,valores);
           end
        end
        %%%%%%%%%
        
        %%%%%%%%%%
        %Funcion para calcular el box-covering con el método overlapping
        %box-covering algorithm
        
        function calculaOBCA(obj,matriz,nodos,combo)
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
            if combo == 0
                obj.draw_box_plot(iteraciones,valores);
            else
                obj.draw_combo_plot(iteraciones,valores);
            end
        end
        %%%%%%%%%%
        
        %%%%%%%%%%
        %Función para calcular el box-covering con el método maximal
        %excluded mass burning.
        function calculaMEMB(obj,matriz,nodos,combo)
            threshold = 0.7;
            matriz(matriz>=threshold)=1;matriz(matriz<threshold)=0;
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
                    [p, p_box] = obj.maxexcmass(U,NCN,distancias,k); %maxexcmass should return node with maximal value for em.
                    C(end+1) = p;
                    NCN = NCN(NCN~=p);
                    %calculate p_box should return nodes 'covered' by p
                    U = setdiff(U,p_box);
                end
                boxes = obj.form_boxes_from_centers(C,nodos,distancias,k);
                
                num_boxes = size(boxes);
                valores(k) = num_boxes(2);
            end
            if combo == 0
                obj.draw_box_plot(iteraciones,valores);
            else
                obj.draw_combo_plot(iteraciones,valores);
            end
        end
        %%%%%%%%%%
        
        %%%%%%%%%
        function boxes = calculaREMCC(obj,matriz,nodos,combo)
            %this function defines k as the distance, not size of the box, much
            %like other algorithms we'll see in this thesis
            tammat = size(matriz);
            threshold = 0.7;
            matriz(matriz>=threshold)=1;matriz(matriz<threshold)=0;
            G = graph(matriz);
            distancias = distances(G);
            distancias(distancias==0)=inf;
            valores = [];
            iteraciones = 1:1:tammat(1);
            sh_paths = [];
            for i=1:tammat(1)
                sh_paths(end+1)= sum(distancias(i,distancias(i,:)~=inf))/(tammat(1)-1);
            end
            for k=1:tammat(1)
                U = nodos; %list of uncovered nodes
                C = []; %list of center nodes
                boxes = {};
                utam = size(U);
                %distances already gives us average shortest distance path between
                %nodes, so now we calculate excluded mass for all nodes,
                while ~isempty(U)
                    masses = obj.excmass(U,nodos,distancias,k); %returns list of masses of U on the first row
                    %and number of node on the second row.
                    if ~any(masses)
                        maxnode = U(randi(length(U),1));
                    else
                        tamm = size(masses);
                        maxf = 0;
                        maxnode = 0;
                        for n=1:tamm(2)
                            f = masses(1,n)*sh_paths(1,n);
                            if f > maxf
                                maxf = f;
                                maxnode = n;
                            end
                        end
                    end
                    
                    C(end+1) = maxnode;
                    maxnodeball = obj.seed_ball(maxnode,U,distancias,k);
                    U = setdiff(U,maxnodeball);
                end
                %box forming of memb with all reached centers
                index = 1:tammat(1)+1:tammat(1)*tammat(2);
                distancias(index)=0;
                boxes = obj.form_boxes_from_centers(C,nodos,distancias,k);
                
                num_boxes = size(boxes);
                valores(k) = num_boxes(2);
            end
            if combo == 0
                obj.draw_box_plot(iteraciones,valores);
            else
                obj.draw_combo_plot(iteraciones,valores);
            end
        end
        %%%%%%%%%
        
        %%%%%%%%%
        %Funcion para calcular el box-covering con el método random sequential
        function calculaRandomSequential(obj,matriz,nodos,combo)
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
                end
                
                num_boxes = size(boxes);
                valores(k) = num_boxes(2);
                fprintf("numcajas en iteracion %d es %d \n",k,num_boxes(2));
            end
            if combo == 0
                obj.draw_box_plot(iteraciones,valores);
            elseif combo == 1
                obj.draw_combo_plot(iteraciones,valores);
            end
        end
        %%%%%%%%%
        
        %%%%%%%%%
        %Función para calcular el box covering con Particle Swarm
        %Optimization
        
        function boxes = calculaPSO(obj,matriz,nodos,combo,g,p,c1,c2)
            %operator (+) is defined as the interger difference of two numbers
            G = graph(matriz);
            distancias = distances(G);
            tammat = size(matriz);
            valores = [];
            iteraciones = 1:1:tammat(1);
            %for k=1:tammat(1)
            for k=1:tammat(1)
                boxes = {};
                P = []; %P stores greedily selected colors for each 'particle', meaning
                %apply greedy algorithm to the initial matriz to have an
                %initial positioning of nodes,for n size of nodes
                gbest = [];
                gbestsize = 1000;
                for i=1:p
                    P(end+1,:) = obj.calculaGreedyColoring(matriz,nodos,0,k); %this returns list of colors of each node
                    nboxes = max(P(end));
                    if ~isempty(gbest) || (nboxes < gbestsize)
                        gbest = P(end,:);
                        gbestsize = nboxes;
                    end
                end
                V = zeros(size(P)); %V stores the velocity of each node (chance to change assigned box
                pbest = P;
                for j=1:g
                    for n=1:p
                        %construct boxes
                        boxes = {};
                        maxcolor = max(P(n,:));
                        for l=1:maxcolor
                            boxes{l} = find(P(n,:)==l);
                        end
                        
                        %leverage bla bla bla
                        for m=1:tammat(1)
                            omega = rand;
                            r1 = rand;
                            r2 = rand;
                            %update speeed
                            V(n,m) = obj.sig(omega*V(n,m)+c1*r1*~isempty(obj.equalbox(pbest(n,m),P(n,m)))+c2*r2*~isempty(obj.equalbox(gbest(m),P(n,m))));
                            
                            if V(n,m) == 1
                                orig = P(n,m);
                                P(n,m) = obj.nbest(m,P(n,m),boxes,distancias,k);
                                
                                if orig ~= P(n,m)
                                    manbox = boxes{orig};
                                    manbox(manbox==m)=[];
                                    boxes{orig} = manbox;
                                    newbox = boxes{P(n,m)};
                                    newbox(end+1) = m;
                                    boxes{P(n,m)} = newbox;
                                end
                            end
                        end
                    end
                end
                num_boxes = size(boxes);
                valores(k) = num_boxes(2);
            end
            if combo == 0
                obj.draw_box_plot(iteraciones,valores);
            else
                obj.draw_combo_plot(iteraciones,valores);
            end
        end


        %%%%%%%%%
        
        %%%%%%%%%
        function draw_combo_plot(obj,iteraciones,valores)
            obj.boxCovering = valores;
            
            x = iteraciones;
            y = obj.boxCovering(iteraciones);
            
            Lb = x;
            Nb = y;
            
            N = log10(Nb)';
            R = log10(Lb)';
            nexttile
            set(gcf,'color',[1 1 1]);
            set(gca, 'FontSize',15);
            p = scatter(R,N,'filled');            
            
            dtt = p.DataTipTemplate;
            dtt.DataTipRows(1).Label = "Lb";
            dtt.DataTipRows(1).Value = Lb;
            dtt.DataTipRows(2).Label = "Nb";
            dtt.DataTipRows(2).Value = Nb;
        end
        %%%%%%%%%
        
        %%%%%%%%%
        function draw_box_plot(obj,iteraciones,valores)
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
        function ball = seed_ball(obj,node,nodes,distancias,k)
            tamn=size(nodes);
            ball = [];
            ball(end+1) = node;
            for i=1:tamn(2)
                if node ~= nodes(1,i) && distancias(node,nodes(1,i)) <= k
                    ball(end+1) = nodes(1,i);
                end
            end
        end
        %%%%%%%%%
        
        %%%%%%%%%
        function masses = excmass(obj,U,nodos,distancias,k)
            tamu = size(U);
            tamn = size(nodos);
            masses = zeros(tamn);
            indexes = U;
            for i=1:tamu(2)
                for j=1:tamu(2)
                    if distancias(U(i),U(j)) <= k && j ~= i
                        masses(U(i)) = masses(U(i)) + 1;
                    end
                end
            end
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
        %Funcion signo para PSO
        function retval = sig(obj,z)
            compare = 1/(1+exp(-z));
            if rand < compare
                retval = 1;
            else
                retval = 0;
            end
        end
        %%%%%%%%%
        
        %%%%%%%%%
        %Funcion para comparar dos cajas (desordenadas)
        function dev = equalbox(obj,b1,b2)
            if isempty(setdiff(b1,b2)) && isempty(setdiff(b2,b1))
                dev = [];
            else
                dev = [setdiff(b1,b2),setdiff(b2,b1)];
            end
        end
        %%%%%%%%%
        
        %%%%%%%%%
        function retval = nbest(obj,n,x,B,distances,lb)
            retval = x;
            while ~isempty(B)
                ind = randi(size(B));
                xk = B{ind};
                B(ind) = [];
                if x == ind
                    continue
                end
                f = arrayfun(@(z) distances(n,z)<lb, xk);
                if ~ismember(0,f)
                    retval = ind;
                    break
                end
            end
        end
        %%%%%%%%%
        
        %%%%%%%%%
        function boxes = form_boxes_from_centers(obj,C,candidates,distancias,k)
            boxes = {};
            tn = size(candidates);
            tc = size(C);
            for i=1:tc(2)
                boxes{end+1} = C(i);
            end
            mindists = [candidates;zeros(1,tn(2));zeros(1,tn(2))+k];
            for i=1:tn(2)
                for j=1:tc(2)
                    if distancias(candidates(i),C(j)) <= mindists(3,i)
                        mindists(3,i) = distancias(candidates(i),C(j));
                        mindists(2,i) = C(j);
                    end
                end
            end
            mindists = sortrows(mindists',3)';
            
            C_bis = mindists(:,mindists(3,:)~=0);
            cbtam = size(C_bis);
            ctam = size(C);
            for i=1:cbtam(2)
                smallers = mindists(:,mindists(3,:)<C_bis(3,i));
                stam = size(smallers);
                smaller = smallers(:,randi(stam(2)));
                for j=1:ctam(2)
                    if ismember(smaller(2,:),boxes{j})
                        boxes{j}(end+1) = C_bis(1,i);
                    end
                end
            end
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
        function [agraph,neighbors] = auxiliary_graph(obj,matriz,lb)
            neighbors = {};
            limit = size(matriz);
            edges = [0 0];
            nlist = [0 0];
            G = graph(matriz);
            distancias = distances(G);
            
            for v1=1:limit(2)
                for v2=1:limit(2)
                    if(distancias(v1,v2)>lb && v1~=v2)
                        if(~any(ismember(edges,[v2 v1],'rows')))
                            edges(end+1,:) = [v1 v2];
                        end
                        nlist(end+1,:) = [v1 v2];
                    end
                end
            end
            edges(1,:) = [];
            nlist(1,:) = [];
            G2 = graph;
            sedge = size(edges);
            for i=1:sedge(1)
                G2 = addedge(G2,edges(i,1),edges(i,2));
            end
            iedges = nlist;
            for i=1:limit(2)
                neighbors{end+1} = iedges(iedges(:,1)==i,2);
            end
            agraph = G2;
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

