%%%%  CARLOS QUESADA GRANJA
%%%%  START: August 10, 2013
%%%%  END: September 5, 2013
%%%%  Laplace's Equation PGD

function [F1, F2] = LaplacesEqPGD06
    %% CONSTANTS
    % Number of nodes
    n = 11;
    % Mesh limits
    Mlim = {[-1 1] [-1 1]};
    % Boundary conditions
    BC = [1 n];
    % Error tolerance
    TOL = 1E-4;
    % Maximum iterations [F-loop R-loop]
    maxIter = [500 251];
    
    %% DEFINITION OF THE SOURCE TERM
    % Example 1
    source = {@(x)cos(2*pi*x); @(y)sin(2*pi*y)};
    % Example 2
    source = {@(x)power(x,2) @(x)ones(1,n); ...
        @(y)ones(1,n) @(y)-power(y,2)};
    % Example 3
%     source = {@(x)2*x.^2 @(x)x @(x)ones(1,n) @(x)ones(1,n) @(x)3*x;
%         @(y)ones(1,n) @(y)ones(1,n) @(y)y.^2 @(y)-0.2*y @(y)y};
    % Example 4
%     source = { ...
%         @(x)(-8*pi^4)*(cos(2*pi^2*x)) @(x)(-8*pi^4)*(sin(pi*x).^2) ...
%         @(x)(8*pi^4)*(cos(2*pi*x));
%         @(y)sin(pi*y).^2 @(y)cos(2*pi*y) @(y)cos(2*pi*y)};

    
    %% INITIALIZATIONS
    % Cell memory preallocation   
    K = cell(2,1);
    M = cell(2,1);
    F = cell(2,1);
    coor = cell(2,1);
    % Axis selection loop
    for xy = 1:2
        % Get coordinates
        coor{xy} = linspace(Mlim{xy}(1), Mlim{xy}(2), n);
        % Definition of stiffness matrices
        K{xy} = zeros(n);
        % Definition of mass matrices
        M{xy} = zeros(n);
        % Definition of F matrix
        F{xy} = [];
    end
    
    %% BOUNDARY CONDITIONS
    % List of nodes
    nonBC = 1:n;
    % Non-boundary conditions
    nonBC(ismember(nonBC, BC)) = [];
    
    %% STIFFNESS AND MASS MATRICES
    % Gauss points and weights
    gPoi = [-3^(-1/2) 3^(-1/2)];    
    gWei = [1 1];
    % Number of Gauss points
    ngp = length(gPoi);
    % Axis selection loop
    for xy = 1:2
        % Integration loop
        for ii = 1:n-1
            % Physical domain of integration of the current element
            gA = coor{xy}(ii);
            gB = coor{xy}(ii+1);
            gL = gB-gA;
            % Translation of gPoi from E-coordinates to x-coordinates
            E2x = gA.*(1-gPoi)./2 + gB.*(gPoi+1)./2;
            % Values of N(x) at x = E2x for ii and ii+1 elements
            % See N in eq. 5.20 of Belytschko & Fish
            Nx = zeros(n,ngp); 
            Nx(ii,:) = (gB-E2x)./gL;
            Nx(ii+1,:) = (E2x-gA)./gL;
            % Values of dN(x) at x = E2x for ii and ii+1 elements
            % See B in eq. 5.20 of Belytschko & Fish
            dNx = zeros(n,ngp);
            dNx(ii,:) = (-1/gL)*ones(1,ngp);
            dNx(ii+1,:) = (1/gL)*ones(1,ngp);
            % Assembly and Gauss weight loop
            for jj = 1:ngp
                % Mass matrix
                M{xy} = M{xy} + gWei(:,jj)*(gL/2)*(Nx(:,jj)*Nx(:,jj)'); 
                % Stiffness matrix
                K{xy} = K{xy} + gWei(:,jj)*(gL/2)*(dNx(:,jj)*dNx(:,jj)');
            end
        end
    end
    % Clear useless variables
    clear dNx E2x gA gB gL gPoi gWei ngp Nx
    
    %% SOURCE TERM
    % Cell memory preallocation
    A = cell(2,1);
    V = cell(2,1);
    Vsize = size(source,2);
    % Axis selection loop
    for xy = 1:2
        % Columns of source loop
        for ii=1:Vsize;
            % Evaluate the source term in all the nodes
            A{xy}(:,ii) = source{xy,ii}(coor{xy});
        end
        % Mass matrix times source matrix
        V{xy} = M{xy}*A{xy};
    end
    
    %% ENRICHMENT AND PROJECTION STAGES
    % Initialization of F error
    Ferr = 1;
    aprt = 0;
    Fiter = 0;
    
    % F loop
    while (Ferr > TOL) && (Fiter < maxIter(1))
        % ENRICHMENT STAGE 
        % Initialization of R
        Riter = 0;
        R{1} = zeros(n,1);
        R{2} = ones(n,1);
        R{2}(BC) = 0;
        % Initialization of R error
        Rerr = 1;
        
        % R loop
        while (Rerr > TOL) && (Riter < maxIter(2))
            % New iteration
            Riter = Riter + 1;
            % Copy of R
            Raux = R;
            
            % Get Rx known Ry 
            % Left term
            D1 = (R{2}'*M{2}*R{2})*K{1} + (R{2}'*K{2}*R{2})*M{1};
            % Source term (V)
            D2 = 0.0;
            for ii = 1:Vsize
                D2 = D2 + V{1}(:,ii)*(R{2}'*V{2}(:,ii));
            end
            % Right term (F)
            for ii = 1:Fiter
                D2 = D2 - (K{1}*F{1}(:,ii))*(R{2}'*M{2}*F{2}(:,ii));
                D2 = D2 - (M{1}*F{1}(:,ii))*(R{2}'*K{2}*F{2}(:,ii));
            end
            % Final assembly with imposition of domain boundaries
            R{1}(nonBC) = D1(nonBC,nonBC)\D2(nonBC);
            % Normalization
            R{1} = R{1}/sqrt(R{1}'*M{1}*R{1});
%             R{1} = R{1}./norm(R{1});

            % Get Ry known Rx 
            % Left term
            D1 = (R{1}'*K{1}*R{1})*M{2} + (R{1}'*M{1}*R{1})*K{2};
            % Source term (V)
            D2 = 0.0;
            for ii = 1:Vsize
                D2 = D2 + V{2}(:,ii)*(R{1}'*V{1}(:,ii));
            end
            % Right term (F)
            for ii = 1:Fiter
                D2 = D2 - (M{2}*F{2}(:,ii))*(R{1}'*K{1}*F{1}(:,ii));
                D2 = D2 - (K{2}*F{2}(:,ii))*(R{1}'*M{1}*F{1}(:,ii));
            end
            % Final assembly with imposition of domain boundaries
            R{2}(nonBC) = D1(nonBC,nonBC)\D2(nonBC);
            % Normalization
            R{2} = R{2}/sqrt(R{2}'*M{2}*R{2});
%             R{2} = R{2}./norm(R{2});
            
            % R error calculation
            Rerr = sqrt(norm(Raux{1}-R{1})+norm(Raux{2}-R{2}));
        end
        
        % New finished iteration
        Fiter = Fiter + 1;
        % Appending R to F
        for xy = 1:2
            F{xy}(:,end+1) = R{xy};
        end
        
        % PROJECTION STAGE
        D1 = zeros(Fiter);
        D2 = zeros(Fiter,1);
        for ii = 1:Fiter
            % Left term (FKF, FMF)
            for jj = 1:Fiter
                D1(ii,jj) = ...
                    (F{1}(:,ii)'*K{1}*F{1}(:,jj)) * ...
                    (F{2}(:,ii)'*M{2}*F{2}(:,jj)) + ...
                    (F{1}(:,ii)'*M{1}*F{1}(:,jj)) * ...
                    (F{2}(:,ii)'*K{2}*F{2}(:,jj));
            end
            % Right term (FV)
            for jj = 1:Vsize
                D2(ii) = D2(ii) + ...
                    (F{1}(:,ii)'*V{1}(:,jj)) * (F{2}(:,ii)'*V{2}(:,jj));
            end
        end
        % Final assembly
        alpha = D1\D2;
        for ii = 1:Fiter
            sqrt_alpha = sqrt(abs(alpha(ii)));
            F{1}(:,ii) = sqrt_alpha.*F{1}(:,ii);
            F{2}(:,ii) = (alpha(ii)/sqrt_alpha).*F{2}(:,ii);
        end
              
        % F error calculation
        Ferr = norm(F{1}(:,Fiter))*norm(F{2}(:,Fiter));
        aprt = max(aprt,sqrt(Ferr));
        Ferr = sqrt(Ferr)/aprt;
    end
    
    %% PLOTTING
    out = F{2}*F{1}';
    figure('Name', mfilename)
    surf(coor{1}, coor{2}, out, 'FaceColor', 'interp');
    F1 = F{1};
    F2 = F{2};
    
end