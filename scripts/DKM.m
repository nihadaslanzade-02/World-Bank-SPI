function [Vdkm,Udkm,Ymdkm, fdkm,indkm]=DKM(X, K, Q, varargin)

% n = number of objects
% J = number of variables
% K = number of ObjectClasses
% m = number of VariableClasses
% X = n x J (objects x variables) matrix
% U = n x K (objects x ObjectClasses) matrix
% V = J x Q (variables x variableClasses) matrix
% Rndstart = number of multistart



% model X = U*Ym*V' + E
%
% problem: min||X-U*Ym*V'||^2

%equivalent to

% max||U*Ym*V'||^2
% subject to
% V binary and row stochastic
% U binary and row stochastic

% initialization

[n,J]=size(X);
un=ones(n,1);
uk=ones(K,1);
um=ones(Q,1);

% 'Stats'    ->    Default value: 'on', print the statistics of the fit of the
%                  model.
%                  If 'off' Statistics are not printed (used for simulation
%                  studies)
% 'Stand'    ->    Default value 'on', standardize variables and therefore compute DFA 
%                  on the correlation matrix.
%                  If 'off' does not standardize variables and therefore
%                  compute DFA on the variance-covariance matrix
% 
% 'Rndst'    ->    an integer values indicating the intital random starts.
%                  Default '20' thus, repeat the anaysis 20 times and retain the
%                  best solution.
% 'MaxIter'  ->    an integer value indicationg the maximum number of
%                  iterations of the algorithm
%                  Default '100'.
% 'ConvToll' ->    an arbitrary samll values indicating the convergence
%                  tollerance of the algorithm, Default '1e-9'.



if nargin < 2
   error('Too few inputs');
end

if ~isempty(X)
    if ~isnumeric(X)
        error('Invalid data matrix');
    end  
    if min(size(X)) == 1
    error(message('Disjoint Factor Analysis:NotEnoughData'));
end

else
    error('Empty input data matrix');
end

if ~isempty(Q)
    if isnumeric(K)
        if Q > J 
              error('The number of latent factors larger that the number of variables');
           end    
    elseif Q < 1 
              error('Invalid number of latent factors');
    end
else
    error('Empty input number of latent factors');
end

% Optional parameters   
pnames = {'Stats' 'Stand' 'Rndst' 'MaxIter' 'ConvToll'};
dflts =  { 'on'    'on'     100       100       1e-9 };
[Stats,Stand,Rndst,MaxIter,ConvToll] = internal.stats.parseArgs(pnames, dflts, varargin{:});


%if ~isempty(eid)
% error(sprintf('Disjoint Factor Analysis; %s',eid), emsg);
%end


if ~isempty(Stats)
    if ischar(Stats)
       StatsNames = {'off', 'on'};
       js = strcmpi(Stats,StatsNames);
           if sum(js) == 0
              error(['Invalid value for the ''Statistics'' parameter: '...
                     'choices are ''on'' or ''off''.']);
           end
       Stats = StatsNames{js}; 
    else  
        error(['Invalid value for the ''Statistics'' parameter: '...
               'choices are ''on'' or ''off''.']);
    end
else 
    error(['Invalid value for the ''Statistics'' parameter: '...
           'choices are ''on'' or ''off''.']);
end


% Standardization
if ~isempty(Stand)
    if ischar(Stand)
       StandNames = {'off', 'on', 'Mahalanobis'};
       js = strcmpi(Stand,StandNames);
           if sum(js) == 0
              error(['Invalid value for the ''Standardization'' parameter: '...
                     'choices are ''on'' or ''off'' or ''Mahalanobis''.']);
           end
       Stand = StandNames{js}; 
       switch Stand
           
       case 'off'
           Xs = X - ones(n,1)*mean(X);            
       case 'on'
           Xs = zscore(X,1);          
       case 'Mahalanobis'
           Jc=eye(n)-(1/n)*ones(n);
           Sx=cov(X,1);
           Xs = Jc*X*Sx^-0.5; 
       end
    else  
        error(['Invalid value for the ''standardization'' parameter: '...
               'choices are ''on'' or ''off'' or ''.']);
    end
else 
    error(['Invalid value for the ''standardization'' parameter: '...
            'choices are ''on'' or ''off'' or ''.']);
end
% end Standardization

% Rndst 
if ~isempty(Rndst)  
    if isnumeric(Rndst)
       if (Rndst < 0) || (Rndst > 1000) 
       error('Rndst must be a value in the interval [0,1000]');
       end
    else
       error('Invalid Number of Random Starts');
    end
else
    error('Invalid Number of Random Starts')
end
% end Rndst

% MAxIter 
if ~isempty(MaxIter)  
    if isnumeric(MaxIter)
       if (MaxIter < 0) || (MaxIter > 1000) 
       error('MaxIter must be a value in the interval [0,1000]');
       end
    else
       error('Invalid Number of Max Iterations');
    end
else
    error('Invalid Number of Max Iterations')
end
% end MaxIter

% ConvToll 
if ~isempty(ConvToll)  
    if isnumeric(ConvToll)
       if (ConvToll < 0) || (ConvToll > 0.1) 
       error('ConvToll must be a value in the interval [0,0.1]');
       end
    else
       error('Invalid Convergence Tollerance');
    end
else
    error('Invalid Convergence Tollerance')
end
% end ConvToll 

% total deviance  
st=sum(sum(Xs.^2));


% Start the algorithm
for loop=1:Rndst
    
    %V=randPU(J,Q); % generate a random partition for variables (col.s)
    %U=randPU(n,K); % generate a random partition for units (rows) 
    %su=sum(U);
    %sv=sum(V);
    %
    [~,U,~,~]=kmeansN(Xs,K,4);
    [~,V,~,~]=kmeansN((Xs'),Q,3);
    su=sum(U);
    sv=sum(V);
 
    % compute initial medoid
    Ym = diag(1./su)*U'*Xs*V*diag(1./sv);
    B=U*Ym*V';
    f0=trace(B'*B)/st;
    %it=0;
    % iteration phase

        
        for it = 1:MaxIter
            
            % given Xm and V update U
            U=zeros(n,K);
            Ymv=Ym*V';
            for i=1:n
                mindif=sum((Xs(i,:)-Ymv(1,:)).^2);
                posmin=1;
                for k=2:K
                    dif=sum((Xs(i,:)-Ymv(k,:)).^2);
                    if dif < mindif
                        mindif=dif;
                        posmin=k;
                    end 
                end
                U(i,posmin)=1;
            end
            su=sum(U);
            while sum(su==0)>0, % if a class is empty split the largest class
                [m,p1]=min(su);
                [m,p2]=max(su);
                ind=find(U(:,p2));
                ind=ind(1:floor(su(p2)/2));
                U(ind,p1)=1;
                U(ind,p2)=0;
                su=sum(U);
            end 
            
            % given U and V updata Ym
      
            su=sum(U);
            Ym = diag(1./su)*U'*Xs*V*diag(1./sv);  

            % diven Xm and U updata V
    
            V=zeros(J,Q);
            Ymu=U*Ym;
            for j=1:J
                mindif=sum((Xs(:,j)-Ymu(:,1)).^2);
                posmin=1;
                for i=2:Q
                    dif=sum((Xs(:,j)-Ymu(:,i)).^2);
                    if dif < mindif
                        mindif=dif;
                        posmin=i;
                    end 
                end
                V(j,posmin)=1;
            end

            sv=sum(V);
            while sum(sv==0)>0, % if a class is empty split the largest class
                [m,p1]=min(sv);
                [m,p2]=max(sv);
                ind=find(V(:,p2));
                ind=ind(1:floor(sv(p2)/2));
                V(ind,p1)=1;
                V(ind,p2)=0;
                sv=sum(V);
            end 

            % given U and V updata Xm
      
            sv=sum(V);
            Ym = diag(1./su)*U'*Xs*V*diag(1./sv);

            B=U*Ym*V';
            f=trace(B'*B)/st;
%   
            fdif = f-f0;
        
            if fdif > ConvToll 
                f0=f; 
            else
                break
            end
        end
        fprintf('DKM: Loop=%g, Explained Variance =%g, iter=%g, fdif=%g\n',loop,f, it,fdif)   
        if loop==1
            Vdkm=V;
            Udkm=U;
            Ymdkm=Ym;
            fdkm=f;
            loopdkm=1;
            indkm=it;
            fdifo=fdif;
        end
        if f > fdkm
            Vdkm=V;
            Udkm=U;
            Ymdkm=Ym;
            fdkm=f;
            loopdkm=loop;
            indkm=it;
            fdifo=fdif;
        end
end
% sort clusters of variables per descending order of cardinality
[~,ic]=sort(diag(Vdkm'*Vdkm), 'descend');
Vdkm=Vdkm(:,ic);
% sort clusters of objects in descending order of cardinality
[~,ic]=sort(diag(Udkm'*Udkm), 'descend');
Udkm=Udkm(:,ic);
fprintf('DKM (Final): Explained Variance =%g, loopdpca=%g, iter=%g, fdif=%g\n',fdkm, loopdkm, indkm,fdifo) 


function [U]=randPU(n,c)

% generates a random partition of n objects in c classes with non-empty
% classes
%
% n = number of objects
% c = number of classes
%
U=zeros(n,c);
U(1:c,:)=eye(c);

U(c+1:n,1)=1;
for i=c+1:n
    U(i,[1:c])=U(i,randperm(c));
end
U(:,:)=U(randperm(n),:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% K-means algorithm %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [loopOtt,UOtt,fOtt,iterOtt]=kmeansN(X,K,Rndstart)
%
% n = number of objects
% J = number of variables
% K = number of clusters of the partition
%
% maxiter=max number of iterations


maxiter=100;
n = size(X,1);
J = size(X,2);
epsilon=0.000001;


% initial partition U0 is given

% best in a fixed number of partitions
%seed=200;             %si può rimettere per ritrovare le soluzioni ottenute
%rand('state',seed)    %
for loop=1:Rndstart
   U0=randPU(n,K);

   su=sum(U0);
   
   % given U compute Xmean (compute centroids)
   Xmean0 = diag(1./su)*U0'*X;

   for iter=1:maxiter
   
   % given Xmean0 assign each units to the closest cluster
   
        U=zeros(n,K);
        for i=1:n
            mindif=sum((X(i,:)-Xmean0(1,:)).^2);
            posmin=1;
            for j=2:K
                dif=sum((X(i,:)-Xmean0(j,:)).^2);
                if dif < mindif
                    mindif=dif;
                    posmin=j;
                end 
            end
            U(i,posmin)=1;
        end
   % given a partition of units 
   % i.e, given U compute Xmean (compute centroids)
   
        su=sum(U);
        while sum(su==0)>0,
            [m,p1]=min(su);
            [m,p2]=max(su);
            ind=find(U(:,p2));
            ind=ind(1:floor(su(p2)/2));
            U(ind,p1)=1;
            U(ind,p2)=0;
            su=sum(U);
        end 
   
   % given U compute Xmean (compute centroids)
        Xmean = diag(1./su)*U'*X;

   
   % compute objective function
   
        BB = U * Xmean - X;
        f = sum(BB(:).^2);
 
%  stopping rule
  
        dif=sum(sum((Xmean-Xmean0).^2));
        if dif > epsilon  
            Xmean0=Xmean;
        else
            break
        end
   end
   if loop==1
        UOtt=U;
        fOtt=f;
        loopOtt=1;
        iterOtt=1;
        XmeanOtt=Xmean;
   end
   if f < fOtt
        UOtt=U;
        fOtt=f;
        loopOtt=loop;
        iterOtt=iter;
        XmeanOtt=Xmean;
   end
end