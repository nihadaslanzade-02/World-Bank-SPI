% %%%%%%%%%%%%%%%%%%%%%%
% Factorial K-means    % 
%                      %
%%%%%%%%%%%%%%%%%%%%%%%%
%
% Nihad Aslanzade July 2024
%
% input
% X (n X J) data matrix
% A (J x Q) loading matrix for dimensionality reduction
% U (n x k) membership matrix for clustering objects
%
% model XA = UYm + E
% 
% problem: min||XA-UYm||^2
% 
% being ||XA||^2 = ||XA-UYm||^2 + ||UYmA||^2
%
% equivalent problem
%
% problem maximize ||UYm||^2
% subject to
%
% U binary and row stochastic
% A Orthonormal A'A=IJ
%
function [Ufkm,Afkm, Yfkm,ffkm,infkm]=FKM(X, K, Q, Rndstart)


%
% initialization
%
maxiter=100;
% convergence tollerance
eps=0.0000000001;

opts.disp=0;
VC=eye(Q);

[n,J]=size(X);

% Standardize data
Xs=zscore(X,1);

st=sum(sum(Xs.^2));

un=ones(n,1);
uk=ones(K,1);
um=ones(Q,1);


for loop=1:Rndstart
    U=randPU(n,K);
    su=sum(U);
    Xmean = diag(1./su)*U'*Xs;
    it=0;
    % update A
    XX=Xs'*U*diag(1./su)*U'*Xs;
    [A,L]=eigs(XX,Q);
    Ymean = Xmean*A;
    Y=Xs*A;
    f0=trace(Ymean'*U'*U*Ymean)/st;
% iteration phase
    fdif=2*eps;
    while fdif > eps | it>=maxiter,
        it=it+1;
      % given Ymean update U
        U=zeros(n,K);
        for i=1:n
            mindif=sum((Y(i,:)-Ymean(1,:)).^2); %factorial k-means

            posmin=1;
            for j=2:K
                dif=sum((Y(i,:)-Ymean(j,:)).^2); %factorial k-means
                if dif < mindif
                    mindif=dif;
                    posmin=j;
                end 
            end
            U(i,posmin)=1;
        end
     %
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
        Xmean = diag(1./su)*U'*Xs;
      
        % given U and Xmean update A
        XX=Xs'*U*diag(1./su)*U'*Xs;
        [A,L]=eigs(XX,Q);
        A=A(:,1:Q);
        Ymean = Xmean*A;
        Y=Xs*A;
        f=trace(Ymean'*U'*U*Ymean)/st;
        fdif = f-f0;
        
        if fdif > eps 
            f0=f; A0=A; 
        else
            break
        end
    end
  fprintf('FKM: Loop=%g, Explained variance=%g, iter=%g, fdif=%g\n',loop,f*100, it,fdif)   
       if loop==1
            Ufkm=U;
            Afkm=A;
            Yfkm=Xs*Afkm;
            ffkm=f;
            loopfkm=1;
            infkm=it;
            fdifo=fdif;
        end
   if f > ffkm
       Ufkm=U;
       ffkm=f;
       Afkm=A;
       Yfkm=Xs*Afkm;
       loopfkm=loop;
       infkm=it;
       fdifo=fdif;
   end
end
% sort components in descend order of variance
% and ritate factors
varYrkm=var(Yfkm,1);
[~,ic]=sort(varYrkm, 'descend');
Afkm=Afkm(:,ic);
Yfkm=Yfkm(:,ic); 
Afkm=rotatefactors(Afkm);
% sort clusters of objects in descending order of cardinality
%
[~,ic]=sort(diag(Ufkm'*Ufkm), 'descend');
Ufkm=Ufkm(:,ic);
fprintf('FKM (Final): Percentage Explained variance=%g, looprkm=%g, iter=%g, fdif=%g\n',ffkm*100, loopfkm, infkm,fdifo)

 function [U]=randPU(n,c)

% generates a random partition of n objects in c classes
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