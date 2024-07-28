% X (n X J) data matrix
% V (J x Q) membership matrix for clustering variables
% U (n x k) membership matrix for clustering objects
%
% model X = UYmV'B + E
% 
% problem: min||X-UYmV'B||^2
% 
% being ||X||^2 = ||X-UYmV'B||^2 + ||UYmV'B||^2
%
% problem maximize ||UYmV'B||^2
% subject to
% V binary and row stochastic
% U binary and row stochastic
% B Diagonal matrix

function [Vcdpca,Ucdpca,Acdpca, Ycdpca,fcdpca,incdpca]=CDPCA(X, K, Q, Rndstart)

% initialization

maxiter=100;
% convergence tollerance
eps=0.0000000001;

opts.disp=0;
VC=eye(Q);

[n,J]=size(X);


% Standardize data
Xs=zscore(X,1);

% compute var-covar matrix 

S=cov(Xs,1);
st=(1./n)*sum(sum(Xs.^2));


for loop=1:Rndstart
    V=randPU(J,Q);
    U=randPU(n,K);
    su=sum(U);
    % update dentroid matrix Xmean 
    Xmean = diag(1./su)*U'*Xs;
    it=0;
    JJ=[1:J]';
    A=zeros(J,Q);
    zJ=zeros(J,1);
    for g=1:Q
        ibCg=V(:,g); 
        JCg=[JJ(ibCg==1)];
        Sg=S(JCg,JCg);
        if sum(ibCg)>1
                [a,c]=eigs(Sg,1,'lm',opts);
                A(JCg,g)=a;
            else
                A(JCg,g)=1;
            end
    end
    Ymean=Xmean*A;
    f0=trace((1./n)*Ymean'*U'*U*Ymean);
    fmax=0;

%   iteration phase
    fdif=2*eps;
    while fdif > eps || it>=maxiter,
      it=it+1;
   
      % given Ymean and A update U
      Y=Xs*A;
      U=zeros(n,K);
        for i=1:n
            mindif=sum((Y(i,:)-Ymean(1,:)).^2);
            posmin=1;
            for j=2:K
                dif=sum((Y(i,:)-Ymean(j,:)).^2);
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
            [~,p1]=min(su);
            [~,p2]=max(su);
            ind=find(U(:,p2));
            ind=ind(1:floor(su(p2)/2));
            U(ind,p1)=1;
            U(ind,p2)=0;
            su=sum(U);
        end 

        % given U and A compute Xmean (compute centroids)
        Xmean = diag(1./su)*U'*Xs;
        
        % given U and Ymean update V and A
        for j=1:J
            posmax=JJ(V(j,:)==1);
            for g=1:Q
                V(j,:)=VC(g,:);
                ibCg=V(:,g);           % new class of V
                ibCpm=V(:,posmax);     % old class of V
                JCg=[JJ(ibCg==1)];
                JCpm=[JJ(ibCpm==1)];
                Sg=S(JCg,JCg);
                S0g=S(JCpm,JCpm);
                if sum(ibCg)>1
                     [a,c]=eigs(Sg,1,'lm',opts);
                      A(:,g)=zJ;                     
                      if sum(a)<0
                          a=-a;
                      end
                      A(JCg,g)=a;
                else
                      A(:,g)=zJ;
                      A(JCg,g)=1;
                end
                if sum(ibCpm)>1
                     [aa,cc]=eigs(S0g,1,'lm',opts);
                     if sum(aa)<0
                        aa=-aa;
                     end

                     A(:,posmax)=zJ;  
                     A(JCpm,posmax)=aa;
                else
                      A(:,posmax)=zJ;  
                      A(JCpm,posmax)=1;
                end  
                Ymean = Xmean*A;
                f=trace((1./n)*Ymean'*(U'*U)*Ymean);
                if f > fmax
                    fmax=f;
                    posmax=g;
                    A0=A;
                else
                    A=A0;
                end
            end
            V(j,:)=VC(posmax,:);
        end
            
        %Y=Xs*A; 
        Ymean = Xmean*A;
        f=trace((1./n)*Ymean'*(U'*U)*Ymean);
        fdif = f-f0;
        
        if fdif > eps 
            f0=f;fmax=f0; A0=A;
        else
            break
        end
    end
  fprintf('CDPCA: Loop=%g, Explained variance=%g, iter=%g, fdif=%g\n',loop,f./st*100, it,fdif)   
       if loop==1
            Vcdpca=V;
            Ucdpca=U;
            Acdpca=A;
            Ycdpca=Xs*Acdpca;
            fcdpca=f;
            loopcdpca=1;
            incdpca=it;
            fdifo=fdif;
        end
   if f > fcdpca
       Vcdpca=V;
       Ucdpca=U;
       fcdpca=f;
       Acdpca=A;
       Ycdpca=Xs*Acdpca;
       loopcdpca=loop;
       incdpca=it;
       fdifo=fdif;
   end
end
% sort clusters of variables in descend order of variance
varYcdpca=var(Ycdpca,1);
[~,ic]=sort(varYcdpca, 'descend');
Acdpca=Acdpca(:,ic);
Vcdpca=Vcdpca(:,ic);
Ycdpca=Ycdpca(:,ic); 
% sort clusters of objects in descending order of variance
%dwc=zeros(K,1);
%for k=1:K
%dwc(k)= trace((Ycdpca-Ucdpca*pinv(Ucdpca)*Ycdpca)'*diag(Ucdpca(:,k))*(Ycdpca-Ucdpca*pinv(Ucdpca)*Ycdpca));
%end
[~,ic]=sort(diag(Ucdpca'*Ucdpca), 'descend');
Ucdpca=Ucdpca(:,ic);
fprintf('CDPCA (Final): Percentage Explained Variance=%g, loopdpca=%g, iter=%g, fdif=%g\n',fcdpca./st*100, loopcdpca, incdpca,fdifo)

 