% Computing pseudo F statistics for clustering

function [pf,Dw,Db] = psF(X,UOtt)
[n,k]=size(UOtt);
if k==1 
    % categorical variable to be coded 
    mx=max(UOtt);
    Ik=eye(mx);
    UOtt=Ik(UOtt,:);
    k=mx;
end
Xm=pinv(UOtt)*X;
Dw=trace((X-UOtt*Xm)'*(X-UOtt*Xm));
Db=trace((UOtt*Xm)'*(UOtt*Xm));
pf=(Db/(k-1))./(Dw/(n-k));
end