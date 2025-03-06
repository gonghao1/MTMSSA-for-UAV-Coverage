function yy=Mutate(x,mu,sigma)
%   x=p.Position; 
    for i=1:size(x,1)-1,
        x1=x{i};
    nVar=numel(x1);
    nMu=ceil(mu*nVar);
    j=randsample(nVar,nMu);   
    y=x1;
    ty=randn(size(j,1));
    y(j)=x1(j)+sigma*ty(1,:);
    x{i}=y;
    end
  yy=x;
end