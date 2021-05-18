function [J Es] =findJacobian(R,x,uc,pc,m,condsel)
uc = repmat({zeros(size(uc{condsel}))},1,2);
X{1} = R.condnames{condsel};
R.condnames = X;
% computes the Jacobian of a function
n=size(x,1);
xvec = x;
% x = zeros(size(x));
plist = 0:0.025:0.3;
for p = 1:length(plist)
    R.IntP.nt = R.IntP.bufferExt + (plist(p)./R.IntP.dt);
    fx = R.IntP.intFx(R,x,uc,pc,m);
    fx = fx{1}(:,end);
    eps=1.e-12;  % could be made better
    xperturb= x;
    for i=1:n
        xperturb(i,:)=xperturb(i,:)+eps;
        fx_st = R.IntP.intFx(R,xperturb,{uc{condsel}+eps},pc,m);
        fx_st = fx_st{1}(:,end);
        J(:,i)=(fx_st-fx)/eps;
        xperturb(i,:)= xvec(i,:);
    end
    e = eig(J);
    Js(:,:,p) = J;
    Es(p) = e(1);
end% J

% J = mean(Js,3);
J = Js(:,:,1);
