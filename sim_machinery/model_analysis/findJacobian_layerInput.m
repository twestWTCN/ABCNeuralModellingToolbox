function [J Es] =findJacobian_layerInput(R,x,uc,pc,m,condsel)
% pertubations across delays
plist = 0:0.025:0.3;


ucm = {};
for rm = 1:size(uc{condsel},2)
    ucm{condsel}{rm} = zeros(R.IntP.bufferExt+ (plist(end)./R.IntP.dt),size(uc{condsel}{rm},2));
end
X{1} = R.condnames{condsel};
R.condnames = X;
% computes the Jacobian of a function
n=size(x,1);
xvec = x;
% x = zeros(size(x));
for p = 1:length(plist)
    R.IntP.nt = R.IntP.bufferExt + (plist(p)./R.IntP.dt);
    fx = R.IntP.intFx(R,x,ucm,pc,m);
    fx = fx{1}(:,end);
    eps=1.e-12;  % could be made better
    xperturb= x;
    for i=1:n
        xperturb(i,:)=xperturb(i,:)+eps;
        fx_st = R.IntP.intFx(R,xperturb,ucm,pc,m);
        fx_st = fx_st{1}(:,end);
        J(:,i)=(fx_st-fx)/eps;
        xperturb(i,:)= xvec(i,:);
    end
    e = eig(J);
    Js(:,:,p) = J;
    Es(p) = e(end);
end% J

J = mean(Js,3);
% J = Js(:,:,1);
