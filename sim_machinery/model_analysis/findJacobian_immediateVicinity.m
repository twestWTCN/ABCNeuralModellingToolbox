function [J Es Js] =findJacobian_immediateVicinity(R,x,uc,pc,m,condsel)
% pertubations across delays


ucm = {};
for rm = 1:size(uc{condsel},2)
    %     ucm{condsel}{rm} = zeros(R.IntP.bufferExt+ (plist(end)./R.IntP.dt),size(uc{condsel}{rm},2));
    ucm{condsel}{rm} = repmat(mean(uc{condsel}{rm}),R.IntP.bufferExt+1,1);
end
% Now set to zero
X{1} = R.condnames{condsel};
R.condnames = X;
% computes the Jacobian of a function
n=size(x,1);
xvec = x;
% x = zeros(size(x));
R.IntP.nt = R.IntP.bufferExt; % remember the first "buffer" isnt actually simulated
fx = R.IntP.intFx(R,x,ucm,pc,m);
fx = fx{1}(:,end);
eps=exp(-8);
xperturb= x;
for i=1:n
        xperturb(i,:)=xperturb(i,:)+eps;
        fx_st = R.IntP.intFx(R,xperturb,ucm,pc,m);
        fx_st = fx_st{1}(:,end);
        J(:,i)=(fx_st-fx)/eps;
        xperturb(i,:)= xvec(i,:);
end
e = eig(J);
Js(:,:) = J;
Es = e(end);

J = Js; %mean(Js,3); % average over delays
% J = Js(:,:,1); % take the first delay
for i = 1:size(Js,3)
    e(:,i) = real(eig(squeeze(Js(:,:,i))));
end
