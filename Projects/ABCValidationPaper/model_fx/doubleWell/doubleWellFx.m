function [xstore,tvec,wflag] = doubleWellFx(R,x,uc,p,m)


% Model priors
pQ.alist = [1 3];
pQ.flist = [14 30];
pQ.theta = [1 -1 1];

alist = pQ.alist.*exp(p.alist);
flist = pQ.flist.*exp(p.flist);
theta = pQ.theta.*exp(p.theta);



wflag = 0;
dt = R.IntP.dt;
t = [0:dt:180];
Wx(:,1:2) = randn(numel(t),2);
x = EM(x(:,1:2),t,dt,Wx,theta);

% subplot(2,1,1)
% plot(t,x(:,1),'r')
% hold on
% plot(t,x(:,2),'b')
% 
% subplot(2,1,2)
% plot(x(:,1),x(:,2))
%     plot(t,x(:,4),'b')

x(:,3) = mixSin(x(:,1),alist,flist,t);
xstore{1} = x';
tvec = t;

function mS = mixSin(X,alist,flist,t)
A(:,1) = alist(1).*sin(2*pi*flist(1)*t + pi/3);
A(:,2) = alist(2).*sin(2*pi*flist(2)*t - pi);

X = minmax(X);
X = [X 1-X];
mS = (X(:,1).*A(:,1)) + (X(:,2).*A(:,2));


function xnar = minmax(X)
xnar = (X-min(X))/(max(X)-min(X));

function x = EM(x,t,dt,Wx,theta)
for i = 1:(numel(t)-1)
    dx = dt*doublewell(t(i),x(i,:),theta);
    x(i+1,:) = (x(i,:)+dx') + (sqrt(dt).*Wx(i,:));
end

function dfx = doublewell(t,x,theta)
dfx(1) = x(2);
dfx(2) = -2*(x(1)-theta(1))*((x(1)-theta(2)).^2)...
    -(2*(x(1)-theta(1)).^2)*(x(1)-theta(2))-(theta(3)*x(2));
dfx = dfx';