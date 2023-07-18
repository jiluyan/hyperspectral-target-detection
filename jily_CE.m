function [CE,u,wu,energy]=jily_CE(Y,d,u,gama)
% clever eye (CE), for hyperspectral/multispectral target detection
%Geng, Xiurui, Luyan Ji, and Kang Sun. 2016. "Clever eye algorithm for target detection of remote sensing imagery." 
%ISPRS Journal of Photogrammetry and Remote Sensing 114:32-9. doi: https://doi.org/10.1016/j.isprsjprs.2015.10.014.
%
%Input: Y: the data matrix, L*N, L is the number of bands, N is the number
%of pixels
%       d:the target vector, L*1
%       u: the initialization for the data origin
%       gama: optimize step
%Output: CE: CE output, 1*L
%        u: the final CE point
%        wu: CE detector
%        energy: the average filter output energy series during the
%        optimization process
% Luyan Ji, jily@mail.ustc.edu.cn
%2023.7.18


m=mean(Y,2);

K1=cov(Y');
[~,nN]=size(Y);
Y1=Y-repmat(m,[1,nN]);
K=Y1*Y1'/nN;

thresu=0.00001;
thresN=100;

deltau=1;
N=0;

energy=zeros(1,thresN);
 uall=u;
while deltau >= thresu && N < thresN 

Y1=Y-repmat(u,[1,nN]);
Ru=Y1*Y1'/nN;
wu=inv(Ru)*(d-u)/((d-u)'*inv(Ru)*(d-u));
Ece=wu'*(Y-repmat(u,[1,nN]));
Ece=Ece*Ece'/nN;
energy(N+1)=Ece;

g11=-2*inv(K)*(d-u);
g12=2*(d-u)'*inv(K)*(m-u)*(1+(m-u)'*inv(K)*(m-u))*inv(K)*(-m+2*u-d)+2*((d-u)'*inv(K)*(m-u))^2*inv(K)*(m-u);
g12=g12/(1+(m-u)'*inv(K)*(m-u))^2;
g1=g11-g12;

u1=u+gama*g1;
deltau=norm(g1);
u=u1;
N=N+1;
 uall=[uall,u];

end

Y1=Y-repmat(u,[1,nN]);
Ru=Y1*Y1'/nN;
wu=inv(Ru)*(d-u)/((d-u)'*inv(Ru)*(d-u));

CE=wu'*Y1;

end