function [ch,wm,wc]=msp(x,P,kappa)
%UKF sigma point calculation (1 and 3 scaling parameters based, basic/recursive)
%function [ch,w]=msp(x,P,kappa)

nx = length(x);


% alpha2 = alpha^2;
% lambda = alpha2 * (nx + kappa) - nx;
% c = nx + lambda;
% pom = sqrt(c)*chol(P)';
% 
% ch(:,1) = x;
% wm(1)=lambda/c;
% wc(1)=lambda/c+(1-alpha^2+beta);
% i = 1;
% for j=2:2:2*(nx)   
%     ch(:,j) = x - pom(:,i);
%     ch(:,j+1) = x + pom(:,i);
%     wm(j)=0.5/(c);
%     wm(j+1)=0.5/(c);
%     wc(j)=0.5/(c);
%     wc(j+1)=0.5/(c);
%     i = i + 1;
% end


pom = chol((nx+kappa)*P)';
ch(:,1) = x;
w(1)=kappa/(nx+kappa);
i = 1;
for j=2:2:2*nx   
    ch(:,j) = x - pom(:,i);
    ch(:,j+1) = x + pom(:,i);
    w(j)=0.5/(nx+kappa);
    w(j+1)=0.5/(nx+kappa);
    i = i + 1;
end
wm = w;
wc = w;


