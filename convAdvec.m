function [p] = convAdvec(h,sigma,mu,q)
%   Input:
%           h: scalar. cellsize
%       sigma: vector. spatially varying coefficient, evaluated at cell centres.
%          mu: scalar. constant coefficient
%           q: vector. Source function, evaluated at nodes.
%  
%   Output: 
%           p: Approximate value of p on the nodes.
%           
%   You must write your discretization as a matrix equation A*p = q.
%   Once you have formed the matrix A, you can solve for p using the 
%   Matlab backslash operator. Make sure that you form A as a sparse
%   matrix. If A is dense the backslash operator will be very inefficient.
%% Defining centers, nodes and the inside centers
n=1/h;
xC = [(h/2):h:(1-(h/2))];
xN = [0:h:1];
xi = xN(2:n);

%% Defining Diffusion Portion
Dnc=1/h*spdiags(ones((n),1)*[-1 1],[0,1],(n)-1,(n));
Dcn = -transpose(Dnc);
diffusion = Dnc*sigma*Dcn;
%% Defining Convection
Avg = (1/2)*spdiags(ones((n),1)*[1 1],[0,1],(n)-1,(n));
convection = mu.*Avg*Dcn;
%% Defining A
A = diffusion + convection;
%% Defining P
p = A\q'
end

%% Test Function (Function produced 2 sinuisoidal functions)

% If h becomes greater than 10^-4 the function will not run
for i = 2:4
h=10.^(-i);
n=1/h;
% x on the nodes
xN = [0:h:1];
% x on the nodes excluding the end points
xi = xN(2:n);
% x on the cell centers
xC = [(h/2):h:(1-(h/2))];
mu = 0.1;
%mu = 10;
% creating sigma points
sig = (1+xC.^2)';
Sig = spdiags(sig, [0], n, n);
% 
diffusion = 4*pi*xi.*cos(2*pi*xi)-4*pi^2*(xi.^2 + 1).*sin(2*pi*xi);
advection = mu*2*pi*cos(2*pi*xi);
q =  diffusion + advection;
Pexp = convAdvec(h,Sig,mu,q);
Ptho = sin(2*pi*xC);
i+1;

figure
hold on
plot(Pexp)
plot(Ptho)
hold off
pause(10);
end
