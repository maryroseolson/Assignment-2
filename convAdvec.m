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

