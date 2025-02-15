% a forward-time centered-space finite difference scheme
% dt is the timestep
% dx is the spatial stepsize
% nx is the number of spatial steps
% alpha is the diffusion coeffiecient
% U0 is the input concentration profile

function y = ftcs(dt,dx,nx,alpha,U0)

r  = alpha*dt/(dx^2);   % scaled diffusion coefficient
r2 = 1-2*r;             % second coefficient

y = zeros(size(U0));    % set up vector for solution

% loop
y(1) = r*U0(2) + r2*U0(1) + r*U0(nx);        % boundary is a ring
for i = 2:nx-1
    y(i) = r*U0(i+1) + r2*U0(i) + r*U0(i-1); % ftcs
end
y(nx) = r*U0(1) + r2*U0(nx) + r*U0(nx-1);    % boundary is a ring

end