function [x, y, u] = elliptic_dir_2D(h, g, f, mu)

%%  solution of the 2D elliptic equation -mu * Δu = f with Dirichlet boundary conditions on the unit square (0,1)x(0,1)
    
% input 
    %   h: space discretization step
    %   g: function handle for the Dirichlet boundary conditions
    %   f: function handle for the source term
    %   mu: diffusion coefficient
% output
    %   x, y: grid points in x and y directions direction
    %   u: numerical solution at the grid points


% initialization of the space grid
M = floor(1/h);
h = 1/M; % adjust h to fit the grid exactly
x = linspace(0,1,M+1)'; 
y = linspace(0,1,M+1)'; 
N = M-1; % number of internal points in each direction

% initialization of the full solution matrix, including boundary points: rows correspond to x, columns correspond to y
u = zeros(M+1, M+1);

% impose Dirichlet boundary conditions
u(1, :)   = g(x', 0);          % bottom boundary, y = 0
u(end, :) = g(x', 1);          % top boundary, y = 1
u(:, 1)   = g(0, y);           % left boundary, x = 0
u(:, end) = g(1, y);           % right boundary, x = 1

% internal grid points
x_int = x(2:end-1);
y_int = y(2:end-1);

% contruction of the matrix A for internal nodes
alpha = 4*mu/h^2; % diagonal entries
beta = -mu/h^2;   % off-diagonal entries

e = ones(N,1);
D = spdiags([beta*e, alpha*e, beta*e], [-1, 0, 1], N, N);
E = beta*speye(N); 

A = kron(speye(N),D) + kron(spdiags([e,e],[-1,1],N,N),E);

% construction of the right-hand side matrix F for internal nodes
F = zeros(N,N); 

% F internal nodes values 
for j = 1:N
    for i = 1:N
        F(j, i) = f(x_int(i), y_int(j));
    end
end

% F corrections for the boundary conditions
F(1, :)   = F(1, :)   + (mu / h^2) * u(1, 2:end-1);       % bottom
F(end, :) = F(end, :) + (mu / h^2) * u(end, 2:end-1);     % top
F(:, 1)   = F(:, 1)   + (mu / h^2) * u(2:end-1, 1);       % left
F(:, end) = F(:, end) + (mu / h^2) * u(2:end-1, end);     % right

% solve the linear system AU = F(:)
U = A\F(:); % solve for internal nodes, F(:) vectorizes F column-wise, output U is a column vector 

% reconstruct the internal solution matrix
u(2:end-1, 2:end-1) = reshape(U, N, N); % reshape U back to a N x N matrix, column-wise

end