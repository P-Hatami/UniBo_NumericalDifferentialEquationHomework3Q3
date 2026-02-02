% Domain
a = 0; b = 1;
c = 0; d = 1;

% Exact solution
u_exact = @(x, y) sin(pi*x).*sin(pi*y);

% Source term
f = @(x, y) (2*pi^2 + y) .* sin(pi*x) .* sin(pi*y);

% Boundary condition
g = u_exact;  % Dirichlet everywhere

% Dummy functions for Neumann/Robin
s = @(x, y) 0; v = @(x, y) 0; r = @(x, y) 0;

% Boundary flags
BC_flags = {'Dirichlet', 'Dirichlet', 'Dirichlet', 'Dirichlet'};

% Grid resolutions to test
resolutions = [20, 40, 80];
errors = zeros(size(resolutions));

for k = 1:length(resolutions)
    Nx = resolutions(k);
    Ny = resolutions(k);

    % Build system
    [A, F] = ellip2D_FD(f, g, s, v, r, a, b, c, d, Nx, Ny, BC_flags);

    % Solve
    U = A \ F;

    % Compare with exact
    x = linspace(a, b, Nx);
    y = linspace(c, d, Ny);
    [X, Y] = meshgrid(x, y);
    U_exact_grid = u_exact(X, Y);
    U_num_grid = reshape(U, Ny, Nx);

    % Compute error (L2 norm)
    errors(k) = sqrt(sum((U_num_grid(:) - U_exact_grid(:)).^2) / numel(U));
end

% Display convergence
disp('Grid size vs. L2 error:');
disp(table(resolutions', errors', 'VariableNames', {'GridSize', 'L2_Error'}));
