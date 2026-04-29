% task 3
mu = 2; 

f = @(x,y) -8 + 0*x + 0*y; 
g = @(x,y) x.^2 + y.^2;
uex = @(x,y) x.^2 + y.^2;

h0 = 1/4;
halv = 5; 

h_vec = zeros(halv+1,1);    
err = zeros(halv+1,1);

for r = 1:halv+1
    h = h0/2^(r-1);
    [x, y, u] = elliptic_dir_2D(h, g, f, mu);
    
    [X, Y] = ndgrid(x,y); 
    Uex = uex(X,Y);      % exact solution on the full grid

    err(r) = max(max(abs(u - Uex)));
    h_vec(r) = h;
end

% estimated convergence orders
p = log2(err(1:end-1)./err(2:end));

disp(table(h_vec, err, ...
    'VariableNames', {'h', 'Error'}));

disp(table(h_vec(2:end), p(:), ...
    'VariableNames', {'h', 'EstimatedOrder'}));

% log-log plot
figure;
loglog(h_vec, err, 'o-', 'LineWidth', 1.5);
grid on;
xlabel('h');
ylabel('E_\infty');
title('Error of the 5-point finite difference scheme');

exportgraphics(gcf, '/Users/stefano/figs3/nme_proj_err.pdf', 'ContentType', 'vector');