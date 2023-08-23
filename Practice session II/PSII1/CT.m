% Parameters
delta = 0.1;
alpha_1 = 0.03;
alpha_2 = 0.05;

Iss=1;
Xss=Iss/delta;
alpha_0= (1+alpha_1/delta-alpha_2)*Iss;

% Varying parameter b_1 range
b_1_range = linspace(0.6, 0.95, 100);

% Initial conditions
X0 = Xss;
I0 = Iss;

% Preallocate arrays to store results
X_vals = zeros(length(b_1_range), 1);
I_vals = zeros(length(b_1_range), 1);

% Preallocate arrays to store eigenvalues
eigenvalues_real = zeros(length(b_1_range), 2);
eigenvalues_imag = zeros(length(b_1_range), 2);

% Loop through varying parameter b_1 values
for i = 1:length(b_1_range)
    b_1 = b_1_range(i);
    
    % Define the system of equations
    dydt = @(t, y) [
        - delta * y(1) + y(2);
       (1/(1-b_1))* ( alpha_0 - alpha_1 * y(1) + (alpha_2 + b_1 - 1) * y(2) )
    ];
    
    % Solve the system using ODE solver
    [t, y] = ode45(dydt, [0 100], [X0; I0]);
    
    % Store the final values for plotting
    X_vals(i) = y(end, 1);
    I_vals(i) = y(end, 2);
    
    % Calculate Jacobian matrix
    J = [
        - delta, 1;
        -alpha_1/(1-b_1), (alpha_2 + b_1 - 1)/(1-b_1)
    ];
    
    % Calculate eigenvalues
    eigvals = eig(J);
    eigenvalues_real(i, :) = real(eigvals);
    eigenvalues_imag(i, :) = imag(eigvals);
end

% Plot results
% figure;
% 
% subplot(3, 1, 1);
% plot(b_1_range, X_vals);
% xlabel('b_1');
% ylabel('X');
% title('X vs. b_1');
% 
% subplot(3, 1, 2);
% plot(b_1_range, I_vals);
% xlabel('b_1');
% ylabel('I');
% title('I vs. b_1');

figure;
plot(b_1_range, eigenvalues_real(:, 1), 'b', b_1_range, eigenvalues_real(:, 2), 'r');
hold on;
plot(b_1_range, eigenvalues_imag(:, 1), 'bo', b_1_range, eigenvalues_imag(:, 2), 'ro');
xlabel('b_1');
ylabel('Eigenvalues');
title('Real and Imaginary Parts of Eigenvalues');
legend('Real (Eigenvalue 1)', 'Real (Eigenvalue 2)', 'Imaginary (Eigenvalue 1)', 'Imaginary (Eigenvalue 2)');
hold off;


