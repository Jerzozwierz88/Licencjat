% Space at which we consider wavelengths and wave numbers
lambdas = logspace(0, 2, 1000); P = 2.*pi ./ lambdas;
% define rheological parameters
[R, N0, N1, U] = deal(0.01, 7, 7, 16);
% figure specifications
f = figure('Renderer', 'painters', 'Position', [10 10 900 600]);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.8, 0.8]);
% base - just q curve imported from Mathematica
semilogx(lambdas, 1-fold(P, R, N0, N1, U), 'k--', 'LineWidth', 1, 'DisplayName','1'); 
hold on;

% shortenings
for s=[0.9, 0.8, 0.7, 0.6]
    tau = log(1/s);
    % ode solver
    exponent = 1/tau .* integral(@(t) 1-fold(P.*exp(2*t), R, N0, N1, U), 0, ...
        tau, "ArrayValued",true);
    % solution to SAS integral in dimensionless time for given shortening
    semilogx(lambdas*s, exponent, 'LineWidth', 0.8, 'DisplayName', num2str(s)); hold on;
end
tilt = append(' n_0=', num2str(round(N0)), ' n_1=', num2str(round(N1)), ...
    ' r=', num2str(R), ' u=', num2str(U));
title({'Efektywny współczynnik wzrostu dla interfejsów w fazie dla'; tilt}); 
xlabel('Końcowa znormalizowana długość fali \lambda(\tau)/H_0'); 
ylabel('$\textit{sgn}(\bar{D}_{xx})(-1+q_{eff})$', Interpreter='latex');
lg = legend('show'); lg.Title.String = 'Wydłużenie L/L_0'; 
hold off;