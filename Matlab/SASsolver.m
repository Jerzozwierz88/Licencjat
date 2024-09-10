clear
global h R U L s N0 N1;
% define rheological parameters
h = 0.5; [R, N0, N1, U] = deal(0.01, 3, 2, 6); 
% figure specifications
f = figure(1);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.8, 0.8]);
clf;
t = tiledlayout(2, 2); t.TileSpacing = 'compact'; t.Padding = 'compact';
nexttile;


% generate initial interface
M = 1000; L=60;
xr = linspace(-L, L, M);
red1 = rednoise(1, M)*1e-6; red1 = red1 - mean(red1);
red2 = rednoise(1, M)*1e-6; red2 = red2 - mean(red2);

% plot initial interface in space (x)
plot(xr, red1, 'k'); ylabel('$\frac{A}{H}$',Interpreter='latex',Rotation=0);
title('Początkowa perturbacja (górny interfejs)');
grid on;

% select shortenings for plotting

for s = [0.9, 0.8, 0.72]
    %evolve interface
    red_sym = (red1+red2)/2;
    red_asym = (red1-red2)/2;
    red_sym = amplification(red_sym, @fold);
    red_asym = amplification(red_asym, @neck);
    
    % compose with assymetric part - which almost does not evolve
    new_interface = red_sym + red_asym;
    new_interface2 = red_sym - red_asym;
    
    nexttile; hold on;
    % plot evolved interface
    plot(xr*s, new_interface*s, 'color', "#77AC30", 'LineWidth', 1.2);
    title(append('Górny interfejs po skróceniu o ', num2str(round((1-s)*100, 2)), '%'));
    grid on; ylabel('$\frac{A}{H}$',Interpreter='latex',Rotation=0);
    hold off;
end
% finish figure specifications
tilt = append('   n_0=', num2str(round(N0)), '   n_1=', num2str(round(N1)), ...
    '   r=', num2str(R), '   u=', num2str(U));
t.Title.String = {'Porównanie wyników FEM i modelu SAS dla'; tilt};

% plot layer
f2 = figure(2);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.3, 0.1, 0.7, 0.8]);
hold on;
% plot interfaces
plot(xr*s, h/s+new_interface, 'k', 'LineWidth', 3);
plot(xr*s, -h/s+new_interface2, 'k', 'LineWidth', 3);
% fill area in between
fill([xr*s, fliplr(xr*s)], [-h/s+new_interface2, fliplr(h/s+new_interface)],[0.4 0.8 0.1]);
xlabel('Położenie');
ylabel('Amplituda'); ylim([-10*h, 10*h]);
title({'Obraz dwóch interfejsów  dla parametrów'; [tilt, sprintf('   s=%0.2f', round(s, 2))]});
hold off;

function [new_interface] = amplification(interface, fun)
    global h R U L s N0 N1;
    % we perform fft to retrieve even modes
    % trick - we mirror interface to get BCs used in folder
    M = length(interface);
    red = fft([interface, interface(end:-1:1)]);
    % cutoff after fft - relevant modes for amplification
    red = red(1:M+1);
    % define space of wave numbers k
    xf = 0:1:M; xf = xf * 2*pi /4/L * (M-1)/M * 2*h;
    
    % amplification
    % perform SAS on whole spectrum of initial wave numbers
    factor = zeros(1, length(xf));
    factor(2:end) = integral(@(t) 1-fun(xf(2:end).*exp(2*t), R, N0, N1, U), ...
        0,log(1/s), "ArrayValued",true);
    % evolve spectrum
    factor(isnan(factor)) = log(1/s);
    amplified_red = red .* exp(factor);
    % compose interface
    % prepare for ifft - add negative wave numbers
    amplified_red = [amplified_red, conj(amplified_red(end-1:-1:2))];
    % real part gives A_r*cos(kx) + A_i *sin(kx)
    new_interface = real(ifft(amplified_red));
    % take only part without mirror reflection
    new_interface = new_interface(1:M);
end