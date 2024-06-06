clc;
clear all;

% Constants and parameters
ro = 19.3; % CZT density (g/cm^3)
Eo = 3; % Initial energy (MeV)
N = 6e3; % Number of particles
M = 500; % Number of steps for energy distribution
L = 5; % Material thickness (cm)

% Energy distribution range
Edisp = logspace(-3, log10(Eo), M)';

% Starting values
Eleak = zeros(M, 1);
leak = 0; pair = 0; comp = 0; foto = 0; unc = 0;

% Install the absorption coefficient data
mu = load('W_mu.txt');

% Monte Carlo simulation
for b = 1:N
    E = Eo;
    u = 2;
    
    % Calculate the initial absorption coefficient
    k = find(mu(:,1) >= E, 1) - 1;
    mutot0 = interp1(mu(k:k+1,1), mu(k:k+1,5), E) * ro;
    l = -log(rand) / mutot0;
    
    if l > L
        unc = unc + 1;
    end
    
    while u > 1
        k = find(mu(:,1) >= E, 1) - 1;
        mutot = interp1(mu(k:k+1,1), mu(k:k+1,5), E) * ro;
        muc = interp1(mu(k:k+1,1), mu(k:k+1,2), E) * ro;
        muf = interp1(mu(k:k+1,1), mu(k:k+1,3), E) * ro;
        mup = interp1(mu(k:k+1,1), mu(k:k+1,4), E) * ro;
        
        d = l - log(rand) / mutot;
        
        if d < L && d > 0
            p = rand;
            pc = muc / (mup + muc + muf);
            pf = muf / (mup + muc + muf);
            
            if p <= pc
                comp = comp + 1;
                E = E / (1 + (E / 0.511) * (1 - rand));
            elseif p <= pc + pf
                foto = foto + 1;
                E = 0;
            else
                pair = pair + 1;
                E = 0.511;
                d2 = d - log(rand) / mutot;
                
                while E > 0
                    k = find(mu(:,1) >= E, 1) - 1;
                    mutot = interp1(mu(k:k+1,1), mu(k:k+1,5), E) * ro;
                    muc = interp1(mu(k:k+1,1), mu(k:k+1,2), E) * ro;
                    muf = interp1(mu(k:k+1,1), mu(k:k+1,3), E) * ro;
                    mup = interp1(mu(k:k+1,1), mu(k:k+1,4), E) * ro;
                    
                    d2 = d2 - log(rand) / mutot;
                    
                    if d2 < L && d2 > 0
                        p = rand;
                        pc = muc / (mup + muc + muf);
                        pf = muf / (mup + muc + muf);
                        
                        if p <= pc
                            comp = comp + 1;
                            E = E / (1 + (E / 0.511) * (1 - rand));
                        elseif p <= pc + pf
                            foto = foto + 1;
                            E = 0;
                        else
                            pair = pair + 1;
                            E = 0.511;
                            d2 = d2 - log(rand) / mutot;
                        end
                    else
                        u = 1;
                        leak = leak + 1;
                        h = find(Edisp >= E, 1);
                        Eleak(h) = Eleak(h) + 1;
                        break;
                    end
                    
                    if E < 1e-3
                        break;
                    end
                end
            end
        else
            u = 1;
            leak = leak + 1;
            h = find(Edisp >= E, 1);
            Eleak(h) = Eleak(h) + 1;
        end
        
        if E < 1e-3
            break;
        end
    end
end

% Calculation and viewing results
Eleak(M) = unc;
D = Eleak .* Edisp;
Dtot = sum(D);
Dunc = D(M);
B_dose = Dtot / Dunc;
B_flux = leak / unc;
uteo = N * exp(-mutot0 * L);
error = abs(uteo - unc) / unc * 100;

disp(['Number of particle: ', num2str(N)]);
disp(['Number of collision: ', num2str(comp)]);
disp(['Photoelectric events: ', num2str(foto)]);
disp(['Pair production: ', num2str(pair)]);
disp(['All leak: ', num2str(leak)]);
disp(['Uncollision leak: ', num2str(unc)]);
disp(['Theoritical uncollision: ', num2str(round(uteo))]);
disp(['Relative error (%): ', num2str(error)]);
disp(['Dose Buildup Factor: ', num2str(B_dose)]);
disp(['Flux Buildup Factor: ', num2str(B_flux)]);

% Grafiği çizdirme
figure;
non_zero_indices = Eleak > 0;
plot(Edisp(non_zero_indices), Eleak(non_zero_indices), 'o');
set(gca, 'YScale', 'log');
xlabel('Energy of \gamma (MeV)','fontsize',20),
ylabel('Number of \gamma','fontsize',20),
title('Linear \gamma Ray (6x10^6 - 3 MeV) Through an Slab (0.5 cm)','fontsize',20)
grid on
disp('===== Program Bitti =====');
