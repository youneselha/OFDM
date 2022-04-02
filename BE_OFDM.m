close all;
clear all;

% Remarque : On considère qu'on travaille avec Fe = 1 HZ (rythme symbole)

% Generation de l'information binaire à transmettre 
%nb_bits = N ;
nb_bits = 12 * 10 ^ 4 ;
bits = randi([0, 1], 1, nb_bits);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Implantation de la chaine de transmission OFDM sans canal %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Nombre de porteuses
N = 16;

% Mapping BPSK
map =  2 * bits - 1;

% Traitement matricielle (chaque ligne est une porteuse / colonne -> temps)
map = reshape(map, N, length(map)/N);

% Construction d'une porteuse, 2 porteuses et 8 porteuses centrees resp.
map1 = [zeros( 4, nb_bits/N); map(5,:); zeros( 11 , nb_bits/N)];
%map1 = [ zeros( N-1, nb_bits/N); map(16,:)];
%map2 = [map(1:2, :); zeros( N-2, nb_bits/N)];
map2 = [zeros( 4, nb_bits/N); map(5,:); zeros( 4 , nb_bits/N); map(10,:); zeros( 6 , nb_bits/N)];
map8centre = [zeros( 4, nb_bits/N); map(5:12, :);zeros( 4, nb_bits/N)];

% Emission OFDM par ifft
signalOFDM = ifft(map);
signalOFDM1 = ifft(map1);
signalOFDM2 = ifft(map2);
signalOFDM8 = ifft(map8centre);

% Calcul de la DSP et visualisation
dsp = abs(fft(reshape(signalOFDM, 1, nb_bits))).^2;
dsp1 = abs(fft(reshape(signalOFDM1, 1, nb_bits))).^2;
dsp2 = abs(fft(reshape(signalOFDM2, 1, nb_bits))).^2;
dsp8 = abs(fft(reshape(signalOFDM8, 1, nb_bits))).^2;

f= [1:N/(length(dsp)-1):N+1];

figure('Name', "Chaine de transmission sans canal")
subplot(3,2,1)
plot(f, dsp1)
xlabel("Numero de la frequence porteuse")
title("DSP du signal pour une porteuse 5")
subplot(3,2,2)
semilogy(f, dsp1)
xlabel("Numero de la frequence porteuse")
title("DSP du signal pour une porteuse 5 -Echelle log-")
subplot(3,2,3)
plot(f, dsp2)
xlabel("Numero de la frequence porteuse")
title("DSP du signal pour 2 porteuses 5 et 10")
subplot(3,2,4)
semilogy(f, dsp2)
xlabel("Numero de la frequence porteuse")
title("DSP du signal pour 2 porteuses 5 et 10 -Echelle log-")
subplot(3,2,5)
plot(f, dsp8)
xlabel("Numero de la frequence porteuse")
title("DSP du signal pour 8 porteuses centrees")
subplot(3,2,6)
semilogy(f, dsp8)
xlabel("Numero de la frequence porteuse")
title("DSP du signal pour 8 porteuses centrees -Echelle log-")



% reception OFDM par fft
signalOFDMrecu = fft(signalOFDM);

% Decision
bitsdecision = sign(signalOFDMrecu);

% Demapping
bitsdemap = (bitsdecision + 1) / 2;
bitsdemap = reshape(bitsdemap, 1, nb_bits);

% Taux d'erreur binaire
TEB = 1 - (length(find(bitsdemap == bits)) / nb_bits);
fprintf("Le TEB en chaine de transmission sans canal : %i\n ", TEB);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Implantation de la chaine de transmission OFDM avec canal multi-trajets sans bruit %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% IG = alpha * (T + IG) avec alpha =  20 à 25 % et T = N * Ts
% Ts = 1 -> N = (IG * ( 1 - alpha) ) / alpha
% y(t) = 0.227x(t) + 0.46x(t −Ts) + 0.688x(t −2Ts) + 0.460x(t −3Ts) + 0.227x(t −4Ts)
% hc = [0.227, 0.46, 0.688, 0.460, 0.227] -> IGmin = 4

% Calcul du nombre de porteuses N pour alpha = 0.25
alpha = 0.25;
N = floor((4 * (1 - alpha)) / alpha);

% Réponse impulsionnelle du canal (Ts = 1)
hc = [0.227, 0.46, 0.688, 0.460, 0.227];

% Réponse fréquentielle et visualisation
Rf = fft(hc, 2^N);
f = [1:N/(length(Rf) - 1):N+1];

figure('Name', "Chaine de transmission avec canal")
subplot(3,1,1)
plot(f, abs(Rf))
xlabel("Numero de la frequence porteuse")
title("Reponse en fréquence (module)")
subplot(3,1,2)
semilogy(f, abs(Rf));
xlabel("Numero de la frequence porteuse")
title("Reponse en fréquence (module) -Echelle log-")
subplot(3,1,3)
plot(f, rad2deg(angle(Rf)))
xlabel("Numero de la frequence porteuse")
title("Reponse en fréquence (phase) en deg")


% Mapping BPSK
    mapcanal =  2 * bits - 1;

% Traitement matricielle et emission OFDM
    mapcanal = reshape(mapcanal, N, length(mapcanal)/N);

    % Emission sans intervalle de garde 
    signalOFDMcanal = ifft(mapcanal);
    % Emission avec intervalle de garde 
    signalOFDMGarde = [zeros(4,nb_bits/N); signalOFDMcanal]; % IG = 4
    % Emission avec prefixe cyclique
    signalOFDMPC = [signalOFDMcanal(end-3:end,:);signalOFDMcanal];
    % Emission avec prefixe cyclique surdimensionne (double)
    signalOFDMPC_surdim = [signalOFDMcanal(end-7:end,:);signalOFDMcanal];

% Passage par canal
    % Sans intervalle de garde
    SignalRecuCanal = filter(hc,1,reshape(signalOFDMcanal, 1,[]));
    % Avec intervalle de garde
    SignalRecuCanalGarde = filter(hc,1,reshape(signalOFDMGarde, 1, [])) ;
    % Avec prefixe cyclique
    SignalRecuPC = filter(hc, 1, reshape(signalOFDMPC, 1, []));
    % Avec prefixe cyclique surdimensionne
    SignalRecuPC_surdim = filter(hc, 1, reshape(signalOFDMPC_surdim, 1, []));

    % Traitement matricielle
    SignalRecuCanal = reshape(SignalRecuCanal, N, []) ;
    SignalRecuCanalGarde = reshape(SignalRecuCanalGarde, N + 4, []) ;
    SignalRecuPC = reshape(SignalRecuPC, N + 4, []) ;
    SignalRecuPC_surdim = reshape(SignalRecuPC_surdim, N + 8, []) ;

% Calcul du DSP du signal avant et apres passage par le canal sans intervalle 
% de garde ni prefixe cyclique et visualisation

dspAvantCanal = abs(fft(reshape(signalOFDMcanal, 1, nb_bits))).^2;
dspCanal = abs(fft(reshape(SignalRecuCanal, 1, nb_bits))).^2;

figure('Name', "Chaine de transmission avec canal")
subplot(2,2,1)
plot(1:N/(length(dspAvantCanal)-1):N+1, dspAvantCanal)
xlabel("Numero de la frequence porteuse")
title("DSP du signal avant passage par canal")
subplot(2,2,2)
semilogy(1:N/(length(dspAvantCanal)-1):N+1, dspAvantCanal)
xlabel("Numero de la frequence porteuse")
title("DSP du signal avant passage par canal -Echelle log-")
subplot(2,2,3)
plot(1:N/(length(dspCanal)-1):N+1, dspCanal)
xlabel("Numero de la frequence porteuse")
title("DSP du signal apres passage par canal")
subplot(2,2,4)
semilogy(1:N/(length(dspCanal)-1):N+1, dspCanal)
xlabel("Numero de la frequence porteuse")
title("DSP du signal apres passage par canal -Echelle log-")

% Reception OFDM (fft)
    % Sans intervalle de garde
    signalOFDMRecuCanal = fft(SignalRecuCanal);
    % Avec intervalle de garde
    signalOFDMRecuCanalGarde = fft(SignalRecuCanalGarde(5:end,:));
    % Avec prefixe cyclique
    signalOFDMRecuPC = fft(SignalRecuPC(5:end,:));

% Reception OFDM (PC avec erreur d'horloge)
    % Attention a la rotation introduite par la non synchronisation
    % exp(2pi*l*Tau/N) *** l : le numéro de porteuse *** tau : le retard *** N : nombre de
    % porteuses (voir cours)
    % donc pour récompenser cette rotation, il faut l'estimer en realite ou 
    % bien faire une estimation globale comme si on a un nouveau canal ! 

    signalOFDMRecuPC_cas1 = fft(SignalRecuPC_surdim(3:2+N,:));
    signalOFDMRecuPC_cas2 = fft(SignalRecuPC_surdim(7:6+N,:));%7 % le cas 6 un peu special
    signalOFDMRecuPC_cas3 = fft([SignalRecuPC_surdim(11:end,:);SignalRecuPC_surdim(1:2,:)]);

% pc et égalisation :
    H = fft(hc,N);
    % Egalisation ZFE
    signalRecuPC_EG = signalOFDMRecuPC ./ (H.');
    % Egalisation ML
    signalRecuPC_EG_ML = signalOFDMRecuPC .* conj(H.');

    % Avec erreur d'horloge
    siganlRecuPC_EG_cas1 = signalOFDMRecuPC_cas1 ./ (H.');
    siganlRecuPC_EG_cas2 = signalOFDMRecuPC_cas2 ./ (H.');
    siganlRecuPC_EG_cas3 = signalOFDMRecuPC_cas3 ./ (H.');

% Constellation en reception sur 2 porteuses
    
    % Sans et avec intervalle de garde 
    figure('Name', "Chaine de transmission avec canal")
    subplot(2,2,1)
    plot(signalOFDMRecuCanal(2, :), '*')
    title('porteuse 2 sans garde')
    subplot(2,2,2)
    plot(signalOFDMRecuCanal(12,:), '*')
    title('porteuse 12 sans garde')
    subplot(2,2,3)
    plot(signalOFDMRecuCanalGarde(2,:), '*')
    title('porteuse 2 avec garde')
    subplot(2,2,4)
    plot(signalOFDMRecuCanalGarde(12,:), '*')
    title('porteuse 12 avec garde')

    % Sans et avec prefixe cyclique
    figure('Name', "Chaine de transmission avec canal")
    subplot(2,2,1)
    plot(signalOFDMRecuCanal(2, :), '*')
    title('porteuse 2 sans pc')
    subplot(2,2,2)
    plot(signalOFDMRecuCanal(12,:), '*')
    title('porteuse 12 sans pc')
    subplot(2,2,3)
    plot(signalOFDMRecuPC(2,:), '*')
    title('porteuse 2 avec pc')
    subplot(2,2,4)
    plot(signalOFDMRecuPC(12,:), '*')
    title('porteuse 12 avec pc')

    % Avec prefixe cyclique et egalisation ZFE ou ML
    % Remarque : on est en 10^-16 --> donc c'est presque le meme point.
    figure('Name', "Chaine de transmission avec canal")
    subplot(3,2,1)
    plot(signalOFDMRecuCanal(2, :), '*')
    title('porteuse 2 sans PC et Egal.')
    subplot(3,2,2)
    plot(signalOFDMRecuCanal(12,:), '*')
    title('porteuse 12 sans PC et Egal.')
    subplot(3,2,3)
    plot(signalRecuPC_EG(2,:), '*')
    title('porteuse 2 avec PC et Egal. ZFE')
    subplot(3,2,4)
    plot(signalRecuPC_EG(12,:), '*')
    title('porteuse 12 avec PC et Egal. ZFE')
    subplot(3,2,5)
    plot(signalRecuPC_EG_ML(2,:), '*')
    title('porteuse 2 avec PC et Egal. ML')
    subplot(3,2,6)
    plot(signalRecuPC_EG_ML(12,:), '*')
    title('porteuse 12 avec PC et Egal. ML')

    % Avec erreur d'horloge (les 3 cas) egalisation ZFE
    figure('Name', "Chaine de transmission avec canal")
    subplot(2,4,1)
    plot(signalOFDMRecuCanal(2, :), '*')
    title('porteuse 2 sans PC et Egal.')
    subplot(2,4,5)
    plot(signalOFDMRecuCanal(12,:), '*')
    title('porteuse 12 sans PC et Egal.')
    
    subplot(2,4,2)
    plot(siganlRecuPC_EG_cas1(2, :), '*')
    title('porteuse 2 PC et Egal. cas1')
    subplot(2,4,6)
    plot(siganlRecuPC_EG_cas1(12,:), '*')
    title('porteuse 12 PC et Egal. cas1')
    
    subplot(2,4,3)
    plot(siganlRecuPC_EG_cas2(2, :), '*')
    title('porteuse 2 PC et Egal. cas2')
    subplot(2,4,7)
    plot(siganlRecuPC_EG_cas2(12,:), '*')
    title('porteuse 12 PC et Egal. cas2')
    
    subplot(2,4,4)
    plot(siganlRecuPC_EG_cas3(2, :), '*')
    title('porteuse 2 PC et Egal. cas 3')
    subplot(2,4,8)
    plot(siganlRecuPC_EG_cas3(12,:), '*')
    title('porteuse 12 PC et Egal. cas 3')

% Decision
signalOFDMRecuCanal= reshape(signalOFDMRecuCanal, 1, []);
bitsdecisioncanal = sign(real(signalOFDMRecuCanal));

signalOFDMRecuCanalGarde = reshape(signalOFDMRecuCanalGarde, 1, []);
bitsdecisioncanalGarde = sign(real(signalOFDMRecuCanalGarde));

signalOFDMRecuPC = reshape(signalOFDMRecuPC, 1, []);
bitsdecisionpc = sign(real(signalOFDMRecuPC));

signalRecuPC_EG = reshape(signalRecuPC_EG, 1, []);
bitsdecisionpc_eg = sign(real(signalRecuPC_EG));

signalRecuPC_EG_ML = reshape(signalRecuPC_EG_ML, 1, []);
bitsdecisionpc_eg_ml = sign(real(signalRecuPC_EG_ML));

% Demapping
bitsdemapcanal = (bitsdecisioncanal + 1) / 2;
bitsdemapcanalGarde = (bitsdecisioncanalGarde + 1) / 2;
bitsdemappc = (bitsdecisionpc + 1) / 2;
bitsdemappc_eg = (bitsdecisionpc_eg + 1) / 2;
bitsdemappc_eg_ml = (bitsdecisionpc_eg_ml + 1) / 2;

%TEB CANAL
TEBCanal = 1 - (length(find(bitsdemapcanal == bits)) / nb_bits);
fprintf("Le TEB en chaine de transmission avec canal multi-trajets " + ...
    "sans intervalle de garde : %i\n ", TEBCanal);
TEBCanalGarde = 1 - (length(find(bitsdemapcanalGarde == bits)) / nb_bits);
fprintf("Le TEB en chaine de transmission avec canal multi-trajets " + ...
    "avec intervalle de garde : %i\n ", TEBCanalGarde);
TEBPC = 1 - (length(find(bitsdemappc == bits)) / nb_bits);
fprintf("Le TEB en chaine de transmission avec canal multi-trajets " + ...
    "avec prefixe cyclique  : %i\n ", TEBPC);
TEB_PC_EG = 1 - (length(find(bitsdemappc_eg == bits)) / nb_bits);
fprintf("Le TEB en chaine de transmission avec canal multi-trajets " + ...
    "avec prefixe cyclique et egalisation ZFE : %i\n ", TEB_PC_EG);
TEB_PC_EG_ML = 1 - (length(find(bitsdemappc_eg_ml == bits)) / nb_bits);
fprintf("Le TEB en chaine de transmission avec canal multi-trajets " + ...
    "avec prefixe cyclique et egalisation ML : %i\n ", TEB_PC_EG_ML);