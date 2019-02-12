function [Y,offset] = MoyenneGlissante(N, speed)
    //N = nombre de point défini de la moyenne glissante et Speed vitesse a moyennée sur intervalle ..)
    h = ones (N,1);
    Y = conv(speed,h)./N;
    offset = N/2+1
    //Y = x(offset+2:length(Y),1)
endfunction
