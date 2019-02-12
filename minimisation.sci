function e = minimisation(c,x,v)
    e = (c(1)*(1-exp(-((x-c(2))/c(3)))) - v);
endfunction 
