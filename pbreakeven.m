function pbe=pbreakeven(s,q)
    pbe=((1-(1-q)^s)^(1/s)-(1-q))/q;
end