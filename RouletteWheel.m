function k=RouletteWheel(P)

P=P./sum(P);
P=cumsum(P);

   k=find(rand<=P,1,'first');





end