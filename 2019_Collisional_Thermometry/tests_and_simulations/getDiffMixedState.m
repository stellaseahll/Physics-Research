function getDiffMixedState(G,gammat,nbar,dnbar,np,isSteady,NP)

z = linspace(0,1,NP);
for i = 1:length(z) 
    fprintf('%.1f%% complete....\n',i/length(z)*100);
    rhop = [z(i) 0; 0 1-z(i)];        
    M = spSinglePassXY(G,gammat,nbar,dnbar,rhop,np,isSteady);
    [F(i,:) Fz(i,:) Fy(i,:)] = M.getAllFish();
    alpha(i,:) = M.getAllExponent(F(i,:));
    alphaz(i,:) = M.getAllExponent(Fz(i,:));
    alphay(i,:) = M.getAllExponent(Fy(i,:));
end

filename = sprintf('DiffMixedState,G%.3fpi,gt%.3f,nbar%.3f,np%d,isSteady%d.mat',G/pi,gammat,nbar,np,isSteady);
save(filename);

end