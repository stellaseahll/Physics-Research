function getDiffnbar(theta,G,gammat,np,isSteady)
rhop = [1+cos(theta*pi) sin(theta*pi); sin(theta*pi) 1-cos(theta*pi)]/2; 
nbar = logspace(-2,2,40);
for i = 1:length(nbar) 
    fprintf('%.1f%% complete....\n',i/length(nbar)*100);
    M = spSinglePassXY(G*pi,gammat,nbar(i),nbar(i)*1e-5,rhop,np,isSteady);
    [F(i,:) Fz(i,:) Fy(i,:)] = M.getAllFish();
    alpha(i,:) = M.getAllExponent(F(i,:));
    alphaz(i,:) = M.getAllExponent(Fz(i,:));
    alphay(i,:) = M.getAllExponent(Fy(i,:));
end

filename = sprintf('DiffNbar,G%.3fpi,gt%.3f,theta%.3fpi,np%d,isSteady%d.mat',G,gammat,theta,np,isSteady);
save(filename);

end