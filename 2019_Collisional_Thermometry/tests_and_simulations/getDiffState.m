function getDiffState(G,gammat,nbar,dnbar,np,isSteady)

theta = (0:0.01:2);
for i = 1:length(theta) 
    fprintf('%.1f%% complete....\n',i/length(theta)*100);
    rhop = [1+cos(theta(i)*pi) sin(theta(i)*pi); sin(theta(i)*pi) 1-cos(theta(i)*pi)]/2;        
    M = spSinglePassXY(G,gammat,nbar,dnbar,rhop,np,isSteady);
    [F(i,:) Fz(i,:) Fy(i,:)] = M.getAllFish();
    alpha(i,:) = M.getAllExponent(F(i,:));
    alphaz(i,:) = M.getAllExponent(Fz(i,:));
    alphay(i,:) = M.getAllExponent(Fy(i,:));
end

filename = sprintf('G%.3fpi,gt%.3f,nbar%.3f,np%d,isSteady%d.mat',G/pi,gammat,nbar,np,isSteady);
save(filename);

end