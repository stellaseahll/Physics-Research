function getGgamma(theta,nbar,dnbar,np,isSteady)

G = 0.025:0.025:1.0;
G(end) = [];
% G(G==0) = [];
gammat = logspace(-3,1,40);
rhop = [1+cos(theta*pi) sin(theta*pi); sin(theta*pi) 1-cos(theta*pi)]/2;

for i = 1:length(G)
    for j = 1:length(gammat)
        fprintf('%.1f%% complete....\n',((i-1)*length(G)+j)/length(G)/length(gammat)*100);
        M = spSinglePassXY(G(i),gammat(j),nbar,dnbar,rhop,np,isSteady);
        [F{i,j} Fz{i,j} Fy{i,j}] = M.getAllFish();
        alpha{i,j}= M.getAllExponent(F{i,j});
        alphaz{i,j} = M.getAllExponent(Fz{i,j});
        alphay{i,j} = M.getAllExponent(Fy{i,j});
    end
end

filename = sprintf('theta%.3f,nbar%.3f,np%d,isSteady%d.mat',theta,nbar,np,isSteady);
save(filename);

end