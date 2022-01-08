function getGgammaMixed(z,nbar,dnbar,np,isSteady)

G = 0.02:0.02:1.0;
G(end) = [];
% G(G==0) = [];
gammat = logspace(-3,1,40);
rhop = [z 0; 0 1-z];

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

filename = sprintf('pz%.3f,nbar%.3f,np%d,isSteady%d.mat',z,nbar,np,isSteady);
save(filename);

end