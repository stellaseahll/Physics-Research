function runDiffMixedState(G,gammat,nbar,np,NP)

    for i = 1:length(G)
        for j = 1:length(gammat)
            getDiffMixedState(G(i)*pi,gammat(j),nbar,0.0001*nbar,np,1,NP);
        end
    end
    
end