function runDiffState(G,nbar,np)

    getDiffState(G*pi,10,nbar,0.0001*nbar,np,1);
    getDiffState(G*pi,10,nbar,0.0001*nbar,np,0);
    getDiffState(G*pi,1,nbar,0.0001*nbar,np,1);
    getDiffState(G*pi,1,nbar,0.0001*nbar,np,0);
    getDiffState(G*pi,0.1,nbar,0.0001*nbar,np,1);
    getDiffState(G*pi,0.1,nbar,0.0001*nbar,np,0);
    getDiffState(G*pi,0.01,nbar,0.0001*nbar,np,1);
    getDiffState(G*pi,0.01,nbar,0.0001*nbar,np,0);
    getDiffState(G*pi,0.001,nbar,0.0001*nbar,np,1);
    getDiffState(G*pi,0.001,nbar,0.0001*nbar,np,0);
    
end