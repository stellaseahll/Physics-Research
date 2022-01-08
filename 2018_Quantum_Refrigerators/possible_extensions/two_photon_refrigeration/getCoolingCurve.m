TH = linspace(1,5,50);
THtmp = interp(TH,10)
TW = linspace(1,21,50);
r1 = [1:-0.1:0; 0:0.1:1];
THzero = zeros(length(r1),length(TW));
for i = 1:length(r1)
    for j = 1:length(TW)
        QCtmp = QC(:,j);
        QCtmp = interp(QCtmp,100);
        THzero(i,j) = THtmp(abs(QCtmp)==min(abs(QCtmp)));
    end
end
    
    
    