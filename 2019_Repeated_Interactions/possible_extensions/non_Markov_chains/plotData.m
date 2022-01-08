figure;
c = colormap(lines);
for j = 1:length(tau)
    subplot(2,2,j)
    for i = 1:6
        plot(tstep(1:10000),real(ez2{j,i}(1:10000)),'color',c(i,:),'linewidth',2); hold on;
        %plot(tstep(1:250:10000),real(ez2{j,i}(1:250:10000)),'color',c(i,:),'Marker','o','MarkerSize',2,'LineStyle','none','LineWidth',3); hold on;
    end
end