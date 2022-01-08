function a = plotFisher(F,doPlot,state)
%takes in the Fisher information, compute alpha based on first two points,
%last two points and best fit
%and second to last points, returns the two alpha
set(0,'defaulttextinterpreter','latex');


np = length(F);
N = 1:np;
a1 = log10(F(2)/F(1))/log10(2);
a2 = log10(F(np)/F(np-1))/log10(np/(np-1));
a3 = polyfit(log10(2:np),log10(F(2:np)),1);
a = [a1 a2 a3(1)];
% a =[a1 a2];
if (~doPlot)
    return;
end
subplot(2,1,1);
switch state
    case 'g'
        loglog(N,F,'bo-');hold on;
        text(1.1,max(N.^a2/np^a2*F(np))*0.8,sprintf('$\\alpha_g = %.3f$',a2),'interpreter','latex','FontSize',14);

    case 'e'
        loglog(N,F,'ro-');hold on;
        text(1.1,max(N.^a2/np^a2*F(np))*0.4,sprintf('$\\alpha_e = %.3f$',a2),'interpreter','latex','FontSize',14);

    case 'p'
        loglog(N,F,'mo-');hold on;
        text(1.1,max(N.^a2/np^a2*F(np))*0.1,sprintf('$\\alpha_p = %.3f$',a2),'interpreter','latex','FontSize',14);

end
% loglog(N,F(1)*N.^a1,'k--');
loglog(N,N.^a2/np^a2*F(np),'k--'); hold on;
% loglog(N,N.^a3(1)*10^(a3(2)),'r--');
% text(1.1,max(F(1)*np.^a1)*0.9,sprintf('$\\alpha_{\\mathrm {fit}} = %.3f$',a3(1)),'interpreter','latex','FontSize',14);
% text(1.1,max(F(1)*np.^a1)*0.5,sprintf('$\\alpha_1 = %.3f,\\alpha_2 = %.3f$',a1,a2),'interpreter','latex','FontSize',14);
% text(1.1,max(F(1)*np.^a1)*0.9,sprintf('$\\alpha_1 = %.3f,\\alpha_2 = %.3f, \\alpha_{\\mathrm {fit}} = %.3f$',a1,a2,a3(1)),'interpreter','latex','FontSize',14);
xlabel('$N$','interpreter','latex','FontSize',14);
ylabel('$F_N$','interpreter','latex','FontSize',14);
subplot(2,1,2);
switch state
    case 'g'
        plot(N,F./N,'b');;hold on;
    case 'e'
        plot(N,F./N,'r');;hold on;
    case 'p'
        plot(N,F./N,'m');;hold on;
end

xlabel('$N$','interpreter','latex','FontSize',14);
ylabel('$F_N/N$','interpreter','latex','FontSize',14);
