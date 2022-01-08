clear; clc;
wp = 1;
ws = 100;
% nc = 1;
kc = logspace(-3,-1,100);
G = logspace(-3,-1,101);
kh = 0.001;
g = 2.5;
% G = 0.1;
nh = 1.0;
Th = ws/log(1/nh+1);
dissType = 'g';
dimp = 30;

nc = 0.1;
Tc = wp/log(1/nc+1);
% xth = sqrt(coth(wp/2/Tc));


%% 

g = 2.5;
x0 = g/wp/sqrt(2); %value according to paper definition. 

%In the function we use D(g/wp) for D^2, which is, according to the paper,
%D(sqrt(2)*x0)

delta = 2*wp*x0*x0;
%sigma = 0.5;
%dx = g/sqrt(2)/wp;
% %s = sigma*dx;
% f=1;
% for i = 1:length(Tc)
%     % Tc = logspace(-1,0.5,10);
%     for j = 1:length(delta)
%         tic
%         fprintf('g=%1.2f, step %i: delta = %1.2f -- ',g,j,delta(j));
%         M = modelCoherentEngine(dimp,ws,wp,g,G,kh,kc,Th,Tc(i),f,dissType,delta(j));
%         M.findSS;
%         [W(j), Qh(j), Qc(j)] = M.getThermoProp;
%         toc
%     end
%     clear M;
%     filename = sprintf('CoherentDataNEW3_Global_wp%.0f_ws%.0f_g%.2f_G%.4f_kh%.4f_kc%.4f_nc%.2f.mat',wp,ws,g,G,kh,kc,nc(i));
%     save(filename);
% end

f = @(x) (1-sign(x))/2;
for i = 1:length(kc)
    
    for j = 1:length(G)
        fprintf('(%d,%d)\n',i,j);
%         tic
%         fprintf('g=%1.2f, step %i: delta = %1.2f -- ',g,j,delta(j));
        M = modelCoherentEngine(dimp,ws,wp,g,G(j),kh,kc(i),Th,Tc,f,dissType,delta);
        M.findSS;
        [W(i,j), Qh(i,j), Qc(i,j)] = M.getThermoProp;
% 
    end

end
clear M;
filename = sprintf('ThesisCoherentData_Global_Step_wp%.0f_ws%.0f_g%.2f_nh%.2f_nc%.2f.mat',wp,ws,g,nh,nc);
save(filename);

%% 

% 
% x0 = 2.5;
% g = sqrt(2)*wp*x0;
% 
% delta = [-4:0.2:-1 linspace(-0.99,2,100) linspace(2.1,4,20)]*g*g/wp;
% 
% f=1;
% for i = 1:length(Tc)
%     % Tc = logspace(-1,0.5,10);
%     for j = 1:length(delta)
%         tic
%         fprintf('g=%1.2f, step %i: delta = %1.2f -- ',g,j,delta(j));
%         M = modelCoherentEngine(dimp,ws,wp,g,G,kh,kc,Th,Tc(i),f,dissType,delta(j));
%         M.findSS;
%         [W(j), Qh(j), Qc(j)] = M.getThermoProp;
%         toc
%     end
%     clear M;
%     filename = sprintf('CoherentDataNEW3_Global_wp%.0f_ws%.0f_g%.2f_G%.4f_kh%.4f_kc%.4f_nc%.2f.mat',wp,ws,g,G,kh,kc,nc(i));
%     save(filename);
% end
% 
% f = @(x) (1-sign(x))/2;
% for i = 1:length(Tc)
%     % Tc = logspace(-1,0.5,10);
%     for j = 1:length(delta)
%         tic
%         fprintf('g=%1.2f, step %i: delta = %1.2f -- ',g,j,delta(j));
%         M = modelCoherentEngine(dimp,ws,wp,g,G,kh,kc,Th,Tc(i),f,dissType,delta(j));
%         M.findSS;
%         [W(j), Qh(j), Qc(j)] = M.getThermoProp;
%         toc
%     end
%     clear M;
%     filename = sprintf('CoherentDataNEW3_Global_Step_wp%.0f_ws%.0f_g%.2f_G%.4f_kh%.4f_kc%.4f_nc%.2f.mat',wp,ws,g,G,kh,kc,nc(i));
%     save(filename);
% end
% 
% 
% %% 
% 
% 
% x0 = 2.0;
% g = sqrt(2)*wp*x0;
% 
% delta = [-4:0.2:-1 linspace(-0.99,2,100) linspace(2.1,4,20)]*g*g/wp;
% 
% f=1;
% for i = 1:length(Tc)
%     % Tc = logspace(-1,0.5,10);
%     for j = 1:length(delta)
%         tic
%         fprintf('g=%1.2f, step %i: delta = %1.2f -- ',g,j,delta(j));
%         M = modelCoherentEngine(dimp,ws,wp,g,G,kh,kc,Th,Tc(i),f,dissType,delta(j));
%         M.findSS;
%         [W(j), Qh(j), Qc(j)] = M.getThermoProp;
%         toc
%     end
%     clear M;
%     filename = sprintf('CoherentDataNEW3_Global_wp%.0f_ws%.0f_g%.2f_G%.4f_kh%.4f_kc%.4f_nc%.2f.mat',wp,ws,g,G,kh,kc,nc(i));
%     save(filename);
% end
% 
% f = @(x) (1-sign(x))/2;
% for i = 1:length(Tc)
%     % Tc = logspace(-1,0.5,10);
%     for j = 1:length(delta)
%         tic
%         fprintf('g=%1.2f, step %i: delta = %1.2f -- ',g,j,delta(j));
%         M = modelCoherentEngine(dimp,ws,wp,g,G,kh,kc,Th,Tc(i),f,dissType,delta(j));
%         M.findSS;
%         [W(j), Qh(j), Qc(j)] = M.getThermoProp;
%         toc
%     end
%     clear M;
%     filename = sprintf('CoherentDataNEW3_Global_Step_wp%.0f_ws%.0f_g%.2f_G%.4f_kh%.4f_kc%.4f_nc%.2f.mat',wp,ws,g,G,kh,kc,nc(i));
%     save(filename);
% end
% 
% 
% 
% %% Do a temperature run
% 
% 
% x0 = 2.5;
% g = sqrt(2)*wp*x0;
% 
% delta = 2*x0^2*wp;
% 
% N = 200;
% 
% %2*x0 must be greater than 1.
% %don't include xth=1, which amounts to Tc=0
% xth = 1 + (1:N)/N*(4*x0-1); %from 1 (Tc -> 0) to 2*x0
% Tc = wp/2./acoth(xth.^2);
% nc = 1./(exp(wp./Tc)-1);
% 
% 
% f=1;
% for i = 1:length(Tc)
%     tic
%     fprintf('g=%1.2f, step %i: nc = %1.2f -- ',g,i,nc(i));
%     M = modelCoherentEngine(dimp,ws,wp,g,G,kh,kc,Th,Tc(i),f,dissType,delta);
%     M.findSS;
%     [W(i), Qh(i), Qc(i)] = M.getThermoProp;
%     toc
%     clear M;
%     filename = sprintf('CoherentDataNEW3_Tc_Global_wp%.0f_ws%.0f_g%.2f_G%.4f_kh%.4f_kc%.4f.mat',wp,ws,g,G,kh,kc);
%     save(filename);
% end
% 
% f = @(x) (1-sign(x))/2;
% for i = 1:length(Tc)
%     tic
%     fprintf('g=%1.2f, step %i: nc = %1.2f -- ',g,i,nc(i));
%     M = modelCoherentEngine(dimp,ws,wp,g,G,kh,kc,Th,Tc(i),f,dissType,delta);
%     M.findSS;
%     [W(i), Qh(i), Qc(i)] = M.getThermoProp;
%     toc
%     clear M;
%     filename = sprintf('CoherentDataNEW3_Tc_Global_Step_wp%.0f_ws%.0f_g%.2f_G%.4f_kh%.4f_kc%.4f.mat',wp,ws,g,G,kh,kc);
%     save(filename);
% end
