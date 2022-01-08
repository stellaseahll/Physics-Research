NT = 




% % spectrum variables give the spectrum that yields the optimum spectrum of
% % a 2 qubit system with energies [0,spectrum_Tmax]. 
% % Degeneracyis not considered here
% % Minimization of the RHS of the inequality in Eq.(23) of notes as
% % discussed
% 
% load('data.mat');
% %Tmin = 1, Tmax = 2.1, alpha = 1
% subplot(2,1,1);
% loglog(x,y_Tmax2);
% title('Tmax=2.1');
% xlabel('Total number of qubits (or log_2 d)');
% ylabel('RHS of bound in Eq.(23)');
% %Tmin = 1, Tmax = 10.1, alpha = 1
% loglog(x,y_Tmax10);
% title('Tmax=10.1');
% xlabel('Total number of qubits (or log_2 d)');
% ylabel('RHS of bound in Eq.(23)');
