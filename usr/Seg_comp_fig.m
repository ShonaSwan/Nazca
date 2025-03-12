
% Segregation Compaction Length
d_sc = sqrt(KD .* zeta);  

% making sure dimensions are the same
wm_in = wm(2:end, 2:end-1);   
wx_in = wx(2:end, 2:end-1);    
Pc_in = Pc(2:end-1, 2:end-1); 

% Horizontal means
w_m_mean = mean(wm_in, 2);     
w_x_mean = mean(wx_in, 2); 
delta_sc_mean = mean(d_sc, 2);  
Pc_mean = mean(Pc_in, 2);       

% Create figure
figure;

% Subplot 1: Segregation Velocity vs. Segregation Compaction Length
subplot(1, 2, 1);
plot(w_x_mean, delta_sc_mean, 'b-', 'LineWidth', 2, 'DisplayName', 'Crystals (wx)');
hold on;
plot(w_m_mean, delta_sc_mean, 'r-', 'LineWidth', 2, 'DisplayName', 'Melt (wm)');
xlabel('Segregation Velocity');
ylabel('Segregation Compaction Length');
legend;
grid on;

% Subplot 2: Compaction Pressure vs. Segregation Compaction Length
subplot(1, 2, 2);
plot(Pc_mean, delta_sc_mean, 'b-', 'LineWidth', 2, 'DisplayName', 'Solid Matrix (Pc)');
xlabel('Compaction Pressure (Pa)');
ylabel('Segregation Compaction Length');
legend;
grid on;

