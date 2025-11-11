function p_rx_dbm = linkBudget_Simplified(p_tx_dbm, g_sat, g_dish, f_c, distance, l_abs, l_total_ptg)
    % Input:
    %   p_tx_dbm   [dBm]     
    %   f_c        [Hz]
    %   distance   [m]             
    %   l_abs (atmospheric absorption loss)     [dB]           
    %   Output:
    %   p_rx_dbm   [dBm]
  
    c       = physconst('LightSpeed');
    % Spreading loss
    l_spr = 10*log10((4*pi*distance*f_c/c)^2);


    p_rx_dbm = p_tx_dbm + g_sat + g_dish - l_spr - l_abs - l_total_ptg;
end

