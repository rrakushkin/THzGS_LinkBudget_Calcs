function p_rx_dbm = linkBudget_Simplified(p_tx_dbm, g_sat, g_dish, f_c, distance, l_abs)
    % LINKBUDGET (insert description)
    %
    % Input:
    %   p_tx_dbm   [dBm]     
    %   f_c        [Hz]
    %   distance   [m]             
    %   l_abs (atmospheric absorption loss)     [dB]           
    %   Output:
    %   p_rx_dbm   [dBm]
  
    
    % Spreading loss
    l_spr = 20*log10(distance * 1e-3) + 20*log10(f_c * 1e-9) + 92.45;


    p_rx_dbm = p_tx_dbm + g_sat + g_dish - l_spr - l_abs;
end

