function p_rx_dbm = linkBudget1(p_tx_dbm, g_sat, directivity_g, gsHPBW, f_c, nElem, sep_dist, distance, l_abs, sat_offAxisAngle, gs_offAxisAngle, gs_polarization, pAngle)
    % LINKBUDGET (insert description)
    %
    % Input:
    %   p_tx_dbm   [dBm]     
    %   f_c        [Hz] 
    %   sep_dist    [m]
    %   distance   [m]             
    %   d = GS dish diameter [m]
    %   offAxisError = pointing error of GS dish [degrees]
    %   gs_polarization = polarization type of GS ('circular' or 'linear')....the satellite is assumed to be linear
    %   pAngle = polarization angle between the tx and rx (degrees)
    %   l_abs (atmospheric absorption loss)     [dB]           
    %   surf_rms    [um] 
    %   Output:
    %   p_rx_dbm   [dBm]
    
    % Constants
     c = physconst('LightSpeed');
     lambda = c ./ f_c;
    
    
    % Spreading loss
    %disp(distance)
    l_spr = 20*log10(distance * 1e-3) + 20*log10(f_c * 1e-9) + 92.45;

    % Pointing loss
    satHPBW = (50.8 .* lambda) ./ (nElem*sep_dist); %taken from payload calculations
    sat_ptg_loss = 12.*((sat_offAxisAngle./ satHPBW).^2); %ITU-R BO.790 (page 3)
    gs_ptg_loss = 12.*((gs_offAxisAngle ./ gsHPBW).^2); %ITU-R BO.790 (page 3) 
    l_ptg = gs_ptg_loss + sat_ptg_loss;
    

    % Polarization loss facotr (PLF)
    switch gs_polarization
        case 'linear'
            l_plf = 10*log10((cosd(pAngle)).^2)
        otherwise
            l_plf = -3
    end
    %disp(l_spr)
    p_rx_dbm = p_tx_dbm + g_sat + directivity_g - l_spr - l_abs - l_ptg + l_plf;

end

