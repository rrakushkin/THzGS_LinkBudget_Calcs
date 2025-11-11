function [snr_db] = receivedSNR(p_rx_dbm, nf_db, bw)
    %LINKCAPACITY Summary of this function goes here
    %
    % Input:
    %   p_rx_dbm   [dBm]    Received power from the link budget
    %   nf_db      [dB]     Receiver noise figure
    %   bw         [Hz]     Bandiwdth
    %
    % Output:
    %   capacity   [bps]    Shannon capacity for the AWGN channel

    % Constants
    kb = physconst('Boltzman');

    % System Noise Temperature
    Ta = 300;                           % [K] Antenna noise temperature when pointing at the Earth (Worst case)
    T0 = 290;                           % [K] Reference noise temperature
    Te = (10^(nf_db/10)-1)*T0;          % [K] Receiver equivalent noise temperature = (F-1)T0
    T_sys = (Ta + Te);                  % [K] System noise temperature

    n_dbm = 10*log10(kb) + 10 * log10(T_sys * bw) + 30;     % [dBm] Noise power
    snr_db = p_rx_dbm - n_dbm;
end

