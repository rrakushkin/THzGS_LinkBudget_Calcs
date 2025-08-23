% Coarsen the grid for a lighter scatter
[EEs, DDs, SSs] = ndgrid(linspace(0.65,0.90,21), linspace(1.5,2.0,21), linspace(20,520,21)*1e-6);
Gcloud = 10*log10( exp(-(4*pi*SSs./lambda).^2).*EEs ) + 20*log10(pi*DDs./lambda);

figure;
scatter3(DDs(:), EEs(:)*100, SSs(:)*1e6, 18, Gcloud(:), 'filled'); % color = gain
xlabel('D (m)'); ylabel('\eta (%)'); zlabel('\sigma (\mum)');
title('Gain colored by dBi'); grid on; colorbar
