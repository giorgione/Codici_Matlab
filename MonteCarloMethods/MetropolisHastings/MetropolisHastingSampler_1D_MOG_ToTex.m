FID = fopen('MetropolisHastingSimulationImage.tex', 'w');
fprintf(FID, '\\begin{figure}\\label{fig: SimulationMetropolisHasting0}\n');

%Generate some statistic
j=0;
for i=1:2:K
if(mod(i,8)==1 && i>1)
    j=j+1;
    fprintf(FID, '\\caption{Simulations %d - %d}\n',j,j*8); 
    fprintf(FID, '\\end{figure}\n');
    fprintf(FID, '\\begin{figure}\\label{fig: SimulationMetropolisHasting%d}\n',j); 
end  
fprintf(FID, '\\begin{tabular}{cc} \n');
fileimg=['MetropolisExample' num2str(i) '.eps'];
fprintf(FID, '\\includegraphics[width=0.5\\textwidth]{ImaginiLatex/%s} &\n',fileimg);
fileimg=['MetropolisExample' num2str(i+1) '.eps'];
fprintf(FID, '\\includegraphics[width=0.5\\textwidth]{ImaginiLatex/%s} \\\\\n',fileimg);
%fileimg=['MetropolisExample' num2str(i+2) '.eps'];
%fprintf(FID, '\\includegraphics[width=0.5\\textwidth]{ImaginiLatex/%s}\\\\ \n',fileimg);
fprintf(FID, '\\textbf{Simulation %d} $\\theta_0=%8.2f$  $\\tau=%8.2f$  & ',i, theta(i,1), tau(i));
fprintf(FID, '\\textbf{Simulation %d} $\\theta_0=%8.2f$  $\\tau=%8.2f$\n',i+1,theta(i+1,1), tau(i+1));
%fprintf(FID, '\\textbf{Simulation %d} $\\theta_0=%8.2f$  $\\sigma=%8.2f$ \n',i+2, theta(i+2,1), sigma(i+2));
fprintf(FID, '\\end{tabular}\n');


end
fprintf(FID, '\\end{figure}\n');
fclose(FID);

movefile('*.eps','/home/giorgio/Desktop/MultiTracking_TESI/ImaginiLatex')
movefile('*.tex','/home/giorgio/Desktop/MultiTracking_TESI/tabelle')
