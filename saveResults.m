
%{
    This script saves results, plots and mystery , output from SDR_receiver.mlx or
    SDR_receiver.m 
    either codes should run first.
%}


path= [pwd '\results\' messageName '\'];
if ~exist(path, 'dir')
   mkdir(path)
end
saveas(figure(1),[path 'IF_BB_signals.png']);
saveas(figure(2),[path 'PLL.png']);
saveas(figure(3),[path 'clock_rec.png']);
saveas(figure(4),[path 'correlation_matching.png']);
saveas(figure(5),[path 'soft_decisions.png']);

fileID = fopen([path messageName '.txt'],'w');
fprintf(fileID,reconstructed_message);