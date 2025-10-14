function [fig,sub] = create_figure(num_groups)
% Figura de referencia donde poner las imágenes. Creo las referencias a los
% ejes para poder dibujar luego
%-----------------------------------
% Superior         Frontal
% Lateral Derecho  Lateral Izquierdo
%------------------------------------

close(figure(50))
sub = zeros(12,1);
fig = figure(50);
set(figure(50),'Color', [ 1 1 1 ],'units','normalized','outerposition',[0 0 1 1]);


switch num_groups
    case 1 % Subplot axis for 2 groups comparations
        for i = 1:2
            sub(i) = subplot(2,2,i,'Position',[(0.3 + (i-1)*0.25) 1.5/3 3/20 8/20]);
            sub(i+2) = subplot(2,2,i+2,'Position', [(0.3 + (i-1)*0.25) 0.2/3 3/20 8/20]);
        end
    case 3 % Subplot axis for 3 groups comparations
        
        sub(1) = subplot(2,6,1,'Position',[(1-1)/6 1.7/3 3/20 8/20]);
        sub(2) = subplot(2,6,2,'Position',[(2-1)/6 1.7/3 3/20 8/20]);
        sub(3) = subplot(2,6,7,'Position', [(1-1)/6 0.2/3 3/20 8/20]);
        sub(4) = subplot(2,6,8,'Position', [(2-1)/6 0.2/3 3/20 8/20]);
        
        sub(5) = subplot(2,6,3,'Position',[(3-1)/6 1.7/3 3/20 8/20]);
        sub(6) = subplot(2,6,4,'Position',[(4-1)/6 1.7/3 3/20 8/20]);
        sub(7) = subplot(2,6,9,'Position', [(3-1)/6 0.2/3 3/20 8/20]);
        sub(8) = subplot(2,6,10,'Position', [(4-1)/6 0.2/3 3/20 8/20]);
        
        sub(12) = subplot(2,6,12,'Position', [(6-1)/6 0.2/3 3/20 8/20]);
        sub(10) = subplot(2,6,6,'Position',[(6-1)/6 1.7/3 3/20 8/20]);
        sub(9) = subplot(2,6,5,'Position',[(5-1)/6 1.7/3 3/20 8/20]);
        sub(11) = subplot(2,6,11,'Position', [(5-1)/6 0.2/3 3/20 8/20]);
        
        annotation('line', [1.95/6 1.95/6] , [0.2 0.9])
        annotation('line', [3.95/6 3.95/6] , [0.2 0.9])

        
end

