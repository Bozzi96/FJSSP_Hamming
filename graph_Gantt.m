function graph_Gantt(sol, G_init, G_j, P, gamma, M0)
    [G,~, M_init, aux, aux_alt] = pre_processing_graph(G_init, P, M0);
    J = length(unique(G_j)); %jobs
    
    mySol = G(gamma  > 0.1,:);
    mySol_init = mySol;
    for i=1:size(mySol,1)
        for j=1:size(mySol,2)
            if(mySol(i,j) > M_init)
                [mySol_init(i,j), ~] = find(mySol(i,j) == aux);
            end
        end
    end
    
%     t= 1:sol.C+1;
%     start_time = sol.s;
%     occupancy = zeros(M, length(t),J);
%     utilz = zeros (M_init,1);
    
    col = maxdistcolor(J, @sRGB_to_OKLab); % From matlab file exchange
    for j=1:J
        current_row = mySol(j,mySol(j,:)~=0); % Remove zero elements
        startDates{j} =sol.s(j,current_row);
        endDates{j}=sol.c(j,current_row);
        varname(j) = "J_{" + num2str(j) + "}"; 
         i=1;
         for i=1:length(startDates{1,j})
             plot([startDates{1,j}(i),endDates{1,j}(i)],[mySol_init(j,i),mySol_init(j,i)],'b','Linewidth',10,'Color',col(j,:), DisplayName=varname(j))
            hold on
        end
    end
   legendUnq();  % From matlab file exchange
   legend('-DynamicLegend', 'NumColumns',5, Location='northeast')
    ylabel('Machine');
    xlabel('Time')
    ylim([0 M_init+1]);
end