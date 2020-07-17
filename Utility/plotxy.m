function  plotxy( variable,x,y1,y2)
%PLOTXY Overlapped x-y plot to make comparison

titlename=variable;
    figure('Name',titlename,'NumberTitle','off');
    plot(x,y1,x,y2,'--');
    legend('NSMOOM','ABAQUS','Location','northeast')
    xlabel('Distance along the path measured from the first node')
    ylabel(variable)
end

