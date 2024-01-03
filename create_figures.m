for n=1151:1300
    filename=sprintf('fourth_variables%d.txt',n);
    if exist(filename,'file')
        data=dlmread(filename);
        
        a=1;
        b=length(data);
        c=1;
        
        name1=sprintf('x1_x2_%d.jpg',n);
        name2=sprintf('y1_y2_%d.jpg',n);
        
        figure('visible','off')
        plot(data(a:c:b,2),data(a:c:b,3),'.','MarkerSize',0.1)
        xlabel('x1')
        ylabel('x2')
        title(filename)
        saveas(gcf,name1)
        close
        
        figure('visible','off')
        plot(data(a:c:b,4),data(a:c:b,5),'.','MarkerSize',0.1)
        xlabel('y1')
        ylabel('y2')
        title(filename)
        saveas(gcf,name2)
        close
    else
        fprintf('File %s does not exist.\n', filename);
    end
end