function map=plotContour(x,y,F)

if all(F(:)<=1)
    red = [192 0 0]/255;
    map = [red; 1 1 1];
    k = [0; 1];
    map = interp1(k,map,linspace(0,1,10));
    map(end,:) = [1 1 1];
    contourf(x,y,F,linspace(0,1,11),'LineColor','none'); hold on;
    caxis([0 1]);
    colorbar;
    colormap(map);
    hold off;
else
 [a,~] = contour(x,y,F,[0 1]);
    a = a(:,2:end);
    for i = 2:length(a)-1
        if (abs(a(1,i)-a(1,i+1))/(a(1,i)+a(1,i+1))> 1 || abs(a(1,i)-a(1,i-1))/(a(1,i)+a(1,i-1))> 1)
            a(:,i)= NaN;
        end
    end
    for i = 2:length(a)-1
        if (abs(a(2,i)-a(2,i+1))/(a(2,i)+a(2,i+1))>1 || abs(a(2,i)-a(2,i-1))/(a(2,i)+a(2,i-1))> 1)
            a(:,i)= NaN;
        end
    end
    idx = find(a(1,:)>x(end) | a(1,:)<x(1));
    a(:,idx)=[];
    idx = find(a(2,:)>y(end) | a(2,:)<y(1));
    a(:,idx)=[];
    blue = [0 82 177]/255;
    red = [192 0 0]/255;
    int = 0:0.1:(max(max(F))+0.1);
    map = [red; 1 1 1];
    map1 = interp1([0;1],map,linspace(0,1,10));
    map = [1 1 1; blue];
    map2 = interp1([0;1],map,linspace(0,1,length(int)-10));
    map2(1,:) = [];
    map = [map1;map2];
    contourf(x,y,F,int,'LineColor','none'); hold on;
    colormap(map);
    plot(a(1,:),a(2,:),'Color','k','Linewidth',1.5,'Linestyle',':');
    colorbar;
    caxis([0 int(end)]);
    hold off;
    
end


