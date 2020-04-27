function xy_meas = rtoxy(r_sim1, r_sim2, r_sim3, T)

    x1 = [-0.05, 1.5];
    x2 = [2, 3];
    x3 = [2.7, 0.05];

    xy_meas = NaN(6,length(r_sim1));
    
    xy_meas(1,1) = 0.5*((x3(2)-x2(2))*(r_sim1(1)^2-r_sim2(1)^2-x1(1)^2+x2(1)^2-x1(2)^2+x2(2)^2)-(x2(2)-x1(2))*(r_sim2(1)^2-r_sim3(1)^2-x2(1)^2+x3(1)^2-x2(2)^2+x3(2)^2));
    xy_meas(1,1) = xy_meas(1,1)/((x3(2)-x2(2))*(x2(1)-x1(1))-(x2(2)-x1(2))*(x3(1)-x2(1)));
    xy_meas(2,1) = 0;xy_meas(3,1) = 0;
    
    xy_meas(4,1) = 0.5*((x3(1)-x2(1))*(r_sim1(1)^2-r_sim2(1)^2-x1(1)^2+x2(1)^2-x1(2)^2+x2(2)^2)-(x2(1)-x1(1))*(r_sim2(1)^2-r_sim3(1)^2-x2(1)^2+x3(1)^2-x2(2)^2+x3(2)^2));
    xy_meas(4,1) = xy_meas(4,1)/((x3(1)-x2(1))*(x2(2)-x1(2))-(x2(1)-x1(1))*(x3(2)-x2(2)));
    xy_meas(5,1) = 0;xy_meas(6,1) = 0;
    
    for j = 2:1:length(r_sim1)
        xy_meas(1,j) = 0.5*((x3(2)-x2(2))*(r_sim1(j)^2-r_sim2(j)^2-x1(1)^2+x2(1)^2-x1(2)^2+x2(2)^2)-(x2(2)-x1(2))*(r_sim2(j)^2-r_sim3(j)^2-x2(1)^2+x3(1)^2-x2(2)^2+x3(2)^2));
        xy_meas(1,j) = xy_meas(1,j)/((x3(2)-x2(2))*(x2(1)-x1(1))-(x2(2)-x1(2))*(x3(1)-x2(1)));
        xy_meas(2,j) = (xy_meas(1,j) - xy_meas(1,j-1))/T;
        xy_meas(3,j) = (xy_meas(2,j) - xy_meas(2,j-1))/T;
        
        xy_meas(4,j) = 0.5*((x3(1)-x2(1))*(r_sim1(j)^2-r_sim2(j)^2-x1(1)^2+x2(1)^2-x1(2)^2+x2(2)^2)-(x2(1)-x1(1))*(r_sim2(j)^2-r_sim3(j)^2-x2(1)^2+x3(1)^2-x2(2)^2+x3(2)^2));
        xy_meas(4,j) = xy_meas(4,j)/((x3(1)-x2(1))*(x2(2)-x1(2))-(x2(1)-x1(1))*(x3(2)-x2(2)));
        xy_meas(5,j) = (xy_meas(4,j) - xy_meas(4,j-1))/T;
        xy_meas(6,j) = (xy_meas(5,j) - xy_meas(5,j-1))/T;
    end

end