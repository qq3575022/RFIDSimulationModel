function [PP, VV, AA, AStates] = groundtruth1D(td)  

p1 = 0;
p2 = 0.64;

ac = 7.35;
vc = 1.2;

ta = vc/ac;
tv = (p2-p1)/vc - vc/ac;

t1 = 0.157;
t2 = t1+ta+tv+ta;

PP = NaN(1,length(td));
VV = NaN(1,length(td));
AA = NaN(1,length(td));

for n = 1:length(td)
    t = td(n);
    if t < t1
        p = p1;
        v = 0;
        a = 0;
    elseif t < t1 + ta
        p = p1+0.5*ac*(t-t1)^2;
        v = ac*(t-t1);
        a = ac;
    elseif t < t2-ta
        p = 0.5*(p1+p2)+vc*(t-0.5*(t1+t2));
        v = vc;
        a = 0;
    elseif t < t2
        p = p2-0.5*ac*(t2-t)^2;
        v = -ac*(t-t2);
        a = -ac;
    else
        p = p2;
        v = 0;
        a = 0;
    end
    PP(n) = p;
    VV(n) = v;
    AA(n) = a;
end
PP = PP + 0.99;
AStates = zeros(3,226);
AStates(1,:) = PP';
AStates(2,:) = VV';
AStates(3,:) = AA';

end
    