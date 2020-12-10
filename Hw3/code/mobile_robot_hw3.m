close all
ds=0;db=0;
CovP = zeros(3,3);

t = 5;



now = 0;
x=[0;0;0];
plotx=[];

ax=subplot(1,1,1);
% axis([0,1000,-1000,1000]);
while now<t
    cla(ax);
%   dx = [now;0;0]
    x = [500*sin(now/t*2*pi/5);-500*cos(now/t*2*pi/5);0];
    
    b = pi/2*now/t;
    ds = 8;%norm(dx);
    
    CovP = UpdateCovP(CovP,x(3),ds,db,b);
    [v,vd] = pcacov(CovP);%eig(CovP);
    for i = 1:3
    v(:,i) = v(:,i)/norm(v(:,i));
    end
    now = now + 0.01;
    plotx=[plotx;x'];
    
    hold on;grid on;
    sc=3;
    p1 = x - v(:,1) * vd(1)*sc/15 ;
    p2 = x - v(:,2) * vd(2)*sc ; 
    p3 = x + v(:,1) * vd(1)*sc/15 ; 
    p4 = x + v(:,2) * vd(2)*sc ; 
    line([p1(1), p3(1)],[p1(2), p3(2)],'color', 'b');
    line([p2(1), p4(1)],[p2(2), p4(2)],'color', 'r');
    xlim([0,2000])
    ylim([-1000,1000])
    plot(plotx(:,1),plotx(:,2));
    
    pause(0.01);

end


function CovX = GetCovX(ds,b,ks,kb)
    CovX = [ks*abs(ds) 0;0 kb];
end
function [Fp,Fsb] = CalculateJocbian(a,ds,db,b)
    L2 = 500;
    t = ds*sin(b+db/2);
    ds
    t1 = a-ds*(cos(b+db/2)/2/L2);
    Fp = [1 0 -t*sin(t1);0 1 t*cos(t1);0 0 1];
    
    t = a-ds*cos(b)/2/L2;
    Fsb = [cos(t)*sin(b)+ds*sin(t)*cos(b)*sin(b)/2/L2, ds*cos(t)*cos(b) - ds^2*sin(t)*sin(b)^2/2/L2 ;
           sin(t)*sin(b)-ds*cos(t)*cos(b)*sin(b)/2/L2, ds*sin(t)*cos(b) + ds^2*cos(t)*sin(b)^2/2/L2 ;
           -cos(b)/L2                                , ds*sin(b)/L2];
end

function CovP = UpdateCovP(CovP,a,ds,db,b)
    Fp=[];Fsb=[];
    [Fp,Fsb] = CalculateJocbian(a,ds,db,b);
    CovX = GetCovX(ds,b,0.02,0.02);
    CovP = Fp*CovP*Fp' + Fsb*CovX*Fsb';
end



