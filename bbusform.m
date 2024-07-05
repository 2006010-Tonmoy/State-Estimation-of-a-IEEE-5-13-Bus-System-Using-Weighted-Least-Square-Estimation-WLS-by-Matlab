% Line Data for B-Bus (Shunt Admittance)Formation.

function bbus = bbusform(num)     % Returns B-bus..

linedata = linedatas(num);
fb = linedata(:,1);
tb = linedata(:,2);
b = linedata(:,5);
nbus = num;    % no. of buses...
nbranch = length(fb);           % no. of branches...
bbus = zeros(nbus,nbus);

 for k=1:nbranch
     bbus(fb(k),tb(k)) = b(k);   %% Only off diagonal are main concerns
     bbus(tb(k),fb(k)) = bbus(fb(k),tb(k));
 end