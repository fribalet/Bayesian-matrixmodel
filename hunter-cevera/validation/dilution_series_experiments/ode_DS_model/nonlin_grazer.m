%Nonlin investigations of ds exp - if grazer numbers increased...
function dy = nonlin_grazer(t,y)
dy=zeros(2,1);
m=0;
psi=0.0;%0.2 %units of zooplankton per phytoplankton
k=0.7;

%Type I ingestion response:
Pstar=500;
cmax=.01;
imax=2.5;

if y(1) <= Pstar
    dy(1)= y(1)*k - y(2)*cmax*y(1) ;
    dy(2)=y(2)*psi*cmax*y(1) -m*y(2);
else
    dy(1)= y(1)*k - y(2)*imax ;
    dy(2)=y(2)*psi*imax -m*y(2);
end
        

end

 






