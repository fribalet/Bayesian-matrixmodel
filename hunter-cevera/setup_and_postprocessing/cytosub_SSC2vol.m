function [ cellvol ] = cytosub_SSC2vol(SSCbu);

%empirically convert SSC in bead units (1 micron red) to cell volume in microliters

%cellvol = 10.^(-0.2028.*log10(SSCbu).^3 + 0.0445.*log10(SSCbu).^2 + 1.3519.*log10(SSCbu) + 0.9592); %12/22 from rob
%cellvol = 10.^(-0.1912.*log10(SSCbu).^3 + 0.0309.*log10(SSCbu).^2 + 1.3486.*log10(SSCbu) + 0.9616);  %FUDGE!!!
%cellvol = 10.^(-0.1968.*log10(SSCbu).^3 + 0.0477.*log10(SSCbu).^2 + 1.3435.*log10(SSCbu) + 0.9545);  %12/23 from rob
cellvol = 10.^(1.1088.*log10(SSCbu) + .9617); %same data as 12/13 from rob

%both of the polynomial fits blow up for SSC smaller than data points
%fplot('-0.1968.*x.^3 + 0.0477.*x.^2 + 1.3435.*x + 0.9545', [-1.5 1.5]) %this one can't work below 
%fplot('-0.2028.*x.^3 + 0.0445.*x.^2 + 1.3519.*x + 0.9592', [-1.5 1.5])

