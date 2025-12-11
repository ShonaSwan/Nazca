%get vertical gradient of bulk SiO2
SiO2ByMgO = c_oxd(:,:,cal.Si)./c_oxd(:,:,cal.Mg);

[~,gradSiO2ByMgO] = gradient(SiO2ByMgO);
gradSiO2ByMgO(gradSiO2ByMgO > 0) = 0;

indLAB    = (T-273.15) < 1200;
[LAB_depth,LAB_iz] = max(ZZ.*indLAB,[],1);

for ix = 1:Nx
gradSiO2ByMgO(ceil(LAB_iz(ix)):end,ix) = 0;%gradSiO2ByMgO(ceil(LAB_iz(ix)/2),ix);
end

[~,MOHO_iz] = max(max(1e-1,-gradSiO2ByMgO),[],1);

% retrieve moho location
MOHO_depth = Zc(MOHO_iz);
