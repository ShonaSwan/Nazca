%get vertical gradient of bulk SiO2
SiO2ByMgO = majx_oxd(:,:,cal.Si)./majx_oxd(:,:,cal.Mg);

indLAB    = (T-273.15) < 1200;
[LAB_depth,LAB_iz] = max(ZZ.*indLAB,[],1);

for ix = 1:Nx
    SiO2ByMgO(LAB_iz(ix):end,ix) = 1;
end

indMOHO              = SiO2ByMgO > 2;
[MOHO_depth,MOHO_iz] = max(ZZ.*indMOHO,[],1);
