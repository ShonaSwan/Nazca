[Soo    , So    ] = deal(So    , S    );
[Cmoo   , Cmo   ] = deal(Cmo   , Cm   );
[Cxoo   , Cxo   ] = deal(Cxo   , Cx   );
[rhooo  , rhoo  ] = deal(rhoo  , rho  );
[rhoxoo , rhoxo ] = deal(rhoxo , rhox );
[Too    , To    ] = deal(To    , T    );
[Tpoo   , Tpo   ] = deal(Tpo   , Tp   );
[chioo  , chio  ] = deal(chio  , chi  );

[dSdtoo     , dSdto    ] = deal(dSdto    , dSdt    );
[dCmdtoo    , dCmdto   ] = deal(dCmdto   , dCmdt   );
[dCxdtoo    , dCxdto   ] = deal(dCxdto   , dCxdt   );
[drhodtoo   , drhodto  ] = deal(drhodto  , drhodt  );

[Gmo,Gxo] = deal(Gm,Gx);

dto = dt;
Pto = Pt;

% reset update history
specrad.S.est   = specrad.S.mean   + 0.*specrad.S.est;
specrad.Cm.est  = specrad.Cm.mean  + 0.*specrad.Cm.est;
specrad.Cx.est  = specrad.Cx.mean  + 0.*specrad.Cx.est;   
specrad.MFD.est = specrad.MFD.mean + 0.*specrad.MFD.est;
specrad.CMP.est = specrad.CMP.mean + 0.*specrad.CMP.est;

GHST.S          = 0.*GHST.S;
GHST.Cm         = 0.*GHST.Cm;   
GHST.Cx         = 0.*GHST.Cx;    
GHST.MFD        = 0.*GHST.MFD;
GHST.CMP        = 0.*GHST.CMP;

FHST.S          = 0.*FHST.S;
FHST.Cm         = 0.*FHST.Cm;
FHST.Cx         = 0.*FHST.Cx;
FHST.MFD        = 0.*FHST.MFD;
FHST.CMP        = 0.*FHST.CMP;