[Soo    , So    ] = deal(So    , S    );
[Coo    , Co    ] = deal(Co    , C    );
[Cmoo   , Cmo   ] = deal(Cmo   , Cm   );
[Cxoo   , Cxo   ] = deal(Cxo   , Cx   );
[Moo    , Mo    ] = deal(Mo    , M    );
[Xoo    , Xo    ] = deal(Xo    , X    );
[rhooo  , rhoo  ] = deal(rhoo  , rho  );
[rhoxoo , rhoxo ] = deal(rhoxo , rhox );
[TRCoo  , TRCo  ] = deal(TRCo  , TRC  );
[TRCmoo , TRCmo ] = deal(TRCmo , TRCm );
[TRCxoo , TRCxo ] = deal(TRCxo , TRCx );
[Too    , To    ] = deal(To    , T    );
[Tpoo   , Tpo   ] = deal(Tpo   , Tp   );
[chioo  , chio  ] = deal(chio  , chi  );

[dSdtoo     , dSdto    ] = deal(dSdto    , dSdt    );
[dCdtoo     , dCdto    ] = deal(dCdto    , dCdt    );
[dCmdtoo    , dCmdto   ] = deal(dCmdto   , dCmdt   );
[dCxdtoo    , dCxdto   ] = deal(dCxdto   , dCxdt   );
[dMdtoo     , dMdto    ] = deal(dMdto    , dMdt    );
[dXdtoo     , dXdto    ] = deal(dXdto    , dXdt    );
[drhodtoo   , drhodto  ] = deal(drhodto  , drhodt  );
[dTRCdtoo   , dTRCdto  ] = deal(dTRCdto  , dTRCdt  );
[dTRCmdtoo  , dTRCmdto ] = deal(dTRCmdto , dTRCmdt );
[dTRCxdtoo  , dTRCxdto ] = deal(dTRCxdto , dTRCxdt );
[dTpdtoo    , dTpdto   ] = deal(dTpdto   , dTpdt   );
[dTdtoo     , dTdto    ] = deal(dTdto    , dTdt    );

[Gmo,Gxo] = deal(Gm,Gx);

dto = dt;
Pto = Pt;

% reset update history
specrad.S.est   = specrad.S.mean   + 0.*specrad.S.est;
specrad.C.est   = specrad.C.mean   + 0.*specrad.C.est;
specrad.Cm.est  = specrad.Cm.mean  + 0.*specrad.Cm.est;
specrad.Cx.est  = specrad.Cx.mean  + 0.*specrad.Cx.est;   
specrad.M.est = specrad.M.mean + 0.*specrad.M.est;
specrad.X.est = specrad.X.mean + 0.*specrad.X.est;
specrad.TRC.est  = specrad.TRC.mean  + 0.*specrad.TRC.est;
specrad.TRCm.est = specrad.TRCm.mean + 0.*specrad.TRCm.est;
specrad.TRCx.est = specrad.TRCx.mean + 0.*specrad.TRCx.est;
specrad.MFD.est = specrad.MFD.mean + 0.*specrad.MFD.est;
specrad.CMP.est = specrad.CMP.mean + 0.*specrad.CMP.est;

GHST.S          = 0.*GHST.S;
GHST.C          = 0.*GHST.C;
GHST.Cm         = 0.*GHST.Cm;   
GHST.Cx         = 0.*GHST.Cx;    
GHST.M          = 0.*GHST.M;
GHST.X          = 0.*GHST.X;
GHST.TRC        = 0.*GHST.TRC;
GHST.TRCm       = 0.*GHST.TRCm;
GHST.TRCx       = 0.*GHST.TRCx;
GHST.MFD        = 0.*GHST.MFD;
GHST.CMP        = 0.*GHST.CMP;

FHST.S          = 0.*FHST.S;
FHST.C          = 0.*FHST.C;
FHST.Cm         = 0.*FHST.Cm;
FHST.Cx         = 0.*FHST.Cx;
FHST.M          = 0.*FHST.M;
FHST.X          = 0.*FHST.X;
FHST.TRC        = 0.*FHST.TRC;
FHST.TRCm       = 0.*FHST.TRCm;
FHST.TRCx       = 0.*FHST.TRCx;
FHST.MFD        = 0.*FHST.MFD;
FHST.CMP        = 0.*FHST.CMP;

MFDCrr = 0.*MFDCrr;
CMPCrr = 0.*CMPCrr;