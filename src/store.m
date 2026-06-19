[Soo    , So    ] = deal(So    , S    );
[Coo    , Co    ] = deal(Co    , C    );
[Moo    , Mo    ] = deal(Mo    , M    );
[Xoo    , Xo    ] = deal(Xo    , X    );
[rhooo  , rhoo  ] = deal(rhoo  , rho  );
[rhoxoo , rhoxo ] = deal(rhoxo  , rhox  );
[TRCmoo  , TRCmo  ] = deal(TRCmo  , TRCm  );
[TRCxoo  , TRCxo  ] = deal(TRCxo  , TRCx  );
[Too    , To    ] = deal(To    , T    );
[Tpoo   , Tpo   ] = deal(Tpo   , Tp   );
[chioo  , chio  ] = deal(chio  , chi   );

[dSdtoo    , dSdto    ] = deal(dSdto    , dSdt    );
[dCdtoo    , dCdto    ] = deal(dCdto    , dCdt    );
[dMdtoo    , dMdto    ] = deal(dMdto    , dMdt    );
[dXdtoo    , dXdto    ] = deal(dXdto    , dXdt    );
[drhodtoo  , drhodto  ] = deal(drhodto  , drhodt  );
[dTRCmdtoo  , dTRCmdto  ] = deal(dTRCmdto  , dTRCmdt  );
[dTRCxdtoo  , dTRCxdto  ] = deal(dTRCxdto  , dTRCxdt  );
[dTpdtoo   , dTpdto   ] = deal(dTpdto   , dTpdt   );
[dTdtoo    , dTdto    ] = deal(dTdto    , dTdt    );

[Gmo,Gxo] = deal(Gm,Gx);

dto = dt;
Pto = Pt;

% reset update history
specrad.S.est   = specrad.S.mean   + 0.*specrad.S.est;
specrad.C.est   = specrad.C.mean   + 0.*specrad.C.est;
specrad.PHS.est = specrad.PHS.mean + 0.*specrad.PHS.est;
specrad.TRCm.est = specrad.TRCm.mean + 0.*specrad.TRCm.est;
specrad.TRCx.est = specrad.TRCx.mean + 0.*specrad.TRCx.est;
specrad.MFD.est = specrad.MFD.mean + 0.*specrad.MFD.est;
specrad.CMP.est = specrad.CMP.mean + 0.*specrad.CMP.est;

GHST.S          = 0.*GHST.S;
GHST.C          = 0.*GHST.C;
GHST.PHS        = 0.*GHST.PHS;
GHST.TRCm       = 0.*GHST.TRCm;
GHST.TRCx       = 0.*GHST.TRCx;
GHST.MFD        = 0.*GHST.MFD;
GHST.CMP        = 0.*GHST.CMP;

FHST.S          = 0.*FHST.S;
FHST.C          = 0.*FHST.C;
FHST.PHS        = 0.*FHST.PHS;
FHST.TRCm       = 0.*FHST.TRCm;
FHST.TRCx       = 0.*FHST.TRCx;
FHST.MFD        = 0.*FHST.MFD;
FHST.CMP        = 0.*FHST.CMP;

% MFDCrr = 0.*MFDCrr;
% CMPCrr = 0.*CMPCrr;