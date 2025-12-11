[Soo    , So    ] = deal(So    , S    );
[Coo    , Co    ] = deal(Co    , C    );
[Moo    , Mo    ] = deal(Mo    , M    );
[Xoo    , Xo    ] = deal(Xo    , X    );
[rhooo  , rhoo  ] = deal(rhoo  , rho  );
[TRCoo  , TRCo  ] = deal(TRCo  , TRC  );
[Too    , To    ] = deal(To    , T    );
[Tpoo   , Tpo   ] = deal(Tpo   , Tp   );
[chioo  , chio  ] = deal(chio  , chi   );

[dSdtoo    , dSdto    ] = deal(dSdto    , dSdt    );
[dCdtoo    , dCdto    ] = deal(dCdto    , dCdt    );
[dMdtoo    , dMdto    ] = deal(dMdto    , dMdt    );
[dXdtoo    , dXdto    ] = deal(dXdto    , dXdt    );
[drhodtoo  , drhodto  ] = deal(drhodto  , drhodt  );
[dTRCdtoo  , dTRCdto  ] = deal(dTRCdto  , dTRCdt  );
[dTpdtoo   , dTpdto   ] = deal(dTpdto   , dTpdt   );
[dTdtoo    , dTdto    ] = deal(dTdto    , dTdt    );

[Gmo,Gxo] = deal(Gm,Gx);

dto = dt;
Pto = Pt;

% reset update history
specrad.S.est   = specrad.S.mean   + 0.*specrad.S.est;
specrad.C.est   = specrad.C.mean   + 0.*specrad.C.est;
specrad.PHS.est = specrad.PHS.mean + 0.*specrad.PHS.est;
specrad.TRC.est = specrad.TRC.mean + 0.*specrad.TRC.est;
specrad.MFD.est = specrad.MFD.mean + 0.*specrad.MFD.est;
specrad.CMP.est = specrad.CMP.mean + 0.*specrad.CMP.est;

GHST.S          = 0.*GHST.S;
GHST.C          = 0.*GHST.C;
GHST.PHS        = 0.*GHST.PHS;
GHST.TRC        = 0.*GHST.TRC;
GHST.MFD        = 0.*GHST.MFD;
GHST.CMP        = 0.*GHST.CMP;

FHST.S          = 0.*FHST.S;
FHST.C          = 0.*FHST.C;
FHST.PHS        = 0.*FHST.PHS;
FHST.TRC        = 0.*FHST.TRC;
FHST.MFD        = 0.*FHST.MFD;
FHST.CMP        = 0.*FHST.CMP;