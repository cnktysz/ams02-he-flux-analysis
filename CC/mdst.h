
enum   EStat1 { BadTrig = 0x0001, BadTRD   = 0x0002, BadRTI   = 0x0004,
		InSAA   = 0x0008, FirstSec = 0x0010, LastSec  = 0x0020,
		IsInL1  = 0x0040, IsInL9   = 0x0080, IsInEcal = 0x0100,
		HasL1   = 0x0200, HasL2    = 0x0400, HasL9    = 0x0800,
		TrdInFS = 0x1000, TrdInTr  = 0x2000, FirstRTI = 0x4000,
		NoPTrk  = 0x8000,
		TrdInL1 = 0x10000, CutTk2nd = 0x20000 };
enum   EStat2 {};

struct Status { ULong64_t status;     // AMSEventR::fStatus
                ULong64_t ustat; };   // Estat1 | (Estat2 << 16)

struct Header { UInt_t run, event, ient, utime, phpat, jmpat;
                Int_t ntrdhit, ntrdseg, ntrhit, ntofcls, 
		      nanticls, nrichhit, necalhit;
                Int_t error; };
struct RTI    { Float_t theta, phi, zenith, cfi, lf, dl1, dl9, cf; };

struct Part   { Float_t mon, emom, mass, emass, beta, ebeta; };
struct Track  { Int_t   bith, bitx;  // Y,X bit pattern
                Int_t   bitr;        // Recovered bit pattern
                Float_t rgt[4], csqx[4], csqy[4]; // 0:In 1:L1 2:L9 3:FS
                Float_t qin, ql1, ql9;
                Float_t p0x, p0y, theta, phi; // At Z=0
                Float_t coox[3], cooy[3];     // 0:L1 1:L2 2:L9
                Float_t resx[2], resy[2];     // 0:L1 1:L9 excl.residual
                Float_t hity[4];              // L1 PG,MD, L9 PG,MD
                Int_t   qsta[2];              // GetQStatus() for 0:L1 and 1:L9
                Float_t qrms;                 // GetInnerQ_all.RMS
                Float_t rgt2[3], csx2, csy2;  // Chikanian,Up,Low
                Int_t   itp, addl;            // Max span, AddLost
                Int_t   rec;                  // GetRecType()
              };
struct TrHit  { Int_t   nhit[6];       // 0:L1xy 1:L1y 2:L9xy 3:L9y 4:Ixy 5:Iy
                Float_t qmax[6];       // 0:L1xy 1:L1y 2:L9xy 3:L9y 4:Ixy 5:Iy
                Float_t qli [7];       // TrTrack::Qlayer (Inner)
                Float_t dx[2], dy[4];  // QmaxHit-Track 0:L1 1:L9 +2:OnlyY
                Int_t   qsi [7];       // GetQStatus()
                Int_t   nh10[7];       // IsTrackPickingUpNoise(10, 0, 0, nh10)
                Float_t rata[7];       // GetTrackerRawSignalRatio(j+2);
                Float_t tkfd[7];       // GetTkFeetDist(j+2)
                Float_t etrx[9], etry[9];
              };
struct TrCls  { Int_t   tkml[9];
                Float_t xcog[9], ycog[9], mdist[9];
                Int_t   flag[9];
                Float_t eta [9], sig [9];
              };

struct Beta   { Int_t pattern; Float_t beta; };
struct BetaH  { Int_t pattern, pbit;                     // Default BetaH
                Float_t beta, chi2t, chi2c, q, ql[4];
                Int_t clsn[4], z[3]; Float_t p[3];
                Int_t type; Float_t qup, qlow;
              };
struct BetaHs { Int_t pattern, pbit;                     // TrTrack independent
                Float_t beta, chi2t, chi2c, q, ql[4]; };
struct TofHit { Int_t   nhit[4];
                Float_t qsum[4], qmax[4];
                Float_t tmin[4], tmax[4]; };

struct Trd    { Float_t coo[3], theta, phi, q;
                Float_t dtrk[3];       // 0:dx 1:dy 2:angle(deg) with TrTrack
                Float_t ql1m, dl1m[2]; // Charge, dx, dy of the closest L1 hit
                Int_t nsegt;
              };
struct TrdK   { Int_t nhits; Float_t llr[3], q, llre[3];
                Int_t nhitt; Float_t llrt[3];
                Int_t nhitd; Float_t llrd[3];
              };
struct TrdG   { Int_t nfit;  Float_t gamma, lkh; };
struct TrdV   {
                Int_t nvtx, ntrk, nhit;
                Float_t chi2, x, y, z;
};
struct TrdHit { Int_t  layer[25];
                Float_t plen[25], amp[25], plmc[25], trfr[25];
	      };
struct Rich   { Int_t status, nhits[3];
                Float_t beta, q, dist, prob, coo[3], npe[3];
                Int_t tile, pmts;
                Float_t nhpe[30];
              };
struct Ecal   { Float_t cog[3], dir[3], s13r, bdt, energy[3]; // 0:D 1:C 2:E
                Float_t dtrk[3];       // 0:dx 1:dy 2:angle(deg) with TrTrack
                Float_t ql9m, dl9m[2]; // Charge, dx, dy of the closest L9 hit
                Int_t catl;
                Float_t enew[2];       // GetCorrectedEnergy
                Float_t dz[2];
              };

struct EcalH  { Int_t apex;
                Float_t csq, rrec, smax, rnt, emip, lmip;
                Float_t edep[10];
              };
struct Estim  { Float_t llre[3], ebdt[2], trklh;
                Int_t catl, clsn[4]; };

struct Mass   { Float_t b1, b2, b3; // TrMass::GetBeta with flag= 1, 10, 100
                Float_t npk, m;     // TrMass::GetNpick, mass(Chikanian,BetaH)
                Float_t v[16], p[16], mql; // TrMass::GetMQL
                Int_t nh;           // Number of inner tracker hits
                Float_t ll[2];      // TrMass::GetLL with (1,1),(1,2)
                Float_t bl2, ru;    // TrMass::GetBL2, TrMass::Rndm
                Int_t zl[4];
                Float_t pz[4], betas;
              };

struct Refit  { Float_t rrgt[4], rcsq[4], resy[4];
                Int_t   slot[2];
                Float_t bcor; };

struct Vertex { Float_t theta, phi, zenith;
                Float_t rgt[2], csqx[2], csqy[2], coo[7], ang[4]; };

struct MCinfo { Float_t coo[3], dir[3], rgt, mcin[2];
                Int_t   part[5]; // 0:L1 1:L2 2:L56 3:L78 4:L9
                Float_t trgt[5]; // 0:L1 1:L2 2:L56 3:L78 4:L9
                Int_t   patmc, npmc;
                Float_t trcx[9], trcy[9], trcz[9];
                Float_t dir1[3], dir9[3];
                Float_t mci2[7], cfw[2];
                Int_t rseed[2];
              };

class DST {
public:
  Status fStatus; Header fHeader; RTI    fRTI;
  Part   fPart;   Track  fTrack;  TrHit  fTrHit;  TrCls  fTrCls;
  Beta   fBeta;   BetaH  fBetaH;  BetaHs fBetaHs; TofHit fTofHit;
  Trd    fTrd;    TrdK   fTrdK;   TrdG   fTrdG;   TrdHit fTrdHit; TrdV fTrdV;
  Rich   fRich;   Ecal   fEcal;   EcalH  fEcalH;  Vertex fVertex; 
  MCinfo fMCinfo; Estim  fEstim;  Refit  fRefit;  Mass   fMass;

  void Clear() { Int_t *ptr = (Int_t *)this;
                 Int_t size = sizeof(DST)/sizeof(Int_t);
		 for (Int_t i = 0; i < size; i++) ptr[i] = 0;
               }

  void Branch(TTree *tree, Int_t mode = 1) {
    tree->Branch("status", &fStatus, "status[2]/i:ustat[2]/i");
    tree->Branch("header", &fHeader, "run/i:event/i:ient/i:"
		                     "utime/i:phpat/i:jmpat/i:"
		                     "ntrdh/I:ntrds/I:ntrh/I:ntofc/I:"
		                     "nantic/I:nrichh/I:necalh/I:error/I");
    tree->Branch("rti",    &fRTI,    "theta/F:phi/F:zenith/F:"
		                     "cfi/F:lf/F:dl1/F:dl9/F:cf/F");
    tree->Branch("part",   &fPart,   "mom/F:emom/F:mass/F:emass/F:"
		                     "beta/F:ebeta/F");
    tree->Branch("track",  &fTrack,  "bith/I:bitx/I:bitr/I:"
		                     "rgt[4]/F:csqx[4]/F:csqy[4]/F:"
		                     "qin/F:ql1/F:ql9/F:"
		                     "p0x/F:p0y/F:theta/F:phi/F:"
		                     "coox[3]/F:cooy[3]/F:"
		                     "resx[2]/F:resy[2]/F:hity[4]/F:"
		                     "qsta[2]/I:qrms/F:"
		                     "rgt2[3]/F:csx2/F:csy2/F:itp/I:addl/I:"
                                     "rect/I");
//  if (mode == 1)
    tree->Branch("trhit",  &fTrHit,  "nhit[6]/I:qmax[6]/F:qli[7]/F:"
		                     "hdx[2]/F:hdy[4]/F:qsi[7]/I:"
		                     "nh10[7]/I:rata[7]/F:tkfd[7]/F:etrx[9]/F:etry[9]/F");
    if (mode != 2)
    tree->Branch("trcls",  &fTrCls,  "tkml[9]/I:"
		                     "xcog[9]/F:ycog[9]/F:mdist[9]/F:"
		                     "flag[9]/I:eta[9]/F:sig[9]/F");
    tree->Branch("beta",   &fBeta,   "pattern/I:beta/F");
    tree->Branch("betah",  &fBetaH,  "pattern/I:pbit/I:beta/F:"
		                     "chi2t/F:chi2c/F:q/F:ql[4]/F:"
		                     "clsn[4]/I:z[3]/I:p[3]/F:"
		                     "type/I:qup/F:qlow/F");
    tree->Branch("betas",  &fBetaHs, "pattern/I:pbit/I:beta/F:"
		                     "chi2t/F:chi2c/F:q/F:ql[4]/F");
  //if (mode == 1)
    tree->Branch("tofhit", &fTofHit, "nhit[4]/I:qsum[4]/F:qmax[4]/F:"
		                     "tmin[4]/F:tmax[4]/F");

    tree->Branch("trd",    &fTrd,    "coo[3]/F:theta/F:phi/F:q/F:"
		                     "dtrk[3]/F:ql1m/F:dl1m[2]/F:ntseg/I");
    tree->Branch("trdk",   &fTrdK,   "nhits/I:llr[3]/F:q/F:llre[3]/F:"
		                     "nhitt/I:llrt[3]/F:"
		                     "nhitd/I:llrd[3]/F");
    tree->Branch("trdg",   &fTrdG,   "nfit/I:gamma/F:lkh/F");
    if (mode == 1)
    tree->Branch("trdhit", &fTrdHit, "layer[25]/I:plen[25]/F:amp[25]/F:"
		                                 "plmc[25]/F:trfr[25]/F");
    tree->Branch("trdv",   &fTrdV, "nvtx/I:nvtrk/I:nvhit/I:chi2/F:"
                                    "trdvx/F:trdvy/F:trdvz/F");

    tree->Branch("rich",   &fRich,   "rstat/I:rnhit[3]/I:"
		                     "beta/F:q/F:dist/F:prob/F:rcoo[3]/F:"
		                     "npe[3]/F:tile/I:pmts/I:nhpe[30]/F");

    tree->Branch("ecal",   &fEcal,   "cog[3]/F:dir[3]/F:s13r/F:bdt/F:"
		                     "energy[3]/F:dtrk[3]/F:ql1m/F:dl1m[2]/F:"
		                     "catl/I:enew[2]/F:dz[2]/F");
    tree->Branch("ecalh",  &fEcalH,  "apex/I:csq/F:rrec/F:smax/F:"
		                     "rnt/F:emip/F:lmip/F:"
		                     "edep[10]/F");

    tree->Branch("vertex", &fVertex, "vth/F:vph/F:vzen/F:vrgt[2]/F:"
		                     "vcsqx[2]/F:vcsqy[2]/F:"
		                     "vcoo[7]/F:vang[4]/F");

    tree->Branch("mcinfo", &fMCinfo, "mcoo[3]/F:mdir[3]/F:mrgt/F:mcin[2]/F:"
		                     "part[5]/I:trgt[5]/F:patmc/I:npmc/I:"
		                     "trcx[9]/F:trcy[9]/F:trcz[9]/F:"
		                     "dir1[3]/F:dir9[3]/F:mci2[7]/F:"
		                     "cfw[2]/F:rseed[2]/I");

    tree->Branch("mass",   &fMass,   "b1/F:b2/F:b3/F:npk/F:m/F:"
		                     "v[16]/F:p[16]/F:mql/F:"
		                     "nh/I:ll[2]/F:bl2/F:ru/F:zl[4]/I:pz[4]/F:betas/F");
/*
    tree->Branch("estim",  &fEstim,  "llre[3]/F:ebdt[2]/F:trklh/F:"
		                     "catl/I:clsn[4]/I");
    tree->Branch("refit",  &fRefit,  "rrgt[4]/F:rcsq[4]/F:resy[4]/F:"
		                     "slot[2]/I:bcor/F");
*/
  }

  void SetAddress(TTree *tree) {
    tree->SetBranchAddress("status", &fStatus);
    tree->SetBranchAddress("header", &fHeader);
    tree->SetBranchAddress("rti",    &fRTI);
if (tree->GetBranch       ("part"))
    tree->SetBranchAddress("part",   &fPart);
    tree->SetBranchAddress("track",  &fTrack);
    tree->SetBranchAddress("trhit",  &fTrHit);
    //tree->SetBranchAddress("trcls",  &fTrCls);
    tree->SetBranchAddress("beta",   &fBeta);
    tree->SetBranchAddress("betah",  &fBetaH);
if (tree->GetBranch       ("betas"))
    tree->SetBranchAddress("betas",  &fBetaHs);
    //tree->SetBranchAddress("tofhit", &fTofHit);
    tree->SetBranchAddress("trd",    &fTrd);
    tree->SetBranchAddress("trdk",   &fTrdK);
if (tree->GetBranch       ("trdg"))
    tree->SetBranchAddress("trdg",   &fTrdG);
if (tree->GetBranch       ("trdv"))
    tree->SetBranchAddress("trdv",   &fTrdV);
  //tree->SetBranchAddress("trdhit", &fTrdHit);
    tree->SetBranchAddress("rich",   &fRich);
    tree->SetBranchAddress("ecal",   &fEcal);
  //tree->SetBranchAddress("ecalh",  &fEcalH);
if (tree->GetBranch       ("vertex"))
    tree->SetBranchAddress("vertex", &fVertex);
    tree->SetBranchAddress("mcinfo", &fMCinfo);
    //tree->SetBranchAddress("mass",   &fMass);
  }
};
