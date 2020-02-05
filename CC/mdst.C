
#include "amschain.h"
#include "root.h"

#include "TrExtAlignDB.h"
#include "Tofrec02_ihep.h"
#include "TrdKCluster.h"
#include "GammaFit.h"
#include "TrCharge.h"
#include "TrRecon.h"
#include "TrMass.h"
#include "EcalH.h"
#include "bcorr.h"
#include "tkdcards.h"


#include "TH2.h"
#include "TF1.h"
#include "TMath.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStopwatch.h"


//#include "tbpos.C"

#include <cstdlib>

#include "mdst.h"

class ParticleP : public ParticleR {
public:
  void setTrack(Int_t itr) { fTrTrack = itr; }
  void setBetaH(Int_t ibh) { fBetaH   = ibh; }
};

class MdSel : public AMSEventR {

public:
  static Int_t   RUN_STD;
  static Int_t   RUN_TTCS1;
  static Int_t   RUN_TTCS2;
  static Int_t   fNtot;
  static TString fFname;
  static TChain *fChain;

  Int_t fNrun;
  Int_t fNevt;
  Int_t fNfil;
  Int_t fIntv;
  Int_t fBrun;
  Int_t fBtim;
  Int_t fBnew;

  Int_t fLevt;
  Int_t fChrg;

  TStopwatch fTimer;

  TFile *fFile;
  TTree *fTree;
  DST    fDst;

  AMSSetupR::RTI fRTI;
  Int_t    fBadRun;
  Int_t    fBadRTI;
  Bool_t   fInSAA;
  UInt_t   fTRTI[2];
  Double_t fPGMD[2];
  Double_t fRcut;
  Double_t fLtf;

  TrTrackR *fTrk;
  Double_t  fBeta;
//EcalHR   *fEcalH;

public:
  MdSel() {
    fNrun = fNevt = fNfil = 0;
    fIntv = 1000; //0;

    fFile = 0;
    fTree = 0;
    fBrun = 0;
    fBtim = 0;
    fBnew = 0;
    fLevt = 0;

    fTRTI[0] = fTRTI[1] = 0;
    fPGMD[0] = fPGMD[1] = 0;
    fBadRun  = fBadRTI  = 0;
    fInSAA   = kFALSE;

  //fEcalH = new EcalHR;

  }
 ~MdSel() { 
   if (fFile) fFile->Write();
   delete fFile;
 }

  bool Notify() { GetBranch(_Tree); fBnew = 1; return true; }

  void UBegin();
  bool UProcessStatus(unsigned long long status);
  bool UProcessCut ();
  void UProcessFill();
  void UTerminate();

  void FillMC   ();
  void FillVtx  ();
  void FillPart ();
  void FillTrack();
  void FillTrHit();
  void FillBeta ();
  void FillTrd  ();
  void FillRich ();
  void FillEcal ();
  void FillMass ();
  void FillHists();

  void  GetClsn(AMSEventR *evt, BetaHR *bh);
  Bool_t CutTk2(AMSEventR *evt, BetaHR *bh);
  Int_t  NtrdSg(AMSEventR *evt);

  void Stat();
  void Fill(const char *hname, Double_t x, Double_t y, Double_t w = 1);

  void RunSel();
  void RTISel();
  void MCinit();
  void MCtrig(const char *hname, Int_t emin = -1, Int_t emax = -1);

  Int_t SpecialCut(AMSPoint pnt, AMSDir dir);

  TrTrackR *GetTrack(void);
};

Int_t   MdSel::RUN_STD   = 1510526122;  //1494600586;1480201390;1464299568;1448664052;
Int_t   MdSel::RUN_TTCS1 = 1412285225;  //1411991496
Int_t   MdSel::RUN_TTCS2 = 1417187198;

Int_t   MdSel::fNtot = 0;
TString MdSel::fFname;
TChain *MdSel::fChain = 0;

extern double memchk(void);

void MdSel::UBegin()
{
  fFile = new TFile(fFname, "recreate");
  fTree = new TTree("tree", "mdst");

  fDst.Branch(fTree);

  TString sbn = "$AMSDataDir/v5.00/phe_bin2.root";
  gSystem->ExpandPathName(sbn);
//TString sbn = "/afs/cern.ch/ams/Offline/AMSDataDir/v5.00/phe_bin2.root";

  TFile *fb = TFile::Open(sbn);
  if (!fb) exit(-1);

  TH1F    *hbin = (TH1F *)fb->Get("hist1");
  Int_t    nbin = hbin->GetNbinsX();
  Double_t *bin = new Double_t[nbin+3];

  for (Int_t i = 0; i <= nbin; i++)
    bin[i] = hbin->GetXaxis()->GetXbins()->fArray[i];

  bin[++nbin] = 10000;
  bin[++nbin] = 20000;

//LxMCcutoff::SetBinning(hbin);

  Double_t pbin[75]
    = { 0.50, 0.65, 0.82, 1.01, 1.22, 1.46, 1.72, 2.00, 2.31, 2.65,
	3.00, 3.36, 3.73, 4.12, 4.54, 5.00, 5.49, 6.00, 6.54, 7.10,
	7.69, 8.30, 8.95, 9.62, 10.3, 11.0, 11.8, 12.6, 13.4, 14.2,
	15.1, 16.1, 17.0, 18.0, 19.0, 20.0, 21.1, 22.2, 23.4, 24.6,
	25.9, 27.2, 28.7, 30.2, 31.8, 33.5, 35.4, 37.3, 39.4, 41.6,
	44.0, 46.6, 49.3, 52.3, 55.6, 59.1, 63.0, 67.3, 72.0, 77.4,
	83.4, 90.2, 98.1, 107., 118., 132., 149., 170., 198., 237.,
	290., 370., 500., 700., 1000 };

  Double_t vbin[68*2+1];
  Double_t ebin[68]
    = { 0.50,  0.65,  0.81,  1.00,  1.21,  1.45,  1.70,  1.97,
	2.28,  2.60,  2.94,  3.30,  3.70,  4.11,  4.54,  5.00,
	5.50,  6.00,  6.56,  7.16,  7.80,  8.50,  9.21,  9.95,
	10.73, 11.54, 12.39, 13.27, 14.19, 15.15, 16.15, 17.18,
	18.25, 19.37, 20.54, 21.76, 23.07, 24.45, 25.87, 27.34,
	28.87, 30.45, 32.10, 33.80, 35.57, 37.40, 40.00, 43.39,
	47.01, 50.87, 54.98, 59.36, 64.03, 69.00, 74.30, 80.00,
	86.00, 92.50, 100.0, 115.1, 132.1, 151.5, 173.5, 206.0,
	260.0, 350.0, 500.0, 1000. };

  for (Int_t i = 0; i < 66; i++) vbin[i]    = -1/ebin[i];
  for (Int_t i = 0; i < 66; i++) vbin[i+67] =  1/ebin[65-i];
                                 vbin[66]   =  0;

  TString sfn = fChain->GetFile()->GetName();
  TString ss0 = sfn(0, sfn.Last('/')-1);
  TString ss  = ss0(1+ ss0.Last('/'), 5);

  cout << "ss= "<<ss.Data() <<endl;
  if (ss.BeginsWith("el.")) {
    nbin = 74;
    bin  = pbin;
    cout << "Using e+/e- binning" << endl;
  }

  fFile->cd();

  new TH1F("hist00", "Texp VS Rcut", nbin, bin);
  new TH1F("hmct00", "Texp VS Rcut", nbin, bin);

  for (Int_t i = 0; i < 4; i++) {
    TString shn = Form("hist%d", i+1);

    new TH2F(shn+"1", "Rcut VS Rgt (p)",     nbin, bin, nbin, bin);
    new TH2F(shn+"2", "Rcut VS Rgt (He)",    nbin, bin, nbin, bin);
    new TH2F(shn+"3", "Rcut VS Rgt (p  ub)", nbin, bin, nbin, bin);
    new TH2F(shn+"4", "Rcut VS Rgt (He ub)", nbin, bin, nbin, bin);

    new TH2F(shn+"5", "Qin  VS Rgt",         nbin, bin, 500, 0, 10);
    new TH2F(shn+"6", "Qin  VS Rgt (ub)",    nbin, bin, 500, 0, 10);
  }

  new TH2F("hist51", "TrdK", 132, vbin, 200,   0,  2);
  new TH2F("hist52", "CFi",  132, vbin,  65,    ebin);
  new TH2F("hist53", "E/p",   65, ebin, 200, -10, 10);


  Int_t pass = 6; //4;
/*if (fChain && fChain->GetListOfFiles() &&
                fChain->GetListOfFiles()->At(0)) {
    TString sfn = fChain->GetListOfFiles()->At(0)->GetTitle();
    if (sfn.Contains("B950/pass6")) pass = 6;
  }
*/
  AMSSetupR::RTI::UseLatest(pass);

 if (!fFname.Contains("tb/") && !fFname.Contains("mdst_1281"))
  TkDBc::UseFinal();

  TRMCFFKEY_DEF ::ReadFromFile = 0;
  TRFITFFKEY_DEF::ReadFromFile = 0;
  TRFITFFKEY.magtemp = 0;

  if (fFname.Contains("tb/"))
    TRFITFFKEY.ErrYL1 = TRFITFFKEY.ErrYL9 = 0;


  fTimer.Start();
}

bool MdSel::UProcessStatus(unsigned long long status)
{
  fNevt++;
  if (fNevt%fIntv == 0 || fNevt == fNtot) Stat();

  if (fNtot == 0 || !fTree) return false;
  return true;
}

bool MdSel::UProcessCut()
{
  MCEventgR *mcg = pMCEventg(0);
  if (mcg) {
    if (fBnew) {
      MCinit();
      fBnew = 0;
    }
    Double_t pmc  = mcg->Momentum;
    Double_t pmin = fPGMD[0];
    Double_t pmax = fPGMD[1];
    if (fPGMD[0] < 0 && fPGMD[1] < 0) { 
      if (!pParticle(0)) return false;

      TrTrackR *trk = pParticle(0)->pTrTrack();
      if (!trk) return false;

      if (trk->HasExtLayers() != 3) return false;
      pmin = -pmin;
      pmax = -pmax;
    }
    if (pmin > 0 && pmax > pmin && (pmc < pmin || pmc > pmax)) return false;
    return true;
  }
  else if (fBnew) {
      if (!mcg && Run() >= RUN_STD) {    // std data
      cout << "Disabing slowcontrol" << endl;
      AMSSetupR::SlowControlR::ReadFromExternalFile = false;
    }

    fBnew = 0;
  }
  return true;
}

void MdSel::UProcessFill()
{
  if (pMCEventg(0)) TRCLFFKEY.UseSensorAlign = 0;

  fDst.Clear();

  fDst.fStatus.status = fStatus;
  fDst.fStatus.ustat  = 0;

  RunSel();
  RTISel();

  if (!pMCEventg(0)) {
    if (fBadRun&4) fDst.fStatus.ustat |= BadTrig;
    if (fBadRun&1) fDst.fStatus.ustat |= BadTRD;
    if (fBadRTI)   fDst.fStatus.ustat |= BadRTI;
    if (fInSAA)    fDst.fStatus.ustat |= InSAA;

    if (UTime() <= fTRTI[0]) fDst.fStatus.ustat |= FirstSec;
    if (UTime() >= fTRTI[1]) fDst.fStatus.ustat |= LastSec;
  }

  fDst.fHeader.run   = Run();
  fDst.fHeader.event = Event();
  fDst.fHeader.ient  = fNevt-1;
  fDst.fHeader.utime = UTime();
  fDst.fHeader.phpat = 0;
  fDst.fHeader.jmpat = 0;

  if (fLevt < Event()) fLevt = Event();

  if (pLevel1(0)) {
    Level1R &lev1 = Level1(0);
    Int_t   sdpat = lev1.JMembPatt;
    Int_t   phpat = lev1.PhysBPatt;
    if (nMCEventg() > 0) lev1.RebuildTrigPatt(sdpat, phpat);
    else                 lev1.RestorePhysBPat();

    fDst.fHeader.phpat = phpat;
    fDst.fHeader.jmpat = sdpat;
  }

  fDst.fHeader.ntrdhit  = nTrdRawHit();
  fDst.fHeader.ntrdseg  = nTrdSegment();
  fDst.fHeader.ntrhit   = nTrRecHit();
  fDst.fHeader.ntofcls  = nTofCluster();
  fDst.fHeader.nanticls = nAntiCluster();
  fDst.fHeader.nrichhit = nRichHit();
  fDst.fHeader.necalhit = nEcalHit();
  fDst.fHeader.error    = EventError();

  if (!pMCEventg(0)) {
    fDst.fRTI.theta  = fRTI.theta;
    fDst.fRTI.phi    = fRTI.phi;
    fDst.fRTI.zenith = fRTI.zenith;
    fDst.fRTI.cfi    = fRTI.cfi[0][1];
    fDst.fRTI.lf     = fLtf;
    fDst.fRTI.dl1    = fPGMD[0];
    fDst.fRTI.dl9    = fPGMD[1];
    fDst.fRTI.cf     = fRTI.cf[0][1];
  }

  // B935 mitigation
  if (pMCEventg(0) && Version() == 935 && pParticle(0) && 
                                          pParticle(0)->pTrdTrack() &&
                                         !pParticle(0)->pTrTrack()) {
    TrdTrackR *trd = pParticle(0)->pTrdTrack();
    AMSPoint ptrd(trd->Coo[0], trd->Coo[1], trd->Coo[2]);
    AMSDir   dtrd(trd->Theta, trd->Phi);

    AMSPoint p, d;
    Int_t tk;
    Int_t is1 = IsInsideTracker(1, ptrd, dtrd, 0, 0, tk, p, d);
    Int_t is9 = IsInsideTracker(9, ptrd, dtrd, 0, 0, tk, p, d);

    float qmax = 0;

    if (is1 && is9) {
      int qopt = TrClusterR::kAsym | TrClusterR::kGain |
	         TrClusterR::kLoss | TrClusterR::kMIP;
      for (Int_t i = 0; i < nTrRecHit(); i++) {
	TrRecHitR *hit = pTrRecHit(i);
	if (!hit || hit->GetLayerJ() != 1) continue;

	Float_t qh = TMath::Sqrt(hit->GetSignalCombination(2, qopt, 1, 1));
	if (qh > qmax) qmax = qh;
      }
    }

    if (qmax > 1.6) {
      TRMCFFKEY_DEF bak = TRMCFFKEY;
      for (Int_t i = 0; i < 4; i++) TRMCFFKEY.OuterSmearing [i/2][i%2] = 
				    TRMCFFKEY.OuterSmearingC[i/2][i%2] = 0;
      TrRecon trec;
      TofRecH::BuildOpt = 0;
      trec.Build(-1111, 1, 0);
      UpdateTrRecon();

      TRMCFFKEY = bak;
    }
  }

  // TTCS-off simulation
  if (0 && pParticle(0)) {
    TrRecon trec;
    TRCLFFKEY.AllowYonlyTracks = 1;
    TRMCFFKEY.TrSim2010_AddNoise[0] = 0.9;
    TRMCFFKEY.TrSim2010_AddNoise[1] = 0.8;
    TRFITFFKEY.MergeExtLimX = TRFITFFKEY.MergeExtLimY = 5;

    VCon *cont = GetVCon()->GetCont("AMSTrTrack");
    cont->eraseC();
    delete cont;

    if (nTrCluster() > 0 && nTrRecHit() < 1000) trec.Build(201111, 1, 1);

    ParticleP *pp = (ParticleP *)pParticle(0);
    if (NTrTrack() >= 1) pp->setTrack(0);
    else pp->setTrack(-1);

    TofRecH tofrec;
    TofRecH::BuildOpt = 0;
    tofrec.ReBuild(0);
  }


  // TTCS-off runs
  if (!pMCEventg(0) &&
       pParticle(0) && RUN_TTCS1 <= Run() && Run() <= RUN_TTCS2) {
    static Int_t first = 1;
    if (first) {
      cout << "TTCS off run" << endl;
      first = 0;
    }

    for (Int_t i = 0; i < nTrRecHit(); i++)
      if (pTrRecHit(i)) pTrRecHit(i)->ClearUsed();

    for (Int_t i = 0; i < nTrCluster(); i++)
      if (pTrCluster(i)) pTrCluster(i)->ClearUsed();

    TrRecon trec;
    TRCLFFKEY.AllowYonlyTracks = 1;

    if (trec.BuildTrTracksSimple(0)) {
      TrTrackR *t1 = pParticle(0)->pTrTrack();
      TrTrackR *t2 = pTrTrack(NTrTrack()-1);

      Int_t chk = (t1 && t2) ? 1 : 0;
      for (Int_t i = 0; chk && t1 && t2 && i < 4; i++) {
	TrRecHitR *h1 = t1->GetHitLJ(i*2+2);
	TrRecHitR *h2 = t2->GetHitLJ(i*2+2);
	if ((h1 && !h2) || (!h1 && h2)) chk = 0;
	if (h1 && h2) {
	  if (h1->GetTkId()       != h2->GetTkId() ||
	      h1->iTrCluster('y') != h2->iTrCluster('y')) chk = 0;
	}
      }

      if (!chk) {
/*	for (Int_t i = 0; i < 2; i++) {
	  TrTrackR *t = (i == 0) ? t1 : t2;
	  cout << Form("TR%d ",i);
	  for (Int_t j = 0; t && j < 4; j++){
	    TrRecHitR *h=t->GetHitLJ(j*2+2);
	    if (h) cout<<Form(" %4d (%2d)",h->GetTkId(),h->iTrCluster('y'));
	    else   cout<<" ---- (--)";
	  }
	  cout<<endl;
	}
	cout<<endl;
*/
	ParticleP *pp = (ParticleP *)pParticle(0);
	pp->setTrack(NTrTrack()-1);
      }
    }
  }
  // TTCS-off runs

  // std data
  if (!pMCEventg(0) && Run() >= RUN_STD  && pParticle(0) &&
                                            pParticle(0)->pTrTrack() &&
      (Run() < RUN_TTCS1 || RUN_TTCS2 < Run())) {

    TrTrackR *trk = pParticle(0)->pTrTrack();

    Int_t itp = TrTrackR::DefaultFitID;
    if (!trk->ParExists(itp)) itp = TrTrackR::kChoutko;
    if (!trk->ParExists(itp)) itp = TrTrackR::kAlcaraz;

    TrRecon trec;
    if (trk->ParExists(itp)) {
      if (trk->GetHitLJ(1)) trk->GetHitLJ(1)->ClearUsed();
      if (trk->GetHitLJ(9)) trk->GetHitLJ(9)->ClearUsed();

      TRFITFFKEY.MergeExtLimX = TRFITFFKEY.MergeExtLimY = 5;
      trec.MergeExtHits(trk, itp);
    }
  }

  FillMC   ();
  FillVtx  ();
  FillPart ();
  FillTrack();
  FillTrHit();
  FillMass ();
  FillBeta ();
  FillTrd  ();
  FillRich ();
  FillEcal ();
  FillHists();

  if (fTree && fTree->Fill() < 0) {
    cout << "Fill failed" << endl;
    fTree = 0;
  }
  fNfil++;
}

void MdSel::FillVtx()
{
//if (pMCEventg(0)) return;
//if (fRTI.zenith < 40 && nTrdTrack() > 0) return;
//if (Run() < 1300000000) return;
//if (Run() >= RUN_STD)   return;  // std data
//if (Version() >= 900)   return;

//TrRecon trec;
//trec.Build(10000, 0, 0);

  VertexR *vtx = 0;
  for (Int_t i = 0; i < NVertex(); i++) {
    VertexR *v = pVertex(i);
    if (!vtx || (v && v->IsPhotonVertex())) vtx = v;
  }
  if (!vtx) return;

  TrTrackR *trk1 = vtx->pTrTrack(0);
  TrTrackR *trk2 = vtx->pTrTrack(1);
  if (!trk1 || !trk2) return;

  Int_t fid = TrTrackR::kVertex;

  AMSDir dirv(vtx->Theta, vtx->Phi);
  if (dirv.z() > 0) dirv = dirv*(-1);

  Double_t vzen = 0;
  if (!pMCEventg(0) &&
      trk1->ParExists(fid) && trk2->ParExists(fid)) {
    Double_t vth = dirv.gettheta();
    Double_t vph = dirv.getphi();

    Int_t res;
    Double_t lng = 0, lat = 0;
    GetGalCoo(res, lng, lat, vth, vph, 1, 4, 3, 0, 3);

    AMSDir dirs(TMath::Pi()/2-fHeader.ThetaS, fHeader.PhiS);
    AMSDir dirg(TMath::DegToRad()*(90-lat), TMath::DegToRad()*lng);
    vzen = TMath::ACos(dirs.prod(dirg))*TMath::RadToDeg();

    fDst.fVertex.rgt [0] = trk1->GetRigidity  (fid);
    fDst.fVertex.rgt [1] = trk2->GetRigidity  (fid);
    fDst.fVertex.csqx[0] = trk1->GetNormChisqX(fid);
    fDst.fVertex.csqx[1] = trk2->GetNormChisqX(fid);
    fDst.fVertex.csqy[0] = trk1->GetNormChisqY(fid);
    fDst.fVertex.csqy[1] = trk2->GetNormChisqY(fid);

    fDst.fVertex.coo [0] = trk1->GetP0x(fid);
    fDst.fVertex.coo [1] = trk1->GetP0y(fid);
    fDst.fVertex.coo [2] = trk2->GetP0x(fid);
    fDst.fVertex.coo [3] = trk2->GetP0y(fid);

    fDst.fVertex.ang [0] = trk1->GetTheta(fid);
    fDst.fVertex.ang [1] = trk1->GetPhi  (fid);
    fDst.fVertex.ang [2] = trk2->GetTheta(fid);
    fDst.fVertex.ang [3] = trk2->GetPhi  (fid);
  }
  else {
    fDst.fVertex.rgt [0] = trk1->GetRigidity  ();
    fDst.fVertex.rgt [1] = trk2->GetRigidity  ();
    fDst.fVertex.csqx[0] = trk1->GetNormChisqX();
    fDst.fVertex.csqx[1] = trk2->GetNormChisqX();
    fDst.fVertex.csqy[0] = trk1->GetNormChisqY();
    fDst.fVertex.csqy[1] = trk2->GetNormChisqY();

    fDst.fVertex.coo [0] = trk1->GetP0x();
    fDst.fVertex.coo [1] = trk1->GetP0y();
    fDst.fVertex.coo [2] = trk2->GetP0x();
    fDst.fVertex.coo [3] = trk2->GetP0y();

    fDst.fVertex.ang [0] = trk1->GetTheta();
    fDst.fVertex.ang [1] = trk1->GetPhi();
    fDst.fVertex.ang [2] = trk2->GetTheta();
    fDst.fVertex.ang [3] = trk2->GetPhi();
  }

  fDst.fVertex.theta   = vtx->Theta;
  fDst.fVertex.phi     = vtx->Phi;
  fDst.fVertex.zenith  = vzen;

  fDst.fVertex.coo [4] = vtx->Vertex[0];
  fDst.fVertex.coo [5] = vtx->Vertex[1];
  fDst.fVertex.coo [6] = vtx->Vertex[2];
}

void MdSel::FillPart()
{
  ParticleR *part = pParticle(0);
  if (!part) return;

  part->Info(0, this);
  fDst.fPart.mon   = part->Momentum;
  fDst.fPart.emom  = part->ErrMomentum;
  fDst.fPart.mass  = part->Mass;
  fDst.fPart.emass = part->ErrMass;
  fDst.fPart.beta  = part->Beta;
  fDst.fPart.ebeta = part->ErrBeta;
}

void MdSel::FillTrack()
{
  fTrk  = 0;
  fBeta = 0;

  ParticleR *part = pParticle(0);
  if (!part) return;

  fTrk = part->pTrTrack();
  if (!fTrk) { fTrk = GetTrack(); if (fTrk) fDst.fStatus.ustat |= NoPTrk; }
  if (!fTrk) return;

  Bool_t ismc = (pMCEventg(0)) ? kTRUE : kFALSE;
  if (ismc) {
    TString str = fFname(fFname.Last('.')-1, 1);
    Int_t rr = str.Atoi();

    SetDefaultMCTuningParameters();

    MCEventgR *mcg = GetPrimaryMC();
  //if (mcg && mcg->Charge == 3) TRMCFFKEY.MCtuneDs[1] = -310.; // Li P6 L1
  //if (mcg && mcg->Charge == 3) TRMCFFKEY.MCtuneDs[1] = -330.; // Li P6 FS
  //if (mcg && mcg->Charge >= 7) TRMCFFKEY.MCtuneDs[1] = -650.; // Use C
    if (mcg && mcg->Charge == 6) TRMCFFKEY.MCtuneDs[1] =    0.;

    static Int_t first = 1;

    if (1) {
      if (first) {
          first = 0;

          cout << "TRMCFFKEY "
                  "MCtuneDmax = "    << TRMCFFKEY.MCtuneDmax  << " "
               << "MCtuneDs = "      << TRMCFFKEY.MCtuneDs[0] << " "
                                     << TRMCFFKEY.MCtuneDs[1] << " "
               << "MCtuneDy9 = "     << TRMCFFKEY.MCtuneDy9   << endl;
          cout << "TRMCFFKEY "
               << "MCscat = "        << TRMCFFKEY.MCscat[0]  << " "
                                     << TRMCFFKEY.MCscat[1]  << " "
                                     << TRMCFFKEY.MCscat[2]  << " "
               << "OuterSmearing = " << TRMCFFKEY.OuterSmearing[0][1] << " "
                                     << TRMCFFKEY.OuterSmearing[1][1] << endl;
        }	

      TRFITFFKEY.Zshift = -1;
      TrExtAlignDB::SmearExtAlign();
      if(Version()>=1106) TRCLFFKEY.ClusterCofGOpt = 1;
    }
    else if (fTrk) {
     Double_t rgt = fTrk->GetRigidity();
     AMSPoint ptr = fTrk->GetP0 ();
     AMSDir   dtr = fTrk->GetDir();
     AMSPoint p, d;
     Int_t tk;
     Int_t is1 = IsInsideTracker(1, ptr, dtr, rgt, 0, tk, p, d);
     Int_t is9 = IsInsideTracker(9, ptr, dtr, rgt, 0, tk, p, d);

     if (is1 || is9 || fTrk->HasExtLayers() == 3) {
      if (TMath::Abs(TRMCFFKEY.noise_fac[1]-1.45) < 0.01) {
	//NearStripMult[0][1][0]=1.3;
	//noise_fac[1]=1.45;
      //TRMCFFKEY.noise_fac[1]        = 1.37;
	TRMCFFKEY.noise_fac[1]        = 1.20;
	TRMCFFKEY.NearStripMult[0][1][0] = 2.;
	TRMCFFKEY.OuterSmearing[0][1] = 8.5e-4;
	TRMCFFKEY.OuterSmearing[1][1] = 9.5e-4;
	cout << "TRMCFFKEY "
	        "noise_fac[1] = "  << TRMCFFKEY.noise_fac[1] << " "
	     << "OuterSmearing = " << TRMCFFKEY.OuterSmearing[0][1] << " "
	                           << TRMCFFKEY.OuterSmearing[1][1] << endl;
      }

      TrRecon trec;
      trec.Build(-1111, 1, 0);
      fTrk = pTrTrack(0);
     }
   }

    if (!fTrk) return;
    //TofRecH tofrec;
    //TofRecH::BuildOpt = 0;
    //tofrec.ReBuild();

  }
  else{ TrLinearEtaDB::SetLinearCluster(); TRFITFFKEY.Zshift = 2; }

  BetaHR *bth = (ismc) ? part->pBetaH() : 0;
  for (Int_t i = 0; !bth && i < NBetaH(); i++)
    if (pBetaH(i) && pBetaH(i)->pTrTrack() && pBetaH(i)->pTrTrack() == fTrk) bth = pBetaH(i);

  fBeta = (bth) ? bth->GetBeta() : part->Beta;
  Double_t qtmp = fTrk->GetInnerQ_all(fBeta).Mean;

  Double_t qthd = (qtmp < 1.5) ? 0 : 1.5;
  fDst.fTrack.addl = fTrk->AddLostHits(1, qthd, 1);

  fDst.fTrack.rec  = fTrk->GetRecType();

  Double_t monz = (qtmp < 1.5) ? TrFit::Mproton : TrFit::Mhelium/2;


  mean_t   qmn  = fTrk->GetInnerQ_all(fBeta, 0, monz);
  Double_t qtrk = qmn.Mean;
  Double_t qrms = qmn.RMS;
  Double_t mass = TrFit::Mproton;
  Double_t chrg = 1;
  Double_t Z = TMath::Floor(qtrk+0.5);
  if ( Z>=2 ) { mass = TrFit::Mhelium*Z/2.; chrg = 2; }

  Float_t bcor = -1;
  if (!ismc && 1300000000 < Run() && Run() < RUN_STD) {
    Float_t bcor1 = 1, bcor2 = 1;
    int bret1 = MagnetVarp::btempcor(bcor1, 0, 1);
    int bret2 = MagnetVarp::btempcor(bcor2, 0, 2);
    if      (bret1 == 0 && bret2 == 0) bcor = (bcor1+bcor2)/2; 
    else if (bret1 != 0 && bret2 == 0) bcor =        bcor2;
  }

  enum { NF = 4+3 };
                                       
  Int_t    alg [NF] = {  1,  1,  1,  1,   3,  1,  1 };
  Int_t    pat [NF] = {  3,  5,  6,  7,   3,  1,  2 };
  Int_t    rfit[NF] = { 23, 23, 23, 23,  23, 23, 23 };
  Int_t    itp [NF] = { -1, -1, -1, -1,  -1, -1, -1 };
  Double_t rgt [NF] = {  0,  0,  0,  0,   0,  0,  0 };
  Double_t csqx[NF] = { -1, -1, -1, -1,  -1, -1, -1 };
  Double_t csqy[NF] = { -1, -1, -1, -1,  -1, -1, -1 };

  if (!ismc) {
    if (Version() >= 933) {
      rfit[0] = 0;
      rfit[1] = rfit[2] = rfit[3] = 21;
      rfit[4] = rfit[5] = rfit[6] = 21;
      bcor = 1;
    }
    if (Version() >= 900 && UTime() > 1305700000) {
      bcor = 1;
      TRFITFFKEY.magtemp = 1;
    }
    else TRFITFFKEY.magtemp = 0;
  }

//if (!ismc && Run() < 1305800000)
//  TkDBc::Head->GetPlane(1)->UpdatePosA()[1] = -0.4e-4;
//  for (Int_t i = 0; i < NF; i++) rfit[i] = 0;

  if(ismc){
    Int_t itp0 = fTrk->iTrTrackPar(1, 0, 23, mass, Z);
    if(itp0<0) return;
    for(Int_t i = 0 ; i < nMCEventg(); i++){
      MCEventgR *mcg = pMCEventg(i);
      Int_t Layer = -1;
      if( mcg->Nskip == -1000 ) Layer = 1;
      if( mcg->Nskip == -1018 ) Layer = 9;
      if( -1013<= mcg->Nskip && mcg->Nskip <= -1007) Layer = -mcg->Nskip-1007+2;
      if( Layer < 0 || Layer>9 || !fTrk->TestHitLayerJ(Layer)) continue;
      fDst.fRefit.tx[Layer-1] = mcg->Coo[0];
      fDst.fRefit.ty[Layer-1] = mcg->Coo[1];
      fDst.fRefit.tz[Layer-1] = mcg->Coo[2];
      fDst.fRefit.fx[Layer-1] = TrTrackR::FitCoo[Layer-1].x();
      fDst.fRefit.fy[Layer-1] = TrTrackR::FitCoo[Layer-1].y();
      fDst.fRefit.fz[Layer-1] = TrTrackR::FitCoo[Layer-1].z();

      TrRecHitR *hit = fTrk->GetHitLJ(Layer);
      Bool_t HasX = hit->GetXCluster()!=NULL;
      fDst.fRefit.nsx [Layer-1] = HasX? hit->GetXCluster()->GetNelem():0;
      fDst.fRefit.nsy [Layer-1] =       hit->GetYCluster()->GetNelem();
    }
  }

  itp[0] = fTrk->iTrTrackPar(alg[0], pat[0], rfit[0], mass, Z);
  if (itp[0] < 0) itp[0] = fTrk->iTrTrackPar(alg[0], pat[0], 0);
  if (itp[0] < 0) return;


  rgt [0] = fTrk->GetRigidity  (itp[0]);
  csqx[0] = fTrk->GetNormChisqX(itp[0]);
  csqy[0] = fTrk->GetNormChisqY(itp[0]);

  if (bcor > 0) rgt[0] *= bcor;      

  fDst.fTrack.p0x   = fTrk->GetP0x  (itp[0]);
  fDst.fTrack.p0y   = fTrk->GetP0y  (itp[0]);
  fDst.fTrack.theta = fTrk->GetTheta(itp[0]);
  fDst.fTrack.phi   = fTrk->GetPhi  (itp[0]);

  AMSPoint pl1 = fTrk->InterpolateLayerJ(1, itp[0]);
  AMSPoint pl2 = fTrk->InterpolateLayerJ(2, itp[0]);
  AMSPoint pl9 = fTrk->InterpolateLayerJ(9, itp[0]);

  fDst.fTrack.coox[0] = pl1.x(); fDst.fTrack.cooy[0] = pl1.y();
  fDst.fTrack.coox[1] = pl2.x(); fDst.fTrack.cooy[1] = pl2.y();
  fDst.fTrack.coox[2] = pl9.x(); fDst.fTrack.cooy[2] = pl9.y();

  Double_t qmin[2] = { 0.7, 0.7 };
  Int_t    mopt[2] = {   8,   9 };
  if (qtrk > 1.5) qmin[0] = 1.5;

  for (Int_t i = 0; i < 2; i++) {
    Int_t lj = (i == 0) ? 1 : 9;

    if (!ismc && !fTrk->TestHitLayerJ(lj))
      for (Int_t j = 0; j < 2; j++) {
        if (fTrk->MergeHits(lj, 100, qmin[j], 0, fBeta, mopt[j])) {
          fDst.fTrack.bitr |= 1 << (i*2+j);
          break;
        }
      }
    TrRecHitR *hit = fTrk->GetHitLJ(lj);
    if (hit) {
      AMSPoint pl  = (i == 0) ? pl1 : pl9;
      AMSPoint coo = hit->GetCoord(-1, 3);  // (PG+MD)/2
      AMSPoint cpg = hit->GetCoord(-1, 1);  //  PG
      AMSPoint cmd = hit->GetCoord(-1, 2);  //  MD
      if (!hit->OnlyY()) fDst.fTrack.resx[i] = coo.x()-pl.x();
      fDst.fTrack.resy[i] = coo.y()-pl.y();
      fDst.fTrack.hity[i*2]   = cpg.y();
      fDst.fTrack.hity[i*2+1] = cmd.y();
      fDst.fTrack.qsta[i] = hit->GetQStatus();
    }
  }

  fDst.fTrack.bith = fTrk->GetBitPatternJ();
  fDst.fTrack.bitx = fTrk->GetBitPatternXYJ();

  for (Int_t i = 1; i < NF; i++) {
    if (!ismc && Run() >= RUN_STD) rfit[i] = 0;  // std data

    itp[i] = fTrk->iTrTrackPar(alg[i], pat[i], rfit[i], mass, chrg);

    if (itp[i] > 0) {
      rgt [i] = fTrk->GetRigidity  (itp[i]);
      csqx[i] = fTrk->GetNormChisqX(itp[i]);
      csqy[i] = fTrk->GetNormChisqY(itp[i]);

      if (bcor > 0) rgt[i] *= bcor;      
    }

  }

  for (Int_t i = 0; i < 4; i++) {
    fDst.fTrack.rgt [i] = rgt [i];
    fDst.fTrack.csqx[i] = csqx[i];
    fDst.fTrack.csqy[i] = csqy[i];
  }

  fDst.fTrack.rgt2[0] = rgt [4];
  fDst.fTrack.csx2    = csqx[4];
  fDst.fTrack.csy2    = csqy[4];

  fDst.fTrack.rgt2[1] = rgt [5];
  fDst.fTrack.rgt2[2] = rgt [6];

  MCEventgR *mc = GetPrimaryMC();
  if (mc && mc->Charge != 0) {
    Double_t rg = mc->Momentum/mc->Charge;
    if (rgt[0] != 0) fDst.fMCinfo.cfw[0] = GetMCCutoffWeight(rg, rgt[0]);
    if (rgt[3] != 0) fDst.fMCinfo.cfw[1] = GetMCCutoffWeight(rg, rgt[3]);
  }

  fDst.fTrack.qin  = qtrk;
  fDst.fTrack.qrms = qrms;
  fDst.fTrack.ql1  = fTrk->GetLayerJQ(1, fBeta);
  fDst.fTrack.ql9  = fTrk->GetLayerJQ(9, fBeta);

  for (Int_t i = 0; i < 7; i++)
    fDst.fTrHit.qli[i] = fTrk->GetLayerJQ(i+2, fBeta);

  AMSPoint ptr = fTrk->GetP0 (itp[0]);
  AMSDir   dtr = fTrk->GetDir(itp[0]);

  AMSPoint p, d;
  Int_t tk;
  Int_t is1 = IsInsideTracker(1, ptr, dtr, rgt[0], 0, tk, p, d);
  Int_t is9 = IsInsideTracker(9, ptr, dtr, rgt[0], 0, tk, p, d);

  if (is1) fDst.fStatus.ustat |= IsInL1;
  if (is9) fDst.fStatus.ustat |= IsInL9;
  if (TMath::Abs(pl9.x()) < 33) fDst.fStatus.ustat |= IsInEcal;

  Bool_t hasl2 = fTrk->TestHitLayerJ(2);
  Bool_t hasl1 = fTrk->TestHitLayerJ(1);
  Bool_t hasl9 = fTrk->TestHitLayerJ(9);

  if (fDst.fTrack.bitr&2) hasl1 = 0;
  if (fDst.fTrack.bitr&8) hasl9 = 0;

  if (hasl2) fDst.fStatus.ustat |= HasL2;
  if (hasl1) fDst.fStatus.ustat |= HasL1;
  if (hasl9) fDst.fStatus.ustat |= HasL9;
}

void MdSel::FillTrHit()
{
  Int_t   ihit[4] = { -1,  -1,  -1, -1 };
  Int_t   nhit[6] = { 0, 0, 0, 0, 0, 0 };
  Float_t qmax[6] = { 0, 0, 0, 0, 0, 0 };

  int  qopt = TrClusterR::kAsym | TrClusterR::kGain |
	      TrClusterR::kLoss | TrClusterR::kMIP  |
              TrClusterR::kBeta | TrClusterR::kRigidity;

  Double_t beta = fBeta;
  Double_t rgt  = (fTrk) ? TMath::Abs(fTrk->GetRigidity()) : 1;

  Int_t nh10[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  if (fTrk) IsTrackPickingUpNoise(10, 0, fTrk, nh10);

  for (Int_t i = 0; i < 7; i++) {
    fDst.fTrHit.nh10[i] = nh10[i+1];
    fDst.fTrHit.rata[i] = GetTrackerRawSignalRatio(i+2);
    fDst.fTrHit.tkfd[i] = GetTkFeetDist(i+2);
  }

  for (Int_t i = 0; i < nTrRecHit(); i++) {
    TrRecHitR *hit = pTrRecHit(i);

    if (Version() < 800 && (hit->GetQStatus()&0x1001FD)) continue; // pass4
    if (Version() > 800 && (hit->GetQStatus()&0x10013D)) continue; // pass6

    Int_t il = (hit->GetLayerJ() == 1) ? 0 : 
              ((hit->GetLayerJ() == 9) ? 2 : 4);
    if (hit->OnlyY()) il += 1;

    nhit[il]++;

    Float_t qhit = TMath::Sqrt(hit->GetSignalCombination(2, qopt, beta, rgt));
    if (qhit > qmax[il]) {
      qmax[il] = qhit;
      if (il < 4) ihit[il] = i;
    }
  }

  for (Int_t i = 0; i < 6; i++) {
    fDst.fTrHit.nhit[i] = nhit[i];
    fDst.fTrHit.qmax[i] = qmax[i];
  }

  for (Int_t i = 0; i < 4; i++) {
    TrRecHitR *hit = (ihit[i] > 0) ? pTrRecHit(ihit[i]) : 0;
    if (!hit) continue;

    Int_t it = (i%2 == 0) ? 0 : 2;
    if (i < 2)
      fDst.fTrHit.dx[i] = hit->GetCoord().x()-fDst.fTrack.coox[it];
      fDst.fTrHit.dy[i] = hit->GetCoord().y()-fDst.fTrack.cooy[it];
  }

  for (Int_t i = 0; fTrk && i < fTrk->GetNhits(); i++) {
    TrRecHitR *hit = fTrk->GetHit(i);
    if (!hit) continue;

    Int_t lj   = hit->GetLayerJ()-1;
    Int_t tkid = hit->GetTkId();
    Int_t imlt = hit->GetResolvedMultiplicity();

  //if (2 <= lj && lj <= 8)
  //  fDst.fTrHit.qsi[lj-2] = hit->GetQStatus();
    if (1 <= lj && lj <= 7)
      fDst.fTrHit.qsi[lj-1] = hit->GetQStatus();   // Bugfix 160407

    TrClusterR *clx = hit->GetXCluster();
    TrClusterR *cly = hit->GetYCluster();
    if (!cly) continue;

    Int_t flag = 0;
    Int_t is = cly->GetSeedIndex();
    if (cly->GetStatus(is-1)) flag += 1;
    if (cly->GetStatus(is)  ) flag += 2;
    if (cly->GetStatus(is+1)) flag += 4;

    flag += cly->GetNelem()*10;

    Double_t eta = cly->GetCofG(2); if (eta < 0) eta += 1;

    fDst.fTrCls.tkml[lj] = TMath::Sign(imlt*1000+TMath::Abs(tkid), tkid);
    fDst.fTrCls.xcog[lj] = (!clx) ? -(hit->GetDummyX()+640)
	                           : clx->GetCofG()+clx->GetSeedAddress();
    fDst.fTrCls.ycog[lj] = cly->GetCofG()+cly->GetSeedAddress();
    fDst.fTrCls.flag[lj] = flag;
    fDst.fTrCls.eta [lj] = eta;
    fDst.fTrCls.sig [lj] = cly->GetTotSignal();
  }

  if (nTrMCCluster() > 0) {
    for (Int_t i = 0; fTrk && i < fTrk->GetNhits(); i++) {
      TrRecHitR *hit = fTrk->GetHit(i);
      if (!hit) continue;

      TrMCClusterR *mc = 0;
      Double_t    dmin = 200e-4;
      for (Int_t j = 0; j < nTrMCCluster(); j++) {
	TrMCClusterR *m = pTrMCCluster(j);
	if (!m || m->GetTkId() != hit->GetTkId()) continue;

	Double_t d = hit->GetCoord().y()-m->GetXgl().y();
	if (TMath::Abs(d) < TMath::Abs(dmin)) {
	  mc   = m;
	  dmin = d;
	}
      }
      if (mc) {
	Int_t lj = hit->GetLayerJ()-1;
	fDst.fMCinfo.trcx[lj] = mc->GetXgl().x();
	fDst.fMCinfo.trcy[lj] = mc->GetXgl().y();
	fDst.fMCinfo.trcz[lj] = mc->GetXgl().z();

	if (lj == 0) {
	  fDst.fMCinfo.dir1[0] = mc->GetDir().x();
	  fDst.fMCinfo.dir1[1] = mc->GetDir().y();
	  fDst.fMCinfo.dir1[2] = mc->GetDir().z();
	}
	if (lj == 8) {
	  fDst.fMCinfo.dir9[0] = mc->GetDir().x();
	  fDst.fMCinfo.dir9[1] = mc->GetDir().y();
	  fDst.fMCinfo.dir9[2] = mc->GetDir().z();
	}
      }
    }
  }

  ParticleR *part = pParticle(0);

  Int_t    imin[9];
  Double_t dmin[9], trky[9];
  for (Int_t i = 0; part && fTrk && i < 9; i++) {
    imin[i] = -1;
    dmin[i] =  0;
    if      (i == 0) trky[i] = fDst.fTrack.cooy[0];
    else if (i == 1) trky[i] = fDst.fTrack.cooy[1];
    else if (i == 8) trky[i] = fDst.fTrack.cooy[2];
    else             trky[i] = part->TrCoo[i-1][1];
  }

  Double_t qmin = (fDst.fTrack.qin > 1.5) ? 1.5 : 0;

  for (Int_t i = 0; fTrk && i< nTrRecHit(); i++) {
    TrRecHitR *hit = pTrRecHit(i);
    if (!hit || hit->OnlyY()) continue;

    Int_t lj = hit->GetLayerJ()-1;
    TrRecHitR *ht = fTrk->GetHitLJ(lj+1);
    if (ht && ht->GetYClusterIndex() == hit->GetYClusterIndex()) continue;

    if (qmin > 0) {
      Float_t qhit = TMath::Sqrt(hit->GetSignalCombination(2, qopt, 1, 1));
      if (qhit < qmin) continue;
    }

    Double_t d = TMath::Abs(hit->GetCoord().y()-trky[lj]);
    if (dmin[lj] == 0 || d < dmin[lj]) {
      dmin[lj] = d;
      imin[lj] = i;
    }
  }

  for (Int_t i = 0; fTrk && i < 9; i++) {
    fDst.fTrCls.mdist[i] = dmin[i];
    if (fDst.fTrCls.tkml[i] == 0 && imin[i] >= 0) {
      TrRecHitR *hit = pTrRecHit(imin[i]);
      if (!hit || hit->GetLayerJ() != i+1) continue;

      TrClusterR *cly = hit->GetYCluster();
      if (!cly) continue;

      Int_t tkid = hit->GetTkId();
      fDst.fTrCls.tkml[i] = TMath::Sign(100000+TMath::Abs(tkid), tkid);
      fDst.fTrCls.xcog[i] = 0;
      fDst.fTrCls.ycog[i] = cly->GetCofG()+cly->GetSeedAddress();
    }
  }
  if(!fTrk) return;
  qopt = TrClusterR::kAsym|TrClusterR::kGain|TrClusterR::kAngle|
    TrClusterR::kLoss|TrClusterR::kMIP;
  if(pMCEventg(0) && Version()>=1106) 
    qopt = TrClusterR::kTotSign2017|TrClusterR::kSimAsym|TrClusterR::kSimSignal|
      TrClusterR::kLoss|TrClusterR::kAngle;
  for(int i=0;i<9;i++){
      TrRecHitR *hit = fTrk->GetHitLJ(i+1);
      if(!hit) continue;
      bool gx = TrCharge::GoodChargeReconHit(hit, 0);
      bool gy = TrCharge::GoodChargeReconHit(hit, 1);
      if(gx) fDst.fTrHit.etrx[i] = hit->GetSignalCombination(0, qopt, 1, 0, 0);
      if(gy) fDst.fTrHit.etry[i] = hit->GetSignalCombination(1, qopt, 1, 0, 0);
  }
}

void MdSel::FillBeta()
{
  ParticleR *part = pParticle(0);
  if (!part) return;

  if (part->pBeta()) {
    fDst.fBeta.pattern = part->pBeta()->Pattern;
    fDst.fBeta.beta    = part->pBeta()->Beta;
  }

  for (Int_t i = 0; i < 4; i++) {
    fDst.fTofHit.nhit[i] = 0;
    fDst.fTofHit.qsum[i] = fDst.fTofHit.qmax[i] = 0;
    fDst.fTofHit.tmin[i] = fDst.fTofHit.tmax[i] = 0;
  }

  for (Int_t i = 0; i < nTofClusterH(); i++) {
    TofClusterHR *tcls = pTofClusterH(i);
    if (!tcls) continue;

    Int_t il = tcls->Layer;
    if (il < 0 || 4 <= il) continue;

    Float_t q = tcls->GetQSignal();
    Float_t t = tcls->Time;

    if (fDst.fTofHit.nhit[il] == 0)
      fDst.fTofHit.tmin[il] = fDst.fTofHit.tmax[il] = t;

    fDst.fTofHit.nhit[il]++;
    fDst.fTofHit.qsum[il] += q;
    if (q > fDst.fTofHit.qmax[il]) fDst.fTofHit.qmax[il] = q;
    if (t > fDst.fTofHit.tmax[il]) fDst.fTofHit.tmax[il] = t;
    if (t < fDst.fTofHit.tmin[il]) fDst.fTofHit.tmin[il] = t;
  }

  Bool_t ismc = (pMCEventg(0)) ? kTRUE : kFALSE;
  BetaHR *bth = (ismc) ? part->pBetaH() : 0;
  for (Int_t i = 0; !bth && i < NBetaH(); i++)
    if (pBetaH(i) && fTrk && pBetaH(i)->pTrTrack() == fTrk) bth = pBetaH(i);

  Double_t rgt = (fTrk) ? TMath::Abs(fTrk->GetRigidity()) : 1;

  if (bth) {
    Int_t   nlay;
    Float_t qrms;
    fDst.fBetaH.pattern = bth->GetBetaPattern();
    fDst.fBetaH.beta    = bth->GetBeta();
    fDst.fBetaH.chi2t   = bth->GetNormChi2T();
    fDst.fBetaH.chi2c   = bth->GetNormChi2C();
    fDst.fBetaH.q       = bth->GetQ(nlay, qrms);

    fDst.fBetaH.pbit    = ((bth->pTrTrack   ()) ? 1 : 0)+
                          ((bth->pTrdTrack  ()) ? 2 : 0)+
                          ((bth->pEcalShower()) ? 4 : 0);

    for (Int_t i = 0; i < 4; i++)
      fDst.fBetaH.ql[i] = bth->GetQL(i);

    //GetClsn(this, bth);
  //GetNTofClustersInTime(bth, fDst.fEstim.clsn);
    GetNTofClustersInTime(bth, fDst.fBetaH.clsn);

    if (CutTk2(this, bth)) fDst.fStatus.ustat |= CutTk2nd;

    fDst.fBetaH.z[0] = bth->GetZ(nlay, fDst.fBetaH.p[0], 0);
    fDst.fBetaH.z[1] = bth->GetZ(nlay, fDst.fBetaH.p[1], 1);
    fDst.fBetaH.z[2] = bth->GetZ(nlay, fDst.fBetaH.p[2], 2);

    fDst.fBetaH.type = bth->GetBuildType();

    Float_t rms;
    Int_t n, opt = TofClusterHR::DefaultQOptIonW;
    fDst.fBetaH.qup  = bth->GetQ(n, rms, 2, opt, 1100, 0, rgt);
    fDst.fBetaH.qlow = bth->GetQ(n, rms, 2, opt,   11, 0, rgt);
  }

  TofRecH::BuildOpt = 1;

  TofRecH tofrec;
  tofrec.ReBuild(0);

  TofRecH::BuildOpt = 0;

  bth = pBetaH(0);
  if (bth) {
    Int_t   nlay;
    Float_t qrms;
    fDst.fBetaHs.pattern = bth->GetBetaPattern();
    fDst.fBetaHs.beta    = bth->GetBeta();
    fDst.fBetaHs.chi2t   = bth->GetNormChi2T();
    fDst.fBetaHs.chi2c   = bth->GetNormChi2C();
    fDst.fBetaHs.q       = bth->GetQ(nlay, qrms);

    fDst.fBetaHs.pbit    = ((bth->pTrTrack   ()) ? 1 : 0)+
                           ((bth->pTrdTrack  ()) ? 2 : 0)+
                           ((bth->pEcalShower()) ? 4 : 0);

    for (Int_t i = 0; i < 4; i++)
      fDst.fBetaHs.ql[i] = bth->GetQL(i);
  }
}

void MdSel::FillTrd()
{
  ParticleR *part = pParticle(0);
  if (!part) return;

  fDst.fTrd.nsegt = NtrdSg(this);

  TrdTrackR *trd = part->pTrdTrack();
  if (trd) {
    fDst.fTrd.coo[0] = trd->Coo[0];
    fDst.fTrd.coo[1] = trd->Coo[1];
    fDst.fTrd.coo[2] = trd->Coo[2];
    fDst.fTrd.theta  = trd->Theta;
    fDst.fTrd.phi    = trd->Phi;
    fDst.fTrd.q      = trd->Q;

    AMSPoint ptrd(trd->Coo[0], trd->Coo[1], trd->Coo[2]);
    AMSDir   dtrd(trd->Theta, trd->Phi);

    AMSPoint p, d;
    Int_t tk;
    Int_t is1 = IsInsideTracker(1, ptrd, dtrd, 0, 0, tk, p, d);
    Int_t is2 = IsInsideTracker(2, ptrd, dtrd, 0, 0, tk, p, d);
    Int_t is5 = IsInsideTracker(5, ptrd, dtrd, 0, 0, tk, p, d);
    Int_t is9 = IsInsideTracker(9, ptrd, dtrd, 0, 0, tk, p, d);
    if (is1)        fDst.fStatus.ustat |= TrdInL1;
    if (is1 && is9) fDst.fStatus.ustat |= TrdInFS;
    if (is2 && is5) fDst.fStatus.ustat |= TrdInTr;

    if (fTrk) {
      AMSPoint pnt;
      AMSDir   dir;
      fTrk->Interpolate(trd->Coo[2], pnt, dir);
      if (dtrd.z()*dir.z() < 0) dtrd = dtrd*(-1);

      fDst.fTrd.dtrk[0] = pnt.x()-trd->Coo[0];
      fDst.fTrd.dtrk[1] = pnt.y()-trd->Coo[1];
      fDst.fTrd.dtrk[2] = TMath::ACos(dir.prod(dtrd))*TMath::RadToDeg();
    }

    Double_t dxt = 10;
    Double_t dyt = 10;
    Int_t   ihit = -1;

    for (Int_t i = 0; i < nTrRecHit(); i++) {
      TrRecHitR *hit = pTrRecHit(i);
      if (!hit || hit->GetLayerJ() != 1 || hit->OnlyY()) continue;

      AMSPoint phit = hit->GetCoord();
      Double_t dy   = ptrd.y()+dtrd.y()/dtrd.z()*(phit.z()-ptrd.z())-phit.y();
      if (TMath::Abs(dy) > TMath::Abs(dyt)) continue;

      Double_t xmin = 4.2;
      for (Int_t k = 0; k < hit->GetMultiplicity(); k++) {
	phit = hit->GetCoord(k);
	Double_t dx = ptrd.x()+dtrd.x()/dtrd.z()*(phit.z()-ptrd.z())-phit.x();
	if (TMath::Abs(dx) < TMath::Abs(xmin)) xmin = dx;
      }
      dyt  = dy;
      dxt  = xmin;
      ihit = i;
    }

    TrRecHitR *hit = (ihit >= 0) ? pTrRecHit(ihit) : 0;
    if (hit) {
      int qopt = TrClusterR::kAsym | TrClusterR::kGain |
	         TrClusterR::kLoss | TrClusterR::kMIP;
      fDst.fTrd.ql1m = TMath::Sqrt(hit->GetSignalCombination(2, qopt, 1, 1));
      fDst.fTrd.dl1m[0] = dxt;
      fDst.fTrd.dl1m[1] = dyt;
    }
  }

//if (!fTrk || pMCEventg(0)) return;
//if (!fTrk) return;                   // 161202

  if (!pMCEventg(0)) {
    if (Run() < 1305800000) TrdKCluster::IsReadGlobalAlignment = 0;
    if (Run() >= RUN_STD) return;
  }

  if (pMCEventg(0)) {
    TrdKCluster::ForceReadAlignment=0;
    TrdKCluster::ForceReadCalibration=0;
    TrdKCluster::ForceReadXePressure=0;
    TrdKCluster::SetDefaultMCXePressure(900);
  }

  // 161202
  if (part->pTrdTrack() && part->pEcalShower()) {
    TrdKCluster tkcl(this, part->pTrdTrack());

    Int_t    nhit = 0;
    Double_t llr[3] = { -1, -1, -1 }, ll[3];

    if (tkcl.IsReadAlignmentOK == 2 && tkcl.IsReadCalibOK == 1) {
      EcalShowerR *ec = part->pEcalShower();
      if (ec) {
	tkcl.GetLikelihoodRatio_TrTrack(15, llr, nhit, ec->EnergyE, ll);
	fDst.fTrdK.nhitt   = nhit;
	fDst.fTrdK.llrt[0] = llr[0];
	fDst.fTrdK.llrt[1] = llr[1];
	fDst.fTrdK.llrt[2] = llr[2];
      }
    }
  }

  if (!fTrk) return;
  // 161202

  TrdKCluster tkcl(this, fTrk, fTrk->Gettrdefaultfit());
  Int_t    nhit = 0;
  Double_t llr[3] = { -1, -1, -1 }, ll[3];

  if (tkcl.IsReadAlignmentOK == 2 && tkcl.IsReadCalibOK == 1) {
    tkcl.GetLikelihoodRatio_TrTrack(15, &llr[0], nhit);
    fDst.fTrdK.nhits  = nhit;
    fDst.fTrdK.llr[0] = llr[0];
    fDst.fTrdK.llr[1] = llr[1];
    fDst.fTrdK.llr[2] = llr[2];

    if (tkcl.CalculateTRDCharge(0, fBeta) == 0) 
      fDst.fTrdK.q = tkcl.GetTRDCharge();

    EcalShowerR *ec = part->pEcalShower();
    if (ec) {
      tkcl.GetLikelihoodRatio_TrTrack(15, &llr[0], nhit, ec->EnergyE, ll);
      fDst.fTrdK.llre[0] = llr[0];
      fDst.fTrdK.llre[1] = llr[1];
      fDst.fTrdK.llre[2] = llr[2];
    }
  }

//if (fDst.fVertex.coo[4] > 0 || fDst.fTrack.qin > 1.5) {
//if (fDst.fVertex.coo[4] != 0 || fDst.fTrack.qin > 1.5) {  // Bugfix 160410

  Int_t tbid[25];
  Double_t de[25], tr[25];

  for (Int_t i = 0; i < 25; i++) tbid[i] = 0;
  for (Int_t i = 0; i < 25; i++) de[i] = tr[i] = 0;

  MCEventgR *mcg = pMCEventg(0);
  if (1 || fDst.fVertex.coo[4] != 0 || fDst.fTrack.qin > 1.5 ||
      (mcg && mcg->Particle == 3)) {                         // 160812

    AMSPoint pnt = tkcl.GetPropogated_TrTrack_P0();
    AMSDir   dir = tkcl.GetPropogated_TrTrack_Dir();
/*
    AMSPoint pmc = pnt;
    AMSDir   dmc = dir;
    if (mcg) {
      pmc.setp(mcg->Coo[0], mcg->Coo[1], mcg->Coo[2]);
      dmc.setp(mcg->Dir[0], mcg->Dir[1], mcg->Dir[2]);
    }
*/
    for (Int_t i = 0, j = 0; i < tkcl.NHits() && j < 25; i++) {
      TrdKHit *hit = tkcl.GetHit(i);
      if (hit) {
	Double_t len = hit->Tube_Track_3DLength(&pnt, &dir);
	Double_t amp = hit->TRDHit_Amp;

	tbid[j] = hit->TRDHit_Layer*10000+
	          hit->TRDHit_Ladder*100+
	          hit->TRDHit_Tube;

	if (len > 0 && amp > 0) {
	  fDst.fTrdHit.layer[j] = hit->TRDHit_Layer;
	  fDst.fTrdHit.plen [j] = len;
	  fDst.fTrdHit.amp  [j] = amp;
//	 if (mcg)
//	  fDst.fTrdHit.plmc [j] = hit->Tube_Track_3DLength(&pmc, &dmc);

	  j++;
	}
      }
    }
  }

  if (mcg) {
    for (Int_t i = 0; i < nTrdMCCluster(); i++) {
      TrdMCClusterR *c = pTrdMCCluster(i);

      for (Int_t j = 0; j < 25; j++)
	if (tbid[j] == c->Layer*10000+c->Ladder*100+c->Tube) {
	  if (c->ProcID == 0 && c->ParticleNo == mcg->Particle) {
	                   de[j] += c->Edep*1e6;
	    fDst.fTrdHit.plmc[j] += c->Step;
	  }
	  if (c->ProcID > 0 && c->ParticleNo == -1 && c->Ekin == 0)
	    tr[j] += c->Edep*1e6;
	  break;
	}
    }
    for (Int_t i = 0; i < 25; i++)
      fDst.fTrdHit.trfr[i] = (de[i]+tr[i] > 0) ? tr[i]/(de[i]+tr[i]) : 0;
  }


  if (fDst.fTrack.rgt[0] > 50 || fDst.fTrack.rgt[0] < -5) {
    GammaFit gf(fTrk);
    Double_t lk = -20;

    fDst.fTrdG.nfit  = gf.GetNhit();
    fDst.fTrdG.gamma = gf.Fit(lk, GammaFit::kProton, GammaFit::kTRDOnly);
    fDst.fTrdG.lkh   = lk;
  }
}

void MdSel::FillRich()
{
  

  ParticleR *part = pParticle(0);
  if (!part) return;

  RichRingR *rich = part->pRichRing();
  if (!rich) return;

  fDst.fRich.status   = rich->Status;
  
  fDst.fRich.nhits[0] = rich->getUsedHits(); //Used Hits in the Ring
  fDst.fRich.nhits[1] = rich->getHits();     //Numbers of hits in the Ring
  fDst.fRich.nhits[2] = rich->getReflectedHits();
  fDst.fRich.beta     = rich->IsClean() ?  rich->getBeta()
              	                        : -rich->getBeta();
  fDst.fRich.q        = rich->getCharge2Estimate();
  fDst.fRich.dist     = rich->DistanceTileBorder();
  fDst.fRich.prob     = rich->getProb();
  fDst.fRich.coo[0]   = rich->getTrackEmissionPoint()[0];
  fDst.fRich.coo[1]   = rich->getTrackEmissionPoint()[1];
  fDst.fRich.coo[2]   = rich->getTrackEmissionPoint()[2];

  fDst.fRich.npe[0]   = rich->getPhotoelectrons();
  fDst.fRich.npe[1]   = rich->getExpectedPhotoElectrons();
  fDst.fRich.npe[2]   = RichHitR::getCollectedPhotoElectrons();

  fDst.fRich.tile     = rich->getTileIndex();
  fDst.fRich.pmts     = rich->getPMTs();

  for(int i=0;i<rich->getHits() && i<30;i++){
    RichHitR *hit = pRichHit(rich->iRichHit(i));
    if(!hit) continue;
    fDst.fRich.nhpe[i] =  hit->Npe;
  }
}

void MdSel::FillEcal()
{
  ParticleR *part = pParticle(0);
  if (!part) return;

  EcalShowerR *ecal = part->pEcalShower();

  if (!ecal) {
    for (Int_t i = 0; i < nEcalShower(); i++) {
      EcalShowerR *ec = pEcalShower(i);
      if (!ecal || ec->EnergyE > ecal->EnergyE) ecal = ec;
    }
  }
  if (!ecal) return;

  fDst.fEcal.cog[0] = ecal->CofG[0];
  fDst.fEcal.cog[1] = ecal->CofG[1];
  fDst.fEcal.cog[2] = ecal->CofG[2];
  fDst.fEcal.dir[0] = ecal->Dir[0];
  fDst.fEcal.dir[1] = ecal->Dir[1];
  fDst.fEcal.dir[2] = ecal->Dir[2];
  fDst.fEcal.s13r   = ecal->S13R;
  fDst.fEcal.bdt    = ecal->GetEcalBDT(4);

  fDst.fEcal.energy[0] = ecal->EnergyD;
  fDst.fEcal.energy[1] = ecal->EnergyC;
  fDst.fEcal.energy[2] = ecal->EnergyE;

  fDst.fEcal.enew[0] = //ecal->GetCorrectedEnergy (2, 2);
  fDst.fEcal.enew[1] = ecal->GetCorrectedEnergy (2, 2);
                     //ecal->GetCorrectedEnergy2(2, 2);

  fDst.fEcal.dz[0] = ecal->DirCR[2];
  fDst.fEcal.dz[1] = ecal->EMDir[2];
  const Int_t CATL = 16384*2*2*2*2*2*2*2*2*2*2*2;
  fDst.fEcal.catl = (ecal->Status&CATL) ? 1 : 0;

  AMSPoint pec(ecal->CofG[0], ecal->CofG[1], ecal->CofG[2]);
  AMSDir   dec(ecal->Dir [0], ecal->Dir [1], ecal->Dir [2]);

  if (fTrk) {
    AMSPoint pnt;
    AMSDir   dir;
    fTrk->Interpolate(ecal->CofG[2], pnt, dir);
    if (dec.z()*dir.z() < 0) dec = dec*(-1);

    fDst.fEcal.dtrk[0] = pnt.x()-ecal->CofG[0];
    fDst.fEcal.dtrk[1] = pnt.y()-ecal->CofG[1];
    fDst.fEcal.dtrk[2] = TMath::ACos(dir.prod(dec))*TMath::RadToDeg();
  }

  Double_t dxt = 10;
  Double_t dyt = 10;
  Int_t   ihit = -1;

  for (Int_t i = 0; i < nTrRecHit(); i++) {
    TrRecHitR *hit = pTrRecHit(i);
    if (!hit || hit->GetLayerJ() != 9 || hit->OnlyY()) continue;

    AMSPoint phit = hit->GetCoord();
    Double_t dy   = pec.y()+dec.y()/dec.z()*(phit.z()-pec.z())-phit.y();
    if (TMath::Abs(dy) > TMath::Abs(dyt)) continue;

    Double_t xmin = 4.2;
    for (Int_t k = 0; k < hit->GetMultiplicity(); k++) {
      phit = hit->GetCoord(k);
      Double_t dx = pec.x()+dec.x()/dec.z()*(phit.z()-pec.z())-phit.x();
      if (TMath::Abs(dx) < TMath::Abs(xmin)) xmin = dx;
    }
    dyt  = dy;
    dxt  = xmin;
    ihit = i;
  }

  TrRecHitR *hit = (ihit >= 0) ? pTrRecHit(ihit) : 0;
  if (hit) {
    int qopt = TrClusterR::kAsym | TrClusterR::kGain |
               TrClusterR::kLoss | TrClusterR::kMIP;
    fDst.fEcal.ql9m = TMath::Sqrt(hit->GetSignalCombination(2, qopt, 1, 1));
    fDst.fEcal.dl9m[0] = dxt;
    fDst.fEcal.dl9m[1] = dyt;
  }

  Int_t ztrk = (fDst.fTrack.qin < 1.5) ? 1 : 2;

  Bool_t ismc = (pMCEventg(0)) ? kTRUE : kFALSE;
  EcalHR etmp;
  EcalHR *ech = pEcalH(0);
  if (!ech) { ech = &etmp; ech->Process(); }

  for (Int_t i = 0; i < ech->Nhit(); i++) {
    Int_t l = ech->Plane(i);
    if (0 <= l && l < 10) fDst.fEcalH.edep[l] += ech->dE(i);
  }

  AMSPoint pp, pl;
  Float_t emip = EcalHR::GetMipEdep(pp, pl);
  Float_t rrec = -1;
  Float_t smax = -10;
  Float_t ecsq = -1;
  if (!ismc) EcalHR::Get(ztrk, rrec, smax, ecsq);

  fDst.fEcalH.csq  = ecsq;
  fDst.fEcalH.rrec = rrec;
  fDst.fEcalH.smax = smax;
  fDst.fEcalH.emip = emip;
  fDst.fEcalH.lmip = TMath::Abs((pp-pl).norm());

//if (fDst.fEcal.bdt > 0 && pEcalH(0)) *fEcalH = *pEcalH(0);
}

void MdSel::FillMC()
{
  MCEventgR *mcg = pMCEventg(0);
  if (!mcg) return;

  fDst.fMCinfo.coo[0] = mcg->Coo[0];
  fDst.fMCinfo.coo[1] = mcg->Coo[1];
  fDst.fMCinfo.coo[2] = mcg->Coo[2];
  fDst.fMCinfo.dir[0] = mcg->Dir[0];
  fDst.fMCinfo.dir[1] = mcg->Dir[1];
  fDst.fMCinfo.dir[2] = mcg->Dir[2];
  fDst.fMCinfo.rgt    = (mcg->Charge != 0) ? mcg->Momentum/mcg->Charge
                                           : mcg->Momentum;

  fDst.fMCinfo.rseed[0] = fHeader.RNDMSeed[0];
  fDst.fMCinfo.rseed[1] = fHeader.RNDMSeed[1];

  vector<double>  z;
  map<Int_t,Double_t> trkid;
  Int_t    pid   = mcg->Particle;
  Int_t    ptkid = mcg->trkID;
  Double_t mom   = mcg->Momentum;
  Double_t mass  = mcg->Mass;
  Double_t charg = mcg->Charge;
  Double_t ekin  = sqrt(mom*mom+mass*mass)-mass;
  Double_t et    = 0;
  Double_t zmax  = -200;
  Double_t zmin  = 200;
  Double_t xmin  = 0;
  Double_t ymin  = 0;
  Double_t emax  = 0;

  for (Int_t k = 1; k < nMCEventg(); k++) {
    MCEventgR &sec = MCEventg(k);
    if (sec.parentID == ptkid && trkid.find(sec.trkID) == trkid.end() &&
	sec.Momentum > 0.05) { // && charg != sec.Charge) {
      trkid[sec.trkID] = sec.Momentum;
      Double_t ek = TMath::Sqrt(sec.Momentum*sec.Momentum+
				sec.Mass*sec.Mass)-sec.Mass;
      if (ek > emax) emax = ek;
      et += ek;
      if (et > 0.05) { //ekin*0.01) {
	if (sec.Coo[2] < zmin) {
	  z.push_back(sec.Coo[2]);
	  zmin = sec.Coo[2];
	  xmin = sec.Coo[0];
	  ymin = sec.Coo[1];
	}
	if (sec.Coo[2] > zmax) zmax = sec.Coo[2];
      }
    }
  }
  fDst.fMCinfo.mcin[0] = et/ekin;
  fDst.fMCinfo.mcin[1] = z.size() ? z[z.size()-1] : -200;
  fDst.fMCinfo.mci2[0] = xmin;
  fDst.fMCinfo.mci2[1] = ymin;
  fDst.fMCinfo.mci2[2] = zmax-zmin;
  fDst.fMCinfo.mci2[3] = emax;

  if (Version() >= 800) {
    MCEventgR *pex = GetPrimaryMC(-1);
    if (pex) {
      fDst.fMCinfo.mci2[4] = pex->Coo[0];
      fDst.fMCinfo.mci2[5] = pex->Coo[1];
      fDst.fMCinfo.mci2[6] = pex->Coo[2];
    }
  }

  enum { NP = 36 };
  Int_t part[NP] = {  1,  2,  3,  5,  6,  7,  8,  9, 11, 12, 13, 14,
		     15, 17, 18, 21, 25, 45, 46, 47, 49, 61, 62, 63,
		     64, 65, 66, 67, 68, 69, 71, 72, 74, 84, 87, 114 };
  Int_t chrg[NP] = {  0,  1, -1,  1, -1,  0,  1, -1,  1, -1,  0,  1,
		     -1,  0,  0, -1,  0,  1,  1,  2,  2,  3,  3,  4,
		      4,  5,  5,  6,  7,  8, 10, 11, 13, 23, 26,  4 };
              // L: 2  3  4  5  6  7  8  1  9
  Int_t ilay[9] = { 1, 2, 2, 2, 2, 3, 3, 0, 4 };

  for (Int_t i = 0; i < nTrMCCluster(); i++) {
    TrMCClusterR *mc = pTrMCCluster(i);
    if (!mc) continue;

    Int_t chg = 0;
    for (Int_t j = 0; j < NP; j++)
      if (TMath::Abs(mc->GetPart()) == part[j]) { chg = chrg[j]; break; }

    if (chg == 0) continue;

    Int_t lay = TMath::Abs(mc->GetTkId())/100;
    if (lay <= 0 || 9 < lay) continue;

    Int_t il = ilay[lay-1];

    Double_t r = mc->GetMomentum()/chg;
    //TODO Add position within 2cm around TrCluster?
    if (TMath::Abs(r) > TMath::Abs(fDst.fMCinfo.trgt[il])) {
      fDst.fMCinfo.trgt[il] = r;
      fDst.fMCinfo.part[il] = mc->GetPart()*((!mc->IsPrimary()) ? -1 : 1);
    }
  }

  // by LD
  double RemoveInteractionFactor = .8;

  fDst.fMCinfo.patmc = 0;
  fDst.fMCinfo.npmc  = 0;
  for (int i = 0; i < nTrMCCluster(); i++) {
    TrMCClusterR *mc = pTrMCCluster(i);
    int Layer = TkDBc::Head->GetJFromLayer(abs(mc->GetTkId())/100);
    if (mc->GetMomentum() > mom*RemoveInteractionFactor)
      fDst.fMCinfo.patmc |= 1<<(Layer-1);
  }
  Int_t n = fDst.fMCinfo.patmc;
  while (n) { fDst.fMCinfo.npmc += (n&1); n >>= 1; }
}

Float_t GetBeta(AMSEventR* ev, Int_t flag, TrTrackR* trk = 0){
  ParticleR *p = ev->pParticle(0);

  if (!trk)  trk = (p) ? p->pTrTrack () : 0;
  TrdTrackR *trd = (p) ? p->pTrdTrack() : 0;
  BetaHR    *bth = (p) ? p->pBetaH()    : 0;
  if(ev->pMCEventg(0) && !bth) bth = ev->pBetaH(0);

  double bm = 0;
  int    nb = 0;

  if (flag%10 && trk) {

    int qopt = TrClusterR::kAsym|TrClusterR::kGain|TrClusterR::kAngle|
      TrClusterR::kLoss|TrClusterR::kMIP;
    if(ev->pMCEventg(0) && ev->Version()>=1106) 
      qopt = TrClusterR::kTotSign2017|TrClusterR::kSimAsym|TrClusterR::kSimSignal|
        TrClusterR::kLoss|TrClusterR::kAngle;

    mean_t mn = TrCharge::GetCombinedMean(
        TrCharge::kInner|TrCharge::kTruncMean|TrCharge::kSqrt, trk, 1,-1,
        qopt , -1, 0);
    //TrClusterR::kBeta|TrClusterR::kRigidity

    bm += 0.94*TMath::Power(1/mn.Mean, 1.1);
    nb++;
  }

  if ((flag/10)%10 && bth) {
    int nlay;
    float qrms;
    double q = bth->GetQ(nlay, qrms, 2, 
        TofRecH::kThetaCor|TofRecH::kBirkCor|
        TofRecH::kReAttCor|TofRecH::kDAWeight|TofRecH::kQ2Q,
        //TofRecH::kBetaCor
        -2, 1, 0);

    bm += 0.99*TMath::Power(1/q, 1.22);
    nb++;
  }

  if ((flag/100)%10 && trd && trk) {
    TrdKCluster tk(ev, trd, trk->GetRigidity());

    AMSPoint p0(trd->Coo);
    AMSDir  dir(trd->Theta, trd->Phi);

    double qs = 0, qmin = 1e9, qmax = 0;
    int    ns = 0;

    for (int j = 0; j < tk.NHits(); j++) {
      TrdKHit *hit = tk.GetHit(j);
      if (!hit) continue;

      double len = hit->Tube_Track_3DLength(&p0, &dir);
      double amp = hit->TRDHit_Amp;
      if (len > 0.3) {
        double q = amp/len;
        qs += q;
        ns++;
        if (q > qmax) qmax = q;
        if (q < qmin) qmin = q;
      }
    }
    if (ns > 5) {
      qs -= qmax+qmin;
      qs /= (ns-2);
      bm += TMath::Sqrt(160/qs);
      nb++;
    }
  }

  return (nb > 0) ? bm/nb : 0;

}
void MdSel::FillMass()
{
  Double_t rgt = fDst.fTrack.rgt2[0]; if (rgt == 0) rgt = fDst.fTrack.rgt[0];
  Double_t ar  = TMath::Abs(rgt);
  Double_t ab  = TMath::Abs(fBeta);

  fDst.fMass.b1  = GetBeta (this,   1, fTrk);
  fDst.fMass.b2  = GetBeta (this,  10);
  fDst.fMass.b3  = GetBeta (this, 100);
  fDst.fMass.npk = TrMass::GetNpick(this, fTrk);

  fDst.fMass.m   = (0 < ab && ab < 1) ? TrMass::GetMass(1, ar, ab) : 0;

  fDst.fMass.mql = TrMass::GetMQL(this, fTrk, 0, fDst.fMass.v, fDst.fMass.p);

  AMSPoint hc[9];
  fDst.fMass.nh = TrMass::GetH(&fDst.fTrCls.tkml[1], &fDst.fTrCls.xcog[1], 
      &fDst.fTrCls.ycog[1], hc,
      fDst.fTrack.theta, fDst.fTrack.phi, 1, 0);


  ParticleR *part = pParticle(0);
  BetaHR *bth = NULL;
  if(part) bth = pMCEventg(0)? part->pBetaH():0;
  for(Int_t i =0;!bth && i<NBetaH(); i++)
    if(pBetaH(i) && fTrk && pBetaH(i)->pTrTrack()==fTrk) bth=pBetaH(i);
  if(!bth) bth = pBetaH(0);
  Int_t nlay; Float_t p;
  
  Int_t bz = bth?bth->GetZ(nlay, p, 0):0;
  if (bz == 1 && ab < 1) {
    Int_t s  = (fDst.fTrack.rgt[0] > 0) ? 1 : -1;
    Int_t a1 = 1, a2 = 2;
    fDst.fMass.ll[0] = TrMass::GetLL(bz, s*a1, ab, hc, fDst.fTrack.theta);
    fDst.fMass.ll[1] = TrMass::GetLL(bz, s*a2, ab, hc, fDst.fTrack.theta);
  }

  AMSPoint ptr(fDst.fTrack.p0x,   fDst.fTrack.p0y, 0);
  AMSDir   dtr(fDst.fTrack.theta, fDst.fTrack.phi);

  fDst.fMass.bl2 = (TMath::Abs(fDst.fTrack.rgt[0]) > 0.01)
                 ? TrMass::GetBL2(ptr, dtr, fDst.fTrack.rgt[0]) : 0; 
  fDst.fMass.ru  = TrMass::Rndm(fDst.fHeader.run);

  if(!bth)return;
  TofChargeHR *tofc = bth->pTofChargeH();
  fDst.fMass.zl[0] = tofc->GetZ(nlay,fDst.fMass.pz[0],0,1000);
  fDst.fMass.zl[1] = tofc->GetZ(nlay,fDst.fMass.pz[1],0,100); 
  fDst.fMass.zl[2] = tofc->GetZ(nlay,fDst.fMass.pz[2],0,10);  
  fDst.fMass.zl[3] = tofc->GetZ(nlay,fDst.fMass.pz[3],0,1);   
  fDst.fMass.betas = bth->GetBetaS();
}


void MdSel::GetClsn(AMSEventR *pev, BetaHR *bh)
{
  for (Int_t i = 0; i < 4; i++) fDst.fEstim.clsn[i] = 0;

      //-->look around(time) used hits:
      for(int il=0;il<4;il++){//layer loop
	int itb=(il<2)?0:1;
	if(bh->TestExistHL(il)){//hit exists
	//barn=bh->GetClusterHL(il)->Bar;//0-9
	  float ltime=bh->GetTime(il);//ns
	  int ntfcls=pev->nTofClusterH();
	  for(int icl=0;icl<ntfcls;icl++){
	    TofClusterHR &tfcl=pev->TofClusterH(icl);
	    if(tfcl.Layer!=il)continue;
	    if(tfcl.NBetaHUsed()>0)continue;//used by BetaHtfcl.Time();
	    if(!tfcl.IsGoodTime())continue;
	    float cltime=tfcl.Time;//ns
	    float ed=tfcl.GetEdep();//mev
	    float dt=cltime-ltime;//later cluster has positive dt
	    //	       if(itb==0)prsh10->Fill(dt,1.);
	    //	       if(itb==1)prsh11->Fill(dt,1.);
	    float tcut(0);
	    if(itb==0)tcut=10;//ns,top
	    else tcut=4;//ns, bot
	    int itm(-1);
	    if(fabs(dt)<=tcut)itm=0;//around-time hits ("in time")
	    if(dt>tcut)itm=1;//later hits ("off time")
	    if(itm>=0){
	      int indx=itm+2*itb;
	      if(indx>=0 && indx<4){
		fDst.fEstim.clsn[indx]++;
	      }
	    }
	  }//-->endof "secondary clust-loop"
	}//-->endof "beta-hit exists"
      }//-->endof "Tof layer loop"
}

Bool_t MdSel::CutTk2(AMSEventR *pev, BetaHR *betah)
{
///--2ndTK Cut
    bool cut2ndtk=0;

    if (!pev->pParticle(0)) return false;
    Int_t nbetah=pev->NBetaH();
    Int_t ibetah=pev->pParticle(0)->iBetaH();
    Int_t itrtrack=pev->pParticle(0)->iTrTrack();

    if(nbetah>1){
      float tk2rig=0;
       float sbeta=betah->GetBeta();
       for(int ibh=0;ibh<nbetah;ibh++){
         if(ibh==ibetah)continue;//not same betah
         BetaHR *betah2=pev->pBetaH(ibh);
        if(betah2->iTrTrack()<0||betah2->iTrTrack()==itrtrack)continue;//not same TK
//---
         float nrig= betah2->pTrTrack()->GetRigidity();
         int ntkhb[2];
         ntkhb[0]=betah2->pTrTrack()->GetBitPatternXYJ();
         ntkhb[1]=betah2->pTrTrack()->GetBitPatternJ();
         int nhit2i[2]={0};
         for(int ilay=0+1;ilay<9-1;ilay++){
             for(int ixy=0;ixy<2;ixy++){
               if((ntkhb[ixy]&(1<<ilay))>0)nhit2i[ixy]++;
             }
         }
//----
         bool tkcuthit=(nhit2i[0]>=3&&nhit2i[1]>=5);
         bool rigcut=(fabs(nrig)>0.5);
         if(tkcuthit&&rigcut&&nrig/sbeta>0){cut2ndtk=1;break;}
//----
      }
   }
    return cut2ndtk;
}

Int_t MdSel::NtrdSg(AMSEventR *pev)
{
    ParticleR *part = pev->pParticle(0);
    if (!part) return 0;

    Int_t nseg = pev->nTrdSegment();
    if (nseg <= 0) return 0;

    Double_t zcen = 120;
    Double_t zmin = 120;
    Double_t zmax = 200;

    AMSPoint pt1; AMSDir dr1;
    AMSPoint pt2; AMSDir dr2;

    TrTrackR *ptrk = fTrk; //part->pTrTrack(); //pev->pTrTrack(fparm.itrk);
    if (!ptrk) return 0;

    ptrk->Interpolate(zmin, pt1, dr1);
    ptrk->Interpolate(zmax, pt2, dr2);

    Double_t tdx = (pt1.x()-pt2.x())/(pt1.z()-pt2.z());
    Double_t tdy = (pt1.y()-pt2.y())/(pt1.z()-pt2.z());
    Double_t tx0 =  pt1.x()-pt1.z()*tdx;
    Double_t ty0 =  pt1.y()-pt1.z()*tdy;
    Double_t txt = tx0+tdx*zcen;
    Double_t tyt = ty0+tdy*zcen;

    Int_t nsgp = 0;

    for (Int_t i = 0; i < nseg; i++) {
        TrdSegmentR *seg = pev->pTrdSegment(i);

        Int_t xy = seg->Orientation;
        if (xy == 0) {
            seg->FitPar[1] *= -1;
            seg->FitPar[0] *= -1;
        }

        Double_t x0 = seg->FitPar[1]+seg->FitPar[0]*zcen;
        Double_t dx = (xy == 1) ? txt-x0 : tyt-x0;

        if (TMath::Abs(dx) < 3) continue;

        Double_t zt = (xy == 1) ? -(tx0-seg->FitPar[1])/(tdx-seg->FitPar[0])
	                        : -(ty0-seg->FitPar[1])/(tdy-seg->FitPar[0]);

        if (zmin < zt && zt < zmax) nsgp++;
    }

    return nsgp;
}

void MdSel::FillHists()
{
  Bool_t ismc = (pMCEventg(0)) ? kTRUE : kFALSE;

  // Inner Tracker cahrge
  Double_t qin  = fDst.fTrack.qin;
  Int_t    qsel = 0;
  if (0.7 < qin && qin < 1.5) qsel = 1;
  if (1.7 < qin && qin < 2.5) qsel = 2;
  if (2.5 < qin)              qsel = 3;

  // TOF charge
  Double_t qtu = fDst.fBetaH.qup;
  Double_t qtl = fDst.fBetaH.qlow;
  if (qsel == 1 && (qtl < 0.5 || 2 < qtl)) qsel = 0;
  if (qsel == 2 &&  qtu < 1.25)            qsel = 0;

  // L1 and L9 charge
  Double_t qlm = (qsel == 1) ? ((ismc) ? 1.05 : 1.0) 
                             : TMath::Floor(qin+0.5);

  Double_t ql1 = fDst.fTrack.ql1;
  Double_t ql9 = fDst.fTrack.ql9;

  Bool_t qs1 = kFALSE;
  Bool_t qs9 = kFALSE;
  if (qlm-0.4 < ql1 && ql1 < qlm+0.9) qs1 = kTRUE;
  if (qlm-0.4 < ql9 && ql9 < qlm+0.9) qs9 = kTRUE;

  Int_t has19 = 0;
  if ((fDst.fTrack.bitx&0x001) && !(fDst.fTrack.bitr&0x3) && qs1) has19 += 1;
  if ((fDst.fTrack.bitx&0x100) && !(fDst.fTrack.bitr&0xc) && qs9) has19 += 2;

  fDst.fTrack.itp = has19;


  // Standard p/He flux selection

  // No TrigRun, Good RTI, No IsInSAA, PGMD diff.
  if (!ismc) {
    if ((fBadRun&4) || fBadRTI || fInSAA || fPGMD[0] > 35
	                                 || fPGMD[1] > 45) return;
  }

  // Beta selection
  if (fDst.fBeta.beta < 0.3 || fDst.fBeta.pattern > 4) return;

  // Hits in L2
  if (!(fDst.fStatus.ustat&HasL2)) return;

  Double_t rcut = (ismc) ? fDst.fMCinfo.rgt : fDst.fRTI.cfi*1.2;

  Int_t phpat = fDst.fHeader.phpat;
  Int_t psel  = (phpat&0x3e) ? 1 : 0;

  for (Int_t i = 0; i < 4; i++) {
    if (qsel == 1 && i >= 2 && !(fDst.fStatus.ustat & IsInEcal)) continue;

    if (i == 1 && !(has19&1)) continue;
    if (i == 2 && !(has19&2)) continue;
    if (i == 3 && has19 != 3) continue;

    if (fDst.fTrack.csqy[i] < 0 || 10 < fDst.fTrack.csqy[i]) continue;

    TString  shn = Form("hist%d", i+1);
    if (qsel == 1) {
      if (psel) Fill(shn+"1", fDst.fTrack.rgt[i], rcut);
      else      Fill(shn+"3", fDst.fTrack.rgt[i], rcut);
    }
    if (qsel == 2) {
      if (psel) Fill(shn+"2", fDst.fTrack.rgt[i], rcut);
      else      Fill(shn+"4", fDst.fTrack.rgt[i], rcut);
    }
    if (ismc || fDst.fTrack.rgt[1] > fRcut) {
      if (psel) Fill(shn+"5", fDst.fTrack.rgt[i], qin);
      else      Fill(shn+"6", fDst.fTrack.rgt[i], qin);
    }
  }

  // e+ fraction

  Int_t ntrk = (fDst.fStatus.status>>13)&3;
  if (ntrk != 1) return;
  
  if (fDst.fBetaH.pattern != 4444 || fDst.fBetaH.beta < 0.8) return;

  Double_t rgt  = fDst.fTrack.rgt [fDst.fTrack.itp];
  Double_t csqx = fDst.fTrack.csqx[fDst.fTrack.itp];
  Double_t csqy = fDst.fTrack.csqy[fDst.fTrack.itp];

  if (rgt == 0 || csqx < 0 || 10 < csqx
               || csqy < 0 || 10 < csqy) return;

  if (TMath::Abs(fDst.fEcal.dtrk[0]) > 2 || 
      TMath::Abs(fDst.fEcal.dtrk[1]) > 5) return;

  Double_t engy = fDst.fEcal.enew[1];

  if (fDst.fEcal.catl)    return;  // ecal->Status&CATLEAK
  if (fDst.fEcal.bdt < 0) return;  // Ecal BDT
  if (engy < 0.5)         return;  // Ecal energy
  
  if (fDst.fTrdK.nhits < 10) return;
  Double_t trdk = fDst.fTrdK.llre[0];

  Double_t ev = (rgt > 0) ? 1/engy : -1/engy;
  Fill("hist51", ev, trdk);
  if (trdk > 0.7) return;

  Fill("hist52", ev, fDst.fRTI.cf*1.2); 
  if (engy < fDst.fRTI.cf*1.2) return;

  Fill("hist53", engy, engy/rgt);
}

void MdSel::RunSel(void)
{
  if (pMCEventg(0) || fBrun == Run()) return;

  if (Run() < 1305800000) {
    TkDBc       ::ForceFromTDV = 0;
    TrExtAlignDB::ForceFromTDV = 0;
    TrExtAlignDB::version      = 2;
    return;
  }


  if (Run() >= RUN_STD) {   // std data
    cout << "std data" << endl;
    fBrun = Run();
    return;
  }


  GetRTIRunTime(Run(), fTRTI);
  fBrun = Run();
  cout << "Run " << Run() << " Time: " << fTRTI[0] << " " << fTRTI[1] << endl;

  fBadRun = 0;

  string stmp;
  if (!pMCEventg(0) && getsetup()->IsBadRun(stmp, UTime(), Run())) {
    TString srb = stmp.c_str();
    if (srb.Contains("TRD_UnusableForAnalysis")) fBadRun |= 1;
    if (srb.Contains("TRD_BadForEPseparation" ) ||
	srb.Contains("TRD_WeakForEPseparation")) fBadRun |= 2;
    if ( Run() == 1306219312 || Run() == 1306219522  || Run() == 1306233745 ||
	(Run() >= 1307125541 && Run() <= 1307218054) || Run() == 1321198167)
      fBadRun |= 4;
  }
}

void MdSel::RTISel(void)
{
  if (pMCEventg(0)) return;
  if (Run() < 1305800000) return;

  Int_t utime = UTime();
  if (utime == fBtim) return;

  fDst.fStatus.ustat |= FirstRTI;

  fBtim = utime;
  fRTI.nev = 0;
  if (getsetup()->getRTI(fRTI, utime) != 0) return;

  fBadRTI  = 0;
  fBadRTI |= (fRTI.nev < 1800)              ? 0 : 0x01;
  fBadRTI |= (fRTI.ntrig/fRTI.nev > 0.98)   ? 0 : 0x02;
  fBadRTI |= (fRTI.npart > 0)               ? 0 : 0x04;
  fBadRTI |= (fRTI.npart/fRTI.ntrig >
	       0.07/1600*fRTI.ntrig &&
	      fRTI.npart/fRTI.ntrig < 0.25) ? 0 : 0x08;
  fBadRTI |= (fRTI.lf > 0.5 && fRTI.lf < 1) ? 0 : 0x10;
  fBadRTI |= (fRTI.nerr >= 0 && 
	      fRTI.nerr/fRTI.nev < 0.1)     ? 0 : 0x20;
  fBadRTI |= (fRTI.zenith < 40)             ? 0 : 0x40;

  fInSAA = IsInSAA();
  fLtf   = fRTI.lf*fRTI.nev/(fRTI.nev+fRTI.nerr);
  fRcut  = TMath::Abs(fRTI.cfi[0][1]); 

  if (fLtf > 1) fLtf = 1;
  if (fLtf < 0) fLtf = 0;

  AMSPoint pn1, pn9, pd1, pd9;
  GetRTIdL1L9(0, pn1, pd1, UTime(), 60);
  GetRTIdL1L9(1, pn9, pd9, UTime(), 60);

  fPGMD[0] = pd1.y();
  fPGMD[1] = pd9.y();

  if ((fBadRun&4) || fBadRTI || fInSAA || fPGMD[0] > 35
                                       || fPGMD[1] > 45) return;

  TH1F *hist = (TH1F *)fFile->Get("hist00");
  if (hist) {
    TAxis *ax = hist->GetXaxis();
    Int_t  ib = ax->FindBin(fRcut*1.2);
    fRcut = ax->GetBinLowEdge(ib+1);
    for (Int_t i = ib; i <= hist->GetNbinsX(); i++)
      hist->Fill(ax->GetBinCenter(i), fLtf);
  }
}

void MdSel::MCinit()
{
  if (!pMCEventg(0)) return;

  if (fLevt > 0) MCtrig("hmct00", 1, fLevt);

  if (!fChain) return;
  fChain->LoadTree(fNevt);

  TFile *file = fChain->GetCurrentFile();
  if (!file) return;

  TObjString *sbj = (TObjString *)file->Get("DataCards");
  if (!sbj) return;

  Int_t  trig = 0;
  TString str = sbj->GetString();

  Double_t pmin = 0;
  Double_t pmax = 0;
  Int_t    flux = 0;
  Int_t    chg  = 1;

  enum { NP = 37 };
  Int_t part[NP] = { 
    1,  2,  3,  5,  6,  7,  8,  9, 11, 12, 13, 14,
    15, 17, 18, 21, 25, 45, 46, 47, 49, 61, 62, 63,
    64, 65, 66, 67, 68, 69, 71, 72, 74, 84, 87, 114,
    145};
  Int_t chrg[NP] = { 
    0,  1, -1,  1, -1,  0,  1, -1,  1, -1,  0,  1,
   -1,  0,  0, -1,  0,  1,  1,  2,  2,  3,  3,  4,
    4,  5,  5,  6,  7,  8, 10, 11, 13, 23, 26,  4,
   -1};


  TObjArray *ar = str.Tokenize("\n");
  for (Int_t i = 0; i < ar->GetEntries(); i++) {
    TString ss = ((TObjString *)ar->At(i))->GetString();
    if (ss.BeginsWith("TRIG=")) {
      ss.ReplaceAll("TRIG=", "");
      trig = ss.Atoi();
    }
    if (ss.Contains("PMIN="))
    {
      TString sss = ss;
      Int_t i = sss.First("=");
      sss.Remove(0,i+1);
      pmin = sss.Atof();
    }
    if (ss.Contains("PMAX="))
    {
      TString sss = ss;
      Int_t i = sss.First("=");
      sss.Remove(0,i+1);
      pmax = sss.Atof();
    }
    if (ss.Contains("AMSFSCRIPT=")) {
      if (ss.Contains("pl1.flux."  ))   flux =   1;
    }
    if(ss.Contains("PART="))
    {
      TString sss = ss;
      Int_t i = sss.First("=");
      sss.Remove(0,i+1);
      Int_t par = sss.Atoi();
      for(int inp=0;inp<NP;inp++)
      {
        if(par==part[inp])
          chg = chrg[inp];
      }
    }
  }
  fPGMD[0] = pmin;
  fPGMD[1] = pmax;
  if (flux) {
    fPGMD[0] = -pmin;
    fPGMD[1] = -pmax;
  }

  cout << Form("TRIG: %d %8d ", Run(), trig)
    << "chg,pmin,pmax= " << chg << " " << pmin << " " << pmax << endl;

  fBrun = Run();
  fChrg = TMath::Abs(chg);
  fLevt = 0;
  MCtrig("hist00");
}

void MdSel::MCtrig(const char *hname, Int_t emin, Int_t emax)
{
  Double_t pmin = fPGMD[0];
  Double_t pmax = fPGMD[1];

  Int_t flux = 0;
  if (pmin < 0 && pmax < 0) {
    flux = 1;
    pmin = -pmin;
    pmax = -pmax;
  }

  Int_t chg  = fChrg;
  Int_t ntrg = emax-emin+1;

  if (emin < 0 && emax < 0) {
    AMSEventR *evt = new AMSEventR;
    TTree *tree = fChain->GetTree();
    Int_t  nent = tree->GetEntries();
    tree->SetBranchAddress("ev.", &evt);
    tree->GetEvent(0);      emin = evt->Event();
    tree->GetEvent(nent-1); emax = evt->Event();

    ntrg = (Int_t)((emax-emin+1.)/nent*(nent+1));
  }

  cout << Form("TRIG: %d %8d ", fBrun, ntrg)  
    << "chg,pmin,pmax= " << chg << " " << pmin << " " << pmax;
  if (flux) cout << " flux";
  cout << endl;

  if (chg != 0) { pmin /= chg; pmax /= chg; }

  TH1F *hist = (TH1F *)fFile->Get(hname);
  if (hist) {
    TString   sct = Form("(%f<x && x<%f)", pmin, pmax);
    TString   sfn = sct+"? 1/x       : 0";
    if (flux) sfn = sct+"? x^(-2.78) : 0";
    TF1 func("func", sfn);

    Double_t hmin = hist->GetXaxis()->GetXmin(); if (hmin < pmin) hmin = pmin;
    Double_t hmax = hist->GetXaxis()->GetXmax(); if (hmax > pmax) hmax = pmax;
    UInt_t   nfil = ntrg*func.Integral(hmin, hmax)
      /func.Integral(pmin, pmax);
    hist->FillRandom("func", nfil);
  }
}

void MdSel::UTerminate()
{
  if (fLevt > 0) MCtrig("hmct00", 1, fLevt);

  if (fFile && LxMCcutoff::GetHead()) {
    fFile->cd();
    //  LxMCcutoff::GetHead()->GetExp()->Write("hexp");
  }
}

void MdSel::Stat()
{
  Double_t rtm = fTimer.RealTime();
  Double_t crm = fTimer.CpuTime();
  Double_t mem = memchk();
  fTimer.Continue();
  if (fNtot <= 0) fNtot = 1;

  cout << Form("%7d %8d (%5.1f%%) %4.0f sec ",
      fNfil, fNevt, 100.*fNevt/fNtot, rtm);

  if (fNevt/rtm < 500) cout << Form("(%5.1f  Hz) ", fNevt/rtm);
  else                 cout << Form("(%5.1f kHz) ", fNevt/rtm*1e-3);

  cout << Form("%4.1f%% %3.0fMB", crm/rtm*100, mem/1024) << endl;
}

void MdSel::Fill(const char *hname, Double_t x, Double_t y, Double_t w)
{
  if (!fFile) return;
  TH2F *hist = (TH2F *)fFile->Get(hname);
  if (hist) hist->Fill(x, y, w);
  else {
    static Int_t nerr = 0;
    if (nerr++ < 20)
      cout << "Hist not exists " << hname << endl;
  }
}

Int_t MdSel::SpecialCut(AMSPoint _coo, AMSDir _dir)
{
  // AMSmceventg::Specualcut

  static AMSPoint cross1;
  static AMSPoint cross9;
  static double par1[4] = { 63.14, 48.4,  158.920, 67.14 };
  static double par9[4] = { 46.62, 34.5, -135.882, 67.14 };
  static bool  initdone = false;
  if (!initdone){
    cross1[0] = par1[0]; cross1[1] = par1[1]; cross1[2] = par1[2];
    cross9[0] = par9[0]; cross9[1] = par9[1]; cross9[2] = par9[2];
    cout << "Specialcut-I-CrossingParametersLayer1: " << cross1 << endl;
    cout << "Specialcut-I-CrossingParametersLayer9: " << cross9 << endl;
    initdone = true;
  }

  bool layer1 = false;
  bool layer9 = false;

  if (_dir[2]) {
    AMSPoint extrap = _coo+_dir*((cross1[2]-_coo[2])/_dir[2]);
    if (fabs(extrap[0]) < fabs(cross1[0]) &&
	fabs(extrap[1]) < fabs(cross1[1]) &&
	sqrt(extrap[0]*extrap[0]+extrap[1]*extrap[1])<par1[3]) layer1 = true;
  }
  if (_dir[2]) {
    AMSPoint extrap = _coo+_dir*((cross9[2]-_coo[2])/_dir[2]);
    if (fabs(extrap[0]) < fabs(cross9[0]) &&
	fabs(extrap[1]) < fabs(cross9[1]) &&
	sqrt(extrap[0]*extrap[0]+extrap[1]*extrap[1])<par9[3]) layer9 = true;
  }

  Int_t isin = 1;
  if (layer1) isin += 2;
  if (layer9) isin += 4;

  return isin;
}

TrTrackR *MdSel::GetTrack(void)
{
  ParticleR *p = pParticle(0); if (!p) return 0;
  if (p->pTrTrack()) return p->pTrTrack();

  BetaR *b = p->pBeta();
  if (b && b->Beta > 0.3 && b->pTrTrack()) return b->pTrTrack();

  BetaHR *bh = p->pBetaH();
  if (bh && bh->GetBeta() > 0.3 && bh->pTrTrack()) return bh->pTrTrack();

  VertexR *v = p->pVertex(); if (!v) return 0;
  TrTrackR *tr = 0;
  for (Int_t i = 0; i < v->NTrTrack(); i++) {
    TrTrackR *t = v->pTrTrack(i);
    if (!tr || TMath::Abs(tr->GetRigidity()) < TMath::Abs(t->GetRigidity()))
      tr = t;
  }
  return tr;
}


Int_t mdsel(AMSChain &ch, const char *oname, Int_t nevt = -1);

Int_t mdsel(const char *flist, Int_t ibeg, Int_t iend, 
	    const char *oname, Int_t nevt = -1)
{
  TFile::SetOpenTimeout(10);

  AMSChain ch;
  if (ch.AddFromFile(flist, ibeg-1, iend) < 0) return -1;

  Int_t ntr = ch.GetNtrees();
  cout << "Nadd= " << iend-ibeg+1 << " " << ntr << endl;

  return mdsel(ch, oname, nevt);
}

Int_t mdsel(const char *fname, const char *oname, Int_t nevt = -1)
{
  AMSChain ch;
  TString sfn = fname;
  if (!sfn.Contains(".root")) sfn += "*.root";
  if (ch.Add(sfn) <= 0) return -1;

  return mdsel(ch, oname, nevt);
}

Int_t mdsel(AMSChain &ch, const char *oname, Int_t nevt)
{
  Int_t ntr  = ch.GetNtrees();
  Int_t nent = ch.GetEntries();
  if (ntr <= 0 || nent <= 0) return -2;

  cout << "Ntr,Nent= " << ntr << " " << nent << endl;

  MdSel::fNtot  = nent;
  MdSel::fFname = oname;
  MdSel::fChain = &ch;

  cout << "New vars" << endl;

  MdSel mdsel;
  if (nevt > 0) ch.Process(&mdsel, "", nevt);
  else          ch.Process(&mdsel);

  return 0;
}

int main(int argc, char *argv[])
{
  if (argc < 3) {
    cout << "mdsel [fname] [oname] (nevt= -1)" << endl;
    cout << "mdsel [fname] [ibeg] [iend] [oname] (nevt= -1)" << endl;
    return 1;
  }

  Int_t ret = -1;
  if (argc == 3) ret = mdsel(argv[1],      argv[2]);
  if (argc == 4) ret = mdsel(argv[1],      argv[2],
			              atoi(argv[3]));
  if (argc == 5) ret = mdsel(argv[1], atoi(argv[2]),
			              atoi(argv[3]), argv[4]);
  if (argc >= 6) ret = mdsel(argv[1], atoi(argv[2]),
			              atoi(argv[3]), argv[4], atoi(argv[5]));

  if (ret  == 0) cout << "Done" << endl;
  return ret;
}
