#include "root.h"
#include "amschain.h"
#include "SpFold.h"

#include "TFile.h"
#include "TKey.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TH3F.h"

#include <iostream>
#include <fstream>
#include <string>
#include <map>

#include "mdst.h"
using namespace std;
Int_t fError = 0;
Int_t LastEventTime = 0;
Int_t iCutoff = 0, TBin = 0;

class DSTFill:public DST{
	public:
		TChain *fCh;
		TFile *fFile;
		TTree *fTree;
		Int_t fEntry;
		Int_t fMode;
		Int_t IsMC;
		Int_t bartel_count;

		Int_t     fBrun;
		Int_t     fEofs;
		Double_t fXtime;

		Int_t    fHasTrd;
		Int_t    fHasEcl;
		Int_t    fHas19;
		Int_t    fMCpl;
		Double_t fMCW;
		Double_t fDMW;
		Double_t fRcut;
		Double_t fRv3;

		Int_t    fBadRTI;
		Int_t    fREFIT;

		TH1F* Counts_TOFEff_sample;
		TH1F* Counts_TOFEff;
		TH1F* Counts_FullSpan;
		TH1F* Counts_FullSpan_unbias;
		TH2F* ExposureIGRF;

		DSTFill(TString fn, TString ofn, Int_t mode,Int_t run_count):fBadRTI(0),fDMW(1),fMCW(1),fMCpl(1){
			fMode = mode%10;
			IsMC  = mode/10;
			bartel_count = run_count;
		
			fCh = new TChain("tree");
			fCh->Add(fn);
			if(fCh->GetEntries()<=0){
				cout<<"Input Error"<<endl;
				exit(0);
			}
			fFile = new TFile(ofn, "recreate");
			if(!fFile || !fFile->IsOpen()){
				cout<<"Output Error"<<endl;
				fFile = new TFile(ofn, "recreate");
				if(!fFile || !fFile->IsOpen()) exit(0);
			}
			SetAddress(fCh);

			TH1F *he = 0;
			TObjArray *ar = fCh->GetListOfFiles();
			for(int i=0;i<ar->GetEntries();i++){
				TString hfn = ar->At(i)->GetTitle();
				TFile *hf = TFile::Open(hfn);
				if(!hf)continue;
				TList *keys = hf->GetListOfKeys();
				for(int j=0;j<keys->GetSize();j++){
					TKey *key = (TKey*)keys->At(j);
					if(!key)continue;
					TString shn = key->GetName();
					if(shn=="hist00"){
						fFile->cd();
						TH1F *htmp = (TH1F*)hf->Get("hist00");
						if(!htmp)continue;
						cout<<"Add hist00 from "<<hfn<<endl;
						if(!he)he =  (TH1F*)htmp->Clone("hexp");
						else he->Add(htmp);
					}
				}
				hf->Close();
				delete hf;
			}
			Init();
			EventLoop();
		}
		void LoadTree(Int_t entry){fEntry = fCh->LoadTree(entry); Clear();}
		void GetEntry(const TString bname){
			TBranch *br = fCh->FindBranch(bname);
			if(br)br->GetEntry(fEntry);
			else{ if(fError++<20)cout<<"Error Branch: "<<bname<<endl;}
		}
		void Init();
		void FillExpTime();
		bool FillTofEff(); 
		bool FillInnTrEff();
		bool FillL1TrEff();
		bool FillL9TrEff();
		void FillFsEff();
		bool FillFlux();
		void EventLoop();
		void Failed(void);
		double Ecal2Rig(Double_t);
		void Fill(TString,Float_t,Float_t,Float_t);
		Double_t GetMCW(Double_t,Double_t[10]);
		Double_t GetMCW_old(Double_t);
};

void DSTFill::Failed(void){
	cout << "Failed" << endl;
	exit(-1);
}
Double_t DSTFill::GetMCW(Double_t rmc,Double_t par[6]){
	if (!IsMC) return 1;

	enum { N = 4 };
	Double_t xn[N]    = {2,5,10,60};
	Double_t bb[2] 	  = {0, 0};

	TF1 *fun = new TF1("spfun", SplFit::SpFunc, 1.92, 60.3, N*2+2);
	for (Int_t i = 0; i < N;   i++) fun->SetParameter(i,   xn[i]);
	for (Int_t i = 0; i < N+2; i++) fun->SetParameter(i+N, par[i]);

	SplFit::fLogX = 1; SplFit::fLogY = 1;
	SplFit::fBlxL = 0; SplFit::fBlxU = 1;
	SplFit::fN    = (fun->GetNpar()-2)/2;

	//cout << "Rigidity: " << rmc << " MCW: " << fun->Eval(rmc) << " Func: " << frt->Eval(rmc) << endl;
	if(rmc<=1.92) return fun->Eval(1.92);
	else if(rmc>=60.3) return fun->Eval(60.3);
	else return fun->Eval(rmc);
}void DSTFill::Fill(TString hname, Float_t x, Float_t y=1, Float_t w = 1){
	TObject *obj = fFile->Get(hname);
	if(!obj){
		cout<<Form("Hist %s not found",hname.Data())<<endl;
		return ;
	}
	TString htype = obj->ClassName();
	if(htype.Contains("TH1")){
		TH1 *h = (TH1 *)obj;
		h->Fill(x,y);
	}
	else{
		TH2 *h = (TH2 *)obj;
		h->Fill(x,y,w);
	}
}

void DSTFill::Init(){
	//Initialization of Histograms
	cout<<"Start of Init"<<endl;
	Int_t nbin = 40;
	Double_t bin[41]=
	{
		1.92,2.15,2.4,
		2.67,2.97,3.29,3.64,4.02,4.43,4.88,5.37,
		5.9,6.47,7.09,7.76,8.48,9.26,10.1,11,12,
		13,14.1,15.3,16.6,18,19.5,21.1,22.8,24.7,
		26.7,28.8,31.1,33.5,36.1,38.9,41.9,45.1,
		48.5,52.2,56.1,60.3
	};
	Double_t time_independent_bin[76]={
		0.8,1,1.16,1.33,1.51,1.71,1.92,2.15,2.4,
		2.67,2.97,3.29,3.64,4.02,4.43,4.88,5.37,
		5.9,6.47,7.09,7.76,8.48,9.26,10.1,11,12,
		13,14.1,15.3,16.6,18,19.5,21.1,22.8,24.7,
		26.7,28.8,31.1,33.5,36.1,38.9,41.9,45.1,
		48.5,52.2,56.1,60.3,64.8,69.7,74.9,80.5,
		86.5,93,100,108,116,125,135,147,160,175,
		192,211,233,259,291,330,379,441,525,643,
		822,1130,1800,3000,8000};

	Int_t necalbin = 26;
	Double_t ecalbin[27]={
		0.8,1,1.16,1.33,1.51,1.71,1.92,2.15,2.4,
		2.67,2.97,3.29,3.64,4.02,4.43,4.88,5.37,
		5.9,6.47,7.09,7.76,8.48,9.26,10.1,147,1800,8000};

	cout << "Use new binning" << endl;
	fFile->cd();

	Int_t     nbrv = nbin*2+2;
	Double_t *vbin = new Double_t[nbrv+1];
	for (Int_t i = 0; i <= nbin; i++) {
		vbin[i]        = -1/bin[i];
		vbin[i+nbin+2] =  1/bin[nbin-i];
	}
	vbin[nbin+1] = 0;

	Double_t lbin[121];
	Double_t rbin[87];
	for (Int_t i = 0; i <= 120; i++) lbin[i] = TMath::Power(10, i*0.05-3);
	for (Int_t i = 0; i <=  86; i++) rbin[i] = TMath::Power(10, i*0.05-0.3);

	//........Creating time binning for days
	float TimeIntegration       = 27.00;      //day time integration
	float TimeBin_width 	    = TimeIntegration*24*3600;
	long unsigned int TimeBegin = 1305417600; //Date2011.05.15;
	long unsigned int TimeEnd   = 1494374400; //Date2017.05.10;
	float NbinTime_f            = (float) (TimeEnd-TimeBegin)/TimeBin_width;
	long int NbinTime           = (long int)NbinTime_f;
	if (NbinTime_f-(float)NbinTime > 0.) NbinTime += 1;
	double *xtimeHist = new double [NbinTime+1];
	for (ULong64_t it=0; it<=NbinTime; it++) xtimeHist[it] = (double)TimeBegin+(it)*(double)TimeBin_width;


	// Histogram Definitions
	new TH2F("ExposureIGRF","Rigidity vs Exposure Time;Time[s];Rigidity [GV]",NbinTime,xtimeHist,nbin,bin); 
	new TH2F("hist11","Selected Events;Rigidity [GV]; CutOff Rigidity [GV]",nbin,bin,nbin,bin); 
	// Flux
	new TH1D("He_generated", "Events; Rigidity [GV];counts",nbin,bin);	
	new TH1D("He_generated_weighted", "Events; Rigidity [GV];counts",nbin,bin);
	new TH2F("He_selected", "Events; Time[s];Rigidity [GV];counts",NbinTime,xtimeHist,nbin,bin);	
	new TH2F("He_selected_weighted", "Events; Time[s];Rigidity [GV];counts",NbinTime,xtimeHist,nbin,bin);
	new TH2F("HeFlux", "Events; Time[s];Rigidity [GV];counts",NbinTime,xtimeHist,nbin,bin);	
	new TH2F("HeFlux_weighted", "Events; Time[s];Rigidity [GV];counts",NbinTime,xtimeHist,nbin,bin);	
	// Trigger Eff.
	new TH2F("Selected_Events_generated_rigidity", "Selected Events; Time[s];Rigidity [GV];counts",NbinTime,xtimeHist,nbin,bin);
	new TH2F("Selected_Events_reweighted_generated_rigidity", "Selected Events; Time[s];Rigidity [GV];counts",NbinTime,xtimeHist,nbin,bin);
	new TH2F("Selected_Events_reconstructed_rigidity", "Selected Events; Time[s];Rigidity [GV];counts",NbinTime,xtimeHist,nbin,bin);
	new TH2F("Selected_Events_reweighted_reconstructed_rigidity", "Selected Events; Time[s];Rigidity [GV];counts",NbinTime,xtimeHist,nbin,bin);
	new TH2F("Counts_Trigger", "Triggered Events; Time[s];Rigidity [GV];counts",NbinTime,xtimeHist,nbin,bin);
	new TH2F("Counts_Trigger_Unbias", "Unbiased Events; Time[s];Rigidity [GV];counts",NbinTime,xtimeHist,nbin,bin);
	new TH2F("Counts_Trigger_weighted", "Triggered Events; Time[s];Rigidity [GV];counts",NbinTime,xtimeHist,nbin,bin);
	new TH2F("Counts_Trigger_Unbias_weighted", "Unbiased Events; Time[s];Rigidity [GV];counts",NbinTime,xtimeHist,nbin,bin);
	new TH1D("Counts_Trigger_1", "Triggered Events; Rigidity [GV];counts",nbin,bin);
	new TH1D("Counts_Trigger_Unbias_1", "Unbiased Events; Rigidity [GV];counts",nbin,bin);
	new TH1D("Counts_Trigger_weighted_1", "Triggered Events; Rigidity [GV];counts",nbin,bin);
	new TH1D("Counts_Trigger_Unbias_weighted_1","Unbiased Events; Rigidity [GV];counts",nbin,bin);
	// TOF Eff.
	new TH2F("Counts_TOFEff_sample" ,"Sample for ToF;Time[s];Rigidity [GV];counts",NbinTime,xtimeHist,nbin,bin);
	new TH2F("Counts_TOFEff" ,"ToF Selection;Time[s];Rigidity [GV];counts",NbinTime,xtimeHist,nbin,bin);
	new TH2F("Counts_TOFEff_sample_weighted" ,"Sample for ToF (weighted);Time[s];Rigidity [GV];counts",NbinTime,xtimeHist,nbin,bin);
	new TH2F("Counts_TOFEff_weighted" ,"ToF Selection (weighted);Time[s];Rigidity [GV];counts",NbinTime,xtimeHist,nbin,bin);
	//Inner Tracker Eff.
	new TH2F("Counts_Tracker_sample","Sample for InnTrack;Time [s];Rigidity [GV];counts",NbinTime,xtimeHist,nbin,bin);
	new TH2F("Counts_Tracker_sample_weighted","Sample for InnTrack;Time [s];Rigidity [GV];counts",NbinTime,xtimeHist,nbin,bin);
	new TH2F("Counts_Tracker_sample_1","Sample for InnTrack;Time [s];Rigidity [GV];counts",NbinTime,xtimeHist,nbin,bin);
	new TH2F("Counts_Tracker_sample_2","Sample for InnTrack;Time [s];Rigidity [GV];counts",NbinTime,xtimeHist,nbin,bin);
	new TH2F("Counts_Tracker_sample_3","Sample for InnTrack;Time [s];Rigidity [GV];counts",NbinTime,xtimeHist,nbin,bin);
	new TH2F("Counts_InnTrack","InnTrack Selection;Time[s];Rigidity [GV];counts",NbinTime,xtimeHist,nbin,bin);
	new TH2F("Counts_InnTrack_weighted","InnTrack Selection;Time[s];Rigidity [GV];counts",NbinTime,xtimeHist,nbin,bin);
	new TH2F("Counts_InnTrack_1","InnTrack Selection;Time[s];Rigidity [GV];counts",NbinTime,xtimeHist,nbin,bin);
	new TH2F("Counts_InnTrack_2","InnTrack Selection;Time[s];Rigidity [GV];counts",NbinTime,xtimeHist,nbin,bin);
	new TH2F("Counts_InnTrack_3","InnTrack Selection;Time[s];Rigidity [GV];counts",NbinTime,xtimeHist,nbin,bin);	
	new TH2F("Cutoff_Tracker","Cutoff Rigidity vs Tracker Rigidity;Rigidity [GV];Rigidity [GV];counts",nbin,bin,nbin,bin);	
	new TH2F("Ecal_Tracker","Ecal Rigidity vs Tracker Rigidity;Rigidity [GV];Rigidity [GV];counts",nbin,bin,nbin,bin);
	new TH2F("Ecal_Tracker_weighted","Ecal Rigidity vs Tracker Rigidity;Rigidity [GV];Rigidity [GV];counts",nbin,bin,nbin,bin);
	//L1 Tracker Eff.
	new TH2F("Counts_L1","L1 Selection;Time[s];Rigidity [GV];counts",NbinTime,xtimeHist,nbin,bin);	
	new TH2F("Counts_L1_weighted","L1 Selection;Time[s];Rigidity [GV];counts",NbinTime,xtimeHist,nbin,bin);
	new TH2F("Counts_L1_sample","Sample for L1 Tracker;Time[s];Rigidity [GV];counts",NbinTime,xtimeHist,nbin,bin);
	new TH2F("Counts_L1_sample_weighted","Sample for L1 Tracker;Time[s];Rigidity [GV];counts",NbinTime,xtimeHist,nbin,bin);	
	//L9 Tracker Eff.
	new TH2F("Counts_L9","L9 Selection;Time[s];Rigidity [GV];counts",NbinTime,xtimeHist,nbin,bin);
	new TH2F("Counts_L9_1","L9 Selection;Time[s];Rigidity [GV];counts",NbinTime,xtimeHist,nbin,bin);
	new TH2F("Counts_L9_2","L9 Selection;Time[s];Rigidity [GV];counts",NbinTime,xtimeHist,nbin,bin);
	new TH2F("Counts_L9_3","L9 Selection;Time[s];Rigidity [GV];counts",NbinTime,xtimeHist,nbin,bin);
	new TH2F("Counts_L9_weighted","L9 Selection;Time[s];Rigidity [GV];counts",NbinTime,xtimeHist,nbin,bin);
	new TH2F("Counts_L9_sample","Sample for L9 Tracker;Time[s];Rigidity [GV];counts",NbinTime,xtimeHist,nbin,bin);
	new TH2F("Counts_L9_sample_weighted","Sample for L9 Tracker;Time[s];Rigidity [GV];counts",NbinTime,xtimeHist,nbin,bin);
	// Migration Matrix
	new TH2F("MigrationMatrix_1","Reconstructed Rigidity vs Generated Rigidity;Generated Rigidity [GV];Reconstructed Rigidity [GV];counts",nbin,bin,nbin,bin);
	new TH2F("MigrationMatrix","Reconstructed Rigidity vs Generated Rigidity;Generated Rigidity [GV];Reconstructed Rigidity [GV];counts",nbin,bin,nbin,bin);
}


void DSTFill::FillExpTime(){
	GetEntry("header");
	Int_t utctime = fHeader.utime;
	static Int_t lastevent = 0;
	GetEntry("rti");
	Double_t rcut = (!IsMC) ? fRTI.cfi*1.2 : fMCinfo.rgt;  
	TH2F *hist = (TH2F *)fFile->Get("ExposureIGRF");
	if (hist && utctime > lastevent)
	{
		TAxis *ax = hist->GetYaxis();
		Int_t  ib = ax->FindBin(rcut);
		for (Int_t i = ib; i <= hist->GetNbinsY(); i++)
			Fill("ExposureIGRF",utctime,ax->GetBinCenter(i),fRTI.lf);
		lastevent = utctime;
	}
}

bool DSTFill::FillL1TrEff(){
	// ToF
	GetEntry("betah");	
	if (fBetaH.beta < 0.4 || fBetaH.pattern != 4444 || fBetaH.q < 1.25) return 0;
	//Physics Trigger
	GetEntry("header");
	Int_t phpat = fHeader.phpat;
	Int_t psel  = (phpat&0x3e) ? 1 : 0;
	if(!psel) return 0;
	Int_t utctime = fHeader.utime;
	// RTI Cutoff
	GetEntry("rti");
	Double_t rcut = (!IsMC) ? fRTI.cfi*1.2 : 0;//fMCinfo.rgt;  
	if(!IsMC)
	{
		if(fRTI.lf < 0.5 || fRTI.zenith > 40) return 0;
		if(fRTI.dl1 > 35 || fRTI.dl9 > 45) return 0;
	}
	// Tracker
	GetEntry("track");
	Double_t rin = fTrack.rgt [0];
	Double_t rfs = fTrack.rgt [3];
	Double_t cinx = fTrack.csqx[0];
	Double_t ciny = fTrack.csqy[0];
	Double_t csq1x = fTrack.csqx[1];
	Double_t csq1y = fTrack.csqy[1];
	Double_t csq9x = fTrack.csqx[2];
	Double_t csq9y = fTrack.csqy[2];
	Double_t cfs = fTrack.csqy[3];
	Double_t qlm = IsMC?2.05:2.0;
	Double_t qin = fTrack.qin;
	Double_t ql1 = fTrack.ql1,
		 ql9 = fTrack.ql9;
	Double_t xl1 = fTrack.coox[0], yl1 = fTrack.cooy[0],
		 xl9 = fTrack.coox[2], yl9 = fTrack.cooy[2];
	Int_t L1XY = fTrack.bitx & 0x001,
	      L9XY = fTrack.bitx & 0x100,
	      InXY = fTrack.bitx & 0x0FE;//254
	Bool_t qsel_l1 = (qlm-0.3)<ql1 && ql1<(qlm+0.5);
	Bool_t qsel_inn = (qlm-0.4)<qin && qin<(qlm+0.9);
	Bool_t qsel_l9= (qlm-0.4)<ql9 && ql9<(qlm+0.9);
	Bool_t csq_inn = ciny < 10 ;
	Bool_t csq_l1 = csq1y < 10;
	Bool_t csq_l9 = csq9y < 10; 

	// L9 & Inner selection
	if (TMath::Abs(yl9) < .5 || TMath::Abs(yl9) > 29.5 || TMath::Abs(xl9) > 33 ) return 0;
	if (1/rin > 1/2.97 || rin < rcut) return 0;
	if(!(qsel_inn && qsel_l9 && (InXY==0x0FE) && (L9XY=0x100) && csq_inn && csq_l9)) return 0;
	// Sample for L1 Tracker
	Fill("Counts_L1_sample",utctime,rin);
	Fill("Counts_L1_sample_weighted",utctime,rin,fMCW);
	// L1 selection
	if(qsel_l1 && L1XY==0x001 && csq_l1)
	{
		Fill("Counts_L1",utctime,rin); 
		Fill("Counts_L1_weighted",utctime,rin,fMCW); 	
		return 1;
	}
	else return 0;
}
bool DSTFill::FillL9TrEff(){
	Int_t stsel = IsInL9;
	if ((fStatus.ustat&stsel) != stsel) return 0;
	//ToF
	GetEntry("betah");	
	if (fBetaH.beta < 0.4 || fBetaH.pattern != 4444 || fBetaH.q < 1.25) return 0;
	//Physics Trigger
	GetEntry("header");
	Int_t phpat = fHeader.phpat;
	Int_t psel  = (phpat&0x3e) ? 1 : 0;
	if(!psel) return 0;
	Int_t utctime = fHeader.utime;
	// RTI Cutoff
	GetEntry("rti");
	Double_t rcut = (!IsMC) ? fRTI.cfi*1.2 : 0;//fMCinfo.rgt;  
	if(!IsMC)
	{	
		if(fRTI.lf < 0.5 || fRTI.zenith > 40) return 0;
		if(fRTI.dl1 > 35 || fRTI.dl9 > 45) return 0;
	}	
	// Tracker
	GetEntry("track");
	Double_t rin = fTrack.rgt [0];
	Double_t rfs = fTrack.rgt [3];
	Double_t cinx = fTrack.csqx[0];
	Double_t ciny = fTrack.csqy[0];
	Double_t csq1x = fTrack.csqx[1];
	Double_t csq1y = fTrack.csqy[1];
	Double_t csq9x = fTrack.csqx[2];
	Double_t csq9y = fTrack.csqy[2];
	Double_t cfs = fTrack.csqy[3];
	Double_t qlm = IsMC?2.05:2.0;
	Double_t qin = fTrack.qin;
	Double_t ql1 = fTrack.ql1,
		 ql9 = fTrack.ql9;
	Double_t xl1 = fTrack.coox[0], yl1 = fTrack.cooy[0],
		 xl9 = fTrack.coox[2], yl9 = fTrack.cooy[2];
	Int_t L1XY = fTrack.bitx & 0x001,
	      L9XY = fTrack.bitx & 0x100,
	      InXY = fTrack.bitx & 0x0FE;//254
	Bool_t qsel_l1 = (qlm-0.4)<ql1 && ql1<(qlm+0.9);
	Bool_t qsel_inn = (qlm-0.4)<qin && qin<(qlm+0.9);
	Bool_t qsel_l9= (qlm-0.3)<ql9 && ql9<(qlm+0.5);
	Bool_t csq_inn = 0 < ciny && ciny < 10 ;
	Bool_t csq_l1 = 0 < csq1y && csq1y < 10;
	Bool_t csq_l9 = 0 < csq9y && csq9y < 10; 
	//Int_t  qstat_l1 = fTrack.qsta[0] & 0x10013D; 
	//Int_t  qstat_l9 = fTrack.qsta[0] & 0x10013D; 

	// L1 & Inner Tracker selection
	if ((xl1*xl1+yl1*yl1) > 62.14*62.14 || TMath::Abs(xl1) > 62.14 || TMath::Abs(yl1) > 47.4 ) return 0;
	if (1/rin > 1/2.97 || rin < rcut)    return 0;
	if(!(qsel_l1 && qsel_inn && (L1XY==0x001) && (InXY==0x0FE) && csq_inn && csq_l1)) return 0;

	// Sample for L9 Tracker
	Fill("Counts_L9_sample",utctime,rin);
	Fill("Counts_L9_sample_weighted",utctime,rin,fMCW);

	// Plot L9 Selection
	if(L9XY) Fill("Counts_L9_1",utctime,rin); 
	if(qsel_l9) Fill("Counts_L9_2",utctime,rin);
	if(csq_l9) Fill("Counts_L9_3",utctime,rin);

	if(qsel_l9 && L9XY==0x100 && csq_l9)
	{
		Fill("Counts_L9",utctime,rin);
		Fill("Counts_L9_weighted",utctime,rin,fMCW);
		return 1;
	}
	else return 0;
}
bool DSTFill::FillInnTrEff(){
	// ToF
	GetEntry("betas");	
	if (fBetaHs.beta < 0.4 || fBetaHs.pattern != 4444 || fBetaHs.q < 1.25) return 0;
	//Physics Trigger
	GetEntry("header");
	Int_t phpat = fHeader.phpat;
	Int_t psel  = (phpat&0x3e) ? 1 : 0;
	if(!psel) return 0;
	Int_t utctime = fHeader.utime;
	// RTI Cutoff
	GetEntry("rti");
	Double_t rcut = (!IsMC) ? fRTI.cfi*1.2 : 0;//fMCinfo.rgt;  
	if(!IsMC){	
		if(fRTI.lf < 0.5 || fRTI.zenith > 40) return 0;
		if(fRTI.dl1 > 35 || fRTI.dl9 > 45) return 0;
	}
	// Tracker
	GetEntry("track");
	Double_t rin = fTrack.rgt [0];
	Double_t rfs = fTrack.rgt [3];
	Double_t cin = fTrack.csqy[0];
	Double_t cfs = fTrack.csqy[3];
	Double_t qlm = IsMC?2.05:2.0;
	Double_t qin = fTrack.qin;
	Double_t ql1 = fTrack.ql1,
		 ql9 = fTrack.ql9;
	Double_t xl1 = fTrack.coox[0], yl1 = fTrack.cooy[0],
		 xl9 = fTrack.coox[2], yl9 = fTrack.cooy[2];
	GetEntry("trhit");
	Int_t L1XY = fTrack.bitx & 0x001,
	      L9XY = fTrack.bitx & 0x100,
	      InXY = fTrack.bitx & 0x0FE;//254
	Int_t inner_hit=0;
	if(InXY & 0x002) inner_hit=inner_hit+1;
	if(InXY & 0x004) inner_hit=inner_hit+1;
	if(InXY & 0x008) inner_hit=inner_hit+1;
	if(InXY & 0x010) inner_hit=inner_hit+1;
	if(InXY & 0x020) inner_hit=inner_hit+1;
	if(InXY & 0x040) inner_hit=inner_hit+1;
	if(InXY & 0x080) inner_hit=inner_hit+1;
	// TRD
	GetEntry("trd");
	Double_t q_trd = fTrd.q;
	GetEntry("trdk");
	Int_t trd_nhit = fTrdK.nhits;
	if(trd_nhit < 15 || q_trd < (qlm-0.4) || q_trd > (qlm+0.9)) return 0;
	// Silicon Tracker
	Bool_t csel = cin < 10;
	Bool_t qsel = (qlm-0.3)<qin && qin<(qlm+0.5);
	Bool_t qsel_l1 = (qlm-0.4)<ql1 && ql1<(qlm+0.9);
	Bool_t qsel_l9 = (qlm-0.4)<ql9 && ql9<(qlm+0.9);
	if(!(qsel_l1 && qsel_l9 && (L1XY==0x001) && (L9XY==0x100))) return 0;

	// New rigidity 
	// We cannot use rigidity from tracker, therefore we obtain from Beta, Cutoff and Ecal respectively. 
	// 0-5.9 GV -> Beta
	// 5.9-19.5 GV -> Geomagnetic Cutoff
	// >19.5 -> Ecal

	Double_t He_mass = 3.727; // He-4 rest mass in GeV/c^2
	Double_t rigidity = 0; //init
	rigidity =  ((He_mass/qlm)*fBetaH.beta)/TMath::Sqrt(1-fBetaH.beta*fBetaH.beta); // Use Beta --> R = (M/Z)(B/sqrt(1-B^2))
	if(rigidity < 19.5 || rigidity > 5.9) rigidity = IsMC ? fMCinfo.rgt : fRTI.cfi;	

	Fill("Counts_Tracker_sample",utctime,rigidity);
	Fill("Counts_Tracker_sample_weighted",utctime,rigidity,fMCW);

	if(qsel && inner_hit > 2 && csel)
	{ 
		Fill("Counts_InnTrack",utctime,rigidity); 	
		Fill("Counts_InnTrack_weighted",utctime,rigidity,fMCW); 	
	}
	return 1;
}
bool DSTFill::FillTofEff()
{
	Int_t stsel = HasL2 | IsInL1 | IsInL9 | HasL1 | HasL9;
	if ((fStatus.ustat&stsel) != stsel) return 0;
	//Physics Trigger
	GetEntry("header");
	Int_t phpat = fHeader.phpat;
	Int_t psel  = (phpat&0x3e) ? 1 : 0;
	Int_t unbias = (phpat&0x1) ? 1 : 0; 
	Int_t utctime = fHeader.utime;
	// RTI Cutoff
	GetEntry("rti");
	Double_t rcut = (!IsMC) ? fRTI.cfi*1.2 : fMCinfo.rgt;  
	if(!IsMC){
		if(fRTI.lf < 0.5 || fRTI.zenith > 40) return 0;
		if(fRTI.dl1 > 35 || fRTI.dl9 > 45) return 0;
	}
	// Tracker
	GetEntry("track");
	Double_t rin = fTrack.rgt [0];
	Double_t rfs = fTrack.rgt [3];
	Double_t cin = fTrack.csqy[0];
	Double_t cfs = fTrack.csqy[3];
	Double_t qlm = IsMC?2.05:2.0;
	Double_t ql1 = fTrack.ql1,
		 ql9 = fTrack.ql9;
	Double_t xl1 = fTrack.coox[0], yl1 = fTrack.cooy[0],
		 xl9 = fTrack.coox[2], yl9 = fTrack.cooy[2];
	// LD's L1 fiducial
	if ((xl1*xl1+yl1*yl1) > 62.14*62.14 || TMath::Abs(xl1) > 62.14 || TMath::Abs(yl1) > 47.4 || TMath::Abs(yl9) < .5 || TMath::Abs(yl9) > 29.5 || TMath::Abs(xl9) > 33 ) return 0;

	if (1/rfs > 1/0.8)    return 0;
	if(rfs < rcut) return 0;

	Bool_t hasx = (fTrack.bitx&0x101) == 0x101;
	Bool_t xfid = fabs(xl9)<=33;
	Bool_t csel = 0 < cfs && cfs < 10  && cin < 10;
	Bool_t qsel = (qlm-0.4)<ql1 && ql1<(qlm+0.9) && (qlm-0.4)<ql9 && ql9<(qlm+0.9);
	if (!(hasx && xfid && csel && psel && qsel)) return 0;	
	// Sample ready for ToF
	Fill("Counts_TOFEff_sample",utctime,rfs);
	if(IsMC) Fill("Counts_TOFEff_sample_weighted",utctime,rfs,fMCW);
	// ToF
	GetEntry("betah");	
	if (fBetaH.beta < 0.4 || fBetaH.pattern != 4444 || fBetaH.q < 1.25) return 0;
	Fill("Counts_TOFEff", utctime,rfs);
	if(IsMC) Fill("Counts_TOFEff_weighted", utctime,rfs,fMCW);
	return 1;
}
bool DSTFill::FillFlux(){
	Int_t stsel = HasL2 | IsInL1 | IsInL9 | HasL1 | HasL9;
	if ((fStatus.ustat&stsel) != stsel) return 0;
	// ToF
	GetEntry("betah");	
	if (fBetaH.beta < 0.4 || fBetaH.pattern != 4444 || fBetaH.q < 1.25) return 0;
	//Physics Trigger
	GetEntry("header");
	Int_t phpat = fHeader.phpat;
	Int_t psel  = (phpat&0x3e) ? 1 : 0;
	Int_t unbias = (phpat&0x1) ? 1 : 0; 
	Int_t utctime = fHeader.utime;
	// RTI Cutoff
	GetEntry("rti");
	GetEntry("mcinfo");
	Double_t rcut = (!IsMC) ? fRTI.cfi*1.2 : 0;//fMCinfo.rgt;  
	if(!IsMC){
		if(fRTI.lf < 0.5 || fRTI.zenith > 40) return 0;
		if(fRTI.dl1 > 35 || fRTI.dl9 > 45) return 0;
	}
	// Tracker
	GetEntry("track");
	Double_t rin = fTrack.rgt [0];
	Double_t rfs = fTrack.rgt [3];
	Double_t cin = fTrack.csqy[0];
	Double_t cfs = fTrack.csqy[3];
	Double_t qlm = IsMC?2.05:2.0;
	Double_t ql1 = fTrack.ql1,
		 ql9 = fTrack.ql9;
	Double_t xl1 = fTrack.coox[0], yl1 = fTrack.cooy[0],
		 xl9 = fTrack.coox[2], yl9 = fTrack.cooy[2];
	GetEntry("trhit");
	Int_t L1XY = fTrack.bitx & 0x001,
	      L9XY = fTrack.bitx & 0x100,
	      InXY = fTrack.bitx & 0x0FE;//254
	Int_t inner_hit=0;
	if(InXY & 0x002) inner_hit=inner_hit+1;
	if(InXY & 0x004) inner_hit=inner_hit+1;
	if(InXY & 0x008) inner_hit=inner_hit+1;
	if(InXY & 0x010) inner_hit=inner_hit+1;
	if(InXY & 0x020) inner_hit=inner_hit+1;
	if(InXY & 0x040) inner_hit=inner_hit+1;
	if(InXY & 0x080) inner_hit=inner_hit+1;


	// LD's L1 fiducial
	if ((xl1*xl1+yl1*yl1) > 62.14*62.14 || TMath::Abs(xl1) > 62.14 || TMath::Abs(yl1) > 47.4 || TMath::Abs(yl9) < .5 || TMath::Abs(yl9) > 29.5 || TMath::Abs(xl9) > 33 ) return 0;

	if (1/rfs > 1/0.8)    return 0;
	if (rfs < rcut) return 0;

	Bool_t hasx = (fTrack.bitx&0x101) == 0x101;
	Bool_t xfid = fabs(xl9)<=33;
	Bool_t csel = 0 < cfs && cfs < 10  && cin < 10;
	Bool_t qsel = (qlm-0.4)<ql1 && ql1<(qlm+0.9) && (qlm-0.4)<ql9 && ql9<(qlm+0.9);
	if (!(hasx && xfid && csel && qsel)) return 0;	
	if(!(inner_hit > 2 && L1XY==0x001 && L9XY==0x100)) return 0;

	// Trigger Efficiency
	if(psel)
	{
		Fill("Counts_Trigger",utctime,rfs);
		Fill("Counts_Trigger_weighted",utctime,rfs,fMCW);
		Fill("hist11",rfs,rcut);
		if(IsMC) Fill("MigrationMatrix",fMCinfo.rgt,rfs);
		return 1;
	}	
	else
	{
		Fill("Counts_Trigger_Unbias",utctime,rfs);
		Fill("Counts_Trigger_Unbias_weighted",utctime,rfs,fMCW);
		return 0;
	}




}
void DSTFill::EventLoop(){
	cout<<"Start of Eventloop"<<endl;
	Int_t Nent = fCh->GetEntries();
	cout<<"Ntr, Nent= "<<fCh->GetNtrees()<<' '<<fCh->GetEntries()<<endl;
	cout<<"IsMC= "<<IsMC<<endl;

	fFile->cd();
	fTree = fCh->CloneTree(0);

	Int_t Load = Nent/10;
	Int_t nfil = 0;
	Int_t fBtime =0;

	TStopwatch timer;
	timer.Start();

	Int_t Nrun = 0;

	Double_t par[6]; // Parameters of the weight fit, will be read from the weight_list.txt		
	if(IsMC)
	{
		// Read MCW from list
		ifstream weight_file ("/afs/cern.ch/work/c/ctuysuz/private/he_flux/root/weight_list.txt");
		string weight;
		int l = 0;
		//cout << "Bartel: " << bartel_count << endl;
		if(bartel_count>48) bartel_count = bartel_count - 3;
		while(l < (bartel_count))
		{
		
			for(int i=0;i<6;i++)	
			{
				getline(weight_file,weight);
				if(l == (bartel_count-1))
				{
					cout << "Parameter[" << i  << "]" << " = " << weight << endl;
					par[i] = stod(weight);
				}
			}	
			l++;
		}

		weight_file.close();
	}

	for(Int_t i=0;i<Nent;i++)
	{
		if( (i+1)%Load==0 || (i+1)==Nent)
		{
			Double_t rtm = timer.RealTime();
			Double_t ctm = timer.CpuTime();
			timer.Continue();
			cout << Form("mdfil_info: %5.1f%% %4.0f sec CPU: %5.1f%%",100.*(i+1)/Nent, rtm, ctm/rtm*100)<<endl;
			/*cout << Form("mdfil_info: %8d %8d (%5.1f%%) %4.0f sec (%5.1f kHz) CPU: %5.1f%%",
			  nfil, i+1, 100.*(i+1)/Nent, rtm, (i+1)/rtm*1e-3, ctm/rtm*100)<<endl;
			  */
		}

		LoadTree(i);


		//Apply Status/RTI Cut here
		GetEntry("status");
		Int_t stcut = BadTrig | InSAA;
		if(fStatus.ustat & stcut )  continue;
		if(fStatus.ustat & BadRTI)  continue;
		// Get MCW
		if(IsMC)
		{
			GetEntry("mcinfo");
			fMCW = GetMCW(fMCinfo.rgt,par);	
		}	
		else fMCW = 1;


		FillExpTime();	// Fill Histograms for Exposure Time
		FillTofEff(); // Fill Histograms for Tof Efficiency
		FillInnTrEff(); //Fill Histograms for Inner Tracker Efficiency
		FillL1TrEff(); //Fill Histograms for L1 Efficiency
		FillL9TrEff(); // Fill Histograms for L9 Efficiency
		FillFlux(); // Fill Trigger and Flux Selections

	}
	fFile->cd();
	fFile->Write();
	fFile->Close();
}
Int_t mdfil(Int_t RID, TString num,Int_t bartel_count)
{
	Int_t mode = 10; // mode=10 for MC, mode = 1 for ISS
	TString sdn,sfn,ofn;

	if(mode == 20) // if MC
	{
		// Read from list
		ifstream myfile ("/afs/cern.ch/work/c/ctuysuz/private/he_flux/run_list/file_list_mc.txt");
		string line;
		int i = 0;
		while(i!=RID+1)
		{
			i++;
			getline(myfile,line);
		}
		myfile.close();

		sdn = "/eos/ams/user/s/selu/mdst/He122/";
		sfn = TString(line);
		ofn = "/afs/cern.ch/work/c/ctuysuz/private/he_flux/LSF/result_mc/he_flux_";	
		ofn += std::to_string(bartel_count);
		ofn += "_";
		ofn +=sfn;	
	}
	else if(mode == 10) // if MC
	{
		// Read from list
		ifstream myfile ("/afs/cern.ch/work/c/ctuysuz/private/he_flux/run_list/file_list_mc_time.txt");
		string line;
		int i = 0;
		while(i!=RID+1)
		{
			i++;
			getline(myfile,line);
		}
		myfile.close();
		sdn = TString(line);
		sfn = "";
		/*
		ofn = "/afs/cern.ch/work/c/ctuysuz/private/he_flux/LSF/result_mc/he_flux_";	
		ofn += std::to_string(bartel_count);
		ofn += "_";
		ofn += std::to_string(RID);	
		ofn += ".root";	
		*/
		ofn = Form("/afs/cern.ch/work/c/ctuysuz/private/he_flux/LSF/result_mc_time/time_%d/he_flux_%d_%d.root",bartel_count,bartel_count,RID);
	}

	else // if ISS
	{
		// Read from list
		ifstream myfile ("/afs/cern.ch/work/c/ctuysuz/private/he_flux/run_list/file_list_iss.txt");
		string line;
		int i = 0;
		while(i!=RID+1)
		{
			i++;
			getline(myfile,line);
		}
		myfile.close();
		sdn = "root://eosuser.cern.ch//eos/user/s/selu/ndst/"; 
		sfn = TString(line);
		ofn = "/afs/cern.ch/work/c/ctuysuz/private/he_flux/LSF/result_iss/he_flux_";	
		ofn +=sfn;

	}

	cout <<"Input Destination: " << sdn+sfn << endl;
	cout <<"Output Destination: " << ofn << endl;
	DSTFill dst(sdn+sfn, ofn, mode, bartel_count);
	return 0;
}
int main(int argc, char **argv)
{
	/*if(argc<2)
	  {
	  cout<<"Not Enough Input"<<endl;
	  return 0;
	  }*/
	Int_t ret = -1;
	//mdfil(atoi(argv[1]),TString(argv[1]),0);
/*	
	for(int i=1;i<=81;i++) 
	{
		if(i==46 || i==47 || i==48) continue; // AMS is not taking data
		mdfil(atoi(argv[1]),TString(argv[1]),i);
	}
*/
	int RID = atoi(argv[1])/81;
	int bartel_count = 1 + atoi(argv[1]) - RID*81;
	cout << "Run ID: " << RID << " Bartel: " << bartel_count << endl; 
	if(bartel_count==46 || bartel_count==47 || bartel_count==48) return 0;		
	mdfil(RID,TString(argv[1]),bartel_count);
	return 0;
}
