
#include "TChain.h"
#include "TROOT.h"
#include "TFile.h"
#include "TKey.h"
#include "TH2.h"
#include "TSystem.h"
#include "TStopwatch.h"
#include "TrMass.h"

#include "amschain.h"

#include "mdst.h"

#include <iostream>
#include <fstream>
using namespace std;

class MdChain : public TChain {
    public:
        Int_t fNent;
        Int_t fBent;

        MdChain(const char* name, const char* title = "")
            : TChain(name, title), fBent(-1), fNent(-1) {}

        Long64_t LoadTree(Long64_t entry) {
           if (fNent > 0 && entry == fNent+1) SetBranchStatus("ecalh", 0);

            Long64_t ent = TChain::LoadTree(entry);
            if (fBent >= 0 && ent < fBent) {
                cout << Form("LoadTree: %9d %4.1f %%",
                        (Int_t)entry, 100.*entry/GetEntriesFast()) << endl;
                TBranch *br = fTree->GetBranch("ecalh");
                if (br) br->SetStatus(0);
                fNent += fTree->GetEntries();
            }
            fBent = ent;
            return ent;
        }
};

Int_t mdsel(MdChain &ch, const char *oname, Int_t mode);

Int_t mdsel(TString fname, const char *oname, Int_t mode=1)
{
    AMSChain ca;

    if(mode/10==0){
        TString sfn = "root://eosams.cern.ch//eos/ams/user/t/tracker/PG/mdst/p6/";
        if(fname.Contains("mdst_")) sfn+=fname;
        if(!fname.Contains(".root")) sfn+="*.root";

        MdChain ch("tree");
        ch.Add(sfn);
        cout<<"Ntr,Nent= "<<ch.GetNtrees()<<' '<<ch.GetEntries()<<endl;
        if(ch.GetNtrees()==0){
            cout<<"No Files Added"<<endl;
            return 0;
        }
        return mdsel(ch, oname, mode);
    }
    else{//MC
        TString sfn = fname;
        if(!fname.Contains(".root")) sfn+="*.root";

        MdChain ch("tree");
        ch.Add(sfn);
        cout<<"Ntr,Nent= "<<ch.GetNtrees()<<' '<<ch.GetEntries()<<endl;
        if(ch.GetNtrees()==0){
            cout<<"No Files Added"<<endl;
            return 0;
        }
        return mdsel(ch, oname, mode);
    }
}

Double_t *getbin(Int_t &nbin)
{
    TString sbn = "$AMSDataDir/v5.00/phe_bin.root";
    gSystem->ExpandPathName(sbn);

    TFile fb(sbn);
    if (!fb.IsOpen()) return 0;

    TH1F *hbin = (TH1F *)fb.Get("hist2");
    nbin = hbin->GetNbinsX();

    Double_t *bin = new Double_t[nbin+1];
    for (Int_t i = 0; i <= nbin; i++)
        bin[i] = hbin->GetXaxis()->GetXbins()->fArray[i];

    return bin;
}

Int_t mdsel(MdChain &ch, const char *oname, Int_t mode)
{
    Int_t ismc = mode/10;
    mode = mode%10;

    cout << "mode,ismc= " << mode << " " << ismc << endl;

    TString ofn = oname;

    Int_t ntr  = ch.GetNtrees();
    Int_t nent = ch.GetEntries();
    if (ntr <= 0 || nent <= 0) return -1;

    cout << "Ntr,Nent= " << ntr << " " << nent << endl;

    TFile of(ofn, "recreate");

    TH1F *hist0 = 0;
    TH1F *hexp0 = 0;
    Int_t nadd  = 0;

    TObjArray *ar = ch.GetListOfFiles();
    for (Int_t i = 0; i < ar->GetEntries(); i++) {
        TFile *f = TFile::Open(ar->At(i)->GetTitle());
        if (!f) continue;

        of.cd();
        TList *kl = f->GetListOfKeys();
        for (Int_t j = 0; j < kl->GetSize(); j++) {
            TKey *key = (TKey *)kl->At(j);
            if (!key) continue;

            TString shn = key->GetName();
            if (shn == "hist00") {
                TH1F *htmp = (TH1F *)f->Get(shn);
                if (!hist0) hist0 = (TH1F *)htmp->Clone();
                else hist0->Add(htmp);
                nadd++;
            }
            if (shn == "hexp") {
                TH1F *htmp = (TH1F *)f->Get(shn);
                if (!hexp0) hexp0 = (TH1F *)htmp->Clone(shn);
                else hexp0->Add(htmp);
            }
        }
        delete f;
    }
    std::cout << "nadd= " << nadd << std::endl;
    if (nadd == 0 || !hist0) return -1;

    ch.fNent = -1;
    ch.fBent = nent;
    TCut cut;
    TCut cutt;
    TCut cutr = "ustat[0]&0x4000";

    if(mode==0){
      TCut cutq = "track.qin>0.8 && track.qin<1.5 && track.csqy[0]>0 && track.csqy[0]<10 &&"
        "ntrh<100";
      TCut cutb = "betah.pattern==4444 && betah.q<1.5 && betah.p[0]>0.99 && betah.z[0]==1 &&"
        "betah.chi2t>0 && betah.chi2t<10 && betah.chi2c>0 && betah.chi2c<10 && "
        "(betah.clsn[0]+betah.clsn[2])<=2";
      TCut cut0 = Form("(ustat[0]&%d)", HasL2);
      TCut cut1 = " (track.rgt[0]>0 && event%1000==0) || (track.rgt[0]*betah.beta<0) ||" 
        " (track.rgt[0]<0 && event%100 ==0 && betah.beta<0)";
      cutt = cutb+cutq+cut0+cut1;
    }
    if(mode==3){ // High Charge
      TCut cutl1 = Form("(qmax[0]>2.5 || qmax[1]>2.5 || ql1>2.5 || (ustat[0]&%d))", IsInL1|HasL1);
      TCut cutl9 = Form("(qmax[2]>2.5 || qmax[3]>2.5 || ql9>2.5 || (ustat[0]&%d))", IsInL9|HasL9);
      TCut cutbh = "(betah.pattern==4444 && betah.chi2t>0 && betah.chi2t<10 && betah.chi2c>0 && betah.chi2c<10 &&"
        "(betah.q>2.5 || betah.z[0]>2.5 || betah.ql[0]>2.5 || betah.ql[1]>2.5))";
      TCut cutbs = "(betas.pattern==4444 && betas.chi2t>0 && betas.chi2t<10 && betas.chi2c>0 && betas.chi2c<10 && "
        "(betas.q>2.5 || betas.ql[0]>2.5 || betas.ql[1]>2.5))";
      TCut cutl2 = Form("(ustat[0]&%d)", HasL2);
      TCut cutin = "track.qin>2.5";
      cutt = (cutbh||cutbs)||(cutl1+cutl2+cutin+cutl9);
      cut = cutt;
    }
    if(mode==2){
      TCut cutq = "1.5< track.qin && track.qin <2.5";
      TCut cutb = "beta.pattern<5";
      TCut cut0 = Form(" (ustat[0]&%d)", HasL2);
      TCut cut1 = Form("(ustat[0]&%d) || ((ustat[0]&%d))", IsInL1, HasL1);
      TCut cut2 = Form("(ustat[0]&%d) || ((ustat[0]&%d))", IsInL9, HasL9);
      TCut cutt = cutb+cutq+cut0+cut1+cut2;

      TCut cuts = "(ustat[0]&0x3000)==0x3000 && (betas.pbit&2) && betas.beta>0.3";
      cuts += "1.5<betas.q && betas.q < 2.5";

      cut = cutt||cuts||cutr;
      ch.SetBranchStatus("vertex", 0);
      ch.SetBranchStatus("trcls",  0);
      ch.SetBranchStatus("trdhit", 0);
      ch.SetBranchStatus("tofhit", 0);
      ch.SetBranchStatus("trdg"  , 0);
      ch.SetBranchStatus("ecalh",  0);
      ch.SetBranchStatus("part", 0);
      ch.SetBranchStatus("mass", 0);
    }
    if(mode==4){ // High p/he
      TCut cut00 = Form("(ustat[0]&%d)", HasL2|TrdInFS|TrdInFS|TrdInTr|IsInEcal);
      TCut cutl1 = Form("ql1>0.7 && ql1<2.5 && (ustat[0]&%d)", HasL1);
      TCut cutl9 = Form("ql1>0.7 && ql9<2.5 && (ustat[0]&%d)", HasL9);
      TCut cutin = Form("track.qin>0.7 && track.qin<2.5 && track.csqy[0]>0 && track.csqy[1]>0 &&track.csqy[1]<20 && track.rgt[1]>100");
      TCut cutbt = Form("betah.pattern==4444 && betah.beta>0.3");
      TCut cutec = Form("ecal.energy[0]>5e3");
      TCut cutph = "(header.phpat&0x3e)!=0";
      cutt = cut00+cutl1+cutl9+cutin+cutbt;
      cut = cutt;
    }
    if(mode==5){
      ch.SetBranchStatus("ecalh", 0);
      ch.SetBranchStatus("betas", 0);
      TCut cut0 = Form("(ustat[0]&%d)", HasL2);
      TCut cuti = Form("track.qin>1.7 && track.qin<2.4 && track.csqy[0]>0 && track.csqy[0]<100 && rgt[0]*rgt2[0]>0 ");
      TCut cutb = Form("betah.pattern==4444 && betah.chi2t<100 && betah.chi2t>0 && betah.chi2c>0 && betah.chi2c<100 && betah.z[0]==2 &&betah.type==1 && betah.q>1.7");
      TCut cutph = "(header.phpat&0x3e)!=0";
      cutt = cut0+cuti+cutb+cutph;

      cut = cutt;
    }

    //ch.SetBranchStatus("trcls",  0);
    //ch.SetBranchStatus("betas" , 0);
    //ch.SetBranchStatus("trdhit", 0);
    //ch.SetBranchStatus("tofhit", 0);
    ch.SetBranchStatus("trdg"  , 0);
    ch.SetBranchStatus("vertex", 0);
    //ch.SetBranchStatus("ecalh",  0);

    TTree *tr = ch.CopyTree(cut);
    cout << "Nsel= " << tr->GetEntries() << endl;

    of.Write();
    of.Close();
    return 0;
}

#include <cstdlib>

int main(int argc, char *argv[])
{
  if (argc < 3) {
    cout << "mdsel [fname] [oname] [mode]"  << endl;
    cout << "mdsel [fname] [oname] [mode] [rid]" << endl;
    cout << "mdsel [fname] [oname] [mode] [if0] [if1]"  << endl;
    cout << "mode: 1:p 2:He 3:Z>2  +10: ismc" << endl;
    return 1;
  }

  Int_t ret = -1;
  if (argc == 3) ret = mdsel(argv[1], argv[2]);
  if (argc == 4) ret = mdsel(argv[1], argv[2], atoi(argv[3]));
  if (ret  == 0) cout << "Done" << endl;
  return ret;
}
