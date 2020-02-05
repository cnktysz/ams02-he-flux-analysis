
TString GetROOTFileName(TString fn, Int_t run){
  return gSystem->GetFromPipe(Form(" cat %s | grep %d | head -n 1", fn.Data(), run));
}
TGraph *GetMean(TH2F *h, int isy = 1){
  TGraph *gr = new TGraph;
  int np = 0;
  int nb = isy? h->GetNbinsX(); h->GetNbinsY();
  for(int i=1;i<nb;i++){
    float x = isy? h->GetXaxis()->GetBinCenter(i):
      h->GetYaxis()->GetBinCenter(i);

    TH1D* hpj= isy? h->ProjectionY("hpj", i,i):
      h->ProjectionX("hpj", i,i);
    float y = hpj->GetMean();
    delete hpj;
    gr->SetPoint(np++, x, y);
  }
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(0.5);
  return gr;
}

void palette(void)
{
    gStyle->SetLabelFont(82,"xyz");
    gStyle->SetLegendFont(82);
    gStyle->SetTitleFont(82,"");
    gStyle->SetTitleFont(82, "xyz");
    gStyle->SetTitleOffset(1.2,"xyz");
    gStyle->SetFillColor(10);

    TColor* Gradient = new TColor();
    Gradient->InitializeColors();
    const Int_t NRGBs = 5;
    const Int_t NCont = 80;

    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red  [NRGBs] = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue [NRGBs] = { 0.51, 1.00, 0.12, 0.00, 0.00 };

    Gradient->CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
}

void rootlogon(void)
{
    gROOT->GetColor(3)->SetRGB(0.0, 0.6, 0.0);
    gROOT->GetColor(5)->SetRGB(0.7, 0.5, 0.0);
    gROOT->GetColor(6)->SetRGB(0.7, 0.0, 0.7);
    gROOT->GetColor(7)->SetRGB(0.0, 0.7, 0.8);

    gSystem->Load("ntuple_slc6_PG.so");
    
    gSystem->SetIncludePath("-I$ROOTSYS/include "
            "-I$AMSSRC/include "
            "-D_PGTRACK_");
    gStyle->SetOptStat(0);
    palette();
}

void setlgd(TLegend *lgd)
{
    lgd->SetFillColor(0);
    lgd->SetBorderSize(0.05);
    lgd->SetTextSize(0.045);
    lgd->SetMargin(0.5);
}
void SetHist(TH1 *h, int sty,int col)
{
    h->SetLineColor(col);
    h->SetMarkerColor(col);
    h->SetMarkerStyle(sty);
    h->SetLineWidth(2);
}
TH2 *geth2(TString fn, TString hn, TString cn="hist_clone", TString htt="",TString xtt="", TString ytt=""){

    TFile f(fn);
    if(!f.IsOpen())return 0;
    TH2 *h = (TH2*)f.Get(hn);
    if(!h){
        cout<<"Hist "<<hn<<" not found in" <<fn<<" . "<<endl;
        return 0;}
    gROOT->cd();
    TH2 *h2 =  (TH2*)h->Clone(cn);
    if(htt.Length()>0)h2->SetTitle(htt);
    else h2->SetTitle(cn);
    if(xtt.Length()>0)h2->SetXTitle(xtt);
    if(ytt.Length()>0)h2->SetYTitle(ytt);
    h2->Sumw2();
    return h2;
}
TH1 *geth(TString fn, TString hn, TString cn="hist_clone", int col=1){

    TFile f(fn);
    if(!f.IsOpen())return 0;
    TH1 *h = (TH1*)f.Get(hn);
    if(!h){
        cout<<"Hist not found. "<<hn<<endl;
        return 0;}

    gROOT->cd();
    TH1 *h2 =  (TH1*)h->Clone(cn);

    h2->SetLineWidth(2);
    return h2;
}
void StripList(TString Input, TString &EOSDir, TString &fname){
    Int_t j = Input.Last('/');
    EOSDir = Input(0,j);
    fname =Input(j+1,Input.Length()-1);
}

