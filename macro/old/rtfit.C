
void rtfit(const char *fname = "he_flux_iss_030419.root", TString shn = "hist",
        Double_t pwr = 2.7)
{
    TFile f(fname);
    if (!f.IsOpen()) return;

    gROOT->cd();

    TH2F *hist_ = (TH2F *)f.Get("ExposureIGRF");
    TH1F *hist0 = (TH1F *)hist_->ProjectionY();
    TH2F *hist1 = (TH2F *)f.Get("Counts_Trigger");
    TH1D   *hrt = fproj(hist1, hist0, "hrt");
    TGraph *grt = SpFold::HtoG(hrt);

    enum { N = 8 };
    Double_t xn[N] = { 2, 5, 10, 20, 50, 100, 200, 500 };
    Double_t bb[2] = { 0, 0 };

    SplFit::fLogX = 1;
    SplFit::fLogY = 1;
    SplFit::fBlxU = 1;
/*
    TF1 *func = SplFit::Fit(grt, N, xn, bb, "q0", 2, 3000);
    TF1  fpw("fpw", "x^[0]"); fpw.SetParameter(0, -pwr);

    SpFold::Scale(hrt, pwr, &fpw);
*/    grt = SpFold::HtoG(hrt);
    grt->SetMarkerStyle(21);
    grt->SetMarkerColor(2);
    grt->SetLineColor(2);

    TCanvas *c1 = new TCanvas;
    c1->SetLogx();
    grt->Draw("ap");
  //  SpFold::Plot(func, pwr)->Draw("l");
/*
    for (Int_t i = 0; i < N+2; i++) 
        cout << Form(" %6.3f,", func->GetParameter(i+N));
    cout << endl;
*/
    //TFile of("rtfit.root", );
 //   gROOT->GetList()->Write();

    c1->SaveAs("png/rtfit.png");
}

TH1D *fproj(TH2F *hist,
        TH1F *hexp, const char *hname = "hprj", Int_t mode = 1,
        Int_t col = 2,
        Int_t sty = 21)
{
    TH1D *hprj = hist->ProjectionY(hname);
    hprj->Reset();

    for (Int_t i = 0; i < hist->GetNbinsY(); i++) {
        Double_t x = TMath::Abs(hist->GetYaxis()->GetBinCenter(i+1));
        Double_t w =            hist->GetYaxis()->GetBinWidth (i+1); 
        Double_t r = (mode == 1) ? x : 1/x;
        Int_t    j = hexp->FindBin(r);
        Double_t t = hexp->GetBinContent(j);
        Double_t c = hist->Integral(i+1, i+1, 1, j-1);
        if (x < 1) continue;

        if (t == 0 && r > 1e3) t = hexp->GetBinContent(hexp->GetNbinsY());
        if (t < 2 && r < 1) { t = 1; c = 1; }
        if (c > 0 && t > 0) {
            hprj->SetBinContent(i+1, c/t/w);
            hprj->SetBinError  (i+1, TMath::Sqrt(c)/t/w);
        }
    }

    hprj->SetMarkerStyle(sty);
    hprj->SetMarkerColor(col);
    hprj->SetLineColor  (col);

    return hprj;
}
