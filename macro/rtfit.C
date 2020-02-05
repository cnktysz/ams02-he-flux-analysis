void rtfit(TString shn = "hist",
		Double_t pwr = 2.7)
{
	TFile *MC_data = new TFile("root/he_flux_mc/he_flux_mc_iter1_190419.root");
	TFile *ISS_data = new TFile("root/he_flux_iss/he_flux_iss_110419_short.root");

	TH1F *hist0 = (TH1F*) ISS_data->Get("hexp")->Clone("hexp");
	TH2F *hist1 = (TH2F*) ISS_data->Get(shn+"11")->Clone(shn+"11");
	TH1D *rate_ratio  = unfold(ISS_data,MC_data);
	/*
	for(int i=1;i<=hist1->GetNbinsX();i++)
		for(int j=1;j<=hist1->GetNbinsY();j++) 
			hist1->SetBinContent(i,j,hist1->GetBinContent(i,j)*rate_ratio->GetBinContent(j)); 
	*/
	TH1D   *hrt = fproj(hist1, hist0, rate_ratio,"hrt");
	TGraph *grt = SpFold::HtoG(hrt);

	enum { N = 8 };
	Double_t xn[N] = {2, 5, 10, 20, 50, 100, 200, 500 };
	Double_t bb[2] = { 0, 0 };

	SplFit::fLogX = 1;
	SplFit::fLogY = 1;
	SplFit::fBlxU = 1;

	TF1 *func = SplFit::Fit(grt, N, xn, bb, "q0", 1, 3000);
	TF1  fpw("fpw", "x^[0]"); fpw.SetParameter(0, -pwr);
	SpFold::Scale(hrt, pwr, &fpw);
	grt = SpFold::HtoG(hrt);

	// Plot New Weights
	TCanvas *c1 = new TCanvas;
	grt->SetMarkerStyle(20);
	grt->SetMarkerColor(2);
	grt->SetLineColor(1);
	c1->SetLogx();
	grt->Draw("ap");
	SpFold::Plot(func, pwr)->Draw("l");
	c1->SaveAs("png/rtfit/ratiofit.png");
	// Print fit result
	for (Int_t i = 0; i < N+2; i++) 
		cout << Form(" %6.3f,", func->GetParameter(i+N));
	cout << endl;

	TFile of("root/rtfit.root", "recreate");
	gROOT->GetList()->Write();

	// Compare with previous fit
	enum { N_ = 9 };
	Double_t xn_ [N_] = { 1, 2, 5, 10, 20, 50, 100, 200, 500 }; // node positions
	Double_t par_[N_+2] = {0.958,  0.838,  0.190, -0.483, -1.254,
		-2.340, -3.172, -4.007, -5.099,  0.592, -2.666}; // fit result
	frt = new TF1("spfun", SplFit::SpFunc, 1, 5000, N_*2+2);
	for (Int_t i = 0; i < N_;   i++) frt->SetParameter(i,   xn_ [i]);
	for (Int_t i = 0; i < N_+2; i++) frt->SetParameter(i+N, par_[i]);
	// Plot old fit and new
	TCanvas *c2 = new TCanvas;
	c2->SetLogx();
	grt->Draw("ap");
	SpFold::Plot(frt, pwr)->Draw("l");
	c2->SaveAs("png/rtfit/old_new_func.png");
}

TH1D *fproj(TH2F *hist,
		TH1F *hexp, TH1D* rate_ratio, const char *hname = "hprj", Int_t mode = 1,
		Int_t col = 2,
		Int_t sty = 21)
{
	TH1D *hprj = hist->ProjectionX(hname);
	hprj->Reset();

	for (Int_t i = 0; i < hist->GetNbinsX(); i++) {
		Double_t x = TMath::Abs(hist->GetXaxis()->GetBinCenter(i+1));
		Double_t w =            hist->GetXaxis()->GetBinWidth (i+1); 
		Double_t r = (mode == 1) ? x : 1/x;
		Int_t    j = hist->GetYaxis()->FindBin(r);
		Double_t t = hexp->GetBinContent(j);
		Double_t c = hist->Integral(i+1, i+1, 1, j-1)*rate_ratio->GetBinContent(i+1);
		if (x < 1) continue;

		if (t == 0 && r > 1e3) t = hexp->GetBinContent(hexp->GetNbinsX());
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

TH1D *unfold(TFile *ISS_data,TFile *MC_data)
{
	bool save_plots 	  = 1;
	bool display_plots        = 0;
	TH2F *ISS_event_rate_time = (TH2F*) ISS_data->Get("Counts_Trigger");
	TH1D *ISS_event_rate      = (TH1D*) ISS_event_rate_time->ProjectionY();
	TH2F *MigrationMatrix     = (TH2F*) MC_data->Get("MigrationMatrix");	
	TH1D *Rec_event_rate      = (TH1D*) MigrationMatrix->ProjectionY();
	TH1D *event_rate_ratio    = (TH1D*) Rec_event_rate->Clone();
	event_rate_ratio->Divide(ISS_event_rate);
	TGraph *gratio = SpFold::HtoG(event_rate_ratio);

	enum { N = 9 };
	Double_t xn[N] = {1, 2, 5, 10, 20, 50, 100, 200, 500 };
	Double_t bb[2] = { 0, 0 };

	SplFit::fLogX = 1;
	SplFit::fLogY = 1;
	SplFit::fBlxU = 1;
	TF1 *func_ratio = SplFit::Fit(gratio, N, xn, bb, "q0", 1, 3000);

	TCanvas *canvas_event_rate = new TCanvas("canvas_event_rate","Event Rate Ratio",600,400);
	gPad->SetLogx();
	event_rate_ratio->GetXaxis()->SetRangeUser(1,3000);
	//event_rate_ratio->GetYaxis()->SetRangeUser(0.5,1.5); 
	event_rate_ratio->GetYaxis()->SetTitleOffset(1.5);	
	event_rate_ratio->GetXaxis()->SetTitleOffset(1.5);
	event_rate_ratio->SetMarkerStyle(20);
	event_rate_ratio->SetMarkerColor(kBlue);
	event_rate_ratio->SetMarkerSize(0.6);
	event_rate_ratio->SetTitle("Event Rate Ratio vs Rigidity");
	event_rate_ratio->GetXaxis()->SetTitle("Rigidity [GV]");
	event_rate_ratio->GetYaxis()->SetTitle("Ratio");	
	event_rate_ratio->Draw("P");
	SpFold::Plot(func_ratio)->Draw("l");
	if(save_plots) canvas_event_rate->SaveAs("png/rtfit/EventRateRatio.png");
	if(canvas_event_rate && !display_plots) canvas_event_rate->Close();

	return event_rate_ratio;
}


