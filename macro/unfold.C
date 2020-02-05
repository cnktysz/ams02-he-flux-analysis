double unfold(Double_t rigidity)
{
	bool save_plots = 0;
	bool display_plots = 0;

	TFile *MC_data = new TFile("root/he_flux_mc/he_flux_mc_220319.root");
	TFile *ISS_data = new TFile("root/he_flux_iss/he_flux_iss_110419.root");

	TH2F *ISS_event_rate_time = (TH2F*) ISS_data->Get("Counts_Trigger")->Clone("Counts_Trigger");
	TH1D *ISS_event_rate      = (TH1D*) ISS_event_rate_time->ProjectionY();
	TH2F *MigrationMatrix     = (TH2F*) MC_data->Get("MigrationMatrix")->Clone("MigrationMatrix");	
	TH1D *Rec_event_rate      = (TH1D*) MigrationMatrix->ProjectionY();
	TH1D *event_rate_ratio    = (TH1D*) Rec_event_rate->Clone();
	event_rate_ratio->Divide(ISS_event_rate);
	TGraph *gratio = SpFold::HtoG(event_rate_ratio);
    	
	enum { N = 8 };
    	Double_t xn[N] = { 2, 5, 10, 20, 50, 100, 200, 500 };
    	Double_t bb[2] = { 0, 0 };

    	SplFit::fLogX = 1;
    	SplFit::fLogY = 1;
    	SplFit::fBlxU = 1;

    	TF1 *func_ratio = SplFit::Fit(gratio, N, xn, bb, "q0", 1, 3000);
	/*        
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
	if(save_plots) canvas_event_rate->SaveAs("png/EventRateRatio.png");
	if(canvas_event_rate && !display_plots) canvas_event_rate->Close();
	*/
	ISS_data->Close();
	MC_data->Close();
	return func_ratio->Eval(rigidity);
}


