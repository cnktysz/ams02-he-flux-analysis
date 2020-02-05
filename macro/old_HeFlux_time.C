gROOT->Reset();
gErrorIgnoreLevel = kWarning; //Turn off Warning //Turn off Warningss
#include <stdio.h>
#include "trigger_eff.C"
#include "ExposureTime.C"
void HeFlux_time()
{
	// Set Plot and Save Settings
	bool save_plots    = 1;
	bool display_plots = 0;
	TString ISS_address = "root/he_flux_iss/he_flux_iss_full_150519.root"; 		  // Set input file locations
	// Read selected events (aka triggered events)
	TFile *ISS_data 	       = new TFile(ISS_address);
	TH2F  *selected_he_events_time = (TH2F*) ISS_data->Get("Counts_Trigger");
	TH1D  *selected_he_events      = (TH1D*) selected_he_events_time->ProjectionY();
	TH2F  *he_flux                 = (TH2F*) selected_he_events_time->Clone();        // Init. Histogram
	TH2F  *he_flux_raw             = (TH2F*) selected_he_events_time->Clone();        // Init. Histogram
	he_flux_raw->SetName("he_flux_raw");
	printf("Total He events: %e\n",selected_he_events->GetEntries());                 // Print total triggered events
	// Plot Selected He Events
	TCanvas *canvas_selected_events = new TCanvas("canvas_selected_events","He Events",600,400);
	gPad->SetLogy();
	selected_he_events_time->Draw("HIST");
	if(save_plots) canvas_selected_events->SaveAs("png/SelectedHeEvents.png");
	if(canvas_selected_events && !display_plots) canvas_selected_events->Close();
	// Choose binning
	Double_t bin_time[41] = {1.92,2.15,2.4,2.67,2.97,3.29,3.64,4.02,4.43,4.88,5.37,5.9,6.47,7.09,7.76,8.48,9.26,10.1,11,12,13,14.1,15.3,16.6,18,19.5,21.1,22.8,24.7,26.7,28.8,31.1,33.5,36.1,38.9,41.9,45.1,48.5,52.2,56.1,60.3};
	// Get results from subdetectors	
	TH2F  *Exposuretime         = (TH2F*) ISS_data->Get("ExposureIGRF");
	TH2F  *trigger_efficiency   = trigger_eff_time_calculate(ISS_data,0,0);	
	TFile *acceptance_file      = new TFile("root/acceptance_time.root");
	TH2F  *effective_acceptance_time = (TH2F*) acceptance_file->Get("effective_acceptance_time");
	
	// CALCULATE FLUX
	TString png_address,plot_title;
	enum{N=4};
	Double_t xn[N] = {3,5,10,30};
	Double_t bb[2] = {0,0};

	cout << "Number of Bartel Rotations: " << he_flux->GetNbinsX() << endl;
	// Loop for Bartel's Rotations
	for(int i=1;i<=he_flux->GetNbinsX();i++)
	{
		// Fit Effective Acceptance
		TH1D *effective_acceptance = (TH1D*) effective_acceptance_time->ProjectionY("",i,i);
		effective_acceptance->SetName("effective_aceptance");
		TGraph *g_acceptance = SpFold::HtoG(effective_acceptance);
		TF1 *f_acceptance    = SplFit::Fit(g_acceptance,N,xn,bb,"q0",1.92,60.3); //fit function
		// Plot Effective Acceptance and Fit
		TCanvas *c1 = new TCanvas;
		gPad->SetLogx();
		effective_acceptance->Draw("P");
		//HIST Settings X-axis
		effective_acceptance->GetXaxis()->SetTitle("Rigidity [GV]");
		effective_acceptance->GetXaxis()->SetRangeUser(1.99,60.1);
		effective_acceptance->GetXaxis()->SetTitleOffset(1.0);
		effective_acceptance->GetXaxis()->SetTitleFont(43);
		effective_acceptance->GetXaxis()->SetTitleSize(15);
		effective_acceptance->GetXaxis()->SetLabelFont(43);
		effective_acceptance->GetXaxis()->SetLabelSize(12);
		//HIST Settings Y-axis
		effective_acceptance->GetYaxis()->SetTitle("Acceptance [m^{2}sr]");	
		effective_acceptance->GetYaxis()->SetRangeUser(0.004,0.018);
		effective_acceptance->GetYaxis()->SetTitleOffset(1.0);
		effective_acceptance->GetYaxis()->SetTitleFont(43);
		effective_acceptance->GetYaxis()->SetTitleSize(15);
		effective_acceptance->GetYaxis()->SetLabelFont(43);
		effective_acceptance->GetYaxis()->SetLabelSize(10);
		//HIST Setings General
		effective_acceptance->SetMarkerStyle(24);
		effective_acceptance->SetMarkerColor(kRed);
		effective_acceptance->SetMarkerSize(0.6);	
		plot_title.Form("Effective Acceptance vs Rigidity [Time Bin: %d]",i);
		effective_acceptance->SetTitle(plot_title);
		effective_acceptance->SetTitleFont(43);
		effective_acceptance->SetTitleSize(10);
		// Fit Settings
		f_acceptance->SetLineColor(1);
		f_acceptance->SetLineWidth(1);
		f_acceptance->SetLineStyle(1);
		f_acceptance->Draw("SAME");
		png_address.Form("png/acceptance_time/acceptance_%d.png",i);
		c1->SaveAs(png_address);
		c1->Close();

		for(int j=1;j<=he_flux->GetNbinsY();j++)
		{
			double event = selected_he_events_time->GetBinContent(i,j);
			double acc = f_acceptance->Eval(effective_acceptance->GetBinCenter(j));
			double trig_eff = trigger_efficiency->GetBinContent(i,j);
			double exposure_time = Exposuretime->GetBinContent(i,j);
			double delta_rigidity = bin_time[j]-bin_time[j-1];
			double den = (acc*trig_eff*exposure_time*delta_rigidity);
			double flux = event/den;
			if(den) he_flux_raw->SetBinContent(i,j,flux);
			else  he_flux_raw->SetBinContent(i,j,0);
			flux = flux * TMath::Power(selected_he_events->GetBinCenter(i),2.7);
			if(den) he_flux->SetBinContent(i,j,flux);
			else he_flux->SetBinContent(i,j,0);
			he_flux->SetBinError(i,j,0);
			he_flux_raw->SetBinError(i,j,0);
		}
	}



	TFile *published_flux_data = new TFile("root/I2_HeliumFlux_pBin_B1036_Hawaii_heV27_Extended.root");
	TH2F  *published_flux      = (TH2F*) published_flux_data->Get("hh_HeliumFlux_TotErr");	
	// Take published bin projections
	TH1D *published_flux_1 = (TH1D*)published_flux->ProjectionX("",7,7);
	published_flux_1->SetName("published_flux_1");
	TH1D *published_flux_2 = published_flux->ProjectionX("",9,9);
	published_flux_2->SetName("published_he_flux_2");
/*	TH1D *he_flux_3 = he_flux_raw->ProjectionX("",6,6);
	he_flux_3->SetName("he_flux_3");
	TH1D *he_flux_4 = he_flux_raw->ProjectionX("",11,11);
	he_flux_4->SetName("he_flux_4");
	TH1D *he_flux_5 = he_flux_raw->ProjectionX("",18,18);
	he_flux_5->SetName("he_flux_5");
	TH1D *he_flux_6 = he_flux_raw->ProjectionX("",27,27);
	he_flux_6->SetName("he_flux_6");
	TH1D *he_flux_7 = he_flux_raw->ProjectionX("",36,36);
	he_flux_7->SetName("he_flux_7");
	TH1D *he_flux_8 = he_flux_raw->ProjectionX("",40,40);
	he_flux_8->SetName("he_flux_8");
*/	

	// Plot Published He Flux 3D
	TCanvas *canvas_published_he_flux = new TCanvas("canvas_published_he_flux","He Flux",600*2,400*2);
	gPad->SetTheta(12); 
  	gPad->SetPhi(240);
  	gPad->Update();
	gPad->SetLogy();
	published_flux->SetTitle("Published He Flux vs Time vs Rigidity");	
	published_flux->SetTitleFont(2);
	published_flux->SetTitleSize(0.2);	
	// Flux axis
	published_flux->GetZaxis()->SetLabelFont(2);
	published_flux->GetZaxis()->SetTitleFont(2);
	published_flux->GetZaxis()->SetTitle("Flux [m^{-2}sr^{-1}s^{-1}GV^{-1}]");
	published_flux->GetZaxis()->SetTitleOffset(1.4);
	// Rigidity Axis
	published_flux->GetYaxis()->SetRangeUser(1.99,10);
	published_flux->GetYaxis()->SetLabelFont(2);
	published_flux->GetYaxis()->SetTitle("Rigidity [GV]");
	published_flux->GetYaxis()->SetTitleFont(2);
	published_flux->GetYaxis()->SetTitleOffset(1.8);
	published_flux->GetYaxis()->SetMoreLogLabels();
	// TIME AXIS Settings
	published_flux->GetXaxis()->SetLabelFont(2);
	published_flux->GetXaxis()->SetLabelOffset(0.03);
	published_flux->GetXaxis()->SetTitle("");
	published_flux->GetXaxis()->SetTimeOffset(0,"UTC");
	published_flux->GetXaxis()->SetTimeDisplay(1);
	published_flux->GetXaxis()->SetNdivisions(8);	
	published_flux->GetXaxis()->SetTimeFormat("#splitline{%b}{%Y}");
	published_flux->Draw("SURF3");
	if(save_plots) canvas_published_he_flux->SaveAs("png/flux_time/HeFlux_published_3D_time.png");
	if(canvas_published_he_flux && !display_plots) canvas_published_he_flux->Close();

	// Plot He Flux 3D
	TCanvas *canvas_he_flux = new TCanvas("canvas_he_flux","He Flux",600*2,400*2);
	gPad->SetTheta(12); 
  	gPad->SetPhi(240);
  	gPad->Update();
	gPad->SetLogy();
	he_flux_raw->SetTitle("He Flux vs Time vs Rigidity");	
	he_flux_raw->SetTitleFont(2);
	he_flux_raw->SetTitleSize(0.2);	
	// Flux axis
	he_flux_raw->GetZaxis()->SetLabelFont(2);
	he_flux_raw->GetZaxis()->SetTitleFont(2);
	he_flux_raw->GetZaxis()->SetTitle("Flux [m^{-2}sr^{-1}s^{-1}GV^{-1}]");
	he_flux_raw->GetZaxis()->SetTitleOffset(1.4);
	// Rigidity Axis
	he_flux_raw->GetYaxis()->SetRangeUser(1.99,10);
	he_flux_raw->GetYaxis()->SetLabelFont(2);
	he_flux_raw->GetYaxis()->SetTitle("Rigidity [GV]");
	he_flux_raw->GetYaxis()->SetTitleFont(2);
	he_flux_raw->GetYaxis()->SetTitleOffset(1.8);
	he_flux_raw->GetYaxis()->SetMoreLogLabels();
	// TIME AXIS Settings
	he_flux_raw->GetXaxis()->SetLabelFont(2);
	he_flux_raw->GetXaxis()->SetLabelOffset(0.03);
	he_flux_raw->GetXaxis()->SetTitle("");
	he_flux_raw->GetXaxis()->SetTimeOffset(0,"UTC");
	he_flux_raw->GetXaxis()->SetTimeDisplay(1);
	he_flux_raw->GetXaxis()->SetNdivisions(8);	
	he_flux_raw->GetXaxis()->SetTimeFormat("#splitline{%b}{%Y}");
	he_flux_raw->Draw("SURF3");
	if(save_plots) canvas_he_flux->SaveAs("png/flux_time/HeFlux3D_time.png");
	if(canvas_he_flux && !display_plots) canvas_he_flux->Close();


	// Take bin projections
	TH1D *he_flux_1 = he_flux_raw->ProjectionX("",1,1);
	he_flux_1->SetName("he_flux_1");
	TH1D *he_flux_2 = he_flux_raw->ProjectionX("",3,3);
	he_flux_2->SetName("he_flux_2");
	TH1D *he_flux_3 = he_flux_raw->ProjectionX("",6,6);
	he_flux_3->SetName("he_flux_3");
	TH1D *he_flux_4 = he_flux_raw->ProjectionX("",11,11);
	he_flux_4->SetName("he_flux_4");
	TH1D *he_flux_5 = he_flux_raw->ProjectionX("",18,18);
	he_flux_5->SetName("he_flux_5");
	TH1D *he_flux_6 = he_flux_raw->ProjectionX("",27,27);
	he_flux_6->SetName("he_flux_6");
	TH1D *he_flux_7 = he_flux_raw->ProjectionX("",36,36);
	he_flux_7->SetName("he_flux_7");
	TH1D *he_flux_8 = he_flux_raw->ProjectionX("",40,40);
	he_flux_8->SetName("he_flux_8");

	// Plot He Flux 1st bin
	TCanvas *canvas_flux_1 = new TCanvas("canvas_flux_1","He Flux",360*4,360*2);
	gStyle->SetOptStat(0);
	TPad *pad1 = new TPad("pad1", "pad1" ,0,0.3,1,1.0);
	pad1->SetBottomMargin(0);
	pad1->SetGridx();
	pad1->SetGridy();
	pad1->Draw();
	pad1->cd();	
	// Plot Settings
	he_flux_1->SetTitleFont(2);
	he_flux_1->SetTitle(Form("He Flux vs Time [%.2f - %.2f] GV",bin_time[0],bin_time[1]));	
	he_flux_1->SetTitleSize(0.06);
	// Y Axis Settings
	he_flux_1->GetYaxis()->SetRangeUser(42,100);
	he_flux_1->GetYaxis()->SetLabelFont(2);	
	he_flux_1->GetYaxis()->SetTitle("Flux [m^{-2}sr^{-1}s^{-1}GV^{-1}]");
	he_flux_1->GetYaxis()->SetTitleFont(2);
	he_flux_1->GetYaxis()->SetTitleSize(0.05);
	he_flux_1->GetYaxis()->SetTitleOffset(0.6);
	// X Axis Settings
	he_flux_1->GetXaxis()->SetLabelFont(2);
	he_flux_1->GetXaxis()->SetLabelOffset(0.03);
	he_flux_1->GetXaxis()->SetTitle("");
	he_flux_1->GetXaxis()->SetTimeOffset(0,"gmt");
	he_flux_1->GetXaxis()->SetTimeDisplay(1);
	he_flux_1->GetXaxis()->SetNdivisions(13);	
	he_flux_1->GetXaxis()->SetTimeFormat("#splitline{%Y}{%b}"); 	
	// Draw Settings 
	he_flux_1->SetMarkerStyle(20);
	he_flux_1->SetMarkerColor(kRed);
	he_flux_1->SetMarkerSize(0.6);
	he_flux_1->Draw("P");
	// Draw Settings	
	published_flux_1->SetMarkerStyle(20);
	published_flux_1->SetMarkerColor(kBlue);
	published_flux_1->SetMarkerSize(0.6);
	published_flux_1->Draw("SAMEP");
	// Set Legend
	TLegend*lgd = new TLegend(0.8, 0.20, 0.9, 0.30);
	lgd->AddEntry(he_flux_1,"METU","lp");
	lgd->AddEntry(published_flux_1,"AMS","lp");
	//lgd->SetMargin(0.35);
	lgd->Draw();
	canvas_flux_1->cd();
	TPad *pad2 = new TPad("pad2","pad2",0,0.05,1,0.3);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.2);
	pad2->SetGridx();
	pad2->Draw();
	pad2->cd();
	// Define ratio plot
	TH1D* ratio1 = (TH1D*)he_flux_1->Clone();
	for(int i=1;i<=ratio1->GetNbinsX();i++)
	{ 
		double div = he_flux_1->GetBinContent(i)/published_flux_1->GetBinContent(i);
		if(div<10) ratio1->SetBinContent(i,div);
		else ratio1->SetBinContent(i,0);
	}
	// Y Axis Settings
	ratio1->SetTitle("");
	ratio1->GetYaxis()->SetTitle("#Phi_{METU} /#Phi_{AMS}");
	ratio1->GetYaxis()->SetRangeUser(0.8,2);
	ratio1->GetYaxis()->SetTitleFont(2);
	ratio1->GetYaxis()->SetTitleSize(0.15);
	ratio1->GetYaxis()->SetTitleOffset(0.2);
	ratio1->GetYaxis()->SetLabelFont(2);
	ratio1->GetYaxis()->SetLabelSize(0.1);
	// X Axis Settings	
	ratio1->GetXaxis()->SetLabelFont(2);
	ratio1->GetXaxis()->SetLabelSize(0.12);
	ratio1->GetXaxis()->SetLabelOffset(0.06);
	ratio1->GetXaxis()->SetTitle("");
	ratio1->GetXaxis()->SetTimeOffset(0,"gmt");
	ratio1->GetXaxis()->SetTimeDisplay(1);
	ratio1->GetXaxis()->SetNdivisions(13);	
	ratio1->GetXaxis()->SetTimeFormat("#splitline{%b}{%Y}"); 	
	// Draw Settings	
	ratio1->SetMarkerColor(kMagenta);
	ratio1->Draw("P");
	// Save settings
	if(save_plots) canvas_flux_1->SaveAs("png/flux_time/HeFlux_time_1.png");
	if(canvas_flux_1 && !display_plots) canvas_flux_1->Close();

	
	// Plot He Flux 2nd bin
	TCanvas *canvas_flux_2 = new TCanvas("canvas_flux_2","He Flux",360*4,360*2);
	gStyle->SetOptStat(0);
	TPad *pad1 = new TPad("pad1", "pad1" ,0,0.3,1,1.0);
	pad1->SetBottomMargin(0);
	pad1->SetGridx();
	pad1->SetGridy();
	pad1->Draw();
	pad1->cd();	
	// Plot Settings
	he_flux_2->SetTitleFont(2);
	he_flux_2->SetTitle(Form("He Flux vs Time [%.2f - %.2f] GV",bin_time[2],bin_time[3]));	
	he_flux_2->SetTitleSize(0.06);
	// Y Axis Settings
	he_flux_2->GetYaxis()->SetRangeUser(35,85);	
	he_flux_2->GetYaxis()->SetLabelFont(2);	
	he_flux_2->GetYaxis()->SetTitle("Flux [m^{-2}sr^{-1}s^{-1}GV^{-1}]");
	he_flux_2->GetYaxis()->SetTitleFont(2);
	he_flux_2->GetYaxis()->SetTitleSize(0.05);
	he_flux_2->GetYaxis()->SetTitleOffset(0.6);
	// X Axis Settings
	he_flux_2->GetXaxis()->SetLabelFont(2);
	he_flux_2->GetXaxis()->SetLabelOffset(0.03);
	he_flux_2->GetXaxis()->SetTitle("");
	he_flux_2->GetXaxis()->SetTimeOffset(0,"gmt");
	he_flux_2->GetXaxis()->SetTimeDisplay(1);
	he_flux_2->GetXaxis()->SetNdivisions(13);	
	he_flux_2->GetXaxis()->SetTimeFormat("#splitline{%Y}{%b}"); 	
	// Draw Settings 
	he_flux_2->SetMarkerStyle(20);
	he_flux_2->SetMarkerColor(kRed);
	he_flux_2->SetMarkerSize(0.6);
	he_flux_2->Draw("P");
	// Draw Settings	
	published_flux_2->SetMarkerStyle(20);
	published_flux_2->SetMarkerColor(kBlue);
	published_flux_2->SetMarkerSize(0.6);
	published_flux_2->Draw("SAMEP");
	// Set Legend
	TLegend*lgd = new TLegend(0.8, 0.20, 0.9, 0.30);
	lgd->AddEntry(he_flux_2,"METU","lp");
	lgd->AddEntry(published_flux_2,"AMS","lp");
	//lgd->SetMargin(0.35);
	lgd->Draw();
	canvas_flux_2->cd();
	TPad *pad2 = new TPad("pad2","pad2",0,0.05,1,0.3);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.2);
	pad2->SetGridx();
	pad2->Draw();
	pad2->cd();
	// Define ratio plot
	TH1D* ratio2 = (TH1D*)he_flux_2->Clone();
	for(int i=1;i<=ratio2->GetNbinsX();i++)
	{ 
		double div = he_flux_2->GetBinContent(i)/published_flux_2->GetBinContent(i);
		if(div<10) ratio2->SetBinContent(i,div);
		else ratio2->SetBinContent(i,0);
	}
	// Y Axis Settings
	ratio2->SetTitle("");
	ratio2->GetYaxis()->SetTitle("#Phi_{METU} /#Phi_{AMS}");
	ratio2->GetYaxis()->SetRangeUser(0.8,2);
	ratio2->GetYaxis()->SetTitleFont(2);
	ratio2->GetYaxis()->SetTitleSize(0.15);
	ratio2->GetYaxis()->SetTitleOffset(0.2);
	ratio2->GetYaxis()->SetLabelFont(2);
	ratio2->GetYaxis()->SetLabelSize(0.1);
	// X Axis Settings	
	ratio2->GetXaxis()->SetLabelFont(2);
	ratio2->GetXaxis()->SetLabelSize(0.12);
	ratio2->GetXaxis()->SetLabelOffset(0.06);
	ratio2->GetXaxis()->SetTitle("");
	ratio2->GetXaxis()->SetTimeOffset(0,"gmt");
	ratio2->GetXaxis()->SetTimeDisplay(1);
	ratio2->GetXaxis()->SetNdivisions(13);	
	ratio2->GetXaxis()->SetTimeFormat("#splitline{%b}{%Y}"); 	
	// Draw Settings	
	ratio2->SetMarkerColor(kMagenta);
	ratio2->Draw("P");
	// Save settings
	if(save_plots) canvas_flux_2->SaveAs("png/flux_time/HeFlux_time_2.png");
	if(canvas_flux_2 && !display_plots) canvas_flux_2->Close();






/*
	// Old Plot He Flux 2nd bin
	TCanvas *canvas_flux_2 = new TCanvas("canvas_flux_2","He Flux",360*4,360*1.5);
	gStyle->SetOptStat(0);
	he_flux_2->SetTitle(Form("He Flux vs Time [%.2f - %.2f] GV",bin_time[2],bin_time[3]));	
	he_flux_2->SetTitleFont(2);
	he_flux_2->SetMarkerStyle(20);
	he_flux_2->SetMarkerColor(kRed);
	he_flux_2->SetMarkerSize(0.6);
	he_flux_2->Draw("P");
	// Y Axis Settings
	he_flux_2->GetYaxis()->SetRangeUser(35,85);	
	he_flux_2->GetYaxis()->SetLabelFont(2);
	he_flux_2->GetYaxis()->SetTitleFont(2);
	he_flux_2->GetYaxis()->SetTitle("Flux [m^{-2}sr^{-1}s^{-1}GV^{-1}]");
	// TIME AXIS Settings
	he_flux_2->GetXaxis()->SetLabelFont(2);
	he_flux_2->GetXaxis()->SetLabelOffset(0.03);
	he_flux_2->GetXaxis()->SetTitle("");
	he_flux_2->GetXaxis()->SetTimeOffset(0,"gmt");
	he_flux_2->GetXaxis()->SetTimeDisplay(1);
	he_flux_2->GetXaxis()->SetNdivisions(20);	
	he_flux_2->GetXaxis()->SetTimeFormat("#splitline{%Y}{%b}");
	if(save_plots) canvas_flux_2->SaveAs("png/flux_time/HeFlux_time_2.png");
	if(canvas_flux_2 && !display_plots) canvas_flux_2->Close();
*/
	// Plot He Flux 3rd bin
	TCanvas *canvas_flux_3 = new TCanvas("canvas_flux_3","He Flux",360*4,360*1.5);
	gStyle->SetOptStat(0);
	he_flux_3->SetTitle(Form("He Flux vs Time [%.2f - %.2f] GV",bin_time[5],bin_time[6]));	
	he_flux_3->SetTitleFont(2);
	he_flux_3->SetMarkerStyle(20);
	he_flux_3->SetMarkerColor(kRed);
	he_flux_3->SetMarkerSize(0.6);
	he_flux_3->Draw("P");
	// Y axis Settings
	he_flux_3->GetYaxis()->SetLabelFont(2);
	he_flux_3->GetYaxis()->SetRangeUser(28,52);
	he_flux_3->GetYaxis()->SetTitleFont(2);
	he_flux_3->GetYaxis()->SetTitle("Flux [m^{-2}sr^{-1}s^{-1}GV^{-1}]");
	// TIME AXIS Settings
	he_flux_3->GetXaxis()->SetLabelFont(2);
	he_flux_3->GetXaxis()->SetLabelOffset(0.03);
	he_flux_3->GetXaxis()->SetTitle("");
	he_flux_3->GetXaxis()->SetTimeOffset(0,"gmt");
	he_flux_3->GetXaxis()->SetTimeDisplay(1);
	he_flux_3->GetXaxis()->SetNdivisions(20);	
	he_flux_3->GetXaxis()->SetTimeFormat("#splitline{%Y}{%b}");
	if(save_plots) canvas_flux_3->SaveAs("png/flux_time/HeFlux_time_3.png");
	if(canvas_flux_3 && !display_plots) canvas_flux_3->Close();

	// Plot He Flux 4th bin
	TCanvas *canvas_flux_4 = new TCanvas("canvas_flux_4","He Flux",360*4,360*1.5);
	gStyle->SetOptStat(0);
	he_flux_4->SetTitle(Form("He Flux vs Time [%.2f - %.2f] GV",bin_time[10],bin_time[11]));	
	he_flux_4->SetTitleFont(2);
	he_flux_4->SetMarkerStyle(20);
	he_flux_4->SetMarkerColor(kRed);
	he_flux_4->SetMarkerSize(0.6);
	he_flux_4->Draw("P");
	// Y axis settings 
	he_flux_4->GetYaxis()->SetLabelFont(2);
	he_flux_4->GetYaxis()->SetRangeUser(12.5,19.5);
	he_flux_4->GetYaxis()->SetTitleFont(2);
	he_flux_4->GetYaxis()->SetTitle("Flux [m^{-2}sr^{-1}s^{-1}GV^{-1}]");
	// TIME AXIS Settings
	he_flux_4->GetXaxis()->SetLabelFont(2);
	he_flux_4->GetXaxis()->SetLabelOffset(0.03);
	he_flux_4->GetXaxis()->SetTitle("");
	he_flux_4->GetXaxis()->SetTimeOffset(0,"gmt");
	he_flux_4->GetXaxis()->SetTimeDisplay(1);
	he_flux_4->GetXaxis()->SetNdivisions(20);	
	he_flux_4->GetXaxis()->SetTimeFormat("#splitline{%Y}{%b}");
	if(save_plots) canvas_flux_4->SaveAs("png/flux_time/HeFlux_time_4.png");
	if(canvas_flux_4 && !display_plots) canvas_flux_4->Close();

	// Plot He Flux 5th bin
	TCanvas *canvas_flux_5 = new TCanvas("canvas_flux_5","He Flux",360*4,360*1.5);
	gStyle->SetOptStat(0);
	he_flux_5->SetTitle(Form("He Flux vs Time [%.2f - %.2f] GV",bin_time[17],bin_time[18]));	
	he_flux_5->SetTitleFont(2);
	he_flux_5->SetMarkerStyle(20);
	he_flux_5->SetMarkerColor(kRed);
	he_flux_5->SetMarkerSize(0.6);
	he_flux_5->Draw("P");
	// Y axis Settings
	he_flux_5->GetYaxis()->SetLabelFont(2);
	he_flux_5->GetYaxis()->SetRangeUser(3.31,4.49);
	he_flux_5->GetYaxis()->SetTitleFont(2);
	he_flux_5->GetYaxis()->SetTitle("Flux [m^{-2}sr^{-1}s^{-1}GV^{-1}]");
	// TIME AXIS Settings
	he_flux_5->GetXaxis()->SetLabelFont(2);
	he_flux_5->GetXaxis()->SetLabelOffset(0.03);
	he_flux_5->GetXaxis()->SetTitle("");
	he_flux_5->GetXaxis()->SetTimeOffset(0,"gmt");
	he_flux_5->GetXaxis()->SetTimeDisplay(1);
	he_flux_5->GetXaxis()->SetNdivisions(20);	
	he_flux_5->GetXaxis()->SetTimeFormat("#splitline{%Y}{%b}");
	if(save_plots) canvas_flux_5->SaveAs("png/flux_time/HeFlux_time_5.png");
	if(canvas_flux_5 && !display_plots) canvas_flux_5->Close();

	// Plot He Flux 6th bin
	TCanvas *canvas_flux_6 = new TCanvas("canvas_flux_6","He Flux",360*4,360*1.5);
	gStyle->SetOptStat(0);
	he_flux_6->SetTitle(Form("He Flux vs Time [%.2f - %.2f] GV",bin_time[26],bin_time[27]));	
	he_flux_6->SetTitleFont(2);
	he_flux_6->SetMarkerStyle(20);
	he_flux_6->SetMarkerColor(kRed);
	he_flux_6->SetMarkerSize(0.6);
	he_flux_6->Draw("P");
	// Y axis Settings
	he_flux_6->GetYaxis()->SetLabelFont(2);
	he_flux_6->GetYaxis()->SetRangeUser(0.53,0.69);
	he_flux_6->GetYaxis()->SetTitleFont(2);
	he_flux_6->GetYaxis()->SetTitle("Flux [m^{-2}sr^{-1}s^{-1}GV^{-1}]");
	// TIME AXIS Settings
	he_flux_6->GetXaxis()->SetLabelFont(2);
	he_flux_6->GetXaxis()->SetLabelOffset(0.03);
	he_flux_6->GetXaxis()->SetTitle("");
	he_flux_6->GetXaxis()->SetTimeOffset(0,"gmt");
	he_flux_6->GetXaxis()->SetTimeDisplay(1);
	he_flux_6->GetXaxis()->SetNdivisions(20);	
	he_flux_6->GetXaxis()->SetTimeFormat("#splitline{%Y}{%b}");
	if(save_plots) canvas_flux_6->SaveAs("png/flux_time/HeFlux_time_6.png");
	if(canvas_flux_6 && !display_plots) canvas_flux_6->Close();

	// Plot He Flux 7th bin
	TCanvas *canvas_flux_7 = new TCanvas("canvas_flux_7","He Flux",360*4,360*1.5);
	gStyle->SetOptStat(0);
	he_flux_7->SetTitle(Form("He Flux vs Time [%.2f - %.2f] GV",bin_time[35],bin_time[36]));	
	he_flux_7->SetTitleFont(2);
	he_flux_7->SetMarkerStyle(20);
	he_flux_7->SetMarkerColor(kRed);
	he_flux_7->SetMarkerSize(0.6);
	he_flux_7->Draw("P");
	// Y axis settings
	he_flux_7->GetYaxis()->SetLabelFont(2);
	he_flux_7->GetYaxis()->SetRangeUser(0.08,0.11);
	he_flux_7->GetYaxis()->SetTitleFont(2);
	he_flux_7->GetYaxis()->SetTitle("Flux [m^{-2}sr^{-1}s^{-1}GV^{-1}]");
	// TIME AXIS Settings
	he_flux_7->GetXaxis()->SetLabelFont(2);
	he_flux_7->GetXaxis()->SetLabelOffset(0.03);
	he_flux_7->GetXaxis()->SetTitle("");
	he_flux_7->GetXaxis()->SetTimeOffset(0,"gmt");
	he_flux_7->GetXaxis()->SetTimeDisplay(1);
	he_flux_7->GetXaxis()->SetNdivisions(20);	
	he_flux_7->GetXaxis()->SetTimeFormat("#splitline{%Y}{%b}");
	if(save_plots) canvas_flux_7->SaveAs("png/flux_time/HeFlux_time_7.png");
	if(canvas_flux_7 && !display_plots) canvas_flux_7->Close();

	// Plot He Flux 8th bin
	TCanvas *canvas_flux_8 = new TCanvas("canvas_flux_8","He Flux",360*4,360*1.5);
	gStyle->SetOptStat(0);
	he_flux_8->SetTitle(Form("He Flux vs Time [%.2f - %.2f] GV",bin_time[39],bin_time[40]));	
	he_flux_8->SetTitleFont(2);
	he_flux_8->SetMarkerStyle(20);
	he_flux_8->SetMarkerColor(kRed);
	he_flux_8->SetMarkerSize(0.6);
	he_flux_8->Draw("P");
	// Y axis Settings
	he_flux_8->GetYaxis()->SetLabelFont(2);
	he_flux_8->GetYaxis()->SetRangeUser(0.0485,0.035);	
	he_flux_8->GetYaxis()->SetTitleFont(2);
	he_flux_8->GetYaxis()->SetTitle("Flux [m^{-2}sr^{-1}s^{-1}GV^{-1}]");
	// TIME AXIS Settings
	he_flux_8->GetXaxis()->SetLabelFont(2);
	he_flux_8->GetXaxis()->SetLabelOffset(0.03);
	he_flux_8->GetXaxis()->SetTitle("");
	he_flux_8->GetXaxis()->SetTimeOffset(0,"gmt");
	he_flux_8->GetXaxis()->SetTimeDisplay(1);
	he_flux_8->GetXaxis()->SetNdivisions(20);	
	he_flux_8->GetXaxis()->SetTimeFormat("#splitline{%Y}{%b}");
	if(save_plots) canvas_flux_8->SaveAs("png/flux_time/HeFlux_time_8.png");
	if(canvas_flux_8 && !display_plots) canvas_flux_8->Close();


	// Save Results
	TFile* file = new TFile("root/result_time.root","RECREATE");
	//selected_he_events->Write("selected_he_events");
	//selected_he_events_time->Write("selected_he_events_time");
	//Exposuretime->Write("Exposuretime");
	//Exposuretime_time->Write("Exposuretime_time");
	//trigger_efficiency->Write("trigger_efficiency");
	//trigger_efficiency_time->Write("trigger_efficiency_time");
	//effective_acceptance->Write("effective_acceptance");
	he_flux->Write("he_flux_2.7");
	//ratio->Write("he_flux_ratio");
	he_flux_raw->Write("he_flux");
	he_flux_1->Write("he_flux_1");
	he_flux_2->Write("he_flux_2");
	he_flux_3->Write("he_flux_3");
	he_flux_4->Write("he_flux_4");
	he_flux_5->Write("he_flux_5");
	he_flux_6->Write("he_flux_6");
	he_flux_7->Write("he_flux_7");
	he_flux_8->Write("he_flux_8");
	//published_flux->Write("published_flux_2.7");
	//published_flux_raw->Write("published_flux");
	file->Close();

	gROOT->Reset();
}
