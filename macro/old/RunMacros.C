gROOT->Reset();
#include "trigger_eff.C"
#include "tof_eff.C"
#include "ExposureTime.C"
#include "acceptance.C"
#include "ext_track_eff.C"
// Global Variables
//-------------------
void HeFlux()
{
	bool save_plots = 1;
	bool display_plots = 0;
	TFile *MC_data = new TFile("mc.root");
	TFile *ISS_data = new TFile("cenk.root");
	TH1D *Exposuretime = ExposureTime(ISS_data,save_plots,display_plots);
	TH1D *trigger_efficiency = trigger_eff_calculate(ISS_data,0,0);
	trigger_eff(ISS_data,MC_data,save_plots,display_plots);
	//TH1D *tof_correction = tof_eff(ISS_data,MC_data,save_plots,display_plots);
	TH1D *effective_acceptance = acceptance(ISS_data,MC_data,save_plots,display_plots);
	//TH1D *ext_track_correction = ext_track_eff(ISS_data,MC_data,save_plots,display_plots);
	return 0;

	TH2F *selected_he_events_time = (TH2F*)ISS_data->Get("Counts_Trigger")->Clone("Counts_Trigger");
	TH1D* selected_he_events = (TH1D*)selected_he_events_time->ProjectionY();
	TH1D* he_flux = (TH1D*) selected_he_events->Clone();

	TCanvas *canvas_selected_events = new TCanvas("canvas_selected_events","He Events",600,400);
	gPad->SetLogx();
	selected_he_events->Draw("HIST");
	if(save_plots) canvas_selected_events->SaveAs("png/SelectedHeEvents.png");
	if(canvas_selected_events && !display_plots) canvas_selected_events->Close();



	Double_t bin[76]={
		0.8,1,1.16,1.33,1.51,1.71,1.92,2.15,2.4,
		2.67,2.97,3.29,3.64,4.02,4.43,4.88,5.37,
		5.9,6.47,7.09,7.76,8.48,9.26,10.1,11,12,
		13,14.1,15.3,16.6,18,19.5,21.1,22.8,24.7,
		26.7,28.8,31.1,33.5,36.1,38.9,41.9,45.1,
		48.5,52.2,56.1,60.3,64.8,69.7,74.9,80.5,
		86.5,93,100,108,116,125,135,147,160,175,
		192,211,233,259,291,330,379,441,525,643,
		822,1130,1800,3000,8000};
	Double_t old_bin[69];
	for(int i=0;i<69;i++) = old_bin[i]=bin[i+6];
	TH1D *published_flux = new TH1D("published_flux","He Flux vs. Rigidity",68,old_bin);
	Double_t published_flux_array[68] = {6.031e+1,5.657e+1,5.174e+1,4.694e+1,4.176e+1,3.650e+1,
				       3.145e+1,2.671e+1,2.250e+1,1.876e+1,1.555e+1,
				       1.282e+1,1.054e+1,8.646e+0,7.081e+0,5.786e+0,4.740e+0,
				       3.866e+0,3.138e+0,2.566e+0,2.100e+0,1.711e+0,1.395e+0,
				       1.134e+0,9.241e-1,7.485e-1,6.100e-1,4.963e-1,4.045e-1,
				       3.295e-1,2.686e-1,2.194e-1,1.797e-1,1.469e-1,1.198e-1,
				       9.824e-2,8.071e-2,6.566e-2,5.415e-2,4.446e-2,3.642e-2,
				       2.984e-2,2.456e-2,2.015e-2,1.654e-2,1.357e-2,1.113e-2,
				       9.068e-3,7.328e-3,5.977e-3,4.905e-3,3.898e-3,3.104e-3,
				       2.433e-3,1.888e-3,1.481e-3,1.121e-3,8.621e-4,6.312e-4,
				       4.614e-4,3.199e-4,2.158e-4,1.397e-4,8.738e-5,4.682e-5,2.288e-5,7.980e-6,2.147e-6};	
	for(int i=1;i<published_flux->GetNbinsX();i++)
	{
	//cout <<"Bin: " << old_bin[i-1] << "-" << old_bin[i] << "P flux : " << published_flux_array[i-1] << endl;
	double flux_ = published_flux_array[i-1] * TMath::Power((old_bin[i-1]+old_bin[i])/2,2.7);
	published_flux->SetBinContent(i,flux_);
	}
	for(int i=1;i<=he_flux->GetNbinsX();i++)
	{
	he_flux->SetBinContent(i,0);
	double event = selected_he_events->GetBinContent(i);
	double acceptance = effective_acceptance->GetBinContent(i);
	double trig_eff = trigger_efficiency->GetBinContent(i);
	double exposure_time = Exposuretime->GetBinContent(i);
	double delta_rigidity = bin[i]-bin[i-1];
	double den = (acceptance*trig_eff*exposure_time*delta_rigidity);
	double flux = event/(acceptance*trig_eff*exposure_time*delta_rigidity);
	flux = flux * TMath::Power((bin[i-1]+bin[i])/2,2.7);
	//cout << "Bin: "<< i  << " -> " << bin[i] <<" Flux: " << flux << endl;
	if(den) he_flux->SetBinContent(i,flux);
	}	
	TCanvas *canvas_HeFlux = new TCanvas("canvas_HeFlux","He Flux",600,400);
	TPad *pad1 = new TPad("pad1", "pad1" ,0,0.3,1,1.0);
	pad1->SetBottomMargin(0);
	pad1->SetGridx();
	pad1->SetLogx();
	pad1->SetGridy();
	pad1->SetLogy();
	pad1->Draw();
	pad1->cd();
	he_flux->SetTitle("He Flux vs Rigidity");	
	he_flux->GetXaxis()->SetRangeUser(1.92,3000);
	he_flux->GetXaxis()->SetTitleOffset(1.5);
	he_flux->GetXaxis()->SetTitle("Rigidity [GV]");
	he_flux->GetYaxis()->SetRangeUser(3e+2,1e+4);
	he_flux->GetYaxis()->SetTitleOffset(1.0);
	he_flux->GetYaxis()->SetTitleSize(15);
	he_flux->GetYaxis()->SetTitleFont(43);
	he_flux->GetYaxis()->SetTitle("Flux R^{2.7} [m^{-2}sr^{-1}s^{-1}GV^{-1}]");
	he_flux->SetMarkerStyle(20);
	he_flux->SetMarkerColor(kRed);
	he_flux->SetMarkerSize(0.6);
	he_flux->Draw("E0P");
	published_flux->SetMarkerStyle(20);
	published_flux->SetMarkerColor(kBlue);
	published_flux->SetMarkerSize(0.6);
	published_flux->Draw("SAMEP");
	TLegend*lgd = new TLegend(0.7, 0.15, 0.90, 0.35);
	lgd->AddEntry(he_flux,"METU","lp");
	lgd->AddEntry(published_flux,"AMS","lp");
	lgd->SetMargin(0.35);
	lgd->Draw();
	canvas_HeFlux->cd();
	TPad *pad2 = new TPad("pad2","pad2",0,0.05,1,0.3);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.2);
	pad2->SetGridx();
	pad2->SetLogx();
	pad2->Draw();
	pad2->cd();
	// Define ratio plot
	TH1D* ratio = (TH1D*)he_flux->Clone();
	double div=1;
	for(int i=1;i<=he_flux->GetNbinsX();i++)
	{
		if(i<=6 || i==75 || i==74) div = 1;
		else div = he_flux->GetBinContent(i)/published_flux->GetBinContent(i-6);
		//cout << "Bin: " << i<<" Ratio: " << div << endl; 
		ratio->SetBinContent(i,div);
	}
	ratio->Draw("P");
	// Histogram Settings
	ratio->SetTitle("");
	ratio->GetYaxis()->SetTitle("#Phi_{METU} /#Phi_{AMS}");
	ratio->GetYaxis()->SetRangeUser(0.8,1.2);
	ratio->GetXaxis()->SetRangeUser(1.92,3000);
	ratio->GetYaxis()->SetTitleSize(15);
	ratio->GetYaxis()->SetTitleFont(43);
	ratio->GetYaxis()->SetTitleOffset(1.0);
	ratio->GetYaxis()->SetLabelFont(43);
	ratio->GetYaxis()->SetLabelSize(8);
	ratio->GetXaxis()->SetTitle("Rigidity [GV]");
	ratio->GetXaxis()->SetTitleSize(12);
	ratio->GetXaxis()->SetTitleFont(43);
	ratio->GetXaxis()->SetTitleOffset(4.);
	ratio->GetXaxis()->SetLabelFont(43);
	ratio->GetXaxis()->SetLabelSize(12);
	ratio->SetMarkerColor(kMagenta);
	//display_plots=1;
	if(save_plots) canvas_HeFlux->SaveAs("png/HeFlux.png");
	if(canvas_HeFlux && !display_plots) canvas_HeFlux->Close();
}


