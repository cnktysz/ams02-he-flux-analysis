#include "trigger_eff.C"
#include "tof_eff.C"
#include "ext_track_eff.C"
#include "inn_track_eff.C"
void acceptance()
{
	bool save_plots     = 1;
	bool display_plots  = 0;
	TString ISS_address = "root/he_flux_iss/he_flux_iss_110419_short.root";
	TString MC_address  = "root/he_flux_mc/he_flux_mc_iter2_250419.root";	
	TFile *MC_data 	    = new TFile(MC_address);	
	TH1D *acceptance_mc = acceptance_MC(MC_data,1);
	TFile *ISS_data     = new TFile(ISS_address);	
	// Get Corrections from other detectors
	TH1D* acceptance_ 	   = (TH1D*)acceptance_mc->Clone(); // Init. effective acceptance
	TH1D* tof_correction 	   = tof_eff(ISS_data,MC_data,save_plots,display_plots);
	TH1D* ext_track_correction = ext_track_eff(ISS_data,MC_data,save_plots,display_plots);
	TH1D* inn_track_correction = inn_track_eff(ISS_data,MC_data,save_plots,display_plots);
	TH1D* total_correction     = (TH1D*)tof_correction->Clone();
	// Compute Effective Acceptance
	double ext,inn,tof,acc,corr;	
	for(int i=1;i<=acceptance_->GetNbinsX();i++)
	{
		if(i<11) ext = 1;
		else ext=ext_track_correction->GetBinContent(i);
		if(i==1) acc=1;
		else acc=acceptance_mc->GetBinContent(i);
		tof  = tof_correction->GetBinContent(i);
		inn  = inn_track_correction->GetBinContent(i);	
		corr = ext*inn*tof;
		total_correction->SetBinContent(i,corr);
		acceptance_->SetBinContent(i,corr*acc);
	}
	// Plot Total Correction
	TCanvas *canvas_correction = new TCanvas("canvas_correction","Correction",600,400);
	total_correction->Draw("HIST");
	//HIST Settings
	gROOT->Reset();
	gPad->SetLogx();
	total_correction->SetMarkerStyle(20);
	total_correction->SetMarkerColor(kBlack);
	total_correction->SetLineColor(kBlack);
	total_correction->SetMarkerSize(0.6);
	total_correction->SetTitle("Total Correction vs Rigidity");
	total_correction->SetTitleFont(43);
	total_correction->SetTitleSize(10);
	//HIST Settings X-axis
	total_correction->GetXaxis()->SetTitle("Rigidity [GV]");
	total_correction->GetXaxis()->SetRangeUser(0.8,3000);
	total_correction->GetXaxis()->SetTitleOffset(1.0);
	total_correction->GetXaxis()->SetTitleFont(43);
	total_correction->GetXaxis()->SetTitleSize(15);
	total_correction->GetXaxis()->SetLabelFont(43);
	total_correction->GetXaxis()->SetLabelSize(12);
	//HIST Settings Y-axis
	total_correction->GetYaxis()->SetTitle("Correction");	
	total_correction->GetYaxis()->SetRangeUser(0.85,1.0);
	total_correction->GetYaxis()->SetTitleOffset(1.0);
	total_correction->GetYaxis()->SetTitleFont(43);
	total_correction->GetYaxis()->SetTitleSize(15);
	total_correction->GetYaxis()->SetLabelFont(43);
	total_correction->GetYaxis()->SetLabelSize(10);
	// Save or Display
	if(save_plots) canvas_correction->SaveAs("png/TotalCorrection.png");
	if(canvas_correction && !display_plots) canvas_correction->Close();
	// Plot Effective Acceptance
	TCanvas *canvas_acceptance = new TCanvas("canvas_acceptance","Acceptance",600,400);
	gPad->SetLogx();
	acceptance_->Draw("P");
	//HIST Settings X-axis
	acceptance_->GetXaxis()->SetTitle("Rigidity [GV]");
	acceptance_->GetXaxis()->SetRangeUser(1,3000);
	acceptance_->GetXaxis()->SetTitleOffset(1.0);
	acceptance_->GetXaxis()->SetTitleFont(43);
	acceptance_->GetXaxis()->SetTitleSize(15);
	acceptance_->GetXaxis()->SetLabelFont(43);
	acceptance_->GetXaxis()->SetLabelSize(12);
	//HIST Settings Y-axis
	acceptance_->GetYaxis()->SetTitle("Acceptance [m^{2}sr]");	
	acceptance_->GetYaxis()->SetRangeUser(0.004,0.018);
	acceptance_->GetYaxis()->SetTitleOffset(1.0);
	acceptance_->GetYaxis()->SetTitleFont(43);
	acceptance_->GetYaxis()->SetTitleSize(15);
	acceptance_->GetYaxis()->SetLabelFont(43);
	acceptance_->GetYaxis()->SetLabelSize(10);
	//HIST Setings General
	acceptance_->SetMarkerStyle(24);
	acceptance_->SetMarkerColor(kRed);
	acceptance_->SetMarkerSize(0.6);
	acceptance_->SetTitle("Effective Acceptance vs Rigidity");
	acceptance_->SetTitleFont(43);
	acceptance_->SetTitleSize(10);
	if(save_plots) canvas_acceptance->SaveAs("png/Acceptance.png");
	if(canvas_acceptance && !display_plots) canvas_acceptance->Close();
	
	TFile *file = new TFile("root/acceptance.root","RECREATE");
	acceptance_->Write("effective_acceptance");
	acceptance_mc->Write("mc_acceptance");	
	file->Close();
}
TH1D* acceptance_MC(TFile* data,bool save_plots)
{
	bool display_plots = 0;
	TH1D* gen = (TH1D*)data->Get("hexp");//hexp
	TH1D* gen_w = (TH1D*)gen->Clone();
	TH2F* selected_time = (TH2F*)data->Get("Counts_Trigger");
	TH1D* selected = selected_time->ProjectionY();
	TH2F* MigrationMatrix = (TH2F*)data->Get("MigrationMatrix");	
	TH1D* acceptance =(TH1D*)selected->Clone();
	for(int i=1;i<=selected->GetNbinsX();i++) acceptance->SetBinContent(i,0); //init. acceptance histogram
	TH1D* MigrationMatrix_x = MigrationMatrix->ProjectionX();
	TH1D* MigrationMatrix_y = MigrationMatrix->ProjectionY();
	TH1D* MMratio =(TH1D*) MigrationMatrix_x->Clone();		
	MMratio->Divide(MigrationMatrix_y);

	// Normalize Migration Matrix
	double total_events = MigrationMatrix->GetEntries();
	for(int i=1;i<=MigrationMatrix->GetNbinsX();i++)
		for(int j=1;j<=MigrationMatrix->GetNbinsY();j++)
			MigrationMatrix->SetBinContent(i,j,MigrationMatrix->GetBinContent(i,j)/total_events);

	TCanvas *canvas_migration_matrix= new TCanvas("canvas_migration_matrix","Migration Matrix",600,400);
	MigrationMatrix->SetTitle("Migration Matrix");
	MigrationMatrix->GetXaxis()->SetTitle("Generated Rigidity [GV]");
	MigrationMatrix->GetYaxis()->SetTitle("Reconstructed Rigidity [GV]");	
	MigrationMatrix->GetXaxis()->SetRangeUser(1,1800);	
	gPad->SetLogx();
	MigrationMatrix->GetYaxis()->SetRangeUser(1,1800);	
	gPad->SetLogy();
	MigrationMatrix->Draw("COLZ");
	if(save_plots) canvas_migration_matrix->SaveAs("png/Acceptance/MigrationMatrix.png");
	canvas_migration_matrix->Close();
	
	TCanvas *canvas_migration_matrix_x= new TCanvas("canvas_migration_matrix_x","Migration Matrix X Proj.",600,400);
	MigrationMatrix_x->SetTitle("Migration Matrix X Proj.");
	MigrationMatrix_x->GetXaxis()->SetTitle("Generated Rigidity [GV]");
	MigrationMatrix_x->GetYaxis()->SetTitle("Events");	
	MigrationMatrix_x->GetXaxis()->SetRangeUser(1,1800);	
	gPad->SetLogx();
	//MigrationMatrix->GetYaxis()->SetRangeUser(1,1800);	
	gPad->SetLogy();
	MigrationMatrix_x->Draw("HIST");
	if(save_plots) canvas_migration_matrix_x->SaveAs("png/Acceptance/MigrationMatrix_x.png");
	canvas_migration_matrix_x->Close();

	TCanvas *canvas_migration_matrix_y= new TCanvas("canvas_migration_matrix_y","Migration Matrix Y Proj.",600,400);
	MigrationMatrix_y->SetTitle("Migration Matrix Y Proj.");
	MigrationMatrix_y->GetXaxis()->SetTitle("Reconstructed Rigidity [GV]");
	MigrationMatrix_y->GetYaxis()->SetTitle("events");	
	MigrationMatrix_y->GetXaxis()->SetRangeUser(1,1800);	
	gPad->SetLogx();
	//MigrationMatrix->GetYaxis()->SetRangeUser(1,1800);	
	gPad->SetLogy();
	MigrationMatrix_y->Draw("HIST");
	if(save_plots) canvas_migration_matrix_y->SaveAs("png/Acceptance/MigrationMatrix_y.png");
	canvas_migration_matrix_y->Close();

	TCanvas *canvas_migration_matrix_ratio= new TCanvas("canvas_migration_matrix_ratio","Migration Matrix Ratio",600,400);
	MMratio->SetTitle("Migration Matrix Ratio");
	MMratio->GetXaxis()->SetTitle("Rigidity [GV]");
	MMratio->GetYaxis()->SetTitle("Ratio");	
	MMratio->GetXaxis()->SetRangeUser(1,3000);	
	//MMratio>SetMarkerStyle(24);
	//MMratio->SetMarkerColor(kRed);
	//MMratio->SetMarkerSize(0.6);
	gPad->SetLogx();
	//MigrationMatrix->GetYaxis()->SetRangeUser(1,1800);	
	gPad->SetLogy();
	MMratio->SetMarkerStyle(24);
	MMratio->SetMarkerColor(kRed);
	MMratio->SetMarkerSize(0.6);
	MMratio->Draw("P");
	if(save_plots) canvas_migration_matrix_ratio->SaveAs("png/Acceptance/MigrationMatrix_ratio.png");
	canvas_migration_matrix_ratio->Close();

	// Plot Generated Events
	TCanvas *canvas_generated = new TCanvas("canvas_generated","Generated Events",600,400);
	gen->SetTitle("Generated Events");
	gen->GetXaxis()->SetTitle("Rigidity [GV]");
	gen->GetYaxis()->SetTitle("Events");	
	gen->GetXaxis()->SetRangeUser(1,3000);	
	gPad->SetLogx();
	gen->Draw("HIST");
	if(save_plots) canvas_generated->SaveAs("png/Acceptance/GeneratedEvents.png");
	canvas_generated->Close();
	// Plot Selected Events
	TCanvas *canvas_selected = new TCanvas("canvas_selected","Selected Events",600,400);
	selected->SetTitle("Selected Events");
	selected->GetXaxis()->SetTitle("Rigidity [GV]");
	selected->GetYaxis()->SetTitle("Events");	
	selected->GetXaxis()->SetRangeUser(1,3000);	
	gPad->SetLogx();
	selected->Draw("HIST");
	if(save_plots) canvas_selected->SaveAs("png/Acceptance/SelectedEvents.png");
	canvas_selected->Close();
	//Compute SEL/GEN ratio
	for(int i=1;i<=selected->GetNbinsX();i++)
	{
		double Ngen = gen->GetBinContent(i);
		double Nsel = selected->GetBinContent(i);
		//cout << "N: " << i << " Nsel: " << Nsel << " Ngen: " << Ngen << endl;//" Result: " << 47.78*Nsel/(Ngen) << endl;
		double corr = MMratio->GetBinContent(i);
		acceptance->SetBinContent(i,corr*47.78*Nsel/Ngen); //47.78
		double error = TMath::Sqrt(Nsel)/Ngen;
		acceptance->SetBinError(i,error);
	}
	// Plot MC Acceptance
	TCanvas *canvas_mc_acceptance = new TCanvas("canvas_mc_acceptance","Acceptance",600,400);
	gPad->SetLogx();
	acceptance->GetXaxis()->SetRangeUser(1,3000);
	//acceptance_mc->GetYaxis()->SetRangeUser(0.004,0.018); //0.004 0.018
	acceptance->GetYaxis()->SetTitleOffset(1.5);
	acceptance->GetXaxis()->SetTitleOffset(1.5);
	acceptance->SetMarkerStyle(24);
	acceptance->SetMarkerColor(kRed);
	acceptance->SetMarkerSize(0.6);
	acceptance->SetTitle("Monte Carlo Acceptance vs Rigidity");
	acceptance->GetXaxis()->SetTitle("Rigidity [GV]");
	acceptance->GetYaxis()->SetTitle("Acceptance [m^{2}sr]");	
	acceptance->Draw("P");
	if(save_plots) canvas_mc_acceptance->SaveAs("png/MC_Acceptance.png");
	if(canvas_mc_acceptance && !display_plots) canvas_mc_acceptance->Close();
	gROOT->Reset();
	return acceptance;
}


