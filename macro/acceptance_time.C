#include "trigger_eff.C"
#include "tof_eff.C"
#include "ext_track_eff.C"
#include "inn_track_eff.C"
void acceptance_time()
{
	bool save_plots     = 1;
	bool display_plots  = 0;
	TString png_address,plot_title;
	TString ISS_address = "root/he_flux_iss/he_flux_iss_full_230519.root";
	TFile *ISS_data     = new TFile(ISS_address);			
	TH2F* acceptance_time      = (TH2F*)ISS_data->Get("Counts_Trigger");
	acceptance_time->SetName("acceptance_time");
	for(int i=1;i<=acceptance_time->GetNbinsX();i++)
	{
		cout << "Bartel: " << i << endl;
		if(i==46 || i==47 || i==48) // Fill No Data time bins with zeros
		{
			for(int j=1;j<=acceptance_time->GetNbinsY();j++)
				acceptance_time->SetBinContent(i,j,0);
		}
		else
		{
			// Read MC data and compute mc acceptance
			TString MC_address  = Form("root/he_flux_mc_time/he_flux_mc_time_%d.root",i);	
			TFile *MC_data 	    = new TFile(MC_address);	
			TH1D *acceptance_mc = acceptance_MC(MC_data,1,i);
			TH1D *acceptance_   = (TH1D*)acceptance_mc->Clone(); // Init. effective acceptance
			acceptance_->SetName("acceptance_");
			// Get Corrections from other detectors
			save_plots = 1;
			TH1D* tof_correction 	   = tof_eff(ISS_data,MC_data,save_plots,display_plots,i);
			TH1D* ext_track_correction = ext_track_eff(ISS_data,MC_data,save_plots,display_plots,i);
			TH1D* inn_track_correction = inn_track_eff(ISS_data,MC_data,save_plots,display_plots,i);
			TH1D* total_correction     = (TH1D*)tof_correction->Clone();
			tof_correction->SetName("tof_correction");
			ext_track_correction->SetName("ext_track_correction");
			inn_track_correction->SetName("inn_track_correction");
			total_correction->SetName("total_correction");
			// Compute Effective Acceptance
			double ext,inn,tof,acc,corr;
			double ext_error,inn_error,tof_error,acc_error,total_error;	
			for(int k=1;k<=acceptance_->GetNbinsX();k++)
			{
				ext  = ext_track_correction->GetBinContent(k);
				acc  = acceptance_mc->GetBinContent(k);
				tof  = tof_correction->GetBinContent(k);
				inn  = inn_track_correction->GetBinContent(k);	
				corr = ext*inn*tof;
				//Calculate Errors
				acc_error   = acceptance_mc->GetBinError(k);
				tof_error   = tof_correction->GetBinError(k);
				inn_error   = inn_track_correction->GetBinError(k);
				ext_error   = ext_track_correction->GetBinError(k);
				ext_error   = 0;
				total_error = corr*acc*TMath::Sqrt((acc_error/acc)**2 + (tof_error/tof)**2 + (inn_error/inn)**2 + (ext_error/ext)**2 ); 
				//cout <<"Bartel: " << i  << " Total Error: "  << total_error <<" Acc. Error: " << acc_error << " ToF Error: " << tof_error;
				//cout << " Inn. Tracker Error: " << inn_error << " Ext. Tracker Error: " << ext_error << endl; 
				total_correction->SetBinContent(k,corr);
				total_correction->SetBinError(k,total_error);
				if((corr*acc) > 1e-5 && (corr*acc) < 1)
				{	
					acceptance_->SetBinContent(k,corr*acc);
					acceptance_time->SetBinContent(i,k,corr*acc);
				}
				else
				{
					acceptance_->SetBinContent(k,0);
					acceptance_time->SetBinContent(i,k,0);	
				}
				acceptance_->SetBinError(k,total_error);
				acceptance_time->SetBinError(i,k,total_error);
			}
			save_plots=1;
			// Plot Effective Acceptance
			TCanvas *canvas_acceptance = new TCanvas("canvas_acceptance","Acceptance",600,400);
			gPad->SetLogx();
			acceptance_->Draw("P");
			//HIST Settings X-axis
			acceptance_->GetXaxis()->SetTitle("Rigidity [GV]");
			acceptance_->GetXaxis()->SetRangeUser(1.99,60.1);
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
			//HIST Settings General
			acceptance_->SetMarkerStyle(24);
			acceptance_->SetMarkerColor(kRed);
			acceptance_->SetMarkerSize(0.6);	
			plot_title.Form("Effective Acceptance vs Rigidity [Time Bin: %d]",i);
			acceptance_->SetTitle(plot_title);
			acceptance_->SetTitleFont(43);
			acceptance_->SetTitleSize(10);
			png_address.Form("png/acceptance_time/acceptance_%d.png",i);
			if(save_plots) canvas_acceptance->SaveAs(png_address);
			if(canvas_acceptance && !display_plots) canvas_acceptance->Close();
			// Close current MC file
			MC_data->Close();
		}
	}
	TFile *file = new TFile("root/acceptance_time.root","RECREATE");
	acceptance_time->Write("effective_acceptance_time");
	file->Close();

}
TH1D* acceptance_MC(TFile* data,bool save_plots,int bartel_count=0)
{
	bool display_plots = 0;
	TH1D* gen = (TH1D*)data->Get("hexp");
	TH1D* gen_w = (TH1D*)gen->Clone();
	TH2F* selected_time = (TH2F*)data->Get("Counts_Trigger");
	TH1D* selected = selected_time->ProjectionY();
	TH2F* MigrationMatrix = (TH2F*)data->Get("MigrationMatrix");	
	TH1D* acceptance =(TH1D*)selected->Clone();
	// Comput Migration Matrix Ratio
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
	//Compute SEL/GEN ratio with corrections
	double Ngen,Nsel,corr,error;
	for(int i=1;i<=selected->GetNbinsX();i++)
	{
		Ngen = gen->GetBinContent(i+6); // Uses different binning
		Nsel = selected->GetBinContent(i);
		corr = MMratio->GetBinContent(i);
		if(Ngen == 0)
		{
			acceptance->SetBinContent(i,0); 
			error = 0;
		}		
		else
		{ 
			acceptance->SetBinContent(i,corr*47.78*Nsel/Ngen); 
			error = TMath::Sqrt(Nsel)/Ngen;
		}	
		acceptance->SetBinError(i,error);
	}
	TString plot_title,png_address;
	// Plot MC Acceptance
	TCanvas *canvas_mc_acceptance = new TCanvas("canvas_mc_acceptance","Acceptance",600,400);
	gPad->SetLogx();
	acceptance->Draw("P");
	//HIST Settings X-axis
	acceptance->GetXaxis()->SetTitle("Rigidity [GV]");
	acceptance->GetXaxis()->SetRangeUser(1.99,60.1);
	acceptance->GetXaxis()->SetTitleOffset(1.0);
	acceptance->GetXaxis()->SetTitleFont(43);
	acceptance->GetXaxis()->SetTitleSize(15);
	acceptance->GetXaxis()->SetLabelFont(43);
	acceptance->GetXaxis()->SetLabelSize(12);
	//HIST Settings Y-axis
	acceptance->GetYaxis()->SetTitle("Acceptance [m^{2}sr]");	
	//acceptance_mc->GetYaxis()->SetRangeUser(0.004,0.018);
	acceptance->GetYaxis()->SetTitleOffset(1.0);
	acceptance->GetYaxis()->SetTitleFont(43);
	acceptance->GetYaxis()->SetTitleSize(15);
	acceptance->GetYaxis()->SetLabelFont(43);
	acceptance->GetYaxis()->SetLabelSize(10);
	//HIST Settings General
	acceptance->SetMarkerStyle(24);
	acceptance->SetMarkerColor(kRed);
	acceptance->SetMarkerSize(0.6);	
	plot_title.Form("Monte Carlo Acceptance vs Rigidity [Time Bin: %d]",bartel_count);
	acceptance->SetTitle(plot_title);
	acceptance->SetTitleFont(43);
	acceptance->SetTitleSize(10);
	png_address.Form("png/mc_acceptance_time/mc_acceptance_%d.png",bartel_count);
	save_plots = 1;
	if(save_plots) canvas_mc_acceptance->SaveAs(png_address);
	canvas_mc_acceptance->Close();

	gROOT->Reset();
	return acceptance;
}


