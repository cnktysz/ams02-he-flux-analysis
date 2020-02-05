TH1D* tof_eff(TFile* ISS_data, TFile* MC_data, bool save_plots,bool display_plots,int time_bin=0)
{
	TH1D *iss_eff = tof_eff_calculate(ISS_data,0,1,time_bin);
	TH1D *mc_eff = tof_eff_calculate(MC_data,1,1,0); 			
	// Plot Eff.
	TCanvas *TofEff = new TCanvas("TofEff","ToF Efficiency vs Rigidity",600,400);
	TPad *pad1 = new TPad("pad1", "pad1" ,0,0.3,1,1.0);
	pad1->SetBottomMargin(0);
	pad1->SetGridx();
	pad1->SetLogx();
	pad1->Draw();
	pad1->cd();
	
	mc_eff->SetMarkerStyle(15);
	mc_eff->SetMarkerColor(1);
	mc_eff->SetLineColor(1);	
	mc_eff->SetMarkerSize(0.6);
	mc_eff->Draw();
	
	iss_eff->SetMarkerStyle(20);
	iss_eff->SetMarkerColor(2);
	iss_eff->SetLineColor(2);	
	iss_eff->SetMarkerSize(0.6);
	iss_eff->Draw("SAME");
	//Add Legend
	TLegend*lgd = new TLegend(0.7, 0.15, 0.90, 0.35);
	lgd->AddEntry(mc_eff,"MC He122","lp");
	lgd->AddEntry(iss_eff,"ISS","lp");
	lgd->SetMargin(0.35);
	lgd->Draw();
	//Histogram Settings
	iss_eff->SetTitle("ToF Efficiency vs Rigidity");
	iss_eff->SetTitleSize(13);
	iss_eff->GetXaxis()->SetTitle("Rigidity [GV]");
	iss_eff->GetYaxis()->SetTitle("ToF Efficiency");
	iss_eff->GetYaxis()->SetRangeUser(0.955,1.01);
	iss_eff->GetXaxis()->SetRangeUser(1.92,60.3);
	iss_eff->GetYaxis()->SetTitleSize(10);
	iss_eff->GetYaxis()->SetTitleFont(43);
	iss_eff->GetYaxis()->SetTitleOffset(1.55);
	iss_eff->GetYaxis()->SetLabelFont(43);
	iss_eff->GetYaxis()->SetLabelSize(10);
	mc_eff->SetTitle("ToF Efficiency vs Rigidity");
	mc_eff->SetTitleSize(13);
	mc_eff->GetXaxis()->SetTitle("Rigidity [GV]");
	mc_eff->GetYaxis()->SetTitle("ToF Efficiency");
	mc_eff->GetXaxis()->SetRangeUser(1.92,60.3);
	mc_eff->GetYaxis()->SetRangeUser(0.955,1.01);	
	mc_eff->GetYaxis()->SetTitleSize(10);
	mc_eff->GetYaxis()->SetTitleFont(43);
	mc_eff->GetYaxis()->SetTitleOffset(1.55);
	mc_eff->GetYaxis()->SetLabelFont(43);
	mc_eff->GetYaxis()->SetLabelSize(10);
	
	//iss_eff->GetYaxis()->SetLineWidth(2);
	//mc_eff->GetYaxis()->SetLineWidth(2);
	// Ratio Plot
	TofEff->cd();
	TPad *pad2 = new TPad("pad2","pad2",0,0.05,1,0.3);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.2);
	pad2->SetGridx();
	pad2->SetLogx();
	pad2->Draw();
	pad2->cd();
	// Define ratio plot
	TH1D* ratio = (TH1D*)iss_eff->Clone();
	ratio->Divide(mc_eff);
	for(int i=1;i<=ratio->GetNbinsX();i++)
	{
		double eff1       = iss_eff->GetBinContent(i);
		double eff2       = mc_eff->GetBinContent(i);
	        double eff1_error = iss_eff->GetBinError(i);	
		double eff2_error = mc_eff->GetBinError(i);
		double error = (eff1/eff2)*TMath::Sqrt((eff1_error/eff1)**2 + (eff2_error/eff2)**2);	
		ratio->SetBinError(i,error);
	}
	ratio->Draw();
	// Histogram Settings
	ratio->SetTitle("");
	ratio->GetXaxis()->SetTitle("Rigidity [GV]");
	ratio->GetYaxis()->SetTitle("ISS/MC");
	ratio->GetYaxis()->SetRangeUser(0.98,1.05);
	ratio->GetXaxis()->SetRangeUser(1.92,60.3);
	ratio->GetYaxis()->SetTitleSize(10);
	ratio->GetYaxis()->SetTitleFont(43);
	ratio->GetYaxis()->SetTitleOffset(1.55);
	ratio->GetYaxis()->SetLabelFont(43);
	ratio->GetYaxis()->SetLabelSize(8);
	ratio->GetXaxis()->SetTitleSize(12);
	ratio->GetXaxis()->SetTitleFont(43);
	ratio->GetXaxis()->SetTitleOffset(4.);
	ratio->GetXaxis()->SetLabelFont(43);
	ratio->GetXaxis()->SetLabelSize(12);

	if(save_plots) TofEff->SaveAs(Form("png/tof_eff_time/TofEfficiency_%d.png",time_bin));
	if(TofEff && !display_plots) TofEff->Close();

	gROOT->Reset();
	return ratio;
}
TH1D *tof_eff_calculate(TFile* data,bool IsMC,bool save_plots,int time_bin)
{
	TH2F* sample_time = (TH2F*)data->Get("Counts_TOFEff_sample")->Clone("Counts_TOFEff_sample");
	TH2F* selection_time = (TH2F*)data->Get("Counts_TOFEff")->Clone("Counts_TOFEff");
	
	if(IsMC)
	{
		TH2F* sampleW_time = (TH2F*)data->Get("Counts_TOFEff_sample_weighted")->Clone("Counts_TOFEff_sample_weighted");
		TH2F* selectionW_time = (TH2F*)data->Get("Counts_TOFEff_weighted")->Clone("Counts_TOFEff_weighted");
		TH1D* sample = sample_time->ProjectionY();
		TH1D* selection = selection_time->ProjectionY();
		TH1D* sampleW = sampleW_time->ProjectionY();
		TH1D* selectionW = selectionW_time->ProjectionY();
	}
	else
	{
		TH2F* sampleW_time = (TH2F*)sample_time->Clone();
		TH2F* selectionW_time = (TH2F*)selection_time->Clone();
		TH1D* sample = sample_time->ProjectionY("",time_bin,time_bin);
		sample->SetName("sample");
		TH1D* selection = selection_time->ProjectionY("",time_bin,time_bin);
		selection->SetName("selection");
		TH1D* sampleW = sampleW_time->ProjectionY("",time_bin,time_bin);
		sampleW->SetName("sampleW");
		TH1D* selectionW = selectionW_time->ProjectionY("",time_bin,time_bin);
		selectionW->SetName("selectionW");
	}

	TH1D* ToF_Eff = (TH1D*)selectionW->Clone();
	ToF_Eff->Divide(sampleW);
	for(int i=1;i<=selectionW->GetNbinsX();i++)
	{	
		double Nsel       = selection->GetBinContent(i);
		double Nsam       = sample->GetBinContent(i);
	        double Nsel_error = selection->GetBinError(i);	
		double Nsam_error = sample->GetBinError(i);
		//double err = TMath::Sqrt(Nsel)/Nsam;
		double error = (ToF_Eff->GetBinContent(i))*TMath::Sqrt((Nsel_error/Nsel)**2 + (Nsam_error/Nsam)**2);	
		ToF_Eff->SetBinError(i,error);
		cout << "Bin: " << i << endl;
		cout <<"Nsel: " << Nsel <<" Nsel_error: " << Nsel_error << " Nsel_error/Nsel: " << Nsel_error/Nsel << endl; 
		cout <<"Nsel: " << Nsel <<" Nsel_error: " << Nsel_error << " Nsel_error/Nsel: " << Nsel_error/Nsel << endl; 
	}
	return ToF_Eff;	
}
