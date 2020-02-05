TH1D *ext_track_eff(TFile* ISS_data, TFile* MC_data, bool save_plots,bool display_plots,int time_bin=0)
{
	//LAYER 1
	TH1D *iss_eff = tracker_eff_calculate(ISS_data,1,0,1,time_bin);
	TH1D *mc_eff = tracker_eff_calculate(MC_data,1,1,1); 		
	//TH1D *iss_eff = (TH1D*)mc_eff->Clone();
	// Plot Eff.
	TCanvas *Canvas_Efficiency = new TCanvas("Canvas_Efficiency","L1 Efficiency vs Rigidity",600,400);
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
	iss_eff->SetTitle("L1 Efficiency vs Rigidity");
	iss_eff->SetTitleSize(18);
	iss_eff->GetXaxis()->SetTitle("Rigidity [GV]");
	iss_eff->GetYaxis()->SetTitle("L1 Efficiency");
	iss_eff->GetYaxis()->SetRangeUser(0.57,0.80);
	iss_eff->GetXaxis()->SetRangeUser(1.92,60.3);
	iss_eff->GetYaxis()->SetTitleSize(10);
	iss_eff->GetYaxis()->SetTitleFont(43);
	iss_eff->GetYaxis()->SetTitleOffset(1.55);
	iss_eff->GetYaxis()->SetLabelFont(43);
	iss_eff->GetYaxis()->SetLabelSize(10);
	mc_eff->SetTitle("L1 Efficiency vs Rigidity");
	mc_eff->SetTitleSize(18);
	mc_eff->GetXaxis()->SetTitle("Rigidity [GV]");
	mc_eff->GetYaxis()->SetTitle("Efficiency");
	mc_eff->GetYaxis()->SetRangeUser(0.57,0.80);	
	mc_eff->GetXaxis()->SetRangeUser(1.92,60.3);
	mc_eff->GetYaxis()->SetTitleSize(10);
	mc_eff->GetYaxis()->SetTitleFont(43);
	mc_eff->GetYaxis()->SetTitleOffset(1.55);
	mc_eff->GetYaxis()->SetLabelFont(43);
	mc_eff->GetYaxis()->SetLabelSize(10);
	// Ratio Plot
	Canvas_Efficiency->cd();
	TPad *pad2 = new TPad("pad2","pad2",0,0.05,1,0.3);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.2);
	pad2->SetGridx();
	pad2->SetLogx();
	pad2->Draw();
	pad2->cd();
	// Define ratio plot (ISS/MC)
	TH1D* ratio = (TH1D*)iss_eff->Clone();
	ratio->Divide(mc_eff);
	for(int i=1;i<=ratio->GetNbinsX();i++)
	{
		double eff1       = iss_eff->GetBinContent(i);
		double eff2       = mc_eff->GetBinContent(i);
	        double eff1_error = iss_eff->GetBinError(i);	
		double eff2_error = mc_eff->GetBinError(i);
		double error = (eff1/eff2)*TMath::Sqrt((eff1_error/eff1)**2 + (eff2_error/eff2)**2);	
		if(ratio->GetBinCenter(i) < 2.97)
		{
			ratio->SetBinContent(i,1);
			ratio->SetBinError(i,0);
		}
		else ratio->SetBinError(i,error);
	}
	ratio->Draw();
	// Histogram Settings
	ratio->SetTitle("");
	ratio->GetXaxis()->SetTitle("Rigidity [GV]");
	ratio->GetYaxis()->SetTitle("ISS/MC");
	ratio->GetXaxis()->SetRangeUser(1.92,60.3);
	ratio->GetYaxis()->SetRangeUser(0.9,1.03);
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

	if(save_plots) Canvas_Efficiency->SaveAs(Form("png/ext_eff_time/L1Efficiency_%d.png",time_bin));
	if(Canvas_Efficiency && !display_plots) Canvas_Efficiency->Close();

	gROOT->Reset();
	TH1D* ext_track_ratio = (TH1D*) ratio->Clone();
	//LAYER 9
	TH1D *iss_eff = tracker_eff_calculate(ISS_data,9,0,1,time_bin);
	TH1D *mc_eff = tracker_eff_calculate(MC_data,9,1,1); 		
	// Plot Eff.
	TCanvas *Canvas_Efficiency = new TCanvas("Canvas_Efficiency","L9 Efficiency vs Rigidity",600,400);
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
	TLegend*lgd = new TLegend(0.7, 0.15, 0.90, 0.3);
	lgd->AddEntry(mc_eff,"MC He122","lp");
	lgd->AddEntry(iss_eff,"ISS","lp");
	lgd->SetMargin(0.35);
	lgd->Draw();
	//Histogram Settings	
	mc_eff->SetTitle("L9 Efficiency vs Rigidity");
	mc_eff->SetTitleFont(43);
	mc_eff->SetTitleSize(0.05);

	mc_eff->GetXaxis()->SetTitle("Rigidity [GV]");
	mc_eff->GetXaxis()->SetRangeUser(1.92,60.3);
	
	mc_eff->GetYaxis()->SetTitle("Efficiency");
	mc_eff->GetYaxis()->SetRangeUser(0.52,0.7);	
	mc_eff->GetYaxis()->SetTitleSize(15);
	mc_eff->GetYaxis()->SetTitleFont(43);
	mc_eff->GetYaxis()->SetTitleOffset(1.0);
	mc_eff->GetYaxis()->SetLabelFont(43);
	mc_eff->GetYaxis()->SetLabelSize(10);	
	// Ratio Plot
	Canvas_Efficiency->cd();
	TPad *pad2 = new TPad("pad2","pad2",0,0.05,1,0.3);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.2);
	pad2->SetGridx();
	pad2->SetLogx();
	pad2->Draw();
	pad2->cd();
	// Define ratio plot (ISS/MC)
	TH1D* ratio = (TH1D*)iss_eff->Clone();
	ratio->Divide(mc_eff);
	for(int i=1;i<=ratio->GetNbinsX();i++)
	{
		double eff1       = iss_eff->GetBinContent(i);
		double eff2       = mc_eff->GetBinContent(i);
	        double eff1_error = iss_eff->GetBinError(i);	
		double eff2_error = mc_eff->GetBinError(i);
		double error = (eff1/eff2)*TMath::Sqrt((eff1_error/eff1)**2 + (eff2_error/eff2)**2);	
		
		if(ratio->GetBinCenter(i) < 2.97)
		{
			ratio->SetBinContent(i,1);
			ratio->SetBinError(i,0);
		}
		else ratio->SetBinError(i,error);
	}
	ratio->Draw();
	// Histogram Settings
	ratio->SetTitle("");
	
	ratio->GetXaxis()->SetTitle("Rigidity [GV]");
	ratio->GetXaxis()->SetRangeUser(1.92,60.3);
	ratio->GetXaxis()->SetTitleSize(12);
	ratio->GetXaxis()->SetTitleFont(43);
	ratio->GetXaxis()->SetTitleOffset(4.);
	ratio->GetXaxis()->SetLabelFont(43);
	ratio->GetXaxis()->SetLabelSize(12);

	ratio->GetYaxis()->SetTitle("#epsilon_{ISS}/#epsilon_{MC}");
	ratio->GetYaxis()->SetRangeUser(0.92,1.0);
	ratio->GetYaxis()->SetTitleSize(15);
	ratio->GetYaxis()->SetTitleFont(43);
	ratio->GetYaxis()->SetTitleOffset(1.0);
	ratio->GetYaxis()->SetLabelFont(43);
	ratio->GetYaxis()->SetLabelSize(8);

	if(save_plots) Canvas_Efficiency->SaveAs(Form("png/ext_eff_time/L9Efficiency_%d.png",time_bin));
	if(Canvas_Efficiency && !display_plots) Canvas_Efficiency->Close();

	gROOT->Reset();	
	ext_track_ratio->Multiply(ratio);
	return ext_track_ratio;
}
TH1D *tracker_eff_calculate(TFile* data,int layer,bool IsMC,bool save_plots,int time_bin=0)
{
	if(layer==1 && IsMC)
	{
		TH2F* sample_time = (TH2F*)data->Get("Counts_L1_sample")->Clone("Counts_L1_sample");
		TH2F* sampleW_time = (TH2F*)data->Get("Counts_L1_sample_weighted")->Clone("Counts_L1_sample_weighted");
		TH2F* selection_time = (TH2F*)data->Get("Counts_L1")->Clone("Counts_L1");
		TH2F* selectionW_time = (TH2F*)data->Get("Counts_L1_weighted")->Clone("Counts_L1_weighted");
	//	cout << "Sample: " << sample_time->GetEntries() << " Selection: " << selection_time->GetEntries() << endl;
	
	}
	else if(layer==1 && !IsMC)
	{
		TH2F* sample_time = (TH2F*)data->Get("Counts_L1_sample")->Clone("Counts_L1_sample");
		TH2F* sampleW_time = (TH2F*)sample_time->Clone();
		TH2F* selection_time = (TH2F*)data->Get("Counts_L1")->Clone("Counts_L1");
		TH2F* selectionW_time = (TH2F*)selection_time->Clone();
	}
	else if(layer==9 && IsMC)
	{
		TH2F* sample_time = (TH2F*)data->Get("Counts_L9_sample")->Clone("Counts_L9_sample");
		TH2F* sampleW_time = (TH2F*)data->Get("Counts_L9_sample_weighted")->Clone("Counts_L9_sample_weighted");
		TH2F* selection_time = (TH2F*)data->Get("Counts_L9")->Clone("Counts_L9");
		TH2F* selectionW_time = (TH2F*)data->Get("Counts_L9_weighted")->Clone("Counts_L9_weighted");
	//	cout << "Sample: " << sample_time->GetEntries() << " Selection: " << selection_time->GetEntries() << endl;
	}
	else if(layer==9 && !IsMC)
	{
		TH2F* sample_time = (TH2F*)data->Get("Counts_L9_sample")->Clone("Counts_L9_sample");
		TH2F* sampleW_time = (TH2F*)sample_time->Clone();
		TH2F* selection_time = (TH2F*)data->Get("Counts_L9")->Clone("Counts_L9");
		TH2F* selectionW_time = (TH2F*)selection_time->Clone();
	}
	else return NULL;

	if(IsMC)
	{
		TH1D* sample = sample_time->ProjectionY();
		TH1D* selection = selection_time->ProjectionY();
		TH1D* sampleW = sampleW_time->ProjectionY();
		TH1D* selectionW = selectionW_time->ProjectionY();
	}
	else
	{
		TH1D* sample = sample_time->ProjectionY("",time_bin,time_bin);
		sample->SetName("sample");
		TH1D* selection = selection_time->ProjectionY("",time_bin,time_bin);
		selection->SetName("selection");
		TH1D* sampleW = sampleW_time->ProjectionY("",time_bin,time_bin);
		sampleW->SetName("sampleW");
		TH1D* selectionW = selectionW_time->ProjectionY("",time_bin,time_bin);
		selectionW->SetName("selectionW");
	}

	TH1D* efficiency = (TH1D*)selectionW->Clone();
	efficiency->Divide(sampleW);
	for(int i=1;i<=selectionW->GetNbinsX();i++)
	{	
		double Nsel       = selection->GetBinContent(i);
		double Nsam       = sample->GetBinContent(i);
	        double Nsel_error = selection->GetBinError(i);	
		double Nsam_error = sample->GetBinError(i);
		double error = (efficiency->GetBinContent(i))*TMath::Sqrt((Nsel_error/Nsel)**2 + (Nsam_error/Nsam)**2);	
		if(selection->GetBinCenter(i)<2.97)
		{
		efficiency->SetBinContent(i,1);
		efficiency->SetBinError(i,0);
		}
		else efficiency->SetBinError(i,error);
	}
	return efficiency;	
}
