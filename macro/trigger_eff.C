void trigger_eff(TFile* ISS_data,TFile* MC_data,bool save_plots=0,bool display_plots=0)
{

	TH1D *iss_eff = trigger_eff_calculate(ISS_data,0,1);
	TH1D *mc_eff = trigger_eff_calculate(MC_data,1,1); 	
	// Plot Eff.
	TCanvas *Trigger_Efficiency = new TCanvas("Trigger_Efficiency","Trigger Efficiency vs Rigidity",600,400);
	TPad *pad1 = new TPad("pad1", "pad1" ,0,0.3,1,1.0);
	pad1->SetBottomMargin(0);
	pad1->SetGridx();
	pad1->SetLogx();
	pad1->Draw();
	pad1->cd();
	mc_eff->SetMarkerStyle(20);
	mc_eff->SetMarkerColor(1);
	mc_eff->SetLineColor(1);	
	mc_eff->SetMarkerSize(0.7);
	mc_eff->Draw();
	iss_eff->SetMarkerStyle(20);
	iss_eff->SetMarkerColor(2);
	iss_eff->SetLineColor(2);	
	iss_eff->SetMarkerSize(0.7);
	iss_eff->Draw("SAME");
	//Add Legend
	TLegend*lgd = new TLegend(0.7, 0.7, 0.90, 0.90);
	lgd->AddEntry(mc_eff,"MC He122","lp");
	lgd->AddEntry(iss_eff,"ISS","lp");
	lgd->SetMargin(0.35);
	lgd->Draw();
	//Histogram Settings
/*
	iss_eff->SetTitle("Trigger Efficiency vs Rigidity");
	iss_eff->SetTitleSize(13);
	iss_eff->GetXaxis()->SetTitle("Rigidity [GV]");
	iss_eff->GetYaxis()->SetTitle("Trigger Efficiency");
	iss_eff->GetYaxis()->SetRangeUser(0,10);
	iss_eff->GetYaxis()->SetTitleSize(10);
	iss_eff->GetYaxis()->SetTitleFont(43);
	iss_eff->GetYaxis()->SetTitleOffset(1.55);
	iss_eff->GetYaxis()->SetLabelFont(43);
	iss_eff->GetYaxis()->SetLabelSize(10);
*/	
	mc_eff->SetTitle("Trigger Efficiency vs Rigidity");
	mc_eff->SetTitleSize(13);
	mc_eff->GetXaxis()->SetTitle("Rigidity [GV]");
	mc_eff->GetYaxis()->SetTitle("Trigger Efficiency");
	mc_eff->GetYaxis()->SetRangeUser(0.905,1.02);	
	mc_eff->GetYaxis()->SetTitleSize(10);
	mc_eff->GetYaxis()->SetTitleFont(43);
	mc_eff->GetYaxis()->SetTitleOffset(1.55);
	mc_eff->GetYaxis()->SetLabelFont(43);
	mc_eff->GetYaxis()->SetLabelSize(10);
	
	//iss_eff->GetYaxis()->SetLineWidth(2);
	//mc_eff->GetYaxis()->SetLineWidth(2);
	// Ratio Plot
	Trigger_Efficiency->cd();
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
	ratio->Draw();
	// Histogram Settings
	ratio->SetTitle("");
	ratio->GetXaxis()->SetTitle("Rigidity [GV]");
	ratio->GetYaxis()->SetTitle("ISS/MC");
	ratio->GetYaxis()->SetRangeUser(0.96,1.04);
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
	
	if(save_plots) Trigger_Efficiency->SaveAs("png/TriggerEfficiency.png");
	if(Trigger_Efficiency && !display_plots) Trigger_Efficiency->Close();

	gROOT->Reset();
	//return(ratio);
}

TH1D* trigger_eff_calculate(TFile* data,bool IsMC,bool save_plots)
{
	gROOT->Reset();
	
	TH2F* trigger_time = (TH2F*)data->Get("Counts_Trigger")->Clone("Counts_Trigger");
	TH2F* unbias_time = (TH2F*)data->Get("Counts_Trigger_Unbias")->Clone("Counts_Trigger_Unbias");
	if(IsMC)
	{
		TH2F* triggerW_time = (TH2F*)data->Get("Counts_Trigger_weighted")->Clone("Counts_Trigger_weighted");
		TH2F* unbiasW_time = (TH2F*)data->Get("Counts_Trigger_Unbias_weighted")->Clone("Counts_Trigger_Unbias_weighted");
	}
	else
	{
		TH2F* triggerW_time = (TH2F*)trigger_time->Clone(); //same with trigger if iss data
		TH2F* unbiasW_time = (TH2F*)unbias_time->Clone(); //same with trigger if iss data	
	}

	TH1D* trigger = trigger_time->ProjectionY();
	TH1D* unbias = unbias_time->ProjectionY();
	TH1D* triggerW = triggerW_time->ProjectionY();	
	TH1D* unbiasW = unbiasW_time->ProjectionY();
	TH1D* Trigger_Eff = (TH1D*)triggerW->Clone(); // init. trigger eff. histogram

	int scl=1;
	if(!IsMC) scl=100; // scale with 100 for real data
	unbiasW->Scale(scl);

	// Calcualte Effs.
	TH1D* hTot =(TH1D*)triggerW->Clone();
	hTot->Add(unbiasW); //obtain total
	Trigger_Eff->Divide(hTot);
	// Calculate errors
	for(int i=1;i<=triggerW->GetNbinsX();i++)
	{	
		double bcu = scl*unbias->GetBinContent(i);
		double bct = trigger->GetBinContent(i);
		double err = TMath::Sqrt(bcu)/(bct);
		Trigger_Eff->SetBinError(i,err);
	}

	bool draw = 0;	
	if(draw) 
	{
		// Plot Triggered events
		TCanvas *C_triggered = new TCanvas("C_triggered","Triggered Events",600,400);
		trigger->Draw();
		C_triggered->SetLogx();
		// Plot Unbiased events
		TCanvas *C_unbias = new TCanvas("C_unbias","Unbias Events",600,400);
		unbias->Draw();
		C_unbias->SetLogx();
		// Plot triggered weighted events
		TCanvas *C_triggeredW = new TCanvas("C_triggeredW","Triggered Events Weighted",600,400); // same with trigger if iss data
		triggerW->Draw();
		C_triggeredW->SetLogx();
		// Plot unbiased weighted events
		TCanvas *C_unbiasW = new TCanvas("C_unbiasW","Unbias Events Weighted",600,400); // same with unbias if iss data
		unbiasW->Draw();
		C_unbiasW->SetLogx();	
		if(save_plots)
		{
			if(IsMC)
			{
				C_triggered->SaveAs("png/TriggerEfficiency_MC_triggered_events.png");			
				C_unbias->SaveAs("png/TriggerEfficiency_MC_unbiased_events.png");	
				C_triggeredW->SaveAs("png/TriggerEfficiency_MC_triggered_events_weighted.png");	
				C_unbiasW->SaveAs("png/TriggerEfficiency_MC_unbiased_events_weighted.png");	
			}
			else
			{
				C_triggered->SaveAs("png/TriggerEfficiency_triggered_events.png");			
				C_unbias->SaveAs("png/TriggerEfficiency_unbiased_events.png");	
			}
			if(C_triggered) C_triggered->Close();
			if(C_unbias) C_unbias->Close();
			if(C_triggeredW) C_triggeredW->Close();
			if(C_unbiasW) C_unbiasW->Close();
		}
	}

	return Trigger_Eff;
}

TH2F* trigger_eff_time_calculate(TFile* data,bool IsMC,bool save_plots)
{
	gROOT->Reset();
	TH2F* trigger_time = (TH2F*)data->Get("Counts_Trigger")->Clone("Counts_Trigger");
	TH2F* unbias_time = (TH2F*)data->Get("Counts_Trigger_Unbias")->Clone("Counts_Trigger_Unbias");
	if(IsMC)
	{
		TH2F* triggerW_time = (TH2F*)data->Get("Counts_Trigger_weighted")->Clone("Counts_Trigger_weighted");
		TH2F* unbiasW_time = (TH2F*)data->Get("Counts_Trigger_Unbias_weighted")->Clone("Counts_Trigger_Unbias_weighted");
	}
	else
	{
		TH2F* triggerW_time = (TH2F*)trigger_time->Clone(); //same with trigger if iss data
		TH2F* unbiasW_time = (TH2F*)unbias_time->Clone(); //same with trigger if iss data	
	}
	TH2F* Trigger_Eff = (TH2F*)triggerW_time->Clone(); // init. trigger eff. histogram
	int scl=1;
	if(!IsMC) scl=100; // scale with 100 for real data
	unbiasW_time->Scale(scl);
	// Calcualte Effs.
	TH2F* hTot =(TH2F*)triggerW_time->Clone();
	hTot->Add(unbiasW_time); //obtain total
	Trigger_Eff->Divide(hTot);
	// Calculate errors
	for(int i=1;i<=triggerW_time->GetNbinsX();i++)
	{	for(int j=1;j<=triggerW_time->GetNbinsY();j++)
		{	
			double bcu = scl*unbias_time->GetBinContent(i,j);
			double bct = trigger_time->GetBinContent(i,j);
			double err = TMath::Sqrt(bcu)/(bct);
			Trigger_Eff->SetBinError(i,j,err);
		}
	}
	return Trigger_Eff;
}

