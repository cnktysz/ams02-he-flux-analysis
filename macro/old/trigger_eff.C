// Takes UTC X axis and converts it to human time
void trigger_eff(){
	TFile *MyFile = new TFile("root/result.root");
	TH2F* biased = (TH2F*)MyFile->Get("Counts_Trigger")->Clone("Counts_Trigger");
	TH2F* unbiased = (TH2F*)MyFile->Get("Counts_Trigger_Unbias")->Clone("Counts_Trigger_Unbias");
	TCanvas *Trigger_Efficiency = new TCanvas("Trigger_Efficiency","Trigger Efficiency vs Rigidity",600,400);

	TH1D* Trigger_Eff = biased->ProjectionY();
	Trigger_Eff->Sumw2();
	TH1D *Trigger_Eff_unbias = unbiased->ProjectionY();
	Trigger_Eff_unbias->Sumw2();
	Trigger_Eff_unbias->Scale(100);
	Trigger_Eff_unbias->Add(Trigger_Eff);
	Trigger_Eff->Divide(Trigger_Eff_unbias);
	Trigger_Eff->SetTitle("Trigger Efficiency vs Rigidity");
	Trigger_Eff->GetXaxis()->SetTitle("Rigidity [GV]");
	Trigger_Eff->GetYaxis()->SetTitle("Trigger Efficiency");
	gPad->SetLogx();
	Trigger_Eff->Draw("E1");
	Trigger_Efficiency->SaveAs("png/TriggerEfficiency.png");
}
