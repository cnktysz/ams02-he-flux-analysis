// Takes UTC X axis and converts it to human time
void tof_eff(){
TFile *MyFile = new TFile("root/result.root");
TH2F* TofEff_sample = (TH2F*)MyFile->Get("Counts_TOFEff_sample")->Clone("Counts_TOFEff_sample");
TH2F* TofEff_selection = (TH2F*)MyFile->Get("Counts_TOFEff")->Clone("Counts_TOFEff");

TCanvas *TofEff = new TCanvas("TofEff","Exposure Time vs Rigidity",600,400);

TH1D* sample = TofEff_sample->ProjectionY();
TH1D* tof = TofEff_selection->ProjectionY();
sample->Sumw2();
tof->Sumw2();
tof->Divide(sample);
tof->SetTitle("ToF Efficiency vs Rigidity");
tof->GetXaxis()->SetTitle("Rigidity [GV]");
tof->GetYaxis()->SetTitle("Efficiency");
TofEff->SetLogx();
tof->Draw();
TofEff->SaveAs("png/TofEff.png");
}
