// Takes UTC X axis and converts it to human time
void ext_track_eff(){
	TFile *MyFile = new TFile("root/result.root");
	TH2F* ExtTracker_sample = (TH2F*)MyFile->Get("Counts_ExtTrack_sample")->Clone("Counts_ExtTrack_sample");
	TH2F* ExtTracker_selection = (TH2F*)MyFile->Get("Counts_ExtTrack")->Clone("Counts_ExtTrack");

	TCanvas *ExtTracker = new TCanvas("ExtTrackerEff","External Tracker Efficiency vs Rigidity",600,400);

	TH1D* sample = ExtTracker_sample->ProjectionY();
	TH1D* ext = ExtTracker_selection->ProjectionY();
	sample->Sumw2();
	ext->Sumw2();
	ext->Divide(sample);
	ext->SetTitle("External Tracker Efficiency vs Rigidity");
	ext->GetXaxis()->SetTitle("Rigidity [GV]");
	ext->GetYaxis()->SetTitle("Efficiency");
	ExtTracker->SetLogx();
	ext->Draw();
	ExtTracker->SaveAs("png/ExtTrackEff.png");
}
