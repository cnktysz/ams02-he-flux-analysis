void inn_track_eff(TFile* MyFile){
	TH2F* Tracker_sample = (TH2F*)MyFile->Get("Counts_Tracker_sample")->Clone("Counts_ExtTrack_sample");
	TH2F* InnTracker_selection = (TH2F*)MyFile->Get("Counts_InnTrack")->Clone("Counts_InnTrack");
	TH2F* Tracker_sample_1 = (TH2F*)MyFile->Get("Counts_Tracker_sample_1")->Clone("Counts_ExtTrack_sample_1");
	TH2F* InnTracker_selection_1 = (TH2F*)MyFile->Get("Counts_InnTrack_1")->Clone("Counts_InnTrack_1");
	TH2F* Tracker_sample_2 = (TH2F*)MyFile->Get("Counts_Tracker_sample_2")->Clone("Counts_ExtTrack_sample_2");
	TH2F* InnTracker_selection_2 = (TH2F*)MyFile->Get("Counts_InnTrack_2")->Clone("Counts_InnTrack_2");
	TH2F* Tracker_sample_3 = (TH2F*)MyFile->Get("Counts_Tracker_sample_3")->Clone("Counts_ExtTrack_sample_3");
	TH2F* InnTracker_selection_3 = (TH2F*)MyFile->Get("Counts_InnTrack_3")->Clone("Counts_InnTrack_3");

	TH1D* sample = Tracker_sample->ProjectionY();
	TH1D* inn = InnTracker_selection->ProjectionY();	
//	TH1D* sample1 = Tracker_sample_1->ProjectionY();
//	TH1D* inn1 = InnTracker_selection_1->ProjectionY();
	TH1D* sample2 = Tracker_sample_2->ProjectionY();
	TH1D* inn2 = InnTracker_selection_2->ProjectionY();	
	TH1D* sample3 = Tracker_sample_3->ProjectionY();
	TH1D* inn3 = InnTracker_selection_3->ProjectionY();	TCanvas *InnTracker = new TCanvas("InnTrackerEff","Inner Tracker Efficiency vs Rigidity",1366,768);

	InnTracker.Divide(2,2);
/*	
	InnTracker.cd(1);
	sample1->Sumw2();
	inn1->Sumw2();
	inn1->Divide(sample1);
	inn1->SetTitle("Inner Tracker Efficiency vs Rigidity (Beta)");
	inn1->GetXaxis()->SetTitle("Rigidity [GV]");
	inn1->GetYaxis()->SetTitle("Efficiency");
	gPad->SetLogx();
	inn1->GetXaxis()->SetRangeUser(0.8,3000);	
	inn1->Draw();
	//InnTracker->SaveAs("png/InnTrackEff.png");
*/
	InnTracker.cd(2);
	sample2->Sumw2();
	inn2->Sumw2();
	inn2->Divide(sample2);
	inn2->SetTitle("Inner Tracker Efficiency vs Rigidity (IGRF)");
	inn2->GetXaxis()->SetTitle("Rigidity [GV]");
	inn2->GetYaxis()->SetTitle("Efficiency");
	gPad->SetLogx();
	inn2->GetXaxis()->SetRangeUser(0.8,3000);	
	inn2->Draw();
	
	InnTracker.cd(3);
	sample3->Sumw2();
	inn3->Sumw2();
	inn3->Divide(sample3);
	inn3->SetTitle("Inner Tracker Efficiency vs Rigidity (ECAL)");
	inn3->GetXaxis()->SetTitle("Rigidity [GV]");
	inn3->GetYaxis()->SetTitle("Efficiency");
	gPad->SetLogx();
	inn3->GetXaxis()->SetRangeUser(0.8,3000);	
	inn3->Draw();
	
	InnTracker.cd(4);
	sample->Sumw2();
	inn->Sumw2();	
	inn->Divide(sample);
	inn->SetTitle("Inner Tracker Efficiency vs Rigidity");
	inn->GetXaxis()->SetTitle("Rigidity [GV]");
	inn->GetYaxis()->SetTitle("Efficiency");
	gPad->SetLogx();
	inn->GetXaxis()->SetRangeUser(0.8,3000);	
	inn->Draw();
	
/*	TCanvas *InnTracker1 = new TCanvas("InnTrackerEff","Inner Tracker Efficiency vs Rigidity",600,400);

	sample2->Add(sample1);
	inn2->Add(inn1);
	sample2->Sumw2();
	inn2->Sumw2();	
	inn2->Divide(sample2);
	inn2->SetTitle("Inner Tracker Efficiency vs Rigidity");
	inn2->GetXaxis()->SetTitle("Rigidity [GV]");
	inn2->GetYaxis()->SetTitle("Efficiency");
	gPad->SetLogx();
	inn2->GetXaxis()->SetRangeUser(0.8,3000);	
	inn2->Draw();
*/


}
