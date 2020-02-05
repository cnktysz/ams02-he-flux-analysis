TH1D *ExposureTime(TFile* MyFile, bool save_plots, bool display_plots){
//	TFile	*MyFile       = new TFile(address);
//	gFile 		      = MyFile;
	TH2F	*ExposureIGRF = (TH2F*)gFile->Get("ExposureIGRF")->Clone("ExposureIGRF");
	TCanvas *ExposureTime = new TCanvas("Exposure Time","Exposure Time vs Rigidity",600,400);
	gPad->SetLogy();
	ExposureIGRF->SetTitle("Exposure Time vs Rigidity");
	ExposureIGRF->GetYaxis()->SetTitle("Rigidity [GV]");
	ExposureIGRF->GetYaxis()->SetRangeUser(0.8,8000);
	ExposureIGRF->GetXaxis()->SetTitle("");
	ExposureIGRF->GetXaxis()->SetTimeOffset(0);
	ExposureIGRF->SetStats(0);
	ExposureIGRF->GetXaxis()->SetTimeDisplay(1);
	ExposureIGRF->GetXaxis()->SetTimeFormat("#splitline{%b}{%Y}");
	ExposureIGRF->Draw("colz");
	if(save_plots) ExposureTime->SaveAs("png/ExposureTime.png");
	if(ExposureTime && !display_plots) ExposureTime->Close();
	// Exposure Time vs Rigidity
	TH1D *reduced_exp = ExposureIGRF->ProjectionY();
	TCanvas *ExposureRigidity = new TCanvas("Exposure Time vs Rigidity","Exposure Time vs Rigidity",600,400);
	TH1D* ExpTime = (TH1D*)reduced_exp->Clone();//ExposureIGRF->ProjectionY("",0,81);
	//cout << ExposureIGRF->GetNbinsX() << endl;
	cout << "Max. Exposure Time: " << ExpTime->GetBinContent(60) << " s"  << endl;
	gPad->SetLogx();
	ExpTime->GetYaxis()->SetTitle("Time [s]");
	ExpTime->Draw();
	if(save_plots) ExposureRigidity->SaveAs("png/ExposureTime_vs_Rigidity.png");
	if(ExposureRigidity && !display_plots) ExposureRigidity->Close();
	gROOT->Reset();
	return ExpTime;
}


