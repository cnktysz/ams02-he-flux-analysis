void RigidityConversion(TFile* MyFile)
{
	TH2F* Ecal_Tracker = (TH2F*)MyFile->Get("Ecal_Tracker")->Clone("Ecal_Tracker");
	int nbin = 75;
	double bin[76]={
		0.8,1,1.16,1.33,1.51,1.71,1.92,2.15,2.4,
		2.67,2.97,3.29,3.64,4.02,4.43,4.88,5.37,
		5.9,6.47,7.09,7.76,8.48,9.26,10.1,11,12,
		13,14.1,15.3,16.6,18,19.5,21.1,22.8,24.7,
		26.7,28.8,31.1,33.5,36.1,38.9,41.9,45.1,
		48.5,52.2,56.1,60.3,64.8,69.7,74.9,80.5,
		86.5,93,100,108,116,125,135,147,160,175,
		192,211,233,259,291,330,379,441,525,643,
		822,1130,1800,3000,8000};
	double coeff[75];	
	for(int n=0; n < nbin ; n++)
	{
	TH1D *slice = Ecal_Tracker->ProjectionX("slice", n,n+1);
	double mean = slice->GetMean();
	coeff[n] = mean/((bin[n]+bin[n+1])/2);
	cout << coeff[n] << ", ";	
	
	}	
	cout << endl;

	/*	
	TCanvas *Ecalslice = new TCanvas("Ecalslice","Ecal2Rig Coeff. vs Ecal Energy",600,400);	
	TH1D *slice_ = Ecal_Tracker->ProjectionX("slice", 23,24);
	//cout << slice_->GetRMS()*3 << endl;
	slice_->Draw();
	*/
	TCanvas *EcalRig = new TCanvas("EcalRig","Ecal2Rig Coeff. vs Ecal Energy",600,400);
	TGraph *Ecal2Rig = new TGraph(nbin,bin,coeff);
	Ecal2Rig->SetTitle("Ecal Energy vs Ecal2Rig Coefficient; Energy [GeV]; Tracker Rigidity / Ecal Energy [GV/GeV]");
	gPad->SetLogx();
	Ecal2Rig->Draw("A*");
	
	TCanvas *EcalTracker = new TCanvas("EcalTracker","Ecal Energy vs  Tracker Rigidity",600,400);
	Ecal_Tracker->SetTitle("Ecal Energy vs. Tracker Full Span Rigidity; Rigidity [GV];Energy [GeV]");	
	
	Ecal_Tracker->Draw();
	TH2F *conversion = new TH2F("conversion","Ecal Energy vs Tracker Rigidity; Energy [GeV]; Rigidity [GV]",nbin,bin,nbin,bin);
	conversion->SetFillColor(kRed);
	conversion->SetMarkerSize(20);
	for(int n=0;n < nbin ; n++) {
	conversion->Fill(((bin[n+1]+bin[n])/2)/coeff[n],((bin[n]+bin[n+1])/2));
	}
	gPad->SetLogx();
	gPad->SetLogy();
	conversion->SetMarkerStyle(32);
	conversion->SetMarkerSize(0.5);		
	conversion->SetMarkerColor(kRed);
	conversion->Draw("SAME");
	
	
}

