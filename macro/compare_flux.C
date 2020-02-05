void compare_flux()
{
	bool save_plots = 1;
	bool display_plots = 0;

	Double_t bin_time[41] = {1.92,2.15,2.4,2.67,2.97,3.29,3.64,4.02,4.43,4.88,5.37,5.9,6.47,7.09,7.76,8.48,9.26,10.1,11,12,13,14.1,15.3,16.6,18,19.5,21.1,22.8,24.7,26.7,28.8,31.1,33.5,36.1,38.9,41.9,45.1,48.5,52.2,56.1,60.3};

	TFile *he_flux_data = new TFile("root/result_time_iter1.root");
	TH2F  *he_flux_raw      = (TH2F*) he_flux_data->Get("he_flux");	
	
	TFile *published_flux_data = new TFile("root/result_time_iter0.root");
	TH2F  *published_flux      = (TH2F*) published_flux_data->Get("he_flux");
	
	// Take published bin projections
	TH1D *published_flux_1 = (TH1D*)published_flux->ProjectionX("",7,7);
	published_flux_1->SetName("published_flux_1");
	TH1D *published_flux_2 = published_flux->ProjectionX("",9,9);
	published_flux_2->SetName("published_he_flux_2");
	TH1D *published_flux_3 = published_flux->ProjectionX("",12,12);
	published_flux_3->SetName("published_flux_3");
	TH1D *published_flux_4 = published_flux->ProjectionX("",17,17);
	published_flux_4->SetName("published_flux_4");
	TH1D *published_flux_5 = published_flux->ProjectionX("",24,24);
	published_flux_5->SetName("published_flux_5");
	TH1D *published_flux_6 = published_flux->ProjectionX("",33,33);
	published_flux_6->SetName("published_flux_6");
	TH1D *published_flux_7 = published_flux->ProjectionX("",42,42);
	published_flux_7->SetName("published_flux_7");
	TH1D *published_flux_8 = published_flux->ProjectionX("",46,46);
	published_flux_8->SetName("published_flux_8");

	// Plot Published He Flux 3D
	TCanvas *canvas_published_he_flux = new TCanvas("canvas_published_he_flux","He Flux",600*2,400*2);
	gPad->SetTheta(12); 
	gPad->SetPhi(240);
	gPad->Update();
	gPad->SetLogy();
	published_flux->SetTitle("Iter0 He Flux vs Time vs Rigidity");	
	published_flux->SetTitleFont(2);
	published_flux->SetTitleSize(0.2);	
	// Flux axis
	published_flux->GetZaxis()->SetLabelFont(2);
	published_flux->GetZaxis()->SetTitleFont(2);
	published_flux->GetZaxis()->SetTitle("Flux [m^{-2}sr^{-1}s^{-1}GV^{-1}]");
	published_flux->GetZaxis()->SetTitleOffset(1.4);
	// Rigidity Axis
	published_flux->GetYaxis()->SetRangeUser(1.99,10);
	published_flux->GetYaxis()->SetLabelFont(2);
	published_flux->GetYaxis()->SetTitle("Rigidity [GV]");
	published_flux->GetYaxis()->SetTitleFont(2);
	published_flux->GetYaxis()->SetTitleOffset(1.8);
	published_flux->GetYaxis()->SetMoreLogLabels();
	// TIME AXIS Settings
	published_flux->GetXaxis()->SetLabelFont(2);
	published_flux->GetXaxis()->SetLabelOffset(0.03);
	published_flux->GetXaxis()->SetTitle("");
	published_flux->GetXaxis()->SetTimeOffset(0,"UTC");
	published_flux->GetXaxis()->SetTimeDisplay(1);
	published_flux->GetXaxis()->SetNdivisions(8);	
	published_flux->GetXaxis()->SetTimeFormat("#splitline{%b}{%Y}");
	published_flux->Draw("SURF3");
	if(save_plots) canvas_published_he_flux->SaveAs("png/flux_time_comparison/HeFlux_iter0_3D_time.png");
	if(canvas_published_he_flux && !display_plots) canvas_published_he_flux->Close();

	// Plot He Flux 3D
	TCanvas *canvas_he_flux = new TCanvas("canvas_he_flux","He Flux",600*2,400*2);
	gPad->SetTheta(12); 
	gPad->SetPhi(240);
	gPad->Update();
	gPad->SetLogy();
	he_flux_raw->SetTitle("Iter1 He Flux vs Time vs Rigidity");	
	he_flux_raw->SetTitleFont(2);
	he_flux_raw->SetTitleSize(0.2);	
	// Flux axis
	he_flux_raw->GetZaxis()->SetLabelFont(2);
	he_flux_raw->GetZaxis()->SetTitleFont(2);
	he_flux_raw->GetZaxis()->SetTitle("Flux [m^{-2}sr^{-1}s^{-1}GV^{-1}]");
	he_flux_raw->GetZaxis()->SetTitleOffset(1.4);
	// Rigidity Axis
	he_flux_raw->GetYaxis()->SetRangeUser(1.99,10);
	he_flux_raw->GetYaxis()->SetLabelFont(2);
	he_flux_raw->GetYaxis()->SetTitle("Rigidity [GV]");
	he_flux_raw->GetYaxis()->SetTitleFont(2);
	he_flux_raw->GetYaxis()->SetTitleOffset(1.8);
	he_flux_raw->GetYaxis()->SetMoreLogLabels();
	// TIME AXIS Settings
	he_flux_raw->GetXaxis()->SetLabelFont(2);
	he_flux_raw->GetXaxis()->SetLabelOffset(0.03);
	he_flux_raw->GetXaxis()->SetTitle("");
	he_flux_raw->GetXaxis()->SetTimeOffset(0,"UTC");
	he_flux_raw->GetXaxis()->SetTimeDisplay(1);
	he_flux_raw->GetXaxis()->SetNdivisions(8);	
	he_flux_raw->GetXaxis()->SetTimeFormat("#splitline{%b}{%Y}");
	he_flux_raw->Draw("SURF3");
	if(save_plots) canvas_he_flux->SaveAs("png/flux_time_comparison/HeFlux3D_iter1_time.png");
	if(canvas_he_flux && !display_plots) canvas_he_flux->Close();


	// Take bin projections
	TH1D *he_flux_1 = he_flux_raw->ProjectionX("",1,1);
	he_flux_1->SetName("he_flux_1");
	TH1D *he_flux_2 = he_flux_raw->ProjectionX("",3,3);
	he_flux_2->SetName("he_flux_2");
	TH1D *he_flux_3 = he_flux_raw->ProjectionX("",6,6);
	he_flux_3->SetName("he_flux_3");
	TH1D *he_flux_4 = he_flux_raw->ProjectionX("",11,11);
	he_flux_4->SetName("he_flux_4");
	TH1D *he_flux_5 = he_flux_raw->ProjectionX("",18,18);
	he_flux_5->SetName("he_flux_5");
	TH1D *he_flux_6 = he_flux_raw->ProjectionX("",27,27);
	he_flux_6->SetName("he_flux_6");
	TH1D *he_flux_7 = he_flux_raw->ProjectionX("",36,36);
	he_flux_7->SetName("he_flux_7");
	TH1D *he_flux_8 = he_flux_raw->ProjectionX("",40,40);
	he_flux_8->SetName("he_flux_8");

	// Plot He Flux 1st bin
	TCanvas *canvas_flux_1 = new TCanvas("canvas_flux_1","He Flux",360*4,360*2);
	gStyle->SetOptStat(0);
	TPad *pad1 = new TPad("pad1", "pad1" ,0,0.3,1,1.0);
	pad1->SetBottomMargin(0);
	pad1->SetGridx();
	pad1->SetGridy();
	pad1->Draw();
	pad1->cd();	
	// Plot Settings
	he_flux_1->SetTitleFont(2);
	he_flux_1->SetTitle(Form("He Flux vs Time [%.2f - %.2f] GV",bin_time[0],bin_time[1]));	
	he_flux_1->SetTitleSize(0.06);
	// Y Axis Settings
	he_flux_1->GetYaxis()->SetRangeUser(42,100);
	he_flux_1->GetYaxis()->SetLabelFont(2);	
	he_flux_1->GetYaxis()->SetTitle("Flux [m^{-2}sr^{-1}s^{-1}GV^{-1}]");
	he_flux_1->GetYaxis()->SetTitleFont(2);
	he_flux_1->GetYaxis()->SetTitleSize(0.05);
	he_flux_1->GetYaxis()->SetTitleOffset(0.6);
	// X Axis Settings
	he_flux_1->GetXaxis()->SetLabelFont(2);
	he_flux_1->GetXaxis()->SetLabelOffset(0.03);
	he_flux_1->GetXaxis()->SetTitle("");
	he_flux_1->GetXaxis()->SetTimeOffset(0,"gmt");
	he_flux_1->GetXaxis()->SetTimeDisplay(1);
	he_flux_1->GetXaxis()->SetNdivisions(13);	
	he_flux_1->GetXaxis()->SetTimeFormat("#splitline{%Y}{%b}"); 	
	// Draw Settings 
	he_flux_1->SetMarkerStyle(20);
	he_flux_1->SetMarkerColor(kRed);
	he_flux_1->SetMarkerSize(0.6);
	he_flux_1->Draw("P");
	// Draw Settings	
	published_flux_1->SetMarkerStyle(20);
	published_flux_1->SetMarkerColor(kBlue);
	published_flux_1->SetMarkerSize(0.6);
	published_flux_1->Draw("SAMEP");
	// Set Legend
	TLegend*lgd = new TLegend(0.8, 0.20, 0.9, 0.30);
	lgd->AddEntry(he_flux_1,"METU","lp");
	lgd->AddEntry(published_flux_1,"AMS","lp");
	//lgd->SetMargin(0.35);
	lgd->Draw();
	canvas_flux_1->cd();
	TPad *pad2 = new TPad("pad2","pad2",0,0.05,1,0.3);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.2);
	pad2->SetGridx();
	pad2->Draw();
	pad2->cd();
	// Define ratio plot
	TH1D* ratio1 = (TH1D*)he_flux_1->Clone();
	for(int i=1;i<=ratio1->GetNbinsX();i++)
	{
		if(i==46 || i==47 || i==48)
		{
			ratio1->SetBinContent(i,0);
			ratio1->SetBinError(i,0);	
		} 
		else
		{
			double my1        = he_flux_1->GetBinContent(i);
			double ams1       = published_flux_1->GetBinContent(i);
			double div1  	  = my1/ams1;
	        	double my1_error  = he_flux_1->GetBinError(i);	
			double ams1_error = published_flux_1->GetBinError(i);
			double error1 = div1*TMath::Sqrt((my1_error/my1)**2 + (ams1_error/ams1)**2);	
			ratio1->SetBinContent(i,div1);
			ratio1->SetBinError(i,error1);
		}
	}
	// Y Axis Settings
	ratio1->SetTitle("");
	ratio1->GetYaxis()->SetTitle("#Phi_{METU} /#Phi_{AMS}");
	ratio1->GetYaxis()->SetRangeUser(0.7,1);
	ratio1->GetYaxis()->SetTitleFont(2);
	ratio1->GetYaxis()->SetTitleSize(0.15);
	ratio1->GetYaxis()->SetTitleOffset(0.2);
	ratio1->GetYaxis()->SetLabelFont(2);
	ratio1->GetYaxis()->SetLabelSize(0.1);
	// X Axis Settings	
	ratio1->GetXaxis()->SetLabelFont(2);
	ratio1->GetXaxis()->SetLabelSize(0.12);
	ratio1->GetXaxis()->SetLabelOffset(0.06);
	ratio1->GetXaxis()->SetTitle("");
	ratio1->GetXaxis()->SetTimeOffset(0,"gmt");
	ratio1->GetXaxis()->SetTimeDisplay(1);
	ratio1->GetXaxis()->SetNdivisions(13);	
	ratio1->GetXaxis()->SetTimeFormat("#splitline{%b}{%Y}"); 	
	// Draw Settings	
	ratio1->SetMarkerColor(kMagenta);
	ratio1->Draw("P");
	// Save settings
	if(save_plots) canvas_flux_1->SaveAs("png/flux_time_comparison/HeFlux_time_1.png");
	if(canvas_flux_1 && !display_plots) canvas_flux_1->Close();


	// Plot He Flux 2nd bin
	TCanvas *canvas_flux_2 = new TCanvas("canvas_flux_2","He Flux",360*4,360*2);
	gStyle->SetOptStat(0);
	TPad *pad1 = new TPad("pad1", "pad1" ,0,0.3,1,1.0);
	pad1->SetBottomMargin(0);
	pad1->SetGridx();
	pad1->SetGridy();
	pad1->Draw();
	pad1->cd();	
	// Plot Settings
	he_flux_2->SetTitleFont(2);
	he_flux_2->SetTitle(Form("He Flux vs Time [%.2f - %.2f] GV",bin_time[2],bin_time[3]));	
	he_flux_2->SetTitleSize(0.06);
	// Y Axis Settings
	he_flux_2->GetYaxis()->SetRangeUser(35,85);	
	he_flux_2->GetYaxis()->SetLabelFont(2);	
	he_flux_2->GetYaxis()->SetTitle("Flux [m^{-2}sr^{-1}s^{-1}GV^{-1}]");
	he_flux_2->GetYaxis()->SetTitleFont(2);
	he_flux_2->GetYaxis()->SetTitleSize(0.05);
	he_flux_2->GetYaxis()->SetTitleOffset(0.6);
	// X Axis Settings
	he_flux_2->GetXaxis()->SetLabelFont(2);
	he_flux_2->GetXaxis()->SetLabelOffset(0.03);
	he_flux_2->GetXaxis()->SetTitle("");
	he_flux_2->GetXaxis()->SetTimeOffset(0,"gmt");
	he_flux_2->GetXaxis()->SetTimeDisplay(1);
	he_flux_2->GetXaxis()->SetNdivisions(13);	
	he_flux_2->GetXaxis()->SetTimeFormat("#splitline{%Y}{%b}"); 	
	// Draw Settings 
	he_flux_2->SetMarkerStyle(20);
	he_flux_2->SetMarkerColor(kRed);
	he_flux_2->SetMarkerSize(0.6);
	he_flux_2->Draw("P");
	// Draw Settings	
	published_flux_2->SetMarkerStyle(20);
	published_flux_2->SetMarkerColor(kBlue);
	published_flux_2->SetMarkerSize(0.6);
	published_flux_2->Draw("SAMEP");
	// Set Legend
	TLegend*lgd = new TLegend(0.8, 0.20, 0.9, 0.30);
	lgd->AddEntry(he_flux_2,"METU","lp");
	lgd->AddEntry(published_flux_2,"AMS","lp");
	//lgd->SetMargin(0.35);
	lgd->Draw();
	canvas_flux_2->cd();
	TPad *pad2 = new TPad("pad2","pad2",0,0.05,1,0.3);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.2);
	pad2->SetGridx();
	pad2->Draw();
	pad2->cd();
	// Define ratio plot
	TH1D* ratio2 = (TH1D*)he_flux_2->Clone();
	for(int i=1;i<=ratio2->GetNbinsX();i++)
	{ 
		if(i==46 || i==47 || i==48)
		{
			ratio2->SetBinContent(i,0);
			ratio2->SetBinError(i,0);	
		} 
		else
		{	
			double my2        = he_flux_2->GetBinContent(i);
			double ams2       = published_flux_2->GetBinContent(i);
			double div2  	  = my2/ams2;
	        	double my2_error  = he_flux_2->GetBinError(i);	
			double ams2_error = published_flux_2->GetBinError(i);
			double error2 = div2*TMath::Sqrt((my2_error/my2)**2 + (ams2_error/ams2)**2);	
			ratio2->SetBinContent(i,div2);
			ratio2->SetBinError(i,error2);
		}
	}

	// Y Axis Settings
	ratio2->SetTitle("");
	ratio2->GetYaxis()->SetTitle("#Phi_{METU} /#Phi_{AMS}");
	ratio2->GetYaxis()->SetRangeUser(0.7,1);
	ratio2->GetYaxis()->SetTitleFont(2);
	ratio2->GetYaxis()->SetTitleSize(0.15);
	ratio2->GetYaxis()->SetTitleOffset(0.2);
	ratio2->GetYaxis()->SetLabelFont(2);
	ratio2->GetYaxis()->SetLabelSize(0.1);
	// X Axis Settings	
	ratio2->GetXaxis()->SetLabelFont(2);
	ratio2->GetXaxis()->SetLabelSize(0.12);
	ratio2->GetXaxis()->SetLabelOffset(0.06);
	ratio2->GetXaxis()->SetTitle("");
	ratio2->GetXaxis()->SetTimeOffset(0,"gmt");
	ratio2->GetXaxis()->SetTimeDisplay(1);
	ratio2->GetXaxis()->SetNdivisions(13);	
	ratio2->GetXaxis()->SetTimeFormat("#splitline{%b}{%Y}"); 	
	// Draw Settings	
	ratio2->SetMarkerColor(kMagenta);
	ratio2->Draw("P");
	// Save settings
	if(save_plots) canvas_flux_2->SaveAs("png/flux_time_comparison/HeFlux_time_2.png");
	if(canvas_flux_2 && !display_plots) canvas_flux_2->Close();

	// Plot He Flux 3nd bin
	TCanvas *canvas_flux_3 = new TCanvas("canvas_flux_3","He Flux",360*4,360*2);
	gStyle->SetOptStat(0);
	TPad *pad1 = new TPad("pad1", "pad1" ,0,0.3,1,1.0);
	pad1->SetBottomMargin(0);
	pad1->SetGridx();
	pad1->SetGridy();
	pad1->Draw();
	pad1->cd();	
	// Plot Settings
	he_flux_3->SetTitleFont(2);
	he_flux_3->SetTitle(Form("He Flux vs Time [%.2f - %.2f] GV",bin_time[5],bin_time[6]));	
	he_flux_3->SetTitleSize(0.06);
	// Y Axis Settings
	he_flux_3->GetYaxis()->SetRangeUser(28,52);
	he_flux_3->GetYaxis()->SetLabelFont(2);	
	he_flux_3->GetYaxis()->SetTitle("Flux [m^{-2}sr^{-1}s^{-1}GV^{-1}]");
	he_flux_3->GetYaxis()->SetTitleFont(2);
	he_flux_3->GetYaxis()->SetTitleSize(0.05);
	he_flux_3->GetYaxis()->SetTitleOffset(0.6);
	// X Axis Settings
	he_flux_3->GetXaxis()->SetLabelFont(2);
	he_flux_3->GetXaxis()->SetLabelOffset(0.03);
	he_flux_3->GetXaxis()->SetTitle("");
	he_flux_3->GetXaxis()->SetTimeOffset(0,"gmt");
	he_flux_3->GetXaxis()->SetTimeDisplay(1);
	he_flux_3->GetXaxis()->SetNdivisions(13);	
	he_flux_3->GetXaxis()->SetTimeFormat("#splitline{%Y}{%b}"); 	
	// Draw Settings 
	he_flux_3->SetMarkerStyle(20);
	he_flux_3->SetMarkerColor(kRed);
	he_flux_3->SetMarkerSize(0.6);
	he_flux_3->Draw("P");
	// Draw Settings	
	published_flux_3->SetMarkerStyle(20);
	published_flux_3->SetMarkerColor(kBlue);
	published_flux_3->SetMarkerSize(0.6);
	published_flux_3->Draw("SAMEP");
	// Set Legend
	TLegend*lgd = new TLegend(0.8, 0.20, 0.9, 0.30);
	lgd->AddEntry(he_flux_3,"METU","lp");
	lgd->AddEntry(published_flux_3,"AMS","lp");
	//lgd->SetMargin(0.35);
	lgd->Draw();
	canvas_flux_3->cd();
	TPad *pad2 = new TPad("pad2","pad2",0,0.05,1,0.3);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.2);
	pad2->SetGridx();
	pad2->Draw();
	pad2->cd();
	// Define ratio plot
	TH1D* ratio3 = (TH1D*)he_flux_3->Clone();
	for(int i=1;i<=ratio3->GetNbinsX();i++)
	{ 
		if(i==46 || i==47 || i==48)
		{
			ratio3->SetBinContent(i,0);
			ratio3->SetBinError(i,0);	
		} 
		else
		{	
			double my3        = he_flux_3->GetBinContent(i);
			double ams3       = published_flux_3->GetBinContent(i);
			double div3  	  = my3/ams3;
	        	double my3_error  = he_flux_3->GetBinError(i);	
			double ams3_error = published_flux_3->GetBinError(i);
			double error3 = div3*TMath::Sqrt((my3_error/my3)**2 + (ams3_error/ams3)**2);	
			ratio3->SetBinContent(i,div3);
			ratio3->SetBinError(i,error3);
		}
	}
	// Y Axis Settings
	ratio3->SetTitle("");
	ratio3->GetYaxis()->SetTitle("#Phi_{METU} /#Phi_{AMS}");
	ratio3->GetYaxis()->SetRangeUser(0.7,1);
	ratio3->GetYaxis()->SetTitleFont(2);
	ratio3->GetYaxis()->SetTitleSize(0.15);
	ratio3->GetYaxis()->SetTitleOffset(0.2);
	ratio3->GetYaxis()->SetLabelFont(2);
	ratio3->GetYaxis()->SetLabelSize(0.1);
	// X Axis Settings	
	ratio3->GetXaxis()->SetLabelFont(2);
	ratio3->GetXaxis()->SetLabelSize(0.12);
	ratio3->GetXaxis()->SetLabelOffset(0.06);
	ratio3->GetXaxis()->SetTitle("");
	ratio3->GetXaxis()->SetTimeOffset(0,"gmt");
	ratio3->GetXaxis()->SetTimeDisplay(1);
	ratio3->GetXaxis()->SetNdivisions(13);	
	ratio3->GetXaxis()->SetTimeFormat("#splitline{%b}{%Y}"); 	
	// Draw Settings	
	ratio3->SetMarkerColor(kMagenta);
	ratio3->Draw("P");
	// Save settings
	if(save_plots) canvas_flux_3->SaveAs("png/flux_time_comparison/HeFlux_time_3.png");
	if(canvas_flux_3 && !display_plots) canvas_flux_3->Close();


	// Plot He Flux 4th bin
	TCanvas *canvas_flux_4 = new TCanvas("canvas_flux_4","He Flux",360*4,360*2);
	gStyle->SetOptStat(0);
	TPad *pad1 = new TPad("pad1", "pad1" ,0,0.3,1,1.0);
	pad1->SetBottomMargin(0);
	pad1->SetGridx();
	pad1->SetGridy();
	pad1->Draw();
	pad1->cd();	
	// Plot Settings
	he_flux_4->SetTitleFont(2);
	he_flux_4->SetTitle(Form("He Flux vs Time [%.2f - %.2f] GV",bin_time[10],bin_time[11]));	
	he_flux_4->SetTitleSize(0.06);
	// Y Axis Settings
	he_flux_4->GetYaxis()->SetRangeUser(12.5,19.5);
	he_flux_4->GetYaxis()->SetLabelFont(2);	
	he_flux_4->GetYaxis()->SetTitle("Flux [m^{-2}sr^{-1}s^{-1}GV^{-1}]");
	he_flux_4->GetYaxis()->SetTitleFont(2);
	he_flux_4->GetYaxis()->SetTitleSize(0.05);
	he_flux_4->GetYaxis()->SetTitleOffset(0.6);
	// X Axis Settings
	he_flux_4->GetXaxis()->SetLabelFont(2);
	he_flux_4->GetXaxis()->SetLabelOffset(0.03);
	he_flux_4->GetXaxis()->SetTitle("");
	he_flux_4->GetXaxis()->SetTimeOffset(0,"gmt");
	he_flux_4->GetXaxis()->SetTimeDisplay(1);
	he_flux_4->GetXaxis()->SetNdivisions(13);	
	he_flux_4->GetXaxis()->SetTimeFormat("#splitline{%Y}{%b}"); 	
	// Draw Settings 
	he_flux_4->SetMarkerStyle(20);
	he_flux_4->SetMarkerColor(kRed);
	he_flux_4->SetMarkerSize(0.6);
	he_flux_4->Draw("P");
	// Draw Settings	
	published_flux_4->SetMarkerStyle(20);
	published_flux_4->SetMarkerColor(kBlue);
	published_flux_4->SetMarkerSize(0.6);
	published_flux_4->Draw("SAMEP");
	// Set Legend
	TLegend*lgd = new TLegend(0.8, 0.20, 0.9, 0.30);
	lgd->AddEntry(he_flux_4,"METU","lp");
	lgd->AddEntry(published_flux_4,"AMS","lp");
	//lgd->SetMargin(0.35);
	lgd->Draw();
	canvas_flux_4->cd();
	TPad *pad2 = new TPad("pad2","pad2",0,0.05,1,0.3);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.2);
	pad2->SetGridx();
	pad2->Draw();
	pad2->cd();
	// Define ratio plot
	TH1D* ratio4 = (TH1D*)he_flux_4->Clone();
	for(int i=1;i<=ratio4->GetNbinsX();i++)
	{ 
		if(i==46 || i==47 || i==48)
		{
			ratio4->SetBinContent(i,0);
			ratio4->SetBinError(i,0);	
		} 
		else
		{	
			double my4        = he_flux_4->GetBinContent(i);
			double ams4       = published_flux_4->GetBinContent(i);
			double div4  	  = my4/ams4;
	        	double my4_error  = he_flux_4->GetBinError(i);	
			double ams4_error = published_flux_4->GetBinError(i);
			double error4 = div4*TMath::Sqrt((my4_error/my4)**2 + (ams4_error/ams4)**2);	
			ratio4->SetBinContent(i,div4);
			ratio4->SetBinError(i,error4);
		}
	}
	// Y Axis Settings
	ratio4->SetTitle("");
	ratio4->GetYaxis()->SetTitle("#Phi_{METU} /#Phi_{AMS}");
	ratio4->GetYaxis()->SetRangeUser(0.8,1.2);
	ratio4->GetYaxis()->SetTitleFont(2);
	ratio4->GetYaxis()->SetTitleSize(0.15);
	ratio4->GetYaxis()->SetTitleOffset(0.2);
	ratio4->GetYaxis()->SetLabelFont(2);
	ratio4->GetYaxis()->SetLabelSize(0.1);
	// X Axis Settings	
	ratio4->GetXaxis()->SetLabelFont(2);
	ratio4->GetXaxis()->SetLabelSize(0.12);
	ratio4->GetXaxis()->SetLabelOffset(0.06);
	ratio4->GetXaxis()->SetTitle("");
	ratio4->GetXaxis()->SetTimeOffset(0,"gmt");
	ratio4->GetXaxis()->SetTimeDisplay(1);
	ratio4->GetXaxis()->SetNdivisions(13);	
	ratio4->GetXaxis()->SetTimeFormat("#splitline{%b}{%Y}"); 	
	// Draw Settings	
	ratio4->SetMarkerColor(kMagenta);
	ratio4->Draw("P");
	// Save settings
	if(save_plots) canvas_flux_4->SaveAs("png/flux_time_comparison/HeFlux_time_4.png");
	if(canvas_flux_4 && !display_plots) canvas_flux_4->Close();

	// Plot He Flux 5th bin
	TCanvas *canvas_flux_5 = new TCanvas("canvas_flux_5","He Flux",360*4,360*2);
	gStyle->SetOptStat(0);
	TPad *pad1 = new TPad("pad1", "pad1" ,0,0.3,1,1.0);
	pad1->SetBottomMargin(0);
	pad1->SetGridx();
	pad1->SetGridy();
	pad1->Draw();
	pad1->cd();	
	// Plot Settings
	he_flux_5->SetTitleFont(2);
	he_flux_5->SetTitle(Form("He Flux vs Time [%.2f - %.2f] GV",bin_time[17],bin_time[18]));	
	he_flux_5->SetTitleSize(0.06);
	// Y Axis Settings
	he_flux_5->GetYaxis()->SetRangeUser(3.31,4.49);
	he_flux_5->GetYaxis()->SetLabelFont(2);	
	he_flux_5->GetYaxis()->SetTitle("Flux [m^{-2}sr^{-1}s^{-1}GV^{-1}]");
	he_flux_5->GetYaxis()->SetTitleFont(2);
	he_flux_5->GetYaxis()->SetTitleSize(0.05);
	he_flux_5->GetYaxis()->SetTitleOffset(0.6);
	// X Axis Settings
	he_flux_5->GetXaxis()->SetLabelFont(2);
	he_flux_5->GetXaxis()->SetLabelOffset(0.03);
	he_flux_5->GetXaxis()->SetTitle("");
	he_flux_5->GetXaxis()->SetTimeOffset(0,"gmt");
	he_flux_5->GetXaxis()->SetTimeDisplay(1);
	he_flux_5->GetXaxis()->SetNdivisions(13);	
	he_flux_5->GetXaxis()->SetTimeFormat("#splitline{%Y}{%b}"); 	
	// Draw Settings 
	he_flux_5->SetMarkerStyle(20);
	he_flux_5->SetMarkerColor(kRed);
	he_flux_5->SetMarkerSize(0.6);
	he_flux_5->Draw("P");
	// Draw Settings	
	published_flux_5->SetMarkerStyle(20);
	published_flux_5->SetMarkerColor(kBlue);
	published_flux_5->SetMarkerSize(0.6);
	published_flux_5->Draw("SAMEP");
	// Set Legend
	TLegend*lgd = new TLegend(0.8, 0.20, 0.9, 0.30);
	lgd->AddEntry(he_flux_5,"METU","lp");
	lgd->AddEntry(published_flux_5,"AMS","lp");
	//lgd->SetMargin(0.35);
	lgd->Draw();
	canvas_flux_5->cd();
	TPad *pad2 = new TPad("pad2","pad2",0,0.05,1,0.3);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.2);
	pad2->SetGridx();
	pad2->Draw();
	pad2->cd();
	// Define ratio plot
	TH1D* ratio5 = (TH1D*)he_flux_5->Clone();
	for(int i=1;i<=ratio5->GetNbinsX();i++)
	{ 
		if(i==46 || i==47 || i==48)
		{
			ratio5->SetBinContent(i,0);
			ratio5->SetBinError(i,0);	
		} 
		else
		{	
			double my5        = he_flux_5->GetBinContent(i);
			double ams5       = published_flux_5->GetBinContent(i);
			double div5  	  = my5/ams5;
	        	double my5_error  = he_flux_5->GetBinError(i);	
			double ams5_error = published_flux_5->GetBinError(i);
			double error5 = div5*TMath::Sqrt((my5_error/my5)**2 + (ams5_error/ams5)**2);	
			ratio5->SetBinContent(i,div5);
			ratio5->SetBinError(i,error5);
		}
	}	
	// Y Axis Settings
	ratio5->SetTitle("");
	ratio5->GetYaxis()->SetTitle("#Phi_{METU} /#Phi_{AMS}");
	ratio5->GetYaxis()->SetRangeUser(0.8,1.2);
	ratio5->GetYaxis()->SetTitleFont(2);
	ratio5->GetYaxis()->SetTitleSize(0.15);
	ratio5->GetYaxis()->SetTitleOffset(0.2);
	ratio5->GetYaxis()->SetLabelFont(2);
	ratio5->GetYaxis()->SetLabelSize(0.1);
	// X Axis Settings	
	ratio5->GetXaxis()->SetLabelFont(2);
	ratio5->GetXaxis()->SetLabelSize(0.12);
	ratio5->GetXaxis()->SetLabelOffset(0.06);
	ratio5->GetXaxis()->SetTitle("");
	ratio5->GetXaxis()->SetTimeOffset(0,"gmt");
	ratio5->GetXaxis()->SetTimeDisplay(1);
	ratio5->GetXaxis()->SetNdivisions(13);	
	ratio5->GetXaxis()->SetTimeFormat("#splitline{%b}{%Y}"); 	
	// Draw Settings	
	ratio5->SetMarkerColor(kMagenta);
	ratio5->Draw("P");
	// Save settings
	if(save_plots) canvas_flux_5->SaveAs("png/flux_time_comparison/HeFlux_time_5.png");
	if(canvas_flux_5 && !display_plots) canvas_flux_5->Close();

	// Plot He Flux 6th bin
	TCanvas *canvas_flux_6 = new TCanvas("canvas_flux_6","He Flux",360*4,360*2);
	gStyle->SetOptStat(0);
	TPad *pad1 = new TPad("pad1", "pad1" ,0,0.3,1,1.0);
	pad1->SetBottomMargin(0);
	pad1->SetGridx();
	pad1->SetGridy();
	pad1->Draw();
	pad1->cd();	
	// Plot Settings
	he_flux_6->SetTitleFont(2);
	he_flux_6->SetTitle(Form("He Flux vs Time [%.2f - %.2f] GV",bin_time[26],bin_time[27]));	
	he_flux_6->SetTitleSize(0.06);
	// Y Axis Settings
	he_flux_6->GetYaxis()->SetRangeUser(0.53,0.69);
	he_flux_6->GetYaxis()->SetLabelFont(2);	
	he_flux_6->GetYaxis()->SetTitle("Flux [m^{-2}sr^{-1}s^{-1}GV^{-1}]");
	he_flux_6->GetYaxis()->SetTitleFont(2);
	he_flux_6->GetYaxis()->SetTitleSize(0.05);
	he_flux_6->GetYaxis()->SetTitleOffset(0.6);
	// X Axis Settings
	he_flux_6->GetXaxis()->SetLabelFont(2);
	he_flux_6->GetXaxis()->SetLabelOffset(0.03);
	he_flux_6->GetXaxis()->SetTitle("");
	he_flux_6->GetXaxis()->SetTimeOffset(0,"gmt");
	he_flux_6->GetXaxis()->SetTimeDisplay(1);
	he_flux_6->GetXaxis()->SetNdivisions(13);	
	he_flux_6->GetXaxis()->SetTimeFormat("#splitline{%Y}{%b}"); 	
	// Draw Settings 
	he_flux_6->SetMarkerStyle(20);
	he_flux_6->SetMarkerColor(kRed);
	he_flux_6->SetMarkerSize(0.6);
	he_flux_6->Draw("P");
	// Draw Settings	
	published_flux_6->SetMarkerStyle(20);
	published_flux_6->SetMarkerColor(kBlue);
	published_flux_6->SetMarkerSize(0.6);
	published_flux_6->Draw("SAMEP");
	// Set Legend
	TLegend*lgd = new TLegend(0.8, 0.20, 0.9, 0.30);
	lgd->AddEntry(he_flux_6,"METU","lp");
	lgd->AddEntry(published_flux_6,"AMS","lp");
	//lgd->SetMargin(0.35);
	lgd->Draw();
	canvas_flux_6->cd();
	TPad *pad2 = new TPad("pad2","pad2",0,0.05,1,0.3);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.2);
	pad2->SetGridx();
	pad2->Draw();
	pad2->cd();
	// Define ratio plot
	TH1D* ratio6 = (TH1D*)he_flux_6->Clone();
	for(int i=1;i<=ratio6->GetNbinsX();i++)
	{ 
		if(i==46 || i==47 || i==48)
		{
			ratio6->SetBinContent(i,0);
			ratio6->SetBinError(i,0);	
		} 
		else
		{	
			double my6        = he_flux_6->GetBinContent(i);
			double ams6       = published_flux_6->GetBinContent(i);
			double div6  	  = my6/ams6;
	        	double my6_error  = he_flux_6->GetBinError(i);	
			double ams6_error = published_flux_6->GetBinError(i);
			double error6     = div6*TMath::Sqrt((my6_error/my6)**2 + (ams6_error/ams6)**2);	
			ratio6->SetBinContent(i,div6);
			ratio6->SetBinError(i,error6);
		}
	}	
	// Y Axis Settings
	ratio6->SetTitle("");
	ratio6->GetYaxis()->SetTitle("#Phi_{METU} /#Phi_{AMS}");
	ratio6->GetYaxis()->SetRangeUser(0.8,1.4);
	ratio6->GetYaxis()->SetTitleFont(2);
	ratio6->GetYaxis()->SetTitleSize(0.15);
	ratio6->GetYaxis()->SetTitleOffset(0.2);
	ratio6->GetYaxis()->SetLabelFont(2);
	ratio6->GetYaxis()->SetLabelSize(0.1);
	// X Axis Settings	
	ratio6->GetXaxis()->SetLabelFont(2);
	ratio6->GetXaxis()->SetLabelSize(0.12);
	ratio6->GetXaxis()->SetLabelOffset(0.06);
	ratio6->GetXaxis()->SetTitle("");
	ratio6->GetXaxis()->SetTimeOffset(0,"gmt");
	ratio6->GetXaxis()->SetTimeDisplay(1);
	ratio6->GetXaxis()->SetNdivisions(13);	
	ratio6->GetXaxis()->SetTimeFormat("#splitline{%b}{%Y}"); 	
	// Draw Settings	
	ratio6->SetMarkerColor(kMagenta);
	ratio6->Draw("P");
	// Save settings
	if(save_plots) canvas_flux_6->SaveAs("png/flux_time_comparison/HeFlux_time_6.png");
	if(canvas_flux_6 && !display_plots) canvas_flux_6->Close();

	// Plot He Flux 7th bin
	TCanvas *canvas_flux_7 = new TCanvas("canvas_flux_7","He Flux",360*4,360*2);
	gStyle->SetOptStat(0);
	TPad *pad1 = new TPad("pad1", "pad1" ,0,0.3,1,1.0);
	pad1->SetBottomMargin(0);
	pad1->SetGridx();
	pad1->SetGridy();
	pad1->Draw();
	pad1->cd();	
	// Plot Settings
	he_flux_7->SetTitleFont(2);
	he_flux_7->SetTitle(Form("He Flux vs Time [%.2f - %.2f] GV",bin_time[35],bin_time[36]));	
	he_flux_7->SetTitleSize(0.06);
	// Y Axis Settings
	he_flux_7->GetYaxis()->SetRangeUser(0.08,0.11);
	he_flux_7->GetYaxis()->SetLabelFont(2);	
	he_flux_7->GetYaxis()->SetTitle("Flux [m^{-2}sr^{-1}s^{-1}GV^{-1}]");
	he_flux_7->GetYaxis()->SetTitleFont(2);
	he_flux_7->GetYaxis()->SetTitleSize(0.05);
	he_flux_7->GetYaxis()->SetTitleOffset(0.6);
	// X Axis Settings
	he_flux_7->GetXaxis()->SetLabelFont(2);
	he_flux_7->GetXaxis()->SetLabelOffset(0.03);
	he_flux_7->GetXaxis()->SetTitle("");
	he_flux_7->GetXaxis()->SetTimeOffset(0,"gmt");
	he_flux_7->GetXaxis()->SetTimeDisplay(1);
	he_flux_7->GetXaxis()->SetNdivisions(13);	
	he_flux_7->GetXaxis()->SetTimeFormat("#splitline{%Y}{%b}"); 	
	// Draw Settings 
	he_flux_7->SetMarkerStyle(20);
	he_flux_7->SetMarkerColor(kRed);
	he_flux_7->SetMarkerSize(0.6);
	he_flux_7->Draw("P");
	// Draw Settings	
	published_flux_7->SetMarkerStyle(20);
	published_flux_7->SetMarkerColor(kBlue);
	published_flux_7->SetMarkerSize(0.6);
	published_flux_7->Draw("SAMEP");
	// Set Legend
	TLegend*lgd = new TLegend(0.8, 0.20, 0.9, 0.30);
	lgd->AddEntry(he_flux_7,"METU","lp");
	lgd->AddEntry(published_flux_7,"AMS","lp");
	//lgd->SetMargin(0.35);
	lgd->Draw();
	canvas_flux_7->cd();
	TPad *pad2 = new TPad("pad2","pad2",0,0.05,1,0.3);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.2);
	pad2->SetGridx();
	pad2->Draw();
	pad2->cd();
	// Define ratio plot
	TH1D* ratio7 = (TH1D*)he_flux_7->Clone();
	for(int i=1;i<=ratio7->GetNbinsX();i++)
	{ 
		if(i==46 || i==47 || i==48)
		{
			ratio7->SetBinContent(i,0);
			ratio7->SetBinError(i,0);	
		} 
		else
		{	
			double my7        = he_flux_7->GetBinContent(i);
			double ams7       = published_flux_7->GetBinContent(i);
			double div7  	  = my7/ams7;
	        	double my7_error  = he_flux_7->GetBinError(i);	
			double ams7_error = published_flux_7->GetBinError(i);
			double error7     = div7*TMath::Sqrt((my7_error/my7)**2 + (ams7_error/ams7)**2);	
			ratio7->SetBinContent(i,div7);
			ratio7->SetBinError(i,error7);
		}
	}		
	// Y Axis Settings
	ratio7->SetTitle("");
	ratio7->GetYaxis()->SetTitle("#Phi_{METU} /#Phi_{AMS}");
	ratio7->GetYaxis()->SetRangeUser(0.8,1.2);
	ratio7->GetYaxis()->SetTitleFont(2);
	ratio7->GetYaxis()->SetTitleSize(0.15);
	ratio7->GetYaxis()->SetTitleOffset(0.2);
	ratio7->GetYaxis()->SetLabelFont(2);
	ratio7->GetYaxis()->SetLabelSize(0.1);
	// X Axis Settings	
	ratio7->GetXaxis()->SetLabelFont(2);
	ratio7->GetXaxis()->SetLabelSize(0.12);
	ratio7->GetXaxis()->SetLabelOffset(0.06);
	ratio7->GetXaxis()->SetTitle("");
	ratio7->GetXaxis()->SetTimeOffset(0,"gmt");
	ratio7->GetXaxis()->SetTimeDisplay(1);
	ratio7->GetXaxis()->SetNdivisions(13);	
	ratio7->GetXaxis()->SetTimeFormat("#splitline{%b}{%Y}"); 	
	// Draw Settings	
	ratio7->SetMarkerColor(kMagenta);
	ratio7->Draw("P");
	// Save settings
	if(save_plots) canvas_flux_7->SaveAs("png/flux_time_comparison/HeFlux_time_7.png");
	if(canvas_flux_7 && !display_plots) canvas_flux_7->Close();

	// Plot He Flux 8th bin
	TCanvas *canvas_flux_8 = new TCanvas("canvas_flux_8","He Flux",360*4,360*2);
	gStyle->SetOptStat(0);
	TPad *pad1 = new TPad("pad1", "pad1" ,0,0.3,1,1.0);
	pad1->SetBottomMargin(0);
	pad1->SetGridx();
	pad1->SetGridy();
	pad1->Draw();
	pad1->cd();	
	// Plot Settings
	he_flux_8->SetTitleFont(2);
	he_flux_8->SetTitle(Form("He Flux vs Time [%.2f - %.2f] GV",bin_time[39],bin_time[40]));	
	he_flux_8->SetTitleSize(0.06);
	// Y Axis Settings
	he_flux_8->GetYaxis()->SetRangeUser(0.0485,0.035);	
	he_flux_8->GetYaxis()->SetLabelFont(2);	
	he_flux_8->GetYaxis()->SetTitle("Flux [m^{-2}sr^{-1}s^{-1}GV^{-1}]");
	he_flux_8->GetYaxis()->SetTitleFont(2);
	he_flux_8->GetYaxis()->SetTitleSize(0.05);
	he_flux_8->GetYaxis()->SetTitleOffset(0.6);
	// X Axis Settings
	he_flux_8->GetXaxis()->SetLabelFont(2);
	he_flux_8->GetXaxis()->SetLabelOffset(0.03);
	he_flux_8->GetXaxis()->SetTitle("");
	he_flux_8->GetXaxis()->SetTimeOffset(0,"gmt");
	he_flux_8->GetXaxis()->SetTimeDisplay(1);
	he_flux_8->GetXaxis()->SetNdivisions(13);	
	he_flux_8->GetXaxis()->SetTimeFormat("#splitline{%Y}{%b}"); 	
	// Draw Settings 
	he_flux_8->SetMarkerStyle(20);
	he_flux_8->SetMarkerColor(kRed);
	he_flux_8->SetMarkerSize(0.6);
	he_flux_8->Draw("P");
	// Draw Settings	
	published_flux_8->SetMarkerStyle(20);
	published_flux_8->SetMarkerColor(kBlue);
	published_flux_8->SetMarkerSize(0.6);
	published_flux_8->Draw("SAMEP");
	// Set Legend
	TLegend*lgd = new TLegend(0.8, 0.20, 0.9, 0.30);
	lgd->AddEntry(he_flux_8,"METU","lp");
	lgd->AddEntry(published_flux_8,"AMS","lp");
	//lgd->SetMargin(0.35);
	lgd->Draw();
	canvas_flux_8->cd();
	TPad *pad2 = new TPad("pad2","pad2",0,0.05,1,0.3);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.2);
	pad2->SetGridx();
	pad2->Draw();
	pad2->cd();
	// Define ratio plot
	TH1D* ratio8 = (TH1D*)he_flux_8->Clone();
	for(int i=1;i<=ratio8->GetNbinsX();i++)
	{ 
		if(i==46 || i==47 || i==48)
		{
			ratio8->SetBinContent(i,0);
			ratio8->SetBinError(i,0);	
		} 
		else
		{	
			double my8        = he_flux_8->GetBinContent(i);
			double ams8       = published_flux_8->GetBinContent(i);
			double div8  	  = my8/ams8;
	        	double my8_error  = he_flux_8->GetBinError(i);	
			double ams8_error = published_flux_8->GetBinError(i);
			double error8     = div8*TMath::Sqrt((my8_error/my8)**2 + (ams8_error/ams8)**2);	
			ratio8->SetBinContent(i,div8);
			ratio8->SetBinError(i,error8);
		}
	}
	// Y Axis Settings
	ratio8->SetTitle("");
	ratio8->GetYaxis()->SetTitle("#Phi_{METU} /#Phi_{AMS}");
	ratio8->GetYaxis()->SetRangeUser(0.8,1.2);
	ratio8->GetYaxis()->SetTitleFont(2);
	ratio8->GetYaxis()->SetTitleSize(0.15);
	ratio8->GetYaxis()->SetTitleOffset(0.2);
	ratio8->GetYaxis()->SetLabelFont(2);
	ratio8->GetYaxis()->SetLabelSize(0.1);
	// X Axis Settings	
	ratio8->GetXaxis()->SetLabelFont(2);
	ratio8->GetXaxis()->SetLabelSize(0.12);
	ratio8->GetXaxis()->SetLabelOffset(0.06);
	ratio8->GetXaxis()->SetTitle("");
	ratio8->GetXaxis()->SetTimeOffset(0,"gmt");
	ratio8->GetXaxis()->SetTimeDisplay(1);
	ratio8->GetXaxis()->SetNdivisions(13);	
	ratio8->GetXaxis()->SetTimeFormat("#splitline{%b}{%Y}"); 	
	// Draw Settings	
	ratio8->SetMarkerColor(kMagenta);
	ratio8->Draw("P");
	// Save settings
	if(save_plots) canvas_flux_8->SaveAs("png/flux_time_comparison/HeFlux_time_8.png");
	if(canvas_flux_8 && !display_plots) canvas_flux_8->Close();
}	
