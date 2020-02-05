#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include <fstream>
using namespace std;

void spenvis()
{
	Double_t ams_bin[48]={1.92,2.15,2.4,
		2.67,2.97,3.29,3.64,4.02,4.43,4.88,5.37,
		5.9,6.47,7.09,7.76,8.48,9.26,10.1,11,12,
		13,14.1,15.3,16.6,18,19.5,21.1,22.8,24.7,
		26.7,28.8,31.1,33.5,36.1,38.9,41.9,45.1,
		48.5,52.2,56.1,60.3,64.8,69.7,74.9,80.5,
		86.5,93,100};	
   	Double_t ams_data[47]=
		{
		6.031e+1,5.657e+1,5.174e+1,4.694e+1,4.176e+1,3.650e+1,
		3.145e+1,2.671e+1,2.250e+1,1.876e+1,1.555e+1,1.282e+1,
		1.054e+1,8.646e+0,7.081e+0,5.786e+0,4.740e+0,3.866e+0,
		3.138e+0,2.566e+0,2.100e+0,1.711e+0,1.395e+0,1.134e+0,
		9.241e-1,7.485e-1,6.100e-1,4.963e-1,4.045e-1,3.295e-1,
		2.686e-1,2.194e-1,1.797e-1,1.469e-1,1.198e-1,9.824e-2,
		8.071e-2,6.566e-2,5.415e-2,4.446e-2,3.642e-2,2.984e-2,
		2.456e-2,2.015e-2,1.654e-2,1.357e-2,1.113e-2
		};
	Double_t spenvis_bin[36]={0.95, 1.05, 1.150, 1.250, 1.550, 1.650,1.950,2.0500,2.350,2.6500,2.9500,3.4500,3.5500,4.45000,4.5500,5.45000,5.5500,7.0500,7.150,8.85000,9.15000,1.0850000E+01,1.1150E+01,1.285000E+01,1.515000E+01,1.685000E+01,1.915000E+01,2.0850000E+01,3.9150000E+01,4.850000E+01,5.9150000E+01,6.850000E+01,7.9150000E+01,8.850000E+01,9.9150000E+01,1.085000E+02};   //35
	Double_t spenvis_data[35]={1.7391E+01, 1.6880E+01, 1.6320E+01, 1.5163E+01, 1.3922E+01, 1.2692E+01,1.1534E+01,1.0481E+01,9.0950E+00,7.9467E+00,6.7355E+00,6.0432E+00,5.1368E+00,4.3146E+00,3.6276E+00,3.0597E+00,2.3521E+00,1.8446E+00,1.4283E+00,1.0973E+00,8.6027E-01,6.8656E-01,5.5709E-01,3.8115E-01,2.7264E-01,2.0226E-01,1.5445E-01,5.3666E-02,2.5047E-02,1.3799E-02,8.4564E-03,5.5819E-03,3.8912E-03,2.8287E-03,2.1253E-03};  //35	


	TH1D *ams_flux = new TH1D("ams_flux","He Flux vs. Rigidity",47,ams_bin); 
	for(int i=0;i<47;i++) ams_flux->SetBinContent(i,ams_data[i]);

	TGraph *g_ams = SpFold::HtoG(ams_flux);

	enum { N = 8 };
	Double_t xn[N] = {2,5,10, 30,35, 40,50 ,100};
	Double_t bb[2] = { 0, 0 };
				
	SplFit::fLogX = 1; SplFit::fLogY = 1;
	SplFit::fBlxL = 0; SplFit::fBlxU = 1;
	SplFit::fN    = (N-2)/2;

	TF1 *func_ams = SplFit::Fit(g_ams, N, xn, bb, "q0", 1.92, 100.0);

	TH1D *spenvis_flux = new TH1D("spenvis_flux","He Flux vs. Rigidity",35,spenvis_bin); 
	for(int i=0;i<35;i++) spenvis_flux->SetBinContent(i,spenvis_data[i]);

	TGraph *g_spenvis = SpFold::HtoG(spenvis_flux);

	TF1 *func_spenvis = SplFit::Fit(g_spenvis, N, xn, bb, "q0", 1.92, 100);
	/*
	TCanvas *c11 = new TCanvas("c11","AMS02 - SPENVIS",600,400);
	gPad->SetLogx();
	gPad->SetLogy();
 	//weights->GetYaxis()->SetRangeUser(1e-6,1e+1);	
	ams_flux->GetXaxis()->SetRangeUser(1.92,100);	
	ams_flux->Draw("P");
	func_ams->Draw("same");
	func_ams->SetLineColor(kBlue);
	func_spenvis->Draw("same");
	func_spenvis->SetLineColor(kRed);
	
	c11->SaveAs("spenvis.pdf");
	c11->Close();
	*/
	TCanvas *canvas_HeFlux = new TCanvas("canvas_HeFlux","He Flux");
	TPad *pad1 = new TPad("pad1", "pad1" ,0,0.3,1,1.0);
	pad1->SetBottomMargin(0);
	pad1->SetGridx();
	pad1->SetLogx();
	pad1->SetGridy();
	pad1->SetLogy();
	pad1->Draw();
	pad1->cd();
	ams_flux->SetTitle("He Flux vs Rigidity");	
	ams_flux->GetXaxis()->SetRangeUser(1.92,100);
	ams_flux->GetXaxis()->SetTitleOffset(1.5);
	ams_flux->GetXaxis()->SetTitle("Rigidity [GV]");
	//ams_flux->GetYaxis()->SetRangeUser(3e+2,1e+4);
	ams_flux->GetYaxis()->SetTitleOffset(1.0);
	ams_flux->GetYaxis()->SetTitleSize(15);
	ams_flux->GetYaxis()->SetTitleFont(43);
	ams_flux->GetYaxis()->SetTitle("Flux [m^{-2}sr^{-1}s^{-1}GV^{-1}]");
	ams_flux->SetMarkerStyle(20);
	ams_flux->SetMarkerColor(kRed);
	ams_flux->SetMarkerSize(0.6);
	ams_flux->Draw("P");
	spenvis_flux->SetMarkerStyle(20);
	spenvis_flux->SetMarkerColor(kBlue);
	spenvis_flux->SetMarkerSize(0.6);
	spenvis_flux->Draw("SAMEP");
	//func_spenvis->Draw("same");
	//func_ams->Draw("same");
	TLegend*lgd = new TLegend(0.7, 0.7, 0.90, 0.9);
	lgd->AddEntry(ams_flux,"AMS-02","lp");
	lgd->AddEntry(spenvis_flux,"SPENVIS","lp");
	lgd->SetMargin(0.35);
	lgd->Draw();
	canvas_HeFlux->cd();
	TPad *pad2 = new TPad("pad2","pad2",0,0.05,1,0.3);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.2);
	pad2->SetGridx();
	pad2->SetLogx();
	pad2->Draw();
	pad2->cd();
	// Define ratio plot
	TH1D* ratio = (TH1D*)ams_flux->Clone();
	for(int i=0;i<ams_flux->GetNbinsX();i++)
	{
		double div = func_ams(ams_flux->GetBinCenter(i)) / func_spenvis(ams_flux->GetBinCenter(i));
		ratio->SetBinContent(i,div);
	}
	ratio->Draw("P");
	// Histogram Settings
	ratio->SetTitle("");
	ratio->GetYaxis()->SetTitle("#Phi_{AMS} /#Phi_{SPENVIS}");
	//ratio->GetYaxis()->SetRangeUser(0.85,1.05);
	ratio->GetXaxis()->SetRangeUser(1.92,100);
	ratio->GetYaxis()->SetTitleSize(15);
	ratio->GetYaxis()->SetTitleFont(43);
	ratio->GetYaxis()->SetTitleOffset(1.0);
	ratio->GetYaxis()->SetLabelFont(43);
	ratio->GetYaxis()->SetLabelSize(8);
	ratio->GetXaxis()->SetTitle("Rigidity [GV]");
	ratio->GetXaxis()->SetTitleSize(12);
	ratio->GetXaxis()->SetTitleFont(43);
	ratio->GetXaxis()->SetTitleOffset(4.);
	ratio->GetXaxis()->SetLabelFont(43);
	ratio->GetXaxis()->SetLabelSize(12);
	ratio->SetMarkerColor(kMagenta);
	canvas_HeFlux->SaveAs("spenvis.pdf");
	




}	
