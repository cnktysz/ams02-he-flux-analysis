#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include <fstream>
using namespace std;

void calculateMCW()
{
	TFile *MC_data 	        = new TFile("root/he_flux_mc/he_flux_mc_iter1.root");
	TFile *Flux_data        = new TFile("root/result.root");
	TH1D  *flux 	        = (TH1D*) Flux_data->Get("he_flux");
	TH1D  *generated_events = (TH1D*) MC_data->Get("hexp");
	TH1D  *weights          = (TH1D*) flux->Clone();
	weights->Reset();
	Double_t norm_flux, norm_events;
	Double_t total_flux = 0;
	Double_t total_events = 0;
	for(int i=1;i<=flux->GetNbinsX();i++)
	{ 
        	total_flux   = total_flux + (double)flux->GetBinContent(i);
		total_events = total_events + (double)generated_events->GetBinContent(i);
	}
	for(int i=1;i<=flux->GetNbinsX();i++)
	{ 
        	norm_flux   = (double)flux->GetBinContent(i)/total_flux;
		norm_events = (double)generated_events->GetBinContent(i)/total_events;
		if(norm_events) weights->SetBinContent(i,norm_flux/norm_events);
	}
	
	TGraph *g_weights = SpFold::HtoG(weights);

	enum { N = 8 };
	Double_t xn[N] = {3, 5,10 ,20,50, 100,300, 1000 };
	Double_t bb[2] = { 0, 0 };
				
	SplFit::fLogX = 1; SplFit::fLogY = 1;
	SplFit::fBlxL = 0; SplFit::fBlxU = 1;
	SplFit::fN    = (N-2)/2;

	TF1 *func = SplFit::Fit(g_weights, N, xn, bb, "q0", 2.97, 2000);

	for (Int_t i = 0; i < N*2+2; i++) 
	{
		if(i==0) cout << "Nodes: "; 
		cout << Form(" %6.3f,", func->GetParameter(i));
		if(i==(N-1)) cout << endl << "Coefs: ";
	}
	cout << endl;
	
	frt = new TF1("spfun", SplFit::SpFunc, 2.97, 2000, N*2+2);
	for (Int_t i = 0; i < N;   i++) frt->SetParameter(i,   xn [i]);
	for (Int_t i = 0; i < N+2; i++) frt->SetParameter(i+N, func->GetParameter(i+N));

	TCanvas *c11 = new TCanvas("c11","Weights",600,400);
	gPad->SetLogx();
	gPad->SetLogy();
 	weights->GetYaxis()->SetRangeUser(1e-6,1e+1);	
	weights->GetXaxis()->SetRangeUser(2.97,2000);	
	weights->Draw("P");
	frt->Draw("same");
	frt->SetLineColor(kBlue);
	c11->SaveAs("png/weights.png");
	c11->Close();
	
	TFile *file = new TFile("root/MCW.root","RECREATE");
	func->Write();
	file->Close();

	ofstream myfile;
	myfile.open("root/weight_list.txt");
	for (Int_t i = 0; i < N+2; i++) myfile << func->GetParameter(i+N) << endl;
	myfile.close();
}	
