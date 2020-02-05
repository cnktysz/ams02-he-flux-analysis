#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include <fstream>
#include "string"
using namespace std;

void calculateMCW_time()
{
	//TFile *Flux_data        = new TFile("/afs/cern.ch/work/c/ctuysuz/private/he_flux/root/I2_HeliumFlux_pBin_B1036_Hawaii_heV27_Extended.root");
	//TH2F  *flux 	        = (TH2F*) Flux_data->Get("hh_HeliumFlux_TotErr");
	TFile *Flux_data        = new TFile("/afs/cern.ch/work/c/ctuysuz/private/he_flux/root/result_time.root");
	TH2F  *flux		= (TH2F*) Flux_data->Get("he_flux");

	int nbin=40;
	Double_t bin[41]=
	{
		1.92,2.15,2.40,2.67,2.97,3.29,3.64,4.02,
		4.43,4.88,5.37,5.90,6.47,7.09,7.76,8.48,
		9.26,10.1,11.0,12.00,13,14.10,15.30,16.6,
		18.0,19.5,21.1,22.8,24.7,26.7,28.8,31.1,
		33.5,36.1,38.9,41.9,45.1,48.5,52.2,56.1,60.3
	};
	Double_t norm_flux, total_flux, norm_events, total_events;
	enum { N = 4 };
	Double_t xn[N] = {2,5,10,60};
	Double_t bb[2] = { 0, 0 };
	SplFit::fLogX = 1; SplFit::fLogY = 1;
	SplFit::fBlxL = 0; SplFit::fBlxU = 1;
	SplFit::fN    = (N-2)/2;
	TString iter;
	ofstream myfile;
	myfile.open("/afs/cern.ch/work/c/ctuysuz/private/he_flux/root/weight_list.txt");	
	cout << "Weights are being printed to root/weights_list.txt..." << endl;
	cout << "Weight plots are in png/weights/ " << endl;

	for(int i=1;i<=flux->GetNbinsX();i++)
	{
		cout << "Iteration: "<< i << endl;
		if(i==46 || i==47 || i==48)
		{ 
			cout << "Skipping iteration " << i << endl;
			continue;
		}
		TFile *MC_data 	        = new TFile(Form("/afs/cern.ch/work/c/ctuysuz/private/he_flux/root/he_flux_mc_time/he_flux_mc_time_%d.root",i));	
		//TFile *MC_data 	        = new TFile(Form("/afs/cern.ch/work/c/ctuysuz/private/he_flux/root/he_flux_mc_time/he_flux_mc_time_%d.root",i));	
		TH1D  *generated_events = (TH1D*) MC_data->Get("hexp");
		TH1D  *weight_		= new TH1D("weight_","weights",nbin,bin);
		TH2F  *weights          = (TH2F*) flux->Clone();
		TH1D  *cenk = (TH1D*)flux->ProjectionY();
		weights->Reset();

		total_flux = 0;
		total_events = 0;
		for(int j=1;j<=flux->GetNbinsX();j++)
		{ 
			total_flux   = total_flux + (double)flux->GetBinContent(i,j);
			total_events = total_events + (double)generated_events->GetBinContent(j+6);
		//	cout << "Bin no: " << j << " Bin Center: " << cenk->GetBinCenter(j+6) << " Flux: " << total_flux << " Events: " << total_events << endl;
		}
		for(int j=1;j<=flux->GetNbinsX();j++)
		{ 
			norm_flux   = (double)flux->GetBinContent(i,j)/total_flux;
			norm_events = (double)generated_events->GetBinContent(j+6)/total_events;
			if(norm_events) weights->SetBinContent(i,j,norm_flux/norm_events);
			if(norm_events) weight_->SetBinContent(j,norm_flux/norm_events);			
		}

		// FIT Weights
		TGraph *g_weights = SpFold::HtoG(weight_);
		TF1 *func = SplFit::Fit(g_weights, N, xn, bb, "q0", 1.92, 60.3);
		// Print results to txt file
		for (int k = 0; k < N+2; k++) myfile << func->GetParameter(k+N) << endl;
		// Plot and save weights
		TCanvas *c11 = new TCanvas("c11","Weights",600,400);
		gPad->SetLogx();
		gPad->SetLogy();
		weight_->GetXaxis()->SetRangeUser(1.92,60.3);	
		weight_->SetMarkerStyle(20);
		weight_->SetMarkerColor(kRed);
		weight_->SetMarkerSize(0.6);
		weight_->Draw("P");
		func->Draw("same");
		func->SetLineColor(kBlue);	
		c11->SaveAs(Form("/afs/cern.ch/work/c/ctuysuz/private/he_flux/png/weights/weights_%d.png",i));
		c11->Close();

	}
	myfile.close();
	cout << "Root files are beng printed to root/MCW_time.root ..." << endl;
	TFile *file = new TFile("/afs/cern.ch/work/c/ctuysuz/private/he_flux/root/MCW_time.root","RECREATE");
	weights->Write("MCW_time");
	file->Close();
	cout << "All weights are calculated.......Done!" << endl;
}	
