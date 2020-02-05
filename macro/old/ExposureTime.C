// Takes UTC X axis and converts it to human time
void ExposureTime(const char* fname){
TFile *MyFile = new TFile(fname);
TH2F* ExposureIGRF = (TH2F*)MyFile->Get("ExposureIGRF")->Clone("ExposureIGRF");
// Exposure Time vs Rigidity vs Time
TCanvas *ExposureTime = new TCanvas("Exposure Time","Exposure Time vs Rigidity",600,400);
//gPad->SetTickx(1);
gPad->SetLogy();
//gPad->SetTicky(1);
ExposureIGRF->SetTitle("Exposure Time vs Rigidity");
ExposureIGRF->GetYaxis()->SetTitle("Rigidity [GV]");
ExposureIGRF->GetXaxis()->SetRangeUser(1,1800);
ExposureIGRF->GetXaxis()->SetTitle("");
ExposureIGRF->GetXaxis()->SetTimeOffset(0);
ExposureIGRF->SetStats(0);
ExposureIGRF->GetXaxis()->SetTimeDisplay(1);
ExposureIGRF->GetXaxis()->SetTimeFormat("#splitline{%Y}{%b,%d}");
ExposureIGRF->GetXaxis()->SetNdivisions(408,kFALSE);
TH1D * projh2Y = ExposureIGRF->ProjectionY();
ExposureIGRF->Draw("colz");
ExposureTime->SaveAs("png/ExposureTime.png");
// Exposure Time vs Rigidity
TCanvas *ExposureRigidity = new TCanvas("Exposure Time vs Rigidity","Exposure Time vs Rigidity",600,400);
TH1D* ExpTime = ExposureIGRF->ProjectionY();
gPad->SetLogx();
ExpTime->GetYaxis()->SetTitle("Time [s]");
ExpTime->Draw();
ExposureRigidity->SaveAs("png/ExposureTime_vs_Rigidity.png");
gROOT->Reset();
}
