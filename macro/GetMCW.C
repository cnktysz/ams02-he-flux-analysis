Double_t GetMCW(Double_t rmc){
	/*enum { N = 8 };
	// Copy node position from calculateMCW.C
	Double_t xn[N] = {3,5,10,20,50,100,300,1000};
	Double_t bb[2] = { 0, 0 };
	// Copy fitting output from calculateMCW.C
	Double_t par[N+2] = {-4.192, -4.458, -5.016, -5.764, -6.782, -7.674, -9.098, -11.058, -0.413, -2.322};	
	
	SplFit::fLogX = 1;
	SplFit::fLogY = 1;
	SplFit::fBlxU = 1;
	
	TF1 *frt = new TF1("spfun", SplFit::SpFunc, 2.97, 2000, N*2+2);
	for (Int_t i = 0; i < N;   i++) frt->SetParameter(i,   xn [i]);
	for (Int_t i = 0; i < N+2; i++) frt->SetParameter(i+N, par[i]);

	for (Int_t i = 0; i < N*2+2; i++) 
		cout << Form(" %6.3f,", frt->GetParameter(i));
	cout << endl;
*/
 	TFile *file = new TFile("root/MCW.root");
	TF1 *frt = (TF1*)file->Get("spfun");
	TCanvas* c1 = new TCanvas;
	frt->Draw("l");
	gPad->SetLogx();
	gPad->SetLogy();
	c1->SaveAs("png/MCW.png");	
	
	if(rmc<=2.97) return frt->Eval(2.97);
	else if(rmc>=2000) return frt->Eval(2000);
	else return frt->Eval(rmc);
}
