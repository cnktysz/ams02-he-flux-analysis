void mcw()
{
	Int_t nbin = 75;
	Double_t bin[76]={
		0.8,1,1.16,1.33,1.51,1.71,1.92,2.15,2.4,
		2.67,2.97,3.29,3.64,4.02,4.43,4.88,5.37,
		5.9,6.47,7.09,7.76,8.48,9.26,10.1,11,12,
		13,14.1,15.3,16.6,18,19.5,21.1,22.8,24.7,
		26.7,28.8,31.1,33.5,36.1,38.9,41.9,45.1,
		48.5,52.2,56.1,60.3,64.8,69.7,74.9,80.5,
		86.5,93,100,108,116,125,135,147,160,175,
		192,211,233,259,291,330,379,441,525,643,
		822,1130,1800,3000,8000};
	


}
Double_t DSTFill::GetMCW(Double_t rmc){
	if (!IsMC) return 1;
	static TF1 *frt = 0;  // Event raw rate
	static TF1 *fuf = 0;  // Raw/Unfold ratio
	if (!frt) {
		enum { N = 9 };
		// Copy node position from rtfit.C
		Double_t xn [N]   = { 1, 2, 5, 10, 20, 50, 100, 200, 500 };

		// Copy fitting output from rtfit.C
		Double_t par[N+2] = {0.958,  0.838,  0.190, -0.483, -1.254,
			-2.340, -3.172, -4.007, -5.099,  0.592, -2.666};

		frt = new TF1("spfun", SplFit::SpFunc, 1, 5000, N*2+2);
		for (Int_t i = 0; i < N;   i++) frt->SetParameter(i,   xn [i]);
		for (Int_t i = 0; i < N+2; i++) frt->SetParameter(i+N, par[i]);
	}
	if (!fuf) {
		enum { N = 12 };
		// Copy node position from rtcmp.C
		Double_t xn[N] = { 1, 1.5, 2, 5, 10, 20, 50,
			100, 200, 500, 1000, 2000 };

		// Taken after 10th iteration
		Double_t par[N+2] = {0};  

		// Read from file : put 0 0 0 0 ... for the first iteration

		TString sfn = "rtcmp.dat";
		ifstream fin(sfn);
		if (fin.good()) {
			cout << "Read from: " << sfn.Data() << endl;
			for (Int_t i = 0; i < N+2; i++) fin >> par[i];
		}


		fuf = new TF1("spfun", SplFit::SpFunc, 1, 5000, N*2+2);
		for (Int_t i = 0; i < N;   i++) fuf->SetParameter(i,   xn [i]);
		for (Int_t i = 0; i < N+2; i++) fuf->SetParameter(i+N, par[i]);
	}
	// Event rate
	SplFit::fLogX = 1; SplFit::fLogY = 1;
	SplFit::fBlxL = 0; SplFit::fBlxU = 1;
	SplFit::fN    = (frt->GetNpar()-2)/2;
	Double_t rt = frt->Eval(rmc);

	// Unfolding
	SplFit::fLogX = 1; SplFit::fLogY = 0;
	SplFit::fBlxL = 1; SplFit::fBlxU = 0;
	SplFit::fN    = (fuf->GetNpar()-2)/2;

	// For the first iteration
	Double_t uf = 4e-6;
	Int_t n = (fuf->GetNpar()-2)/2;
	if (fuf->GetParameter(n)   > 0 &&
			fuf->GetParameter(n+1) > 0) uf = 4e-6*fuf->Eval(rmc);

	// Unfolded event rate
	Double_t mcw = rt*uf*rmc;

	return mcw;
}

