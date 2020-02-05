TH1D* migrationmatrix(bool save_plots)
{
	TH2F* MigrationMatrix = (TH2F*)file->Get("MigrationMatrix")->Clone("MigrationMatrix");	
	

	TH1D* MigrationMatrix_x = MigrationMatrix->ProjectionX();
	TH1D* MigrationMatrix_y = MigrationMatrix->ProjectionY();
	TH1D* MMratio =(TH1D*) MigrationMatrix_x->Clone();		
	MMratio->Divide(MigrationMatrix_y);

	// Normalize Migration Matrix
	double total_events = MigrationMatrix->GetEntries();
	for(int i=1;i<=MigrationMatrix->GetNbinsX();i++)
		for(int j=1;j<=MigrationMatrix->GetNbinsY();j++)
			MigrationMatrix->SetBinContent(i,j,MigrationMatrix->GetBinContent(i,j)/total_events);

	TCanvas *canvas_migration_matrix= new TCanvas("canvas_migration_matrix","Migration Matrix",600,400);
	MigrationMatrix->SetTitle("Migration Matrix");
	MigrationMatrix->GetXaxis()->SetTitle("Generated Rigidity [GV]");
	MigrationMatrix->GetYaxis()->SetTitle("Reconstructed Rigidity [GV]");	
	MigrationMatrix->GetXaxis()->SetRangeUser(1,1800);	
	gPad->SetLogx();
	MigrationMatrix->GetYaxis()->SetRangeUser(1,1800);	
	gPad->SetLogy();
	MigrationMatrix->Draw("COLZ");
	if(save_plots) canvas_migration_matrix->SaveAs("png/Acceptance/MigrationMatrix.png");
	canvas_migration_matrix->Close();
	
	TCanvas *canvas_migration_matrix_x= new TCanvas("canvas_migration_matrix_x","Migration Matrix X Proj.",600,400);
	MigrationMatrix_x->SetTitle("Migration Matrix X Proj.");
	MigrationMatrix_x->GetXaxis()->SetTitle("Generated Rigidity [GV]");
	MigrationMatrix_x->GetYaxis()->SetTitle("Events");	
	MigrationMatrix_x->GetXaxis()->SetRangeUser(1,1800);	
	gPad->SetLogx();
	//MigrationMatrix->GetYaxis()->SetRangeUser(1,1800);	
	gPad->SetLogy();
	MigrationMatrix_x->Draw("HIST");
	if(save_plots) canvas_migration_matrix_x->SaveAs("png/Acceptance/MigrationMatrix_x.png");
	canvas_migration_matrix_x->Close();

	TCanvas *canvas_migration_matrix_y= new TCanvas("canvas_migration_matrix_y","Migration Matrix Y Proj.",600,400);
	MigrationMatrix_y->SetTitle("Migration Matrix Y Proj.");
	MigrationMatrix_y->GetXaxis()->SetTitle("Reconstructed Rigidity [GV]");
	MigrationMatrix_y->GetYaxis()->SetTitle("events");	
	MigrationMatrix_y->GetXaxis()->SetRangeUser(1,1800);	
	gPad->SetLogx();
	//MigrationMatrix->GetYaxis()->SetRangeUser(1,1800);	
	gPad->SetLogy();
	MigrationMatrix_y->Draw("HIST");
	if(save_plots) canvas_migration_matrix_y->SaveAs("png/Acceptance/MigrationMatrix_y.png");
	canvas_migration_matrix_y->Close();

	TCanvas *canvas_migration_matrix_ratio= new TCanvas("canvas_migration_matrix_ratio","Migration Matrix Ratio",600,400);
	MMratio->SetTitle("Migration Matrix Ratio");
	MMratio->GetXaxis()->SetTitle("Rigidity [GV]");
	MMratio->GetYaxis()->SetTitle("Ratio");	
	MMratio->GetXaxis()->SetRangeUser(1,3000);	
	//MMratio>SetMarkerStyle(24);
	//MMratio->SetMarkerColor(kRed);
	//MMratio->SetMarkerSize(0.6);
	gPad->SetLogx();
	//MigrationMatrix->GetYaxis()->SetRangeUser(1,1800);	
	gPad->SetLogy();
	MMratio->SetMarkerStyle(24);
	MMratio->SetMarkerColor(kRed);
	MMratio->SetMarkerSize(0.6);
	MMratio->Draw("P");
	if(save_plots) canvas_migration_matrix_ratio->SaveAs("png/Acceptance/MigrationMatrix_ratio.png");
	canvas_migration_matrix_ratio->Close();

	return MMratio;
}
