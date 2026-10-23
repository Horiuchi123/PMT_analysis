#include "declare.h"
//#include "MetaData/MetaData.h"
#include "ExtractCharacteristics.h"
#include "TFitter.h"
#include "mQFitter.h"
#include "Color.C"
#include <tuple>
#include <cstdlib>
#include <algorithm>




class h_settings{
public:
	TString varexp;
	TString selection = "";
	TString option  = "";

	TString h_name  = "";
	TString x_title = "";
	TString y_title = "";
	TString h_title = "";

	Double_t bin_per_unit = 1;
	Double_t quantile = 0.01;
	Double_t content_thr = 0.01;

	Int_t h_min =  1;
	Int_t h_max = -1;//h_max < h_min (defalt)

};



class save_can{
	public:
		inline static TString output_dir = "./";
		inline static TString suffix = "";
		static void vas(TString canvas_name, TString root_file_name, TString pdf_name){
			// root file
			pdf_name = output_dir + pdf_name;
			root_file_name = output_dir + root_file_name;
			TFile* out_f = TFile::Open(root_file_name,"update");
			TCanvas* canvas = (TCanvas*) gROOT->FindObject(canvas_name);
			canvas->Write();
			out_f->Close();
	
			// pdf file
			pdf_name += suffix;
			canvas->Print(pdf_name);
	
			return;
		}
};



int hist_draft(TTree *&tree, h_settings &s){
	if (s.h_min > s.h_max){
	Long64_t n_events = tree->Draw(s.varexp,s.selection);//Returns -1 in case of error or number of selected events in case of success.
	TH1D *htemp = (TH1D*) gROOT->FindObject("htemp");
	if (!htemp) { puts("htemp could not be drawn!"); return -1;}
	const Double_t p[2] = {s.quantile,1-s.quantile};
	Double_t xp[2];
	Int_t two = htemp->GetQuantiles(2,xp,p);
	s.h_min = TMath::FloorNint(xp[0]);
	s.h_max = TMath::CeilNint(xp[1]);
	}
	

	delete dynamic_cast<TH1D*>(gROOT->FindObject(s.h_name));// Even if the object is not found, ```delete nullptr``` is safe.
	TH1D *h = new TH1D(s.h_name, s.h_title+";"+s.x_title+";"+s.y_title, s.bin_per_unit*(s.h_max-s.h_min), s.h_min, s.h_max);
	TString varexp_with_h = s.varexp + ">>" + s.h_name;
	Long64_t n_events = tree->Draw(varexp_with_h,s.selection,s.option);


	if (s.content_thr>0 && s.content_thr<1){
		Double_t ymax = h->GetBinContent(h->GetMaximumBin());
		Int_t top = TMath::Nint(h->GetXaxis()->GetBinLowEdge(h->GetMaximumBin()));
		Int_t above1 = TMath::Nint(h->GetXaxis()->GetBinLowEdge(h->FindFirstBinAbove(0.2*ymax)));
		Int_t above2 = TMath::Nint(h->GetXaxis()->GetBinLowEdge(h->FindLastBinAbove(0.2*ymax)+1));
		s.h_min = top - (top - above1)*15;
		s.h_max = top + (above2 - top)*15;
		delete h;
		h = new TH1D(s.h_name, s.h_title+";"+s.x_title+";"+s.y_title, s.bin_per_unit*(s.h_max-s.h_min), s.h_min, s.h_max  );
		n_events = tree->Draw(varexp_with_h,s.selection,s.option);
	}

	return 0;
}



int one_event(TTree* tree, FitResultSPE &res){

	enum Hist { T, PED, MPED, Q, MQ, MT, COUNT};
	static std::vector<TCanvas*> canvases = {new TCanvas("T"   ,"Timing Distribution"     ,800,600),
						 new TCanvas("PED" ,"Pedestal"                ,800,600),
						 new TCanvas("MPED","MPPC Pedestal"           ,800,600),
						 new TCanvas("Q"   ,"Charge Distribution"     ,800,600),
						 new TCanvas("MQ"  ,"MPPC Charge Distribution",800,600),
						 new TCanvas("MT"  ,"MPPC Timing Distribution",1400,600)};
	std::vector<TString> c_name(COUNT);
	std::vector<TString> pdf_name = {"timing","pedestal","Mpedestal","charge","Mcharge","Mtiming"};
	for (int i = 0; i < COUNT; i++){
		c_name[i] = Form("c%d_%06d", i, res.run);
		canvases[i]->SetName(c_name[i]);
		pdf_name[i] = pdf_name[i]+".pdf";
	}
	TString canvas_root_file = "canvas.root";
	//std::vector<TH1D*>       h(Hist::COUNT);
	TFitResultPtr frp;
	int ret_n;


	//#0 Timing Distribution Histogram
	h_settings setT;
	setT.varexp       = "(NormalHits.time - TriggerHits.time)*1e9";
	setT.selection    = "channel==1&&type==0";
	setT.h_name       = Form("h%d_%d",T,res.run);
	setT.x_title      = "time/ ns";
	setT.h_title      = "Detection Timing Distribution";
	setT.bin_per_unit = 4;//3.99;
	setT.quantile     = 0.001;
	setT.content_thr  = 0.01;
	setT.y_title      = Form("Number of events/bin (%05.4f ns)",1/setT.bin_per_unit);


	canvases[T]->cd();
	ret_n = hist_draft(tree,setT);

	if (ret_n){ return -1;}
	TH1D* h_t = (TH1D*) canvases[T]->GetPrimitive(setT.h_name);	
	frp = TFit::ter(h_t);
        std::tie(res.tmu,res.tsigma,res.ttau) = std::make_tuple(frp->Parameter(1),frp->Parameter(2),frp->Parameter(3));

	if (gROOT->FindObject("ft")){delete (TF1*) gROOT->FindObject("ft");}
	TF1* f_t = new TF1("ft", EMGt, setT.h_min, setT.h_max, 5);
        ExtractCharacteristics(f_t,frp,res.tdc);
	res.pmtftn = f_t->Integral(setT.h_min,setT.h_max)*setT.bin_per_unit;



	//#1 Pedestal Distribution Histogram
	h_settings set;
	set.varexp       = "hgain";
	set.selection    = "channel==1&&type==1";
	set.h_name       = Form("h%d_%d",PED,res.run);
	set.x_title      = "ADC count";
	set.y_title      = "Number of events/bin (1 ADC count)";
	set.h_title      = "pedestal distribution";
	set.bin_per_unit = 1;
	set.content_thr  = -1;
	set.option       = "";// Error if option = "goff".
	
	canvases[PED]->cd();
	ret_n = hist_draft(tree,set);
	if (ret_n) { return -1;}
	TH1D* h_ped = (TH1D*) canvases[PED]->GetPrimitive(set.h_name) ;
	Double_t ped_ave = h_ped->GetMean();



	//#2 MPPC Pedestal Charge Distribution Histogram
	set.selection    = "channel==2&&type==1";
	set.h_name       = Form("h%d_%d",MPED,res.run);
	set.h_min        = 1;
	set.h_max        = -1;
	canvases[MPED]->cd();
	ret_n = hist_draft(tree,set);
	if (ret_n) { return -1;}
	TH1D* h_mped = (TH1D*) canvases[MPED]->GetPrimitive(set.h_name) ;
	Double_t mppc_ped_ave = h_mped->GetMean();
	set.varexp       += Form("-%lf",mppc_ped_ave);
	set.h_min        = 1;
	set.h_max        = -1;
	ret_n = hist_draft(tree,set);
	if (ret_n) { return -1;}
	h_mped = (TH1D*) canvases[MPED]->GetPrimitive(set.h_name) ;



	//#2 Charge Distribution Histogram
	double l_bound, u_bound;
	TFit::GetBound(l_bound,u_bound);
	//std::cout << "lower and upper bound: (" <<  l_bound << ", " << u_bound << ")" << std::endl;
	set.varexp     = Form("hgain-%lf",ped_ave);
	//set.selection  = "channel==1&&type==0";
	set.selection  = Form("channel==1&&type==0&&%s>%f&&%s<%f",setT.varexp.Data(),l_bound,setT.varexp.Data(),u_bound);
	//TString::Format("%.3f", l_bound);
	set.h_name     = Form("h%d_%d",Q,res.run);
	set.h_title    = "SPE distribution";
	set.quantile   = 0.001;
	set.h_min        = 1;
	set.h_max        = -1;

	canvases[Q]->cd();
	ret_n = hist_draft(tree,set);
	if (ret_n) { return -1;}
	TH1D* h_q = (TH1D*) canvases[Q]->GetPrimitive(set.h_name) ;
	res.pmthn = h_q->Integral();

	frp = QFitter(h_q, &canvases[Q]);

	if (gROOT->FindObject("fq")){delete (TF1*) gROOT->FindObject("fq");}
	TF1* f_q = new TF1("fq", SPE, set.h_min, set.h_max, 8);
        ExtractCharacteristics( f_q, frp, res.qdc);
	res.pmtfqn = f_q->Integral(set.h_min,set.h_max)*set.bin_per_unit;


	//#4 MPPC Charge
	set.varexp       = Form("hgain-%lf",mppc_ped_ave);
	set.selection    = "channel==2";
	set.h_name       = Form("h%d_%d", MQ, res.run);
	set.h_title      = "MPPC Charge distribution";
	set.quantile     = 0.0001;
	set.h_min        = 1;
	set.h_max        = -1;
	canvases[MQ]->cd();
	ret_n = hist_draft(tree,set);
	if (ret_n) { return -1;}
	TH1D* h_mq = (TH1D*) canvases[MQ]->GetPrimitive(set.h_name) ;
	h_mped->SetLineColor(2);
	h_mped->Draw("same");
	


	//#5 MPPC timing
	double m_maximumx = h_mq->GetXaxis()->GetBinCenter(h_mq->GetMaximumBin());
	setT.selection    = "channel==2&&type==0";
	setT.selection   += Form("&&hgain-%lf>%lf&&hgain-%lf<%lf",mppc_ped_ave,m_maximumx*0.5,mppc_ped_ave,m_maximumx*1.5);
	setT.h_name       = Form("h%d_1pe_%d",MT,res.run);
	setT.bin_per_unit = 4;//3.99;
	setT.quantile     = 0.001;
	setT.content_thr  = 0.01;
	setT.y_title      = Form("Number of events/bin (%05.4f ns)",1/setT.bin_per_unit);
	setT.h_min        = 1;
	setT.h_max        = -1;
	static bool first_loop = true;
	if (first_loop) {canvases[MT]->Divide(2,1);first_loop = false;}
	canvases[MT]->cd(1);
	ret_n = hist_draft(tree,setT);
	if (ret_n) { return -1;}
	TH1D* h_mt1 = (TH1D*) gPad->GetPrimitive(setT.h_name) ;
	frp = TFit::ter(h_mt1);

	TFit::GetBound(l_bound,u_bound);
	res.mppcn[0] = h_mt1->Integral(h_mt1->FindBin(l_bound),h_mt1->FindBin(u_bound));

	//#5-2 MPPPC 2pe timing
	setT.selection    = "channel==2&&type==0";
	setT.selection   += Form("&&hgain-%lf>%lf&&hgain-%lf<%lf",mppc_ped_ave,m_maximumx*1.5,mppc_ped_ave,m_maximumx*2.5);
	setT.h_name       = Form("h%d_2pe_%d",MT,res.run);
	setT.content_thr  = -1;
	canvases[MT]->cd(2);
	ret_n = hist_draft(tree,setT);
	if (ret_n) { return -1;}
	TH1D* h_mt2 = (TH1D*) gPad->GetPrimitive(setT.h_name) ;
	frp = TFit::ter(h_mt2);
	TFit::GetBound(l_bound,u_bound);
	res.mppcn[1] = h_mt2->Integral(h_mt1->FindBin(l_bound),h_mt1->FindBin(u_bound));

	res.mppcnpe = res.mppcn[0]+res.mppcn[1]*2;

	// Save Canvases
	for (int i = 0; i<COUNT; i++){
		save_can::vas(c_name[i], canvas_root_file, pdf_name[i]);
	}

	return 0;
	
}




void draw_hist(TString out_dir = "../../result/", TString out_f_name = "trash.root"){
	// Output
	save_can::output_dir = out_dir;
	//save_can::output_dir = "/media/hkpd/KIOXIA/06_scan_result_xy/phi_000/";

	// Input
	std::vector<TString> input = {"/media/hkpd/KIOXIA/06_scan/","run","_eventtree_0_500e-9.root"};
	//std::vector<TString> input = {"/media/hkpd/KIOXIA/test/","test","_eventtree_0_500e-9.root"};
	TString tree_name = "event";

	// Run Number
	Int_t FirstRun = 100001;
	Int_t LastRun  = FirstRun+1386;//19*73=1387
	


	
	// Global Parameters settings
	Color_def();
  	gStyle->SetOptFit(1111);
	//gSystem->Load("MetaData/libMetaData.so");
	


	for (int i_run = FirstRun; i_run <= LastRun; i_run++){
		if ((i_run == LastRun) == (i_run == FirstRun)){
			save_can::suffix = "";
	       	} else if ( i_run == FirstRun ){
			save_can::suffix = "(";
		} else if ( i_run == LastRun ){
			save_can::suffix = ")";
		}
		TString in_f_name = input[0]+input[1]+Form("%06d",i_run)+input[2];
		static TFile* in_f = nullptr;
		if (in_f != nullptr ){
		  if (in_f->IsOpen()){
	            in_f->Close();
		  }
		}
		in_f = TFile::Open(in_f_name);
		if (!in_f){std::cout<< "Error opening file\n input file address:" << in_f << std::endl;exit(-1);}
		TTree *tree = (TTree*) in_f->Get(tree_name);
		MetaData* metadata = (MetaData*)in_f->Get("metadata");
		FitResultSPE res;
		res.run    = i_run;
		res.cond[gHV][VAL]    = std::stod((*metadata)["HV"]);
		res.cond[gBX][VAL]    = std::stod((*metadata)["Bx"]);
		res.cond[gBY][VAL]    = std::stod((*metadata)["By"]);
		res.cond[gBZ][VAL]    = std::stod((*metadata)["Bz"]);
		res.cond[gLD][VAL]    = std::stod((*metadata)["LD_current/mA"]);
		res.cond[gTHETA][VAL] = std::stod((*metadata)["theta"]);
		res.cond[gPHI][VAL]   = std::stod((*metadata)["phi"]);
		res.cond[gALPHA][VAL] = std::stod((*metadata)["alpha"]);
		int phi_int = std::stoi((*metadata)["phi"]);
		if ( phi_int >= 90 || phi_int < 0 ){continue;}
		//&& phi_int != "90") && ((*metadata)["phi"] != "180" && (*metadata)["phi"] != "270")){ continue;}
		  std::cout << "Run number: " << i_run << std::endl;
		  std::cout << "theta: " << res.cond[gTHETA][VAL] << std::endl;
		  std::cout << "phi  : " << res.cond[gPHI][VAL] << std::endl;
		

		int ret_m = one_event(tree, res);
		if (ret_m) {continue;}
		result_fill(save_can::output_dir+out_f_name,res);
		
		if (0){
		  std::cout << "###  Summery ###" << std::endl;
		  std::cout << "Condition" << std::endl;
		  std::cout << "Run number: " << i_run << std::endl;
		  std::cout << "theta: " << res.cond[gTHETA][VAL] << std::endl;
		  std::cout << "phi  : " << res.cond[gPHI][VAL] << std::endl;
		  std::cout << "\nCharge Distribution" << std::endl;
		  std::cout << "res.qdc[MPV][VAL]  : " << res.qdc[MPV][VAL] << std::endl;
		  std::cout << "res.qdc[MPV][ERR]  : " << res.qdc[MPV][ERR] << std::endl;
		  std::cout << "\nTiming Distribution" << std::endl;
                  std::cout << "MPV = "  << res.tdc[MPV][VAL]  << "+/-" << res.tdc[MPV][ERR]  << std::endl;
                  std::cout << "FWHM = " << res.tdc[FWHM][VAL] << "+/-" << res.tdc[FWHM][ERR] << std::endl;
		  std::cout << "\nMPPC" << std::endl;
                  std::cout << "single pe = "  << res.mppcn[0]  << std::endl;
                  std::cout << "double pe = "  << res.mppcn[1]  << std::endl;
		  std::cout << "\nPMT" << std::endl;
                  std::cout << "number of pe by histogram = "  << res.pmthn  << std::endl;
                  std::cout << "number of pe by timing function  = "  << res.pmtftn  << std::endl;
                  std::cout << "number of pe by charge function  = "  << res.pmtfqn  << std::endl;

		}
	}
}
