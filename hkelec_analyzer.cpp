#include "declare.h"
//#include "MetaData/MetaData.h"
#include "ExtractCharacteristics.h"
#include "TFitter.h"
#include "mQFitter.h"
#include "Color.C"
#include <tuple>
#include <cstdlib>
#include <algorithm>



enum hltot {TRIGT, TOT, HGAIN, LGAIN, hltotCOUNT};
using Responses = std::array<double,3>;
struct Channel_Signals{
	std::vector<Responses> sig;
	std::vector<Responses> ped;
};
using All_Signals = std::array<Channel_Signals,12>;


class h_settings{
public:
	TString varexp;
	TString selection = "";
	TString option  = "";

	TString h_name  = "";
	TString x_title = "";
	TString y_title = "";
	TString h_title = "";

	hltot tq;

	Double_t bin_width = 1;
	Double_t quantile = 0.01;
	Double_t content_thr = 0.01;

	//T h_min =  1;
	//T h_max = -1;//h_max < h_min (defalt)

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


//int vec_to_TH3(const std::vector<Responses>& content, h_settings &s){
//	//using Iter = std::vector<Responses>::iterator;
//	using PairIter = std::pair<typename std::vector<T>::const_iterator, typename std::vector<T>::const_iterator>;
//	
//	if (s.h_min > s.h_max){
//	  //std::pair<Iter, Iter> minmax_ij = std::minmax_element(content.begin(), content.end(), [s](const auto& a, const auto& b){ return a[s.tq] < b[s.tq]; });
//	  PairIter p = std::minmax_element(vec.begin(),vec.end());
//	  s.h_min = *p.first;
//	  s.h_max = *p.second;
//	  int nbinsx = (s.h_max - s.h_min)/s.bin_width + 1 ;
//	  TH1D* h0 = new TH1D(s.h_name,s.h_title+";"+s.x_title+";"+s.y_title,nbinsx,s.h_min-s.bin_width*0.5,s.h_max+s.bin_width*0.5);
//	  for(std::vector<T> v: content){h0->Fill(v);}
//	  h0->Draw();
////	  const Double_t p[2] = {s.quantile,1-s.quantile};
////	  Double_t xp[2];
////	  Int_t two = h0->GetQuantiles(2,xp,p);
////	  s.h_min = TMath::FloorNint(xp[0]);
////	  s.h_max = TMath::CeilNint(xp[1]);
////	  h0->Draw();
//	  //delete h0;
//	}
//	return 0;
//	
//	//delete dynamic_cast<TH1D*>(gROOT->FindObject(s.h_name));// Even if the object is not found, ```delete nullptr``` is safe.
////	TH1D *h = new TH1D(s.h_name, s.h_title+";"+s.x_title+";"+s.y_title, s.bin_width*(s.h_max-s.h_min), s.h_min, s.h_max);
////	for(std::array<long,4> v: content){h->Fill(v[s.tq]);}
////	//h->Draw();
////	return 0;
//
//
//	//if (s.content_thr>0 && s.content_thr<1){
//	//	Double_t ymax = h->GetBinContent(h->GetMaximumBin());
//	//	Int_t top = TMath::Nint(h->GetXaxis()->GetBinLowEdge(h->GetMaximumBin()));
//	//	Int_t above1 = TMath::Nint(h->GetXaxis()->GetBinLowEdge(h->FindFirstBinAbove(0.2*ymax)));
//	//	Int_t above2 = TMath::Nint(h->GetXaxis()->GetBinLowEdge(h->FindLastBinAbove(0.2*ymax)+1));
//	//	s.h_min = top - (top - above1)*15;
//	//	s.h_max = top + (above2 - top)*15;
//	//	delete h;
//	//	h = new TH1D(s.h_name, s.h_title+";"+s.x_title+";"+s.y_title, s.bin_width*(s.h_max-s.h_min), s.h_min, s.h_max  );
//	//	n_events = tree->Draw(varexp_with_h,s.selection,s.option);
//	//}
//
//}





All_Signals read_eventtree(TTree* tree){

	std::vector<Hit> *v_hit = nullptr;
	std::vector<Hit> *v_trig = nullptr;
	tree->SetBranchAddress("NormalHits", &v_hit);
	tree->SetBranchAddress("TriggerHits", &v_trig);
	
	All_Signals allsig;
	Responses pons;
	/* Process events */
	long nEvents = tree->GetEntries();
	for (long iEvent=0; iEvent<nEvents; ++iEvent) {
	  tree->GetEntry(iEvent);

	  for (const Hit &hit : *v_hit) {
	    if (v_trig->size() != 1){/*std::cout << "This event has " << v_trig->size() << " of triggers! This event is skipped."<< std::endl;*/continue;}
	    Hit trig = v_trig->at(0);
	    if (!std::isfinite(trig.time)){std::cout << " Trigger time is not finite! " << std::endl;continue;}
	    if (!std::isfinite(hit.hgain)||!std::isfinite(hit.lgain)||!isfinite(hit.tdc)){puts("No value is finite!");continue;}
	    pons[TRIGT] = hit.tdc>=trig.tdc?(hit.tdc-trig.tdc)/4.0:-((trig.tdc-hit.tdc)/4.0);//implicit type conversion from unsinged long to double
	    //pons[TRIGT] = (hit.time-trig.time)*1e9;//implicit type conversion from unsinged long to double
	    pons[TOT]   = hit.tot/4.0;// UShort/double
	    pons[HGAIN] = static_cast<double>(hit.hgain);
            if (hit.type == 0){
	      allsig[hit.channel].sig.emplace_back(pons);
            } else if (hit.type == 1){
	      allsig[hit.channel].ped.emplace_back(pons);
	    }
	  }
	}
	return allsig;
}


int PMT_run(Channel_Signals chsig, Channel_Signals mppcsig, FitResultSPE &res){

	enum Hist { T,TQ,TTOT,QTOT, Q, MT, PED, MPED,  MQ, histCOUNT};
	static std::vector<TCanvas*> canvases = {new TCanvas("T"   ,"Timing Distribution"     ,800,600),
						 new TCanvas("TQ"  ,"Timing vs hgain"         ,800,600),
						 new TCanvas("TTOT","Timing vs tot"           ,800,600),
						 new TCanvas("QTOT","hgain vs tot"            ,800,600),
						 new TCanvas("Q"   ,"Charge Distribution"     ,800,600),
						 new TCanvas("MT"  ,"MPPC Timing Distribution",1400,600),
						 new TCanvas("PED" ,"Pedestal"                ,800,600),
						 new TCanvas("MPED","MPPC Pedestal"           ,800,600),
						 new TCanvas("MQ"  ,"MPPC Charge Distribution",800,600)
	                                         };
	std::vector<TString> c_name(histCOUNT);
	std::vector<TString> pdf_name = {"timing","TQmap","TTOT","QTOT","charge","Mtiming","pedestal","Mpedestal","Mcharge"};
	for (int i = 0; i < histCOUNT; i++){
		c_name[i] = Form("c%d_%06d", i, res.run);
		canvases[i]->SetName(c_name[i]);
		pdf_name[i] = pdf_name[i]+".pdf";
	}
	TString canvas_root_file = "canvas.root";
	//std::vector<TH1D*>       h(Hist::COUNT);
//	TFitResultPtr frp;
	int ret_n;


	//#0 Timing Distribution Histogram
	h_settings setT;
	setT.tq           = hltot::TRIGT;
	setT.h_name       = Form("h%d_%d",T,res.run);
	setT.x_title      = "time [ns]";
	setT.h_title      = "Detection Timing Distribution";
	setT.bin_width    = 4;
	setT.quantile     = 0.001;
	setT.content_thr  = 0.01;
	setT.y_title      = "Number of Events / tdc";


	canvases[T]->cd();
//	ret_n = hist_draft(chsig.sig, setT);
//	using Iter = std::vector<Responses>::iterator;
//	std::pair<Iter, Iter> minmax_ij = std::minmax_element(chsig.sig.begin(), chsig.sig.end(), [setT](const auto& a, const auto& b){ return a[setT.tq] < b[setT.tq]; });
//	double minimum = (*minmax_ij.first)[setT.tq];
//	double maximum = (*minmax_ij.second)[setT.tq];
	//TH1D* t_hist = new TH1D("t_hist", setT.h_title, (maximum - minimum)+1,minimum,maximum);
	double mini = 25.125;
	double maxi = 275.125;
	unsigned int nbinx = 1000;
	//TH1D* t_hist = new TH1D("t_hist", setT.h_title, nbinx, mini, maxi);
	std::cout << "size: " <<chsig.sig.size()<< std::endl;
	//for (Responses v: chsig.sig){
	//  auto va = v[setT.tq];
	//  if (va >= mini && va < maxi){
	//    t_hist->Fill(va);
	//  }
	//}
	//t_hist->Draw();

        
	TH3D *h_3d = new TH3D("h3d","",nbinx,mini,maxi,400,-0.125,99.125,300,400,700);
	for (Responses v: chsig.sig){
	  h_3d->Fill(v[0],v[1],v[2]);
	}
	TH2D *h_tq = (TH2D*) h_3d->Project3D("xz");
	canvases[TQ]->cd();
	h_tq->Draw();
	TH2D *h_ttot = (TH2D*) h_3d->Project3D("xy");
	canvases[TTOT]->cd();
	h_ttot->Draw();
	TH2D *h_qtot = (TH2D*) h_3d->Project3D("yz");
	canvases[QTOT]->cd();
	h_qtot->Draw();
	TH1D *h_t = (TH1D*) h_3d->Project3D("x");
	canvases[T]->cd();
	h_t->Draw();
	TH1D *h_q = (TH1D*) h_3d->Project3D("z");
	canvases[Q]->cd();
	h_q->Draw();

	res.pmthn = h_t->GetEntries();
	std::cout << "entry: " << res.pmthn << std::endl;
	




	//MPPC
	double mmini = -100.125;
	double mmaxi = 150.125;
	unsigned int mnbinx = 251;
	TH3D *m_hist = new TH3D("m_hist","",mnbinx,mmini,mmaxi,800,-0.125, 199.125,600,400,1000);
	for (Responses v: mppcsig.sig){
	  if (v[2] > 530){
	    m_hist->Fill(v[0],v[1],v[2]);
	  }
	}
	canvases[MT]->Divide(3,2);
	TH2D *mh_tq = (TH2D*) m_hist->Project3D("xz");
	canvases[MT]->cd(1);
	mh_tq->Draw();
	TH2D *mh_ttot = (TH2D*) m_hist->Project3D("xy");
	canvases[MT]->cd(2);
	mh_ttot->Draw();
	TH2D *mh_qtot = (TH2D*) m_hist->Project3D("yz");
	canvases[MT]->cd(3);
	mh_qtot->Draw();
	canvases[MT]->cd(4);
	//TH3D *m_hist2 = new TH3D("m_hist2","",mnbinx,mmini,mmaxi,800,-0.125, 199.125,600,400,1000);
	//for (Responses v: mppcsig.sig){
	//  if (v[2] > 530){
	//    m_hist2->Fill(v[0],v[1],v[2]);
	//  }
	//}
	//TH2D *mh_tq2 = (TH2D*) m_hist2->Project3D("xz");
	//mh_tq2->Draw();
//	if (ret_n){ return -1;}
//	TH1D* h_t = (TH1D*) canvases[T]->GetPrimitive(setT.h_name);	
//	frp = TFit::ter(h_t);
//        std::tie(res.tmu,res.tsigma,res.ttau) = std::make_tuple(frp->Parameter(1),frp->Parameter(2),frp->Parameter(3));
//
//	if (gROOT->FindObject("ft")){delete (TF1*) gROOT->FindObject("ft");}
//	TF1* f_t = new TF1("ft", EMGt, setT.h_min, setT.h_max, 5);
//        ExtractCharacteristics(f_t,frp,res.tdc);
//	res.pmtftn = f_t->Integral(setT.h_min,setT.h_max)*setT.bin_width;
//
//
//
//	//#1 Pedestal Distribution Histogram
//	h_settings set;
//	set.varexp       = "NormalHits.hgain";
//	set.selection    = "NormalHits.channel==1&&type==1";
//	set.h_name       = Form("h%d_%d",PED,res.run);
//	set.x_title      = "ADC count";
//	set.y_title      = "Number of events/bin (1 ADC count)";
//	set.h_title      = "pedestal distribution";
//	set.bin_width = 1;
//	set.content_thr  = -1;
//	set.option       = "";// Error if option = "goff".
//	
//	canvases[PED]->cd();
//	ret_n = hist_draft(tree,set);
//	if (ret_n) { return -1;}
//	TH1D* h_ped = (TH1D*) canvases[PED]->GetPrimitive(set.h_name) ;
//	Double_t ped_ave = h_ped->GetMean();
//
//
//
//	//#2 MPPC Pedestal Charge Distribution Histogram
//	set.selection    = "NormalHits.channel==2&&type==1";
//	set.h_name       = Form("h%d_%d",MPED,res.run);
//	set.h_min        = 1;
//	set.h_max        = -1;
//	canvases[MPED]->cd();
//	ret_n = hist_draft(tree,set);
//	if (ret_n) { return -1;}
//	TH1D* h_mped = (TH1D*) canvases[MPED]->GetPrimitive(set.h_name) ;
//	Double_t mppc_ped_ave = h_mped->GetMean();
//	set.varexp       += Form("-%lf",mppc_ped_ave);
//	set.h_min        = 1;
//	set.h_max        = -1;
//	ret_n = hist_draft(tree,set);
//	if (ret_n) { return -1;}
//	h_mped = (TH1D*) canvases[MPED]->GetPrimitive(set.h_name) ;
//
//
//
//	//#2 Charge Distribution Histogram
//	double l_bound, u_bound;
//	TFit::GetBound(l_bound,u_bound);
//	std::cout << "lower and upper bound: (" <<  l_bound << ", " << u_bound << ")" << std::endl;
//	set.varexp     = Form("NormalHits.hgain-%lf",ped_ave);
//	//set.selection  = "NormalHits.channel==1&&NormalHits.type==0";
//	set.selection  = Form("NormalHits.channel==1&&NormalHits.type==0&&%s>%f&&%s<%f",setT.varexp.Data(),l_bound,setT.varexp.Data(),u_bound);
//	std::cout << Form("NormalHits.channel==1&&NormalHits.type==0&&%s>%f&&%s<%f",setT.varexp.Data(),l_bound,setT.varexp.Data(),u_bound) << std::endl;
//	//TString::Format("%.3f", l_bound);
//	set.h_name     = Form("h%d_%d",Q,res.run);
//	set.h_title    = "SPE distribution";
//	set.quantile   = 0.001;
//	set.h_min        = 1;
//	set.h_max        = -1;
//
//	canvases[Q]->cd();
//	ret_n = hist_draft(tree,set);
//	if (ret_n) { return -1;}
//	TH1D* h_q = (TH1D*) canvases[Q]->GetPrimitive(set.h_name) ;
//	res.pmthn = h_q->Integral();
//
//	frp = QFitter(h_q, &canvases[Q]);
//
//	if (gROOT->FindObject("fq")){delete (TF1*) gROOT->FindObject("fq");}
//	TF1* f_q = new TF1("fq", SPE, set.h_min, set.h_max, 8);
//        ExtractCharacteristics( f_q, frp, res.qdc);
//	res.pmtfqn = f_q->Integral(set.h_min,set.h_max)*set.bin_width;
//
//
//	//#4 MPPC Charge
//	set.varexp       = Form("NormalHits.hgain-%lf",mppc_ped_ave);
//	set.selection    = "NormalHits.channel==2";
//	set.h_name       = Form("h%d_%d", MQ, res.run);
//	set.h_title      = "MPPC Charge distribution";
//	set.quantile     = 0.0001;
//	set.h_min        = 1;
//	set.h_max        = -1;
//	canvases[MQ]->cd();
//	ret_n = hist_draft(tree,set);
//	if (ret_n) { return -1;}
//	TH1D* h_mq = (TH1D*) canvases[MQ]->GetPrimitive(set.h_name) ;
//	h_mped->SetLineColor(2);
//	h_mped->Draw("same");
//	
//
//
//	//#5 MPPC timing
//	double m_maximumx = h_mq->GetXaxis()->GetBinCenter(h_mq->GetMaximumBin());
//	setT.selection    = "NormalHits.channel==2&&NormalHits.type==0";
//	setT.selection   += Form("&&NormalHits.hgain-%lf>%lf&&NormalHits.hgain-%lf<%lf",mppc_ped_ave,m_maximumx*0.5,mppc_ped_ave,m_maximumx*1.5);
//	setT.h_name       = Form("h%d_1pe_%d",MT,res.run);
//	setT.bin_width = 4;//3.99;
//	setT.quantile     = 0.001;
//	setT.content_thr  = 0.01;
//	setT.y_title      = Form("Number of events/bin (%05.4f ns)",1/setT.bin_width);
//	setT.h_min        = 1;
//	setT.h_max        = -1;
//	static bool first_loop = true;
//	if (first_loop) {canvases[MT]->Divide(2,1);first_loop = false;}
//	canvases[MT]->cd(1);
//	ret_n = hist_draft(tree,setT);
//	if (ret_n) { return -1;}
//	TH1D* h_mt1 = (TH1D*) gPad->GetPrimitive(setT.h_name) ;
//	frp = TFit::ter(h_mt1);
//
//	TFit::GetBound(l_bound,u_bound);
//	res.mppcn[0] = h_mt1->Integral(h_mt1->FindBin(l_bound),h_mt1->FindBin(u_bound));
//
//	//#5-2 MPPPC 2pe timing
//	setT.selection    = "NormalHits.channel==2&&NormalHits.type==0";
//	setT.selection   += Form("&&NormalHits.hgain-%lf>%lf&&NormalHits.hgain-%lf<%lf",mppc_ped_ave,m_maximumx*1.5,mppc_ped_ave,m_maximumx*2.5);
//	setT.h_name       = Form("h%d_2pe_%d",MT,res.run);
//	setT.content_thr  = -1;
//	canvases[MT]->cd(2);
//	ret_n = hist_draft(tree,setT);
//	if (ret_n) { return -1;}
//	TH1D* h_mt2 = (TH1D*) gPad->GetPrimitive(setT.h_name) ;
//	frp = TFit::ter(h_mt2);
//	TFit::GetBound(l_bound,u_bound);
//	res.mppcn[1] = h_mt2->Integral(h_mt1->FindBin(l_bound),h_mt1->FindBin(u_bound));
//
//	res.mppcnpe = res.mppcn[0]+res.mppcn[1]*2;
//
	// Save Canvases
	for (int i = 0; i<2/*<histCOUNT*/; i++){
		save_can::vas(c_name[i], canvas_root_file, pdf_name[i]);
	}

	return 0;
      
}




void hkelec_analyzer(TString out_f_name = "trash.root"){
	// Input
	//std::vector<TString> input = {"/mnt/f/06_scan_11_2000V_5s/","run","_eventtree_0_500e-9.root"};
	//std::vector<TString> input = {"/mnt/f/06_scan/","run","_eventtree_0_500e-9.root"};
	//std::vector<TString> input = {"/home/keio/horiuchi/hkdigitizer/scripts/","trash","_eventtree.root"};
	std::vector<TString> input = {"/home/hkpd/individuals/horiuchi/data/IO-data/test/","run","_eventtree_-200e-9_425e-9.root"};
	//std::vector<TString> input = {"/media/hkpd/KIOXIA/06_scan/","run","_eventtree_0_500e-9.root"};
	//std::vector<TString> input = {"/media/hkpd/KIOXIA/test/","test","_eventtree_0_500e-9.root"};
	TString tree_name = "event";
	int PMTch = 5;
	int MPPCch = 7;

	// Run Number
	Int_t FirstRun = 9;
	//Int_t FirstRun = 116644;
	Int_t LastRun  = 9;//1387-18;//FirstRun+0;//19*73=1387


	// Output
	//save_can::output_dir = "/media/hkpd/KIOXIA/06_scan_result_xy/phi_000/";
	TString out_dir = input[0]+"_result/";
	save_can::output_dir = out_dir+"figures/";
	gSystem->Exec(TString("mkdir -p ")+save_can::output_dir);


	
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
		FitResultSPE res;
		res.run    = i_run;
		MetaData* metadata = (MetaData*)in_f->Get("metadata");
		if (metadata->HasKey("HV"))res.cond[gHV][VAL]    = std::stod((*metadata)["HV"]);
		if (metadata->HasKey("Bx"))  res.cond[gBX][VAL]    = std::stod((*metadata)["Bx"]);
		if (metadata->HasKey("By"))res.cond[gBY][VAL]    = std::stod((*metadata)["By"]);
		if (metadata->HasKey("Bz"))  res.cond[gBZ][VAL]    = std::stod((*metadata)["Bz"]);
		if (metadata->HasKey("LDmA"))res.cond[gLD][VAL]    = std::stod((*metadata)["LDmA"]);
		if (metadata->HasKey("theta")) res.cond[gTHETA][VAL] = std::stod((*metadata)["theta"]);
		if (metadata->HasKey("phi"))   res.cond[gPHI][VAL]   = std::stod((*metadata)["phi"]);
		if (metadata->HasKey("alpha")) res.cond[gALPHA][VAL] = std::stod((*metadata)["alpha"]);
		//int phi_int = std::stoi((*metadata)["phi"]);
		if ( res.cond[gTHETA][VAL] > 1 ){continue;}
		std::cout << "theta: " << res.cond[gTHETA][VAL] << std::endl;
		std::cout << "save_can::suffix: " << save_can::suffix << std::endl;
		//&& phi_int != "90") && ((*metadata)["phi"] != "180" && (*metadata)["phi"] != "270")){ continue;}
		
	        All_Signals allsig = read_eventtree(tree);
              
		int ret_m = PMT_run(allsig[PMTch], allsig[MPPCch], res);
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
