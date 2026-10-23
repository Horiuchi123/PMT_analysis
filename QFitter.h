// This code requires root6. root5 can't be used. 
Double_t SPE(Double_t *x, Double_t *par){
   double mah0 = (x[0]-par[2])/par[3]/TMath::Sqrt2();//Mahalanobis distance/√2
   double mah1 = (x[0]-par[2])/par[3]/par[4]/TMath::Sqrt2();
   double mah2 = (x[0]-par[2]*par[5])/par[6]/TMath::Sqrt2();
   return par[0] *( par[1] /(par[2]*(1-par[5])) /2. *( TMath::Erfc(mah1) - TMath::Erfc(mah2) ) + (1-par[1])/(TMath::Sqrt(2*TMath::Pi())*par[3]) * exp(- mah0*mah0) );
}
Double_t SPE2(Double_t *x, Double_t *par){
   double mah0 = (x[0]-par[2])/par[3]/TMath::Sqrt2();//Mahalanobis distance/√2
   double mah1 = (x[0]-par[2])/par[4]/TMath::Sqrt2();
   double mah2 = (x[0]-par[5])/par[6]/TMath::Sqrt2();
   return par[0] *( par[1] /(par[2]-par[5]) /2. *( TMath::Erfc(mah1) - TMath::Erfc(mah2) ) + (1-par[1])/(TMath::Sqrt(2*TMath::Pi())*par[3]) * exp(- mah0*mah0) );
}


struct DPE {
   using Func = double (*)(double*, double*);
   DPE(Func fFunc, Double_t xmin, Double_t xmax)  : f(fFunc), a(xmin), b(xmax) {}
   double operator()(double *x, double *par) const {
      int npoint = 100;
      double doub = 0;
      double weight;
      double it[1], ixt[1];
      for (int i = 0; i <= npoint; ++i){
        it[0]  = a + i * (b - a) / npoint;
        ixt[0] = x[0] - it[0];
        weight = (i == 0 || i == npoint) ? 0.5 : 1.0;
        doub += weight * f(it, par) * f(ixt, par) / par[0];
      }
      doub *= (b - a) / npoint;
      return (1 - par[7]) * f(x, par) + par[7] * doub;
   }
   Func f;
   Double_t a;
   Double_t b;
};

 

TFitResultPtr QFitter(TH1D*& h, TCanvas**c = nullptr, Bool_t del = false){

  TH1D *h1 = (TH1D*)h->Clone();
  TH1D *h2 = (TH1D*)h->Clone(); 
  h1->SetName("h2r");
  TCanvas *c1 = new TCanvas("c1r","single photoelectron",1400,800);

  c1->Divide(2,2);
  c1->cd(1);
  double xmin = h1->GetBinLowEdge(1);
  double xmax = h1->GetBinLowEdge(h1->GetNbinsX() + 1);


  h1->Fit("gaus","Q");
  TF1 *gaus = h1->GetFunction("gaus");

  TF1* spe[3];
  for ( int i_fnc = 0; i_fnc < 3; i_fnc++){
  spe[i_fnc] = new TF1(Form("spe[%d]",i_fnc),SPE,xmin,xmax,7);
  spe[i_fnc]->SetNpx(1000);
  spe[i_fnc]->SetLineColor(i_fnc+2);
  }
  spe[0]->SetParNames("N","erf_ratio","#mu","#sigma","#sigma_ratio","trig_thr","#sigma_thr");
  spe[0]->SetParameter(0,gaus->GetParameter(0)*TMath::Sqrt(2*TMath::Pi())*gaus->GetParameter(2));
  spe[0]->SetParLimits(0,gaus->GetParameter(0)*TMath::Sqrt(2*TMath::Pi())*gaus->GetParameter(2)*0.1,gaus->GetParameter(0)*TMath::Sqrt(2*TMath::Pi())*gaus->GetParameter(2)*2.);
  spe[0]->SetParameter(1,0.5);
  spe[0]->SetParLimits(1,0.01,0.99);
  spe[0]->SetParameter(2,gaus->GetParameter(1));
  spe[0]->SetParameter(3,gaus->GetParameter(2));
  spe[0]->FixParameter(4,1/TMath::Sqrt2());
  //spe[0]->SetParLimits(5,xmin,gaus->GetParameter(1));
  spe[0]->SetParameter(5,0.1);
  spe[0]->SetParLimits(5,0.01,0.99);
  spe[0]->SetParameter(6,0);
  spe[0]->SetParLimits(6,0   ,gaus->GetParameter(2));
  TFitResultPtr frp = h1->Fit("spe[0]","SQ");

  spe[1]->SetParameter(0,frp->Parameter(0)*(1-frp->Parameter(1)));
  spe[2]->SetParameter(0,frp->Parameter(0)* frp->Parameter(1));
  spe[1]->SetParameter(1,0);
  spe[2]->SetParameter(1,1);
  
  for (int i_fnc = 1; i_fnc < 3; i_fnc++){
    for (int i_par = 2; i_par < 7; i_par++){
      spe[i_fnc]->SetParameter(i_par,frp->Parameter(i_par));
    }
    spe[i_fnc]->Draw("same");
  }


  // Second Fit
  c1->cd(2);
  DPE*  dpe_by_spe = new DPE(SPE,-10,90);
  TF1* dpe[4];
  for ( int i = 0; i<4; i++){
    dpe[i] = new TF1(Form("dpe[%d]",i),dpe_by_spe,xmin,xmax,8);
    //dpe[i] = new TF1(Form("dpe[%d]",i),DPE,xmin,xmax,8);
    dpe[i]->SetNpx(1000);
    dpe[i]->SetLineColor(2+i);
  }
  double spe_mean = spe[0]->Mean(xmin,xmax);
  double his_mean = h2->GetMean();
  double dpe_ratio_est = TMath::Min(0.1, TMath::Abs(1 - spe_mean/his_mean)*0.5);
  dpe[0]->SetParNames("N","erf_ratio","#mu","#sigma","#sigma_ratio","trig_thr","#sigma_thr","DPE_ratio");
  dpe[0]->SetParameters(frp->Parameter(0),frp->Parameter(1),frp->Parameter(2),frp->Parameter(3),frp->Parameter(4),frp->Parameter(5),frp->Parameter(6),dpe_ratio_est );
  //dpe[0]->SetParLimits(0,frp->GetParameter(0)*0.1,gaus->GetParameter(0)*2.);
  dpe[0]->SetParLimits(1,0.01,0.99);
  dpe[0]->SetParLimits(5,0.01,0.99);
  dpe[0]->SetParLimits(7,0   ,0.1);
  TFitResultPtr frp2 = h2->Fit("dpe[0]","S","",0,80);


  dpe[1]->SetParameter(0,frp2->Parameter(0)*(1-frp2->Parameter(7))*(1-frp2->Parameter(1)) );
  dpe[2]->SetParameter(0,frp2->Parameter(0)*(1-frp2->Parameter(7))* frp2->Parameter(1));
  dpe[3]->SetParameter(0,frp2->Parameter(0)*frp2->Parameter(7) );
  dpe[1]->SetParameter(1,0);
  dpe[2]->SetParameter(1,1);
  dpe[3]->SetParameter(1,frp2->Parameter(1));
  dpe[1]->SetParameter(7,0);
  dpe[2]->SetParameter(7,0);
  dpe[3]->SetParameter(7,1);
  for (int i_fnc = 1; i_fnc < 4; i_fnc++){
    for (int i_par = 2; i_par < 7; i_par++){
      dpe[i_fnc]->SetParameter(i_par,frp2->Parameter(i_par));
    }
  dpe[i_fnc]->Draw("same");
  }


  // Third Fit
  c1->cd(3);
  DPE*  dpe_by_spe2 = new DPE(SPE2,-10,90);
  const int n_f = 4;
  TF1*  dpe2[n_f];
  for ( int i = 0; i<n_f; i++){
    dpe2[i] = new TF1(Form("dpe2[%d]",i),dpe_by_spe2,xmin,xmax,8);
    dpe2[i]->SetNpx(1000);
    dpe2[i]->SetLineColor(2+i);
  }
  dpe2[0]->SetParNames("N","erf_ratio","#mu","#sigma","#sigma_erf","trig_thr","#sigma_thr","DPE_ratio");
  dpe2[0]->SetParameters(frp2->Parameter(0),frp2->Parameter(1),frp2->Parameter(2),frp2->Parameter(3),frp2->Parameter(3)*frp2->Parameter(4),frp2->Parameter(2)*frp2->Parameter(5),frp2->Parameter(6),frp2->Parameter(7) );
  TFitResultPtr frp3 = h->Fit("dpe2[0]","S");

  dpe2[1]->SetParameter(0,frp3->Parameter(0)*(1-frp3->Parameter(7))*(1-frp3->Parameter(1)) );
  dpe2[2]->SetParameter(0,frp3->Parameter(0)*(1-frp3->Parameter(7))* frp3->Parameter(1));
  dpe2[3]->SetParameter(0,frp3->Parameter(0)*frp3->Parameter(7) );
  dpe2[1]->SetParameter(1,0);
  dpe2[2]->SetParameter(1,1);
  dpe2[3]->SetParameter(1,frp3->Parameter(1));
  dpe2[1]->SetParameter(7,0);
  dpe2[2]->SetParameter(7,0);
  dpe2[3]->SetParameter(7,1);
  for (int i_fnc = n_f - 1; i_fnc >= 0 ; i_fnc--){
    for (int i_par = 2; i_par < 7; i_par++){
      dpe2[i_fnc]->SetParameter(i_par,frp3->Parameter(i_par));
    }
    dpe2[i_fnc]->Draw("same");
    if ( c!= nullptr){
      (*c)->cd();
      dpe2[i_fnc]->Draw("same");
      c1->cd(3);
    }
  }



  if (del){
    delete c1;
    delete dpe_by_spe;
    delete dpe_by_spe2;
    delete spe[0]; delete spe[1]; delete spe[2];
    delete dpe[0]; delete dpe[1]; delete dpe[2]; delete dpe[3];
  }
  return frp3;


}
/*
void Mean(TF1* func, TMatrixDSym cov,double *mean, double *variance, double *stddev, double *resolution){
  double xmin,xmax;  func->GetRange(xmin,xmax);
  unsigned int Npar = func->GetNpar();
  
  TF1* xfunc  = new TF1("xfunc",[func](double *x, double *par){ return x[0]*func->EvalPar(x,par); }, xmin, xmax, Npar);
  TF1* xxfunc = new TF1("xxfunc",[func](double *x, double *par){ return x[0]*x[0]*func->EvalPar(x,par); }, xmin, xmax, Npar);
  TF1* Fmean  = new TF1("Fmean",[xfunc,func,xmin,xmax](double *x, double *par){ xfunc->SetParameters(par);func->SetParameters(par); return xfunc->Integral(xmin,xmax)/func->Integral(xmin,xmax);},0,1,Npar);//This is a constant function: ∂f/∂x = 0
  TF1* Fvari = new TF1("Fvari",[xxfunc,func,Fmean,xmin,xmax](double *x, double *par){ xxfunc->SetParameters(par);func->SetParameters(par);Fmean->SetParameters(par); return xxfunc->Integral(xmin,xmax)/func->Integral(xmin,xmax) - TMath::Sq(Fmean->Eval(0));},0,1,Npar);
  TF1* Freso = new TF1("Freso",[Fmean,Fvari,xmin,xmax](double *x, double *par){return TMath::Sqrt(Fvari->EvalPar(x,par))/Fmean->EvalPar(x,par);},0,1,Npar);
  Fmean->SetParameters(func->GetParameters());
  Fvari->SetParameters(func->GetParameters());
  Freso->SetParameters(func->GetParameters());
  double* dfdp_arr_mean = new double[Npar];
  double* dfdp_arr_vari = new double[Npar];
  double* dfdp_arr_reso = new double[Npar];
  double nonsense = 0;
  Fmean->GradientPar(&nonsense,dfdp_arr_mean,1e-6);
  Fvari->GradientPar(&nonsense,dfdp_arr_vari,1e-6);
  Freso->GradientPar(&nonsense,dfdp_arr_reso,1e-6);
  TVectorT<double> dfdp_m(Npar,dfdp_arr_mean);
  TVectorT<double> dfdp_v(Npar,dfdp_arr_vari);
  TVectorT<double> dfdp_r(Npar,dfdp_arr_reso);
  mean[1] = TMath::Sqrt(cov.Similarity(dfdp_m));
  variance[1] = TMath::Sqrt(cov.Similarity(dfdp_v));
  resolution[1] = TMath::Sqrt(cov.Similarity(dfdp_r));

  mean[0] = Fmean->Eval(0);
  variance[0] = Fvari->Eval(0);
  resolution[0] = Freso->Eval(0);
  stddev[0] = TMath::Sqrt(variance[0]);
  stddev[1] = 0.5*variance[1]/TMath::Sqrt(variance[0]);

  delete[] dfdp_arr_mean; delete[] dfdp_arr_vari; delete[] dfdp_arr_reso;

}
 *
  xfunc->SetParameters(func->GetParameters());
  xfunc->SetParErrors(func->GetParErrors());
  double mean2 = xfunc->Integral(xmin,xmax)/area;
  double mean2_err = xfunc->IntegralError(xmin,xmax)/area;
  std::cout << "mean2: " << mean2 << "\tmean2_err: " << mean2_err << std::endl;
  for (int i = 0; i<Npar;i++){
    std::cout << "dfdp["<<i<<"] = " << dfdp_arr[i] << std::endl;
  }
  cov.Print();
  std::cout << "area:  " << area  << "\tarea_err:  " << area_err  << std::endl;
  double area     = func->Integral(xmin,xmax);
  double area_err = func->IntegralError(xmin,xmax);
  double mean     = func->Mean(xmin,xmax);// = func->Moment(1,xmin,xmax);
  double vari = func->Variance(xmin,xmax);// = func->CentralMoment(2,xmin,xmax);
Double_t SPE_mod_pre(Double_t *x, Double_t *par){
   double mah0 = (x[0]-par[1])/par[2]/TMath::Sqrt2();//Mahalanobis distance/√2
   double mah1 = (x[0]-par[1])/par[2]/par[4]/TMath::Sqrt2();
   double mah2 = (x[0]-par[5])/par[6]/TMath::Sqrt2();
   return par[0]/4. *TMath::Erfc(mah1) *TMath::Erfc(-mah2) *(par[3] + (1-par[3])*exp(- mah0*mah0));
}
Double_t DPE_old_version(Double_t *x, Double_t *par){
   double tmin = -10;
   double tmax = 90;
   int npoint = 100;
   double doub = 0;
   double weight;
   double it[1],ixt[1];
   for (int i = 0; i<=npoint; ++i){
     it[0] = tmin + i*(tmax-tmin)/npoint;
     ixt[0] = x[0] - it[0];
     weight = ((i == 0 || i == npoint)?0.5:1.);
     doub += weight*SPE(it,par)*SPE(ixt,par)/par[0];
   }
   doub *= (tmax-tmin)/npoint;
   return (1-par[7])*SPE(x,par) + par[7]*doub;
}
  */
  //TMatrixDSym cov2 = frp2->GetCovarianceMatrix();
  //Mean(dpe[4],cov2,res.f2_mean,res.f2_variance,res.f2_stddev,res.f2_resolution);


  // Summary
  //std::cout << "mean: " << res.f_mean[0] << "+/-" << res.f_mean[1]  << std::endl;
  //std::cout << "mean_qdc: " << res.qdc[MEAN][VAL] << "+/-" << res.qdc[MEAN][ERR]  << std::endl;
  //std::cout << "mean2: " << res.f2_mean[0] << "+/-" << res.f2_mean[1]  << std::endl;
  //std::cout << "stddev: " << res.f_stddev[0] << "+/-" << res.f_stddev[1]  << std::endl;
  //std::cout << "stddev2: " << res.f2_stddev[0] << "+/-" << res.f2_stddev[1]  << std::endl;
  //std::cout << "stddev2_qdc: " << res.qdc[STDDEV][VAL] << "+/-" << res.qdc[STDDEV][ERR]  << std::endl;
  //std::cout << "resolution2: " << res.f2_resolution[0] << "+/-" << res.f2_resolution[1]  << std::endl;
  //std::cout << "resolution2_qdc: " << res.qdc[RESO][VAL] << "+/-" << res.qdc[RESO][ERR]  << std::endl;
  //std::cout << "mean_diff:  " << res.f2_mean[0] - dpe[4]->Mean(xmin,xmax)  << std::endl;
  //std::cout << "stdd_diff:  " << res.f2_stddev[0] - res.f_stddev[0]  << std::endl;
  //std::cout << "vari_diff:  " << res.f2_variance[0] - res.f_variance[0]  << std::endl;
  //std::cout << "vari_diff:  " << res.f2_variance[0] - dpe[4]->Variance(xmin,xmax)  << std::endl;
  //
  //std::cout << "MPV: " << res.qdc[MPV][VAL] << "+/-" << res.qdc[MPV][ERR]  << std::endl;
  //std::cout << "FWHM: " << res.qdc[FWHM][VAL] << "+/-" << res.qdc[FWHM][ERR]  << std::endl;
//void Mean(TF1* func, TMatrixDSym cov,double *mean, double *variance, double *stddev, double *resolution);
//#include "declare.h"
//#ifndef CHARA
//#define CHARA
//#include "ExtractCharacteristics.h"
////DistChara ExtractCharacteristics(TF1* func, TFitResultPtr frp)
//#endif
//#include <functional>
