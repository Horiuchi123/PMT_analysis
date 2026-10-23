// This code requires root6. root5 can't be used. 
#include <boost/math/special_functions/gamma.hpp>
struct Multi_Photoelectron{
  Multi_Photoelectron(int number_of_photoelectron): m(number_of_photoelectron) {}
  double operator() (double *x, double *par) const {
    if (par[6]<=0||par[5]<=0){return 0;}
    double ret = 0;	    
    for (int i = 1; i <= m; i++){
      ret += NthPE(i,x,par);	   
    }
    return par[0] * ret * (0.5)* TMath::Erfc( -(x[0]-par[5])/(TMath::Sqrt2()*par[6]) );
  }


  int m;
  double NthPE(int n, double *x, double *par) const {
    if (par[2] <= 0||par[7]<=0){return 0;}
    if (par[1] <  0||par[1] > 1){return 0;}
    auto f = [&](int k, double es, double mu){
      if (es <= 0){return 0.0;}
      double z = (x[0]-mu)*(x[0]-mu)/(2.0*es);
      double norm = 1/TMath::Sqrt(2*TMath::Pi()*es);
      if (!k){
	return norm * TMath::Exp(-z);
      } else {
        double ret = 0;
        for (int i = 0; i <= k-1; i++){
	  double ii = (i + 1)*0.5;
          double coeff = TMath::Power(x[0]-mu,k-1-i)/(TMath::Factorial(i)*TMath::Factorial(k-1-i)) *TMath::Power(2*es,ii);
	  double gamma= (x[0]>mu && (i%2 == 0)) ? std::tgamma(ii) + boost::math::tgamma_lower(ii, z) : boost::math::tgamma(ii, z);
	  ret += coeff*gamma;
	}
	return ret*norm/2.0;
      }
    };

    double val = 0;
    for (int k = 0; k <= n; k++){
      double sum = 0; 
      for ( int ell = 0; ell <= k; ell++){
	double ef;
	if (n==k&&ell==0){ef = ((x[0] > 0 ? 1.0:0.0)*TMath::Power(x[0] - 0,k-1)/TMath::Factorial(k-1));}
	else {ef = f(k,(n-k)*par[3]+ell*par[4], (n-k+ell)*par[2] ); }
        sum += TMath::Binomial(k,ell)*(ell%2?-1:1)* ef ; //( (n==k&&ell==0)?((x[0]>mu?1.0:0.0)*TMath::Power(x[0]-mu,k-1)/TMath::Factorial(k-1)):f(k,(n-k)*par[3]+ell*par[4], (n-k+ell)*par[2] ) );
      }
      val += TMath::Binomial(n,k) * TMath::Power(1-par[1],n-k) * TMath::Power(par[1]/par[2],k) * sum;
    }
    return ROOT::Math::poisson_pdf(n,par[7]) * val / (1 - ROOT::Math::poisson_pdf(0,par[7]));
  }
};



TFitResultPtr mQFitter(TH1D*& h, TCanvas**c = nullptr){

  TH1D *h1 = (TH1D*)h->Clone();
  TH1D *h2 = (TH1D*)h->Clone(); 
  h1->SetName("h2r");
  static TCanvas *c1 = new TCanvas("c1r","single photoelectron",1400,800);
  static bool first_call = true;
  if (first_call){
  	c1->Divide(2,2);
	first_call = false;
  }
  c1->cd(1);
  double xmin = h1->GetBinLowEdge(1);
  double xmax = h1->GetBinLowEdge(h1->GetNbinsX() + 1);


  h1->Fit("gaus","Q");
  TF1 *gaus = h1->GetFunction("gaus");


  static Multi_Photoelectron * sPE = new Multi_Photoelectron(1);
  static Multi_Photoelectron * dPE = new Multi_Photoelectron(2);
  static Multi_Photoelectron * tPE = new Multi_Photoelectron(3);
  static Multi_Photoelectron * qPE = new Multi_Photoelectron(4);
  TF1* mPE[3];
  for ( int i_fnc = 0; i_fnc < 3; i_fnc++){
  if (gROOT->FindObject( Form("mPE[%d]",i_fnc) )){ delete (TF1*) gROOT->FindObject(Form("mPE[%d]",i_fnc));}
  mPE[i_fnc] = new TF1(Form("mPE[%d]",i_fnc),tPE,xmin,xmax,8);
  mPE[i_fnc]->SetNpx(1000);
  mPE[i_fnc]->SetLineColor(i_fnc+2);
  }
  mPE[0]->SetParNames("N","erf_ratio","#mu","#sigma^{2}","#sigma_{erf}^{2}","trig_thr","#sigma_{thr}^{2}");
  mPE[0]->SetParameter(0,gaus->GetParameter(0)*TMath::Sqrt(2*TMath::Pi())*gaus->GetParameter(2));
  mPE[0]->SetParLimits(0,gaus->GetParameter(0)*TMath::Sqrt(2*TMath::Pi())*gaus->GetParameter(2)*0.1,gaus->GetParameter(0)*TMath::Sqrt(2*TMath::Pi())*gaus->GetParameter(2)*2.);
  mPE[0]->SetParameter(1,0.5);
  mPE[0]->SetParLimits(1,0.01,0.99);
  mPE[0]->SetParameter(2,gaus->GetParameter(1));
  mPE[0]->SetParameter(3,TMath::Sq(gaus->GetParameter(2)));
  mPE[0]->FixParameter(4,TMath::Sq(gaus->GetParameter(2))/2.0);
  //mPE[0]->SetParLimits(5,xmin,gaus->GetParameter(1));
  mPE[0]->SetParameter(5,gaus->GetParameter(1)*0.1);
  //mPE[0]->SetParLimits(5,0.01,0.99);
  mPE[0]->SetParameter(6,gaus->GetParameter(2)*0.1);
  //mPE[0]->SetParLimits(6,0   ,gaus->GetParameter(2));
  mPE[0]->SetParameter(7,0.1);
  TFitResultPtr frp = h1->Fit("mPE[0]","SQ");

  mPE[1]->SetParameter(0,frp->Parameter(0)*(1-frp->Parameter(1)));
  mPE[2]->SetParameter(0,frp->Parameter(0)* frp->Parameter(1));
  mPE[1]->SetParameter(1,0);
  mPE[2]->SetParameter(1,1);
  
  for (int i_fnc = 1; i_fnc < 3; i_fnc++){
    for (int i_par = 2; i_par < 7; i_par++){
      mPE[i_fnc]->SetParameter(i_par,frp->Parameter(i_par));
    }
    mPE[i_fnc]->Draw("same");
  }
  return frp;

}







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

 

TFitResultPtr QFitter(TH1D*& h, TCanvas**c = nullptr){

  TH1D *h1 = (TH1D*)h->Clone();
  TH1D *h2 = (TH1D*)h->Clone(); 
  TH1D *h3 = (TH1D*)h->Clone(); 
  h1->SetName("h1");
  h2->SetName("h2");
  static TCanvas *c1 = new TCanvas("c1","single photoelectron",1400,800);
  static bool first_call = true;
  if (first_call){
  	c1->Divide(2,2);
	first_call = false;
  }
  c1->cd(1);
  double xmin = h1->GetBinLowEdge(1);
  double xmax = h1->GetBinLowEdge(h1->GetNbinsX() + 1);

  h1->Fit("gaus","Q");
  TF1 *gaus = h1->GetFunction("gaus");

  TF1* spe[3];
  for ( int i_fnc = 0; i_fnc < 3; i_fnc++){
  if (gROOT->FindObject( Form("spe[%d]",i_fnc) )){ delete (TF1*) gROOT->FindObject(Form("spe[%d]",i_fnc));}
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
  static DPE*  dpe_by_spe = new DPE(SPE,-10,90);
  const int n_f = 4;
  TF1* dpe[n_f];
  for ( int i = 0; i<n_f; i++){
    if ( gROOT->FindObject( Form("dpe[%d]",i) ) ){ delete (TF1*) gROOT->FindObject(Form("dpe[%d]",i));}
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
  for (int i_fnc = 1; i_fnc < n_f; i_fnc++){
    for (int i_par = 2; i_par < 7; i_par++){
      dpe[i_fnc]->SetParameter(i_par,frp2->Parameter(i_par));
    }
  dpe[i_fnc]->Draw("same");
  }


  // Third Fit
  c1->cd(3);
  static DPE*  dpe_by_spe2 = new DPE(SPE2,-10,90);
  TF1*  dpe2[n_f];
  for ( int i = 0; i<n_f; i++){
    if ( gROOT->FindObject( Form("dpe2[%d]",i) ) ){ delete (TF1*) gROOT->FindObject(Form("dpe2[%d]",i));}
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
  return frp3;

  // Final Fit
  //c1->cd(4);
  //static Multi_Photoelectron * TPE = new Multi_Photoelectron(3);
  //TF1* tpe;
  //tpe = new TF1("tpe",TPE,xmin,xmax,8);
  //tpe->SetNpx(1000);
  //tpe->SetParNames("N","erf_ratio","#mu","#sigma^{2}","#sigma_{erf}^{2}","trig_thr","#sigma_{thr}","poisson_mean");
  //tpe->SetParameters(frp3->Parameter(0),frp3->Parameter(1),frp3->Parameter(2),TMath::Sq(frp3->Parameter(3)),TMath::Sq(frp3->Parameter(4)),frp3->Parameter(5),frp2->Parameter(6),0.01 );
  //for (int i = 0; i<7; i++){    tpe->FixParameter(i,tpe->GetParameter(i));  }    
  //tpe->FixParameter(7,0.01);
  //TFitResultPtr frp4 = h3->Fit("tpe","S");
  


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
