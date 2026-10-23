
Double_t EMGt(Double_t *x, Double_t *par){
   return par[0]*0.5/par[3]*exp( (2*par[1] + par[2]*par[2]/par[3] - 2.*x[0]) *0.5/par[3])
          *TMath::Erfc((par[1] +par[2]*par[2]/par[3] -x[0])/(sqrt(2.)*par[2]));
   // par[0] = γ, par[1] = μ, par[2] = σ, par[3] = τ, x[0] = x 
   // f(x) = γ/2τ * e^((2μ + σσ/τ - 2x))/2τ * erfc ( (μ + σσ/τ - x) / 2σ)
}
Double_t EMGtC(Double_t *x, Double_t *par){
   return EMGt(x,par) + par[4];
}
Double_t EMGkC(Double_t *x, Double_t *par){
   return par[4] + par[0]*0.5/par[2]/par[3]*exp((par[1]-x[0])/par[3]/par[2] + 0.5/par[3]/par[3])
          *TMath::Erfc(((par[1]-x[0])/par[2] + par[3] )/sqrt(2.));
   // par[0] = γ, par[1] = μ, par[2] = σ, par[3] = τ/σ = K, x[0] = x 
}


class TFit{
  public:
    inline static TF1 *emgtc;
    inline static double xmin, xmax;

    static TFitResultPtr ter(TH1D*& h1){
      xmin = h1->GetBinLowEdge(1);
      xmax = h1->GetBinLowEdge(h1->GetNbinsX() + 1);
      Double_t ymax = h1->GetBinContent(h1->GetMaximumBin());
      Int_t top = TMath::Nint(h1->GetXaxis()->GetBinLowEdge(h1->GetMaximumBin()));
      Int_t above1 = TMath::Nint(h1->GetXaxis()->GetBinLowEdge(h1->FindFirstBinAbove(0.5*ymax)));
      Int_t above2 = TMath::Nint(h1->GetXaxis()->GetBinLowEdge(h1->FindLastBinAbove(0.5*ymax)+1));
      Double_t rough_area = h1->GetEntries()*h1->GetXaxis()->GetBinWidth(1);
      // preliminary fitting
      TF1 *emgkc = new TF1("emgkc", EMGkC, -500, 500, 5);
      emgkc->SetParameters(rough_area, top, above2-above1, 1, h1->GetEntries()/h1->GetNbinsX());
      emgkc->SetParNames("N","#mu","#sigma","K=#tau/#sigma","C");
      emgkc->SetParLimits(0,rough_area*0.1,rough_area*2);
      emgkc->SetParLimits(1,xmin,xmax);
      emgkc->SetParLimits(2,0,xmax-xmin);
      emgkc->FixParameter(3,1);
      //emgkc->SetParLimits(2,0.,1e15);// This made fitting fail.
      //emgkc->SetParLimits(3,0.01,100);
      h1->Fit(emgkc,"Q","",h1->GetBinCenter(0),h1->GetBinCenter(h1->GetNbinsX()));
      // Fitting
      if (gROOT->FindObject("emgtc")){delete (TF1*) gROOT->FindObject("emgtc");}
      emgtc = new TF1("emgtc", EMGtC, xmin, xmax, 5);
      emgtc->SetNpx(1000);
      emgtc->SetParNames("#gamma","#mu","#sigma","#tau","C");
      emgtc->SetParameters(emgkc->GetParameter(0),emgkc->GetParameter(1),emgkc->GetParameter(2),emgkc->GetParameter(3)*emgkc->GetParameter(2),emgkc->GetParameter(4));
    
      delete emgkc;
    
      return h1->Fit(emgtc,"S","",h1->GetBinCenter(0),h1->GetBinCenter(h1->GetNbinsX()));//Reference: https://root.cern/manual/fitting/;
    }


  static void GetBound(double &lb, double &ub){
    double maximumx = emgtc->GetMaximumX(xmin, xmax);
    double bound_val = emgtc->GetParameter("C") + 1;
    lb = emgtc->GetX(bound_val, xmin, maximumx);
    ub = emgtc->GetX(bound_val, maximumx, xmax);
    return;
  }
};





// Below is trash
// Below is trash
// Below is trash
// Below is trash
// Below is trash




//#ifndef CHARA
//#define CHARA
//#include "ExtractCharacteristics.h"
////DistChara ExtractCharacteristics(TF1* func, TFitResultPtr frp)
//#endif



Double_t dfdx(Double_t *x, Double_t *par){
  Double_t f = EMGt(x,par);
  return 1/par[3] *( - f + par[0]*TMath::Gaus(x[0],par[1],par[2],kTRUE) ); 
}
Double_t EMGwidth(Double_t *y, Double_t *par){
  TF1* f = new TF1("f",EMGt,par[4],par[5],4);
  f->SetParameters(par[0],par[1],par[2],par[3]);
  double mpv       =   f->GetMaximumX();
  return f->GetX(y[0],mpv,par[5]) - f->GetX(y[0],par[4],mpv);
}

void error_calculation(TFitResultPtr &frp);

void error_calculation(TFitResultPtr &frp){
  TMatrixDSym cov  = frp->GetCovarianceMatrix();
  TMatrixDSym cor  = frp->GetCorrelationMatrix();
  double gamma     = frp->Parameter(0);
  double mu        = frp->Parameter(1);
  double sigma     = frp->Parameter(2);
  double tau       = frp->Parameter(3);
  double gamma_err = frp->ParError(0);
  double mu_err    = frp->ParError(1);
  double sigma_err = frp->ParError(2);
  double tau_err   = frp->ParError(3);
  double xmax =  mu +100;
  double xmin =  mu -100;
  TF1* emgt   = new TF1("emgt",EMGt,xmin,xmax,4);
  TF1* demgdx = new TF1("dEMGdx",dfdx,xmin,xmax,4);
  TF1* fw     = new TF1("fw",EMGwidth,xmin,xmax,6);
  emgt  ->SetParameters(gamma,mu,sigma,tau);
  demgdx->SetParameters(gamma,mu,sigma,tau);
  fw    ->SetParameters(gamma,mu,sigma,tau,xmin,xmax);
  double mpv      = emgt->GetMaximumX();
  double c_m      = - demgdx->GradientPar(1,&mpv)/demgdx->Derivative(mpv);
  double c_s      = - demgdx->GradientPar(2,&mpv)/demgdx->Derivative(mpv);
  double c_t      = - demgdx->GradientPar(3,&mpv)/demgdx->Derivative(mpv);
  double mpv_err  = TMath::Sqrt( TMath::Sq(c_m*mu_err) + TMath::Sq(c_s*sigma_err) + TMath::Sq(c_t*tau_err) + 2*c_m*c_s*cov[1][2] + 2*c_m*c_t*cov[1][3] + 2*c_s*c_t*cov[2][3] );
  double hm       = emgt->Eval(mpv)/2.;// Half Maximum
  double fwhm     = fw->Eval(hm);      // = emgt->GetX(emgt->Eval(mpv) /2.,mpv,xmax) - emgt->GetX(emgt->Eval(mpv) /2.,xmin,mpv);
  double dfwdy    = fw->Derivative(hm);// = 1/demgdx->Eval(emgt->GetX(hm,mpv,xmax)) - 1/demgdx->Eval(emgt->GetX(hm,xmin,mpv));
  double coef_m   = dfwdy*0.5*( emgt->GradientPar(1,&mpv) + demgdx->Eval(mpv)* c_m ) + fw->GradientPar(1,&hm);// = 0
  double coef_s   = dfwdy*0.5*( emgt->GradientPar(2,&mpv) + demgdx->Eval(mpv)* c_s ) + fw->GradientPar(2,&hm);
  double coef_t   = dfwdy*0.5*( emgt->GradientPar(3,&mpv) + demgdx->Eval(mpv)* c_t ) + fw->GradientPar(3,&hm);
  double fwhm_err = TMath::Sqrt( TMath::Sq(coef_s * sigma_err) + TMath::Sq(coef_t * tau_err) + 2.*coef_s*coef_t*cov[2][3] );

  std::cout << "mpv=" << mpv << "±" <<  mpv_err << "\n" << "fwhm =" << fwhm << "±" << fwhm_err << std::endl;
  
}
