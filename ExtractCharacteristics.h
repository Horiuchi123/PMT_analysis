// If you encounter the error below, try ```sudo apt-get install libgsl-dev```.
// Warning in <ROOT::Math::IntegratorOneDim::CreateIntegrator>: Error loading one dimensional GSL integrator - use Gauss integrator
#include "declare.h"

void ExtractCharacteristics(TF1* func, TFitResultPtr frp, Double_t chara[n_chara][2]){
  double xmin,xmax;  func->GetRange(xmin,xmax);
  Int_t Npar = func->GetNpar();
  const std::vector<double> Parameters = frp->Parameters();
  const double* ParErrors = frp->GetErrors();
  TMatrixDSym cov = frp->GetCovarianceMatrix();

  
  TF1* xfunc  = new TF1("xfunc", [func](double *x, double *par){ return x[0]*func->EvalPar(x,par); },      xmin, xmax, Npar);
  TF1* xxfunc = new TF1("xxfunc",[func](double *x, double *par){ return x[0]*x[0]*func->EvalPar(x,par); }, xmin, xmax, Npar);


  int grp1[] = {MEAN,VARI,RESO,MPV,FWHM};
  std::vector<TF1*>                F(n_chara);
  std::vector<TVectorT<double>> dfdp(n_chara);

  //Functions below is a constant function: ∂f/∂x = 0
  F[MEAN] = new TF1("Fmean",[xfunc,func,xmin,xmax](double *x, double *par){ xfunc->SetParameters(par);func->SetParameters(par); return xfunc->Integral(xmin,xmax)/func->Integral(xmin,xmax);},0,1,Npar);
  F[VARI] = new TF1("Fvari",[xxfunc,func,F,xmin,xmax](double *x, double *par){ xxfunc->SetParameters(par);func->SetParameters(par);F[MEAN]->SetParameters(par); return xxfunc->Integral(xmin,xmax)/func->Integral(xmin,xmax) - TMath::Sq(F[MEAN]->Eval(0));},0,1,Npar);
  F[RESO] = new TF1("Freso",[F,xmin,xmax](double *x, double *par){return TMath::Sqrt(F[VARI]->EvalPar(x,par))/F[MEAN]->EvalPar(x,par);},0,1,Npar);
  F[MPV]  = new TF1("Fmpv",[func,xmin,xmax](double *x, double *par){ func->SetParameters(par); return func->GetMaximumX(xmin,xmax);},0,1,Npar);
  F[FWHM] = new TF1("Ffwhm",[func,xmin,xmax](double *x, double *par){ func->SetParameters(par); double max_px = func->GetMaximumX(xmin,xmax); double max_py = func->GetMaximum(xmin,xmax); return func->GetX(0.5*max_py,max_px,xmax) - func->GetX(0.5*max_py,xmin,max_px) ;},0,1,Npar);

  double nonsense = 0;
  for (int i : grp1 ){
    F[i]->SetParameters(Parameters.data());
    F[i]->SetParErrors(ParErrors);
    dfdp[i].ResizeTo(Npar);
    F[i]->GradientPar(&nonsense,dfdp[i].GetMatrixArray());
    chara[i][ERR] = TMath::Sqrt(cov.Similarity(dfdp[i]));
    chara[i][VAL] = F[i]->Eval(0);
  }


  chara[STDDEV][ERR] = 0.5*chara[VARI][ERR]/TMath::Sqrt(chara[VARI][VAL]);
  chara[STDDEV][VAL] = TMath::Sqrt(chara[VARI][VAL]);


  for (int i: grp1){
    delete F[i];
  }
  delete xfunc; delete xxfunc;

}
