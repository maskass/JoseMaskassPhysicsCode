#include <TMath.h>

using namespace TMath;

Double_t CB(Double_t *x,Double_t *par) {
  Double_t t=(x[0]-par[1])/par[2];
  Double_t absAlpha=fabs((Double_t)par[3]);

  if (t > -par[3]) {
    return par[0]*exp(-0.5*t*t);
  }
  else {
    Double_t a =  TMath::Power(par[4]/absAlpha,par[4])*exp(-0.5*absAlpha*absAlpha);
    Double_t b= (par[4]/absAlpha)-absAlpha;

    return par[0]*a*TMath::Power(b-t, -par[4]);
  }
}

Double_t Matulevich(Double_t *x,Double_t *par) {
  Double_t G=exp((-4*TMath::Log2((x[0]-par[1])*(x[0]-par[1])))/(par[2]*par[2]));

  if (x[0] >= par[1]) {
    return par[0]*G;
  }
  else {
    return par[0]*(G+exp((x[0]-par[1])/par[3])*(1-G));
  }
}

Double_t Novosibirsk(Double_t *x,Double_t *par) {
  Double_t Lambda=SinH(par[3]*Sqrt(Log(4)))/(par[2]*par[3]*Sqrt(Log(4)));
 
  return par[0]*Exp(Power(par[3],2)-0.5*Power(Log(1+Lambda*par[3]*Power(x[0]-par[1],2)),2)/Power(par[3],2));
 
}

Double_t Novo(Double_t *x,Double_t *par) {
double qa=0,qb=0,qc=0,qx=0,qy=0;

if(TMath::Abs(par[3]) < 1.e-7) 
     qc = 0.5*TMath::Power(((x[0]-par[1])/par[2]),2);
     else {
  qa = par[3]*sqrt(log(4.));
  qb = sinh(qa)/qa;
  qx = (x[0]-par[1])/par[2]*qb;
  qy = 1.+par[3]*qx;
  
  //---- Cutting curve from right side

  if( qy > 1.E-7) 
    qc = 0.5*(TMath::Power((log(qy)/par[3]),2) + par[3]*par[3]);
  else
    qc = 15.0;
     }

//---- Normalize the result

return par[0]*exp(-qc);

}
Double_t AsymGaus(Double_t *x,Double_t *par) {

  Double_t C=par[0];
  Double_t mean=par[1];
  Double_t sigmaL=par[2];
  Double_t sigmaR=par[3];

  Double_t arg=x[0]-mean;
  Double_t coef;
  
  if (arg < 0.0) 
    coef = -0.5/(sigmaL*sigmaL);
  else
    coef = -0.5/(sigmaR*sigmaR);

  return C*TMath::Exp(coef*arg*arg);
}
