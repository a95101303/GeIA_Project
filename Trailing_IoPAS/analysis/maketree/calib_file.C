
Double_t linear(Double_t *x,Double_t *par)
{
  return par[0]+par[1]*x[0];
} 

void calib_file()
{


//GET FIT PARAMETER
  Double_t mean[2]={0};
  Double_t sigma[2]={0};
  
  mean[0]=0.1;  
  sigma[0]=10.0;
  mean[1]=3.09965e+04;  
  sigma[1]=5.91525e+01;
 
  Double_t real_energy_of_xray[4]={0};
   real_energy_of_xray[0]=0.0;
   real_energy_of_xray[1]=59.5;
   
  for(int ii =0; ii<2;ii++)
    printf("%f  %f   %f \n",real_energy_of_xray[ii],mean[ii],sigma[ii]);  
  
  TF1 *linefitfun=new TF1("linefitfun",linear,0,31000,2);
 
  linefitfun->SetLineColor(2);
  linefitfun->FixParameter(0,59.5/3.09965e+05);  
  linefitfun->FixParameter(1,3.15988322e-02);
  TCanvas *lucky_cali = new TCanvas("lucky_cali");
  gStyle->SetOptFit(1);
  TGraphErrors *cali = new TGraphErrors(2,mean,real_energy_of_xray,sigma);
  cali->SetMarkerColor(2);
  cali->SetMarkerStyle(21);
  //cali->Fit("linefitfun","","",0,30000);  
  cali->Draw("ap"); 
  
  TF1 *fitfun=new TF1("fitfun","1.9195e-03*x + 3.15988322e-02",0,31000);
  
  fitfun->Draw("same");
  //printf("%.8e\t%.8e\n",linefitfun->GetParameter(0),linefitfun->GetParameter(1));
 

}
