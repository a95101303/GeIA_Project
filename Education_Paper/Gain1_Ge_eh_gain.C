#include <stdio.h>
#include <stdlib.h>
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLine.h"
#include "TGraphErrors.h"
#include "TGraph2D.h"
#include "TGraph.h"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include "TRandom.h"
#include "TRandom3.h"
#include <math.h>

#include "TCutG.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TAxis.h"

int Option=0; //0 is the electron and 1 is the hole
//print 'Option: '+str(Option)
//(2) Temperature is the one of the important parameters in this studies
double Temperature=77;
//print 'temperature: '+str(Temperature)
//(3) Electric field is another one of the important parameters in this studies also
double Electric_field=100000; //1000-100000(Lower-bound:15)
//#print 'Electric_field: '+str(Electric_field*100)
//#(4) Charge Constant
double CC = 1.6 * TMath::Power(10,-19);
//#(5) Boltzman COnstantKB = 8.617 * TMath::Power(10,-5)
double KB = 8.617 * TMath::Power(10,-5);
//(6) Electron mass with the unit of eV/(cm2/s2)

//==============================================================Avalanche================================================================
//==============================================================Mean_free_path=================================================================
double EMF(int Option, double T)//Effective_mass_factor, 0 n type 1 p type
{
    if(Option==0 and T>30) return (0.27/0.263);
    if(Option==0 and T<30) return (0.27/0.261);
    if(Option==1 and T>30) return (0.37/0.37);
    if(Option==1 and T<30) return (0.27/0.36);
}
//Temperature Dependence of Silicon Carrier E ective Masses with Application to FemtosecondRe ectivity Measurements(Paper)
double SV(int Option, double T)//Saturation_Velocity: 0 for electrons and 1 for holes
{
    double V_Sat_300K_e=0.7e7;double V_Sat_300K_h=0.63e7;
    double A_v_e=0.55;double A_v_h=0.61;
    
    if(Option==0) return ( V_Sat_300K_e/( 1 - A_v_e + A_v_e*(T/300.)) );
    if(Option==1) return ( V_Sat_300K_h/( 1 - A_v_h + A_v_h*(T/300.)) );
}
//==========================================
double Mobility(int Option, double Electric_field_1,double T)//Get the Mobility under the different electric fields and temperatures, 0 for electrons and 1 for holes
{
    //E_Saturate = 50 + (450/63)*(Temperature-4)
        //E_Saturate = 50 + (450/63)*(Temperature-4)
        double E_Saturate;double Drift_Velocity;double Mobility_0; double Mobility;double Mobility_Original;
        if(T==300)
        {
            if(Option==0)E_Saturate = 1.5e4;Mobility_0=(3.9e3);
            if(Option==1)E_Saturate = 4.0e4;Mobility_0=(1.9e3);
        }
        if(T==77)
        {
            if(Option==0)E_Saturate = 2.0e3;Mobility_0=(3.6e4);
            if(Option==1)E_Saturate = 2.0e3;Mobility_0=(4.2e4);
        }
        if(T==4)
        {
            if(Option==0)E_Saturate = 1.0e3;Mobility_0=(1.0e6);
            if(Option==1)E_Saturate = 1.0e3;Mobility_0=(1.0e6);
        }

        cout << "Mobility_0: " << Mobility_0 << endl;
        cout << "SV(T): " << SV(Option,T) << endl;

        //Mobility = Mobility_0 / (1 + (Mobility_0*Electric_field_1)/(SV(Option,T)) );
        cout << "Mobility: " << Mobility << endl;
        return Mobility_0;

}
double EFV(double Mobility, double E)//Emperical_formula_V
{
    return Mobility*E;
}
//==========================================1
double RT(double Mobility,double Effective_mass_factor)//Relaxation_time:vThe time that the electron will bump into other electrons for the first time
{
 // return (Mobility * 9.1 * Effective_mass_factor * TMath::Power(10,-31)/(1.6*TMath::Power(10,-19)*10000));
   return ((1./9.) * 1e-16 * Mobility * Effective_mass_factor * 0.51 * 1e6 * 1e-4);
}
//============================================================================================================
double MFP(int Option, double T)//Mean_free_path:The length that the particle will bump into another electrons for the first time.
{
    if(T==77 and Option==0)return 4.2e-5;//cm
    if(T==77 and Option==1)return 6.5e-5;//cm
    if(T==4 and Option==0)return 2.65e-4;//cm
    if(T==4 and Option==1)return 3.53e-4;//cm
     // return 1e-2*V*(Relaxation_time);//MFP=tau*V
    //return 6.5e-5;//cm

}
//============================================================================================================
double Ionization_rate(int Option, double ionization_energy,double E_x,double T)//0 n type 1 p type
{
    double MFP_E = MFP(0,T);
    double MFP_H = MFP(1,T);

    double A_s; double B_s; double B_n; double B_p; double Z_E;
    if(Option==0)A_s = (1/MFP_E);
    if(Option==1)A_s = (1/MFP_H);
    B_s = ionization_energy * A_s;
    //cout << "A_s: " << A_s << endl;
    //cout << "ionization_energy: " << ionization_energy << endl;
   // cout << "B_s: " << B_s << endl;
    B_n = (ionization_energy)*(1/(MFP_E));
    B_p = (ionization_energy)*(1/(MFP_H));
    Z_E = 1.0 + (B_n/E_x)*TMath::Power(2.718,-B_n/E_x)+ (B_p/E_x)*TMath::Power(2.718,-B_p/E_x);
    //cout << "(A_s/Z_E): " << (A_s/Z_E) << endl;
    //cout << "B_s: " << B_s << endl;cout << "Z_E: " << Z_E << endl;
    cout << "(B_n/E_x): " << (B_n/E_x) << endl;
    cout << "(B_p/E_x): " << (B_p/E_x) << endl;
return (A_s/Z_E)*TMath::Power(2.718,(-B_s/E_x));
}

double Ionization_rate_300K(int Option, double E_x)//0 n type 1 p type
{
    if(Option==0)return 3.8e6*TMath::Power(2.718,(-1.75e6/E_x));
    if(Option==1)return 2.25e7*TMath::Power(2.718,(-3.26e6/E_x));
}


//============================================================================================================
void Gain1_Ge_eh_gain()
{
double xarray[13];
double yarray_4K_h[13];
double yarray_4K_e[13];
double yarray_77K_h[13];
double yarray_77K_e[13];
double xarrayerror[13];
double yarrayerror[13];

double MFP_Ele_4K;double MFP_Hol_4K;
double MFP_Ele_77K;double MFP_Hol_77K;



for(int i=1; i<13; i++)
{
    int element=i-1;
    xarrayerror[element] = 0;
    yarrayerror[element] = 0;
    xarray[element]      = (i)*(1e4);
    
    cout << "KKKK!!! " << endl;
    yarray_4K_e[element] = (Ionization_rate(0,3,xarray[element],4));
    //yarray_4K_h[element] = TMath::Log10(Ionization_rate(1,0.72,xarray[element],4));
    yarray_4K_h[element] = (Ionization_rate(1,3,xarray[element],4));
    //cout << "KKKKK!!! " << endl;
    cout << "xarray[element]: " << xarray[element] << endl;
    cout << "yarray_4K_e[element]: " << yarray_4K_e[element] << endl;
    cout << "=====================" << endl;
    yarray_77K_e[element] = (Ionization_rate(0,3,xarray[element],77));
    yarray_77K_h[element] = (Ionization_rate(1,3,xarray[element],77));
    cout << "yarray_77K_e[element]: " << yarray_77K_e[element] << endl;

}
 
    
    
    TGraphErrors *gr_4K_e  = new TGraphErrors(12,xarray,yarray_4K_e,xarrayerror,yarrayerror);
    TGraphErrors *gr_4K_h  = new TGraphErrors(12,xarray,yarray_4K_h,xarrayerror,yarrayerror);
    TGraphErrors *gr_77K_e = new TGraphErrors(12,xarray,yarray_77K_e,xarrayerror,yarrayerror);
    TGraphErrors *gr_77K_h = new TGraphErrors(12,xarray,yarray_77K_h,xarrayerror,yarrayerror);

 
    TCanvas * c = new TCanvas("c", "c", 0,0,1000,1000);
    gStyle->SetOptStat(0);
    gStyle->SetTitleSize(0.035,"XY");
    gStyle->SetTitleFont(62,"XY");
    gStyle->SetLegendFont(62);




 //Ionization Rate
gStyle->SetOptFit();

TF1 *gr_4K_e_fa = new TF1("gr_4K_e","Ionization_rate(0,0.72,x,4)",1e5,1e6);
TF1 *gr_4K_h_fa = new TF1("gr_4K_h","Ionization_rate(1,0.72,x,4)",1e5,1e6);
TF1 *gr_77K_e_fa = new TF1("gr_77K_e","Ionization_rate(0,0.72,x,77)",1e5,1e6);
TF1 *gr_77K_h_fa = new TF1("gr_77K_h","Ionization_rate(1,0.72,x,77)",1e5,1e6);

gr_4K_h->SetLineColor(0);
    
gr_4K_e_fa->SetLineColor(2);
gr_4K_e_fa->SetLineWidth(3);
gr_4K_e_fa->SetLineStyle(1);
    
gr_4K_h_fa->SetLineColor(2);
gr_4K_h_fa->SetLineWidth(3);
gr_4K_h_fa->SetLineStyle(10);
    
gr_77K_e_fa->SetLineColor(3);
gr_77K_e_fa->SetLineWidth(3);
gr_77K_e_fa->SetLineStyle(1);

gr_77K_h_fa->SetLineColor(3);
gr_77K_h_fa->SetLineWidth(3);
gr_77K_h_fa->SetLineStyle(10);

gr_4K_h->SetTitle(";E(V/cm) ;Ionization Rate(/cm) ");
gr_4K_h->GetHistogram()->SetMaximum(1e5);
gr_4K_h->GetHistogram()->SetMinimum(0);
gr_4K_h->GetXaxis()->SetRangeUser(1e5,1e7);
gr_4K_h_fa->GetXaxis()->CenterTitle();
gr_4K_h_fa->GetYaxis()->SetTitle("Ionization Rate(/cm)");
gr_4K_h_fa->GetXaxis()->SetTitle("E(V/cm)");
gr_4K_h_fa->GetYaxis()->CenterTitle();
//================================
gr_4K_h_fa->GetYaxis()->SetTitleSize(0.02);
gr_4K_h_fa->GetXaxis()->SetTitleSize(0.02);
gr_4K_h_fa->GetXaxis()->SetLabelSize(0.02);
gr_4K_h_fa->GetYaxis()->SetLabelSize(0.02);
gr_4K_h_fa->GetXaxis()->SetLabelFont(22);
gr_4K_h_fa->GetYaxis()->SetLabelFont(22);
gr_4K_h_fa->GetXaxis()->SetTitleColor(1);
gr_4K_h_fa->GetYaxis()->SetTitleColor(1);

TLegend *leg1= new TLegend(0.1,0.6,0.6,0.9);
leg1->SetFillColor(0);
leg1->SetFillStyle(0);
leg1->SetTextSize(0.05);
leg1->SetBorderSize(0);
leg1->SetTextFont(22);
leg1->AddEntry("","GeIA group","");
leg1->AddEntry(gr_4K_e_fa,"Ge-Electron(4K)","lp");
leg1->AddEntry(gr_4K_h_fa,"Ge-Hole(4K)","lp");
leg1->AddEntry(gr_77K_e_fa,"Ge-Electron(77K)","lp");
leg1->AddEntry(gr_77K_h_fa,"Ge-Hole(77K)","lp");

gr_4K_h->Draw("AL");
gr_4K_h_fa->Draw("Lsame");
gr_4K_e_fa->Draw("Lsame");
gr_77K_h_fa->Draw("Lsame");
gr_77K_e_fa->Draw("Lsame");

    
c->SetLogx();
//c.SetLogy()
c->Draw();
leg1->Draw();
c->Print("IR_77K_4K.pdf");
 
}

















































