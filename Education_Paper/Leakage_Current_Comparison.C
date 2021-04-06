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

double Coulomb_Constant = 1.6*1e-19;
double Permittivity = 8.8*1e-12;// F/m = C/V.m
double Permittivity_cm = 8.8*1e-14;// F/cm = C/V.cm

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
double Electron_mass = (0.51*TMath::Power(10,6)) / (9*TMath::Power(10,20));

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
double Mobility(double Electric_field_1,double T)//Get the Mobility under the different electric fields and temperatures
{
    //E_Saturate = 50 + (450/63)*(Temperature-4)
        //E_Saturate = 50 + (450/63)*(Temperature-4)
        double E_Saturate;double Drift_Velocity;double Mobility_0; double Mobility;double Mobility_Original;
        if(T==77){E_Saturate = 2000;Mobility_0=(3.6e4);}
        if(T==4) {E_Saturate = 1000;Mobility_0=(1e6);}
        //cout << "SV(T): " << SV(1,T) << endl;

        Mobility = Mobility_0 / (1 + (Mobility_0*Electric_field_1)/(SV(1,T)) );
        //cout << "Mobility: " << Mobility << endl;
        return Mobility;

}
//==========================================1
double Relaxation_time(double Mobility,double Effective_mass_factor)//The time that the electron will bump into other electrons for the first time
{
 // return (Mobility * 9.1 * Effective_mass_factor * TMath::Power(10,-31)/(1.6*TMath::Power(10,-19)*10000));
   return ((1./9.) * 1e-16 * Mobility * Effective_mass_factor * 0.51 * 1e6 * 1e-4);
}
//============================================================================================================
double Mean_free_path(double Relaxation_time)//The length that the particle will bump into another electrons for the first time.
{
    //cout << "Relaxation_time: " << Relaxation_time << endl;
    //cout << "MEAN FREE PATH: " << (TMath::Power(10,7))*(Relaxation_time) << endl;
 //return (TMath::Power(10,7))*(Relaxation_time);//MFP=tau*V
    return SV(1,4)*(Relaxation_time);//cm

}
//============================================================================================================
double Ionization_rate(int Option, double ionization_energy,double E_x,double T)//0 n type 1 p type
{
    double MFP_E = Mean_free_path(Relaxation_time(Mobility(E_x,T),(0.12)));
    double MFP_H = Mean_free_path(Relaxation_time(Mobility(E_x,T),(0.21)));
    double A_s; double B_s; double B_n; double B_p; double Z_E;
    if(Option==0)A_s = (1/MFP_E);
    if(Option==1)A_s = (1/MFP_H);
    B_s = ionization_energy * A_s;
   // cout << "A_s: " << A_s << endl;
   // cout << "ionization_energy: " << ionization_energy << endl;
   // cout << "B_s: " << B_s << endl;
    B_n = (ionization_energy)*(1/(MFP_E));
    B_p = (ionization_energy)*(1/(MFP_H));
    Z_E = 1.0 + (B_n/E_x)*TMath::Power(2.718,-B_n/E_x)+ (B_p/E_x)*TMath::Power(2.718,-B_p/E_x);
    //cout << "(A_s/Z_E): " << (A_s/Z_E) << endl;
    //cout << "B_s: " << B_s << endl;cout << "Z_E: " << Z_E << endl;
   // cout << "(B_n/E_x): " << (B_n/E_x) << endl;
   // cout << "(B_p/E_x): " << (B_p/E_x) << endl;
return (A_s/Z_E)*TMath::Power(2.718,(-B_s/E_x));
}

double BC(int Option)//Boltzmann_constant(1)
{
    double Boltzmann_constant_Unit;
    if(Option==0)Boltzmann_constant_Unit=1.3*1e-23;//(J.K^-1)
    if(Option==1)Boltzmann_constant_Unit=8.6*1e-5; //(eV.K^-1)
    if(Option==2)Boltzmann_constant_Unit=1.3*1e-16;//(erg.K^-1)
    return Boltzmann_constant_Unit;
}
double Bulk_Leakage_Current(double Impurity_Level, double T)
{//Impurity_Level(count/cm^3),T = Temperature(K)
    double Factor = TMath::Exp(-0.68/(2*BC(1)*T));
    return Impurity_Level*Factor*Coulomb_Constant;
}
double Contact_Leakage_Current(double T)
{//W = 0.28eV(Barrier Potential)
    double Initial = (130*130*TMath::Exp(-0.28/(BC(1)*130)));
    double Comparison = (T*T*TMath::Exp(-0.28/(BC(1)*T)));
    double Ratio = Initial/Comparison;

    return (5*1e-7/(Ratio));
}
double Conductivity_for_Ge(double T)
{
    cout << "T: " << T << endl;
    double TM_factor       = 1/(TMath::Power(T,0.25));
    cout << "TM_factor: " << TM_factor  << endl;
    double Ln_Conductivity = -234.2*TM_factor + 51.5;
    cout << "Ln_Conductivity: " << Ln_Conductivity << endl;
    double Conductivity    = TMath::Power(2.71828,Ln_Conductivity);
    cout << "Conductivity: " << Conductivity << endl;
    cout << "Resistivity: " << 1/(Conductivity) << endl;

    return Conductivity;
}
double Surface_Leakage_Current(double T)
{//J = conductivity * E
    double Applied_electric_field = 5*1e-11/(Conductivity_for_Ge(77));
    cout << "Applied_electric_field: " << Applied_electric_field << endl;
    
    return Conductivity_for_Ge(T)*Applied_electric_field;
}
//
double Cross_Section_electronics(double Radius)
{
    return TMath::Pi()*Radius*Radius;
}
double Ohm_Cu(double T, double Length, double Radius)//m^2,m
{
    return 1.68*1e-8*(1+(0.00404)*(T-300))*(Length/Cross_Section_electronics(Radius));
}
double Johnson_Nyquist_noise_Current(double T, double Bandwidth, double Length, double Radius)
{
    //cout << "Ohm_Cu(T,Length,Radius): " << Ohm_Cu(T,Length,Radius) << endl;
    cout << "sqrt(4*BC(1)*T*Bandwidth*Coulomb_Constant/50): " << sqrt(4*BC(1)*T*Bandwidth*Coulomb_Constant/Ohm_Cu(T,Length,Radius)) << endl;
    return sqrt(4*BC(1)*T*Bandwidth*Coulomb_Constant/Ohm_Cu(T,Length,Radius));
}
double Response_time_to_eV(double Current, double Voltage, double Time_Duration)
{
    return Current*Voltage;
}
double Depletion_Voltage(double Depletion_Length, double Impurity_Level, double T)
{
    //return sqrt(   (Depletion_Length*Depletion_Length*Coulomb_Constant*Impurity_Level*TMath::Exp(-0.01/(BC(1)*T)))/(2*Permittivity_cm*16)  );
    return (   (Depletion_Length*Depletion_Length*Coulomb_Constant*Impurity_Level*TMath::Exp(-0.01/(BC(1)*T)))/(2*Permittivity_cm*16)  );
}

double Ex(double V, double r, double Impurity_Level, double R2, double R1)//R2 outer, R1 inner
{
    
    double Ex1      = (Impurity_Level*Coulomb_Constant*r)/(2*Permittivity*16);//
    double Ex2_Up   = V + (Impurity_Level*Coulomb_Constant/(4*Permittivity*16))*(R2*R2-R1*R1);
    double Ex2_Down = r*log(R2/R1);
    return (Ex1) - (Ex2_Up/Ex2_Down);
}
double Ex_PPC1(double r)//T=77K (E>1e4V/cm) as (r<0.1cm), from Prof. Dongming's version
{
    return -16019.51 + 66380.48*TMath::Exp(-9.636358*r);
}
double Ex_PPC2(double r)//T=77K
{
    return 1380.7+1.008e5*TMath::Exp(-34*r);
}

double IRF(int T, double E)//Ionization_Rate_fitting
{
    //Hole for 4K =>y = 1.545584e-14*E^3.532872(/cm)
    //Hole for 77K=>y = 7.549454e-21*E^4.624588(/cm)
    if(T==77)return 7.549454e-21*TMath::Power(E,4.624588);
    if(T==4) return 1.545584e-14*TMath::Power(E,3.532872);
}

double Ionization_rate_300K(int Option, double E_x)//0 n type 1 p type
{
    if(Option==0)return 3.8e6*TMath::Power(2.718,(-1.75e6/E_x));
    if(Option==1)return 2.25e7*TMath::Power(2.718,(-3.26e6/E_x));
}

double Breakdown_ds(double Gain, double T)
{
    return 0.6*1e6*TMath::Power((Gain*T*TMath::Power(2,15))/(16*300),0.33);//(V/cm)
}
//
void Leakage_Current_Comparison()
{
    
    TCanvas * c1 = new TCanvas("c", "c", 0,0,1000,1000);
    gStyle->SetOptStat(0);
    gStyle->SetTitleSize(0.035,"XY");
    gStyle->SetTitleFont(62,"XY");
    gStyle->SetLegendFont(62);
     /*
    const int Number_of_temperature=9;
    double Temperature_Array[Number_of_temperature]={4,10,20,30,40,50,60,70,77};//K
    */
    /*
    double r_array[20];double r_mean_array[20];double r_difference_array[20];
    double Ex_array[80];
    
    
    for(int kkk=1; kkk<20; kkk++)
    {
        if(kkk<10)r_array[kkk-1] = 1e-3*(kkk);
        if(kkk>=10)r_array[kkk-1] = 1e-2*(kkk-9);
        cout << "r_array[kkk]: " << r_array[kkk-1] << endl;
    }
    
    double Voltage_for_Gain[80];
    double Gain_77K_voltage[80];double Gain_4K_voltage[80];double Breakdown_voltage[80];

    for(int jjj=0; jjj<80; jjj++)
    {
            Voltage_for_Gain[jjj] = 10000*(jjj+1);
            //Ex_array[jjj]= Ex(-Voltage_for_Gain[jjj],1e-1,1e10*TMath::Exp(-0.01/(BC(1)*77)),7.5,0.039);//Real case from Akash
            double Strip_Ex_V1 = Ex(-3500,0.05,1e10*TMath::Exp(-0.01/(BC(1)*300)),7.5,0.039);//Real case from Akash
            double Strip_Ex_V2 = Ex(-Voltage_for_Gain[jjj],0.05,1e10*TMath::Exp(-0.01/(BC(1)*300)),7.5,0.039);
            Ex_array[jjj] = (Ex_PPC2(0.004)*Strip_Ex_V2)/(Strip_Ex_V1);
            cout << "Voltage_for_Gain[jjj]: " << Voltage_for_Gain[jjj] << endl;
            cout << "Ex_array[jjj]: " << Ex_array[jjj] << endl;
    }

    
    TGraph *Ex_77= new TGraph(80,Voltage_for_Gain,Ex_array);
    Ex_77->SetLineColor(2);
    Ex_77->SetTitle("");
    Ex_77->GetXaxis()->SetTitle("Voltage(V)");
    Ex_77->GetYaxis()->SetTitle("Max E(V/cm)");
    Ex_77->GetXaxis()->SetRangeUser(0,1e6);
    Ex_77->GetYaxis()->SetRangeUser(0,100);
    Ex_77->SetMarkerColor(2);
    Ex_77->SetMarkerStyle(8);

    TLine *l  =new TLine(0,Breakdown_ds(1,300),1e6,Breakdown_ds(1,300));
    l->SetLineColor(3); l->SetLineWidth(3);
    TLine *l1  =new TLine(0,Breakdown_ds(1,77),1e6,Breakdown_ds(1,77));
    l1->SetLineColor(4); l1->SetLineWidth(3);
    TLine *l2  =new TLine(0,Breakdown_ds(1,4),1e6,Breakdown_ds(1,4));
    l2->SetLineColor(5); l2->SetLineWidth(3);
    
    
    TLegend *leg= new TLegend(0.7,0.1,0.9,0.3);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.06);
    leg->SetBorderSize(0);
    leg->SetTextFont(20);

    leg->AddEntry(l,"300K","L");
    leg->AddEntry(l1,"77K","L");
    leg->AddEntry(l2,"4K","L");


    Ex_77->Draw("ALP");
    l->Draw("same");l1->Draw("same");l2->Draw("same");
    
    leg->Draw();
    c1->SetLogy();
    c1->SetLogx();
    c1->Print("Ex_PPC_Max_check.pdf");
    */
    
    
    double r_array[20];double r_mean_array[20];double r_difference_array[20];
    double Ex_array_77K[20];double Ex_array_4K[20];
    
    for(int kkk=1; kkk<20; kkk++)
    {
        if(kkk<10)r_array[kkk-1] = 1e-3*(kkk);
        if(kkk>=10)r_array[kkk-1] = 1e-2*(kkk-9);
        cout << "r_array[kkk]: " << r_array[kkk-1] << endl;
    }
    
    double Voltage_for_Gain[200];
    double Gain_77K_voltage[200];double Gain_4K_voltage[200];

    for(int jjj=7; jjj<9; jjj++)
    {
        Voltage_for_Gain[jjj] = 500*(jjj+1);
        for(int kkk=0; kkk<18; kkk++)
        {
            r_mean_array[kkk] = (r_array[kkk]+r_array[kkk+1])/2;
            r_difference_array[kkk] = (r_array[kkk+1] - r_array[kkk]);

            double Strip_Ex_V1 = Ex(-4000,r_mean_array[kkk],1e10*TMath::Exp(-0.01/(BC(1)*77)),1,0.039);//Real case from Akash
            double Strip_Ex_V2 = Ex(-Voltage_for_Gain[jjj],r_mean_array[kkk],1e10*TMath::Exp(-0.01/(BC(1)*77)),1,0.039);
            Ex_array_77K[kkk] = (Ex_PPC2(r_mean_array[kkk])*Strip_Ex_V2)/(Strip_Ex_V1);
            
            cout << "r_mean_array[kkk]: " << r_mean_array[kkk] << endl;
            cout << "Ex_array_77K[kkk]: " << Ex_array_77K[kkk] << endl;
                        
        }
        
        double Gain_77K=0; double Gain_4K=0;
        for(int kkk=0; kkk<18 ; kkk++)
        {
            if(r_mean_array[kkk]>=0.004){
                if(Ex_array_77K[kkk]>9e4){
                cout << "=========================================================" << endl;
                cout << "r_mean_array[kkk]: " << r_mean_array[kkk] << endl;
                cout << "r_difference_array[kkk]: " << r_difference_array[kkk] << endl;


            Gain_77K = Gain_77K + r_difference_array[kkk]*Ionization_rate_300K(1,Ex_array_77K[kkk]);
            Gain_4K  = Gain_4K  + r_difference_array[kkk]*Ionization_rate_300K(1,Ex_array_77K[kkk]);
                cout << "Gain_77K: " << Gain_77K << endl;
                cout << "Gain_4K: " << Gain_4K << endl;
                cout << "=========================================================" << endl;
                }
            }
        }
        //Gain_77K = Gain_77K + r_mean_array[0]*IRF(77,Ex_array_77K[0]);
        //Gain_4K = Gain_4K + r_mean_array[0]*IRF(4,Ex_array_77K[0]);

        Gain_77K_voltage[jjj] = (1+Gain_77K);
        Gain_4K_voltage[jjj]  = (1+Gain_4K);
        
        cout << "Voltage_for_Gain[jjj]: " << Voltage_for_Gain[jjj] << endl;
        cout << "Gain_77K: " << Gain_77K_voltage[jjj] << endl;
        cout << "Gain_4K: "  << Gain_4K_voltage[jjj]  << endl;
    }
    TGraph *Gain_77K_Comparison= new TGraph(200,Voltage_for_Gain,Gain_77K_voltage);
    Gain_77K_Comparison->SetLineColor(2);
    Gain_77K_Comparison->SetTitle("");
    Gain_77K_Comparison->GetXaxis()->SetTitle("Voltage(V)");
    Gain_77K_Comparison->GetYaxis()->SetTitle("Gain");
    Gain_77K_Comparison->GetXaxis()->SetRangeUser(1,1e5);
    Gain_77K_Comparison->GetYaxis()->SetRangeUser(1,1e4);
    Gain_77K_Comparison->SetMarkerColor(2);
    Gain_77K_Comparison->SetMarkerStyle(8);

    TGraph *Gain_4K_Comparison= new TGraph(200,Voltage_for_Gain,Gain_4K_voltage);
    Gain_4K_Comparison->SetMarkerColor(3);
    Gain_4K_Comparison->SetMarkerStyle(8);

    TLegend *leg= new TLegend(0.1,0.6,0.4,0.9);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.06);
    leg->SetBorderSize(0);
    leg->SetTextFont(20);

    //
    TLine *l=new TLine(0,1000,6e4,1000);
    TLine *l_77=new TLine(5.2e4,0,5.2e4,1000);
    TLine *l_4=new TLine(4e4,0,4e4,1000);

    TLine *l1=new TLine(0,100,6e4,100);
    TLine *l1_77=new TLine(3.2e4,0,3.2e4,100);
    TLine *l1_4=new TLine(2e4,0,2e4,100);

    TLine *l2=new TLine(0,10,6e4,10);
    TLine *l2_77=new TLine(1.9e4,0,1.9e4,10);
    TLine *l2_4=new TLine(1e4,0,1e4,10);
     //
    TLine *l3=new TLine(0,1,6e4,1);
    
    cout << "4K: " << (4e4/6e4) << "," << (2e4/6e4) << "," << (1e4/6e4) << endl;
    cout << "77K: " << (5.2e4/2e5) << "," << (3.2e4/2e5) << "," << (1.9e4/2e5) << endl;
    leg->AddEntry(Gain_77K_Comparison,"77K","lp");
    leg->AddEntry(Gain_4K_Comparison,"4K","lp");

    Gain_77K_Comparison->Draw("ALP");
    Gain_4K_Comparison->Draw("LPsame");
    //
    l->Draw("same");l_4->Draw("same");l_77->Draw("same");
    l1->Draw("same");l1_4->Draw("same");l1_77->Draw("same");
    l2->Draw("same");l2_4->Draw("same");l2_77->Draw("same");
    //
    l3->Draw("same");
    leg->Draw();
    c1->SetLogy();
    c1->SetLogx();
    c1->Print("Gain_PCC_2.pdf");
    
    
    //Depletion_Voltage
    /*
    double Depletion_cm[Number_of_temperature];
    for(int kkk=0; kkk<9; kkk++)
    {
        Depletion_cm[kkk] = Depletion_Voltage(2.37,1e10,Temperature_Array[kkk]);
    }

    
    TGraph *Depletion_Length_Calculated = new TGraph(Number_of_temperature,Temperature_Array,Depletion_cm);
    Depletion_Length_Calculated->SetLineColor(2);
    Depletion_Length_Calculated->SetTitle("");
    Depletion_Length_Calculated->GetXaxis()->SetTitle("T(K)");
    Depletion_Length_Calculated->GetYaxis()->SetTitle("Voltage(V)");
    Depletion_Length_Calculated->GetXaxis()->SetRangeUser(0,80);
    Depletion_Length_Calculated->GetYaxis()->SetRangeUser(1e-10,1e3);
    Depletion_Length_Calculated->SetMarkerColor(2);
    Depletion_Length_Calculated->SetMarkerStyle(8);

    TLegend *leg= new TLegend(0.6,0.1,0.9,0.4);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.06);
    leg->SetBorderSize(0);
    leg->SetTextFont(20);

    leg->AddEntry(Depletion_Length_Calculated,"","");

    Depletion_Length_Calculated->Draw("ALP");
    //leg->Draw();
    c1->SetLogy();
    c1->Print("Depletion_Length.pdf");
    */

    /*
const int Number_of_temperature=21;
double Temperature_Array[Number_of_temperature]={4,10,20,30,40,50,60,70,77,80,90,100,110,120,130,140,150,160,170,180,300};//K
double BLK[Number_of_temperature];double CLK[Number_of_temperature];double SLK[Number_of_temperature];double TLK[Number_of_temperature];
double Number_of_electron[Number_of_temperature];double Number_of_electron_Log[Number_of_temperature];double Number_of_electron_G[Number_of_temperature];
    for(int kkk=0; kkk<21; kkk++)
    {
        double BLK_constant=Bulk_Leakage_Current(1e23,Temperature_Array[kkk]);
        double CLK_constant=Contact_Leakage_Current(Temperature_Array[kkk]);
        double SLK_constant=Surface_Leakage_Current(Temperature_Array[kkk]);

        
        BLK[kkk] = TMath::Log10(BLK_constant);
        CLK[kkk] = TMath::Log10(CLK_constant);
        SLK[kkk] = TMath::Log10(SLK_constant);
        TLK[kkk] = TMath::Log10(BLK_constant+CLK_constant+SLK_constant);
         
        Number_of_electron[kkk]       = (3*3*sqrt((BLK_constant+CLK_constant+SLK_constant)*1e-6/(1.6e-19)));
        Number_of_electron_Log[kkk]   = TMath::Log10(Number_of_electron[kkk]);
        Number_of_electron_G[kkk]     = TMath::Log10(Number_of_electron[kkk]/1e3);
        
        cout << "Temperature_Array[kkk]: " << Temperature_Array[kkk] << endl;
        cout << "Number_of_electron[kkk]: " << Number_of_electron[kkk] << endl;
        
        if(BLK[kkk]<-100) BLK[kkk]=-300;if(CLK[kkk]<-100) CLK[kkk]=-300;
        
        //cout << "BLK[kkk]: " << BLK[kkk]  << endl;
        //cout << "CLK[kkk]: " << CLK[kkk]  << endl;
        //cout << "SLK[kkk]: " << SLK[kkk]  << endl;
        
    }
    
    //double A_limited = TMath::Log10(3*sqrt((1.6e-19)*1e-6/(1.6e-19)));
    double A_limited = TMath::Log10(68./100.);

    TLine *l3=new TLine(77,-20,77,10);
    TLine *l4=new TLine(4,-20,4,10);
    TLine *l5=new TLine(0,A_limited,300,A_limited);

    cout << "A_Limited: " << A_limited << endl;
    TGraph *Electron_number = new TGraph(Number_of_temperature,Temperature_Array,Number_of_electron_Log);
    Electron_number->SetLineColor(1);
    Electron_number->SetTitle("");
    Electron_number->GetXaxis()->SetTitle("T(K)");
    Electron_number->GetXaxis()->CenterTitle();
    Electron_number->GetYaxis()->SetTitle("Log(Threshold)[eV]");
    Electron_number->GetYaxis()->CenterTitle();
    Electron_number->GetXaxis()->SetRangeUser(0,300);
    Electron_number->GetYaxis()->SetRangeUser(-20,10);
    Electron_number->SetMarkerColor(1);
    //BLK_Graph->SetMarkerStyle(8);
    Electron_number->SetLineWidth(3);

    TGraph *TGraph_Number_of_electron = new TGraph(Number_of_temperature,Temperature_Array,Number_of_electron_G);
    TGraph_Number_of_electron->SetLineColor(2);
    TGraph_Number_of_electron->SetLineStyle(5);
    TGraph_Number_of_electron->SetMarkerColor(2);
    //TLK_Graph->SetMarkerStyle(8);
    TGraph_Number_of_electron->SetLineWidth(2);

    Electron_number->Draw("ALP");
    TGraph_Number_of_electron->Draw("LPsame");
    l3->Draw("LPsame");
    l4->Draw("LPsame");
    l5->Draw("LPsame");

    //c1->SetLogy();
    c1->Print("eh_pair_number_Summary.pdf");
    */
}
/*
TGraph *BLK_Graph = new TGraph(Number_of_temperature,Temperature_Array,BLK);
BLK_Graph->SetLineColor(1);
BLK_Graph->SetTitle("");
BLK_Graph->GetXaxis()->SetTitle("T(K)");
BLK_Graph->GetYaxis()->SetTitle("Log_{10}I(A/cm^{2})");
BLK_Graph->GetXaxis()->SetRangeUser(0,300);
BLK_Graph->GetYaxis()->SetRangeUser(-100,20);
BLK_Graph->SetMarkerColor(1);
//BLK_Graph->SetMarkerStyle(8);
BLK_Graph->SetLineWidth(3);

TGraph *CLK_Graph = new TGraph(Number_of_temperature,Temperature_Array,CLK);
CLK_Graph->SetLineColor(4);
CLK_Graph->SetMarkerColor(4);
//CLK_Graph->SetMarkerStyle(8);
CLK_Graph->SetLineWidth(3);

TGraph *SLK_Graph = new TGraph(Number_of_temperature,Temperature_Array,SLK);
SLK_Graph->SetLineColor(8);
SLK_Graph->SetMarkerColor(8);
//SLK_Graph->SetMarkerStyle(8);
SLK_Graph->SetLineWidth(3);

TGraph *TLK_Graph = new TGraph(Number_of_temperature,Temperature_Array,TLK);
TLK_Graph->SetLineColor(2);
TLK_Graph->SetLineStyle(5);
TLK_Graph->SetMarkerColor(2);
//TLK_Graph->SetMarkerStyle(8);
TLK_Graph->SetLineWidth(2);

TLine *l3=new TLine(0,TMath::Log10(1.6e-12),300,TMath::Log10(1.6e-12));

BLK_Graph->Draw("ALP");
CLK_Graph->Draw("LPsame");
SLK_Graph->Draw("LPsame");
TLK_Graph->Draw("LPsame");
l3->Draw("LPsame");

TLegend *leg= new TLegend(0.6,0.1,0.9,0.4);
leg->SetFillColor(0);
leg->SetFillStyle(0);
leg->SetTextSize(0.06);
leg->SetBorderSize(0);
leg->SetTextFont(20);

leg->AddEntry(BLK_Graph,"Bulk","LP");
leg->AddEntry(CLK_Graph,"Contact","LP");
leg->AddEntry(SLK_Graph,"Surface","LP");
leg->AddEntry(TLK_Graph,"Total","LP");

leg->Draw();
//c1->SetLogy();
c1->Print("Leakage_Current_Summary.pdf");
 */
