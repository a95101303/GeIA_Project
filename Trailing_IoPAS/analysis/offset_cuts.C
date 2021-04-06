#include "TMath.h"
#include <TH2F.h>
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLine.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TFile.h"
#include "TTree.h"
#include "TCut.h"
#include "TCutG.h"
#include "TH1F.h"
#include "TH1.h"
#include "TGraph.h"
#include "TTree.h"


void offset_cuts()
{  

    cout << "1" << endl;

  const int data_set = 1;
  char rootfile_name[data_set][100];
  int file_date[data_set];
    cout << "4" << endl;

  FILE *lhs;
  lhs = fopen("rootfile_names.txt","r");
  for(int jj= 0; jj<data_set; jj++)
    {
      fscanf(lhs,"%s  %d\n",&rootfile_name[jj],&file_date[jj]);
    }
    cout << "7" << endl;

  for (Int_t n=0; n<data_set; n++)
    { 
      char fname[100];
      sprintf(fname,"%s.root",rootfile_name[n]);
      printf("\nopen....%s date %d\n",rootfile_name[n],file_date[n]);
      TFile *fin=new TFile(fname);
    
      TTree *tr=(TTree*)fin->Get("tr");
 
        cout << "1" << endl;
      TCut gran = Form("random_trig_on_off==1");  
      TCut garan = Form("random_trig_on_off==0"); 
        cout << "3" << endl;

      TCut os60m_ch0_cut = Form("max_60m[0]<32750");

        cout << "5" << endl;
        TCut q60m_ch0_cut =  Form("");
        TCut min60m_ch0_cut =  Form("");
        TCut ped60m_ch0_cut =  Form("");
        TCut pedt60m_ch0_cut =  Form("");
        TCut ped_diff60m_ch0_cut = Form("");

        /*
        TCut q60m_ch0_cut =  Form("q_60m[0]/max_60m[0]>-500.0 && q_60m[0]/max_60m[0]<-300.0");
        TCut min60m_ch0_cut =  Form("min_60m[0]>-500.0 && min_60m[0]<-300.0");
        TCut ped60m_ch0_cut =  Form("ped_60m[0]>-230.0 && ped_60m[0]<-160.0");
        TCut pedt60m_ch0_cut =  Form("pedt_60m[0]>-200 && pedt_60m[0]<0.0");
        TCut ped_diff60m_ch0_cut = Form("(ped_60m[0]-pedt_60m[0])>-200&&(ped_60m[0]-pedt_60m[0])<0");
         */
        /*
        TCut min60m_ch0_cut =  Form("min_60m[0]>-300.0 && min_60m[0]<-200.0");
        TCut ped60m_ch0_cut =  Form("ped_60m[0]>-220.0 && ped_60m[0]<-160.0");
        TCut pedt60m_ch0_cut =  Form("pedt_60m[0]>-220 && pedt_60m[0]<0.0");
        TCut ped_diff60m_ch0_cut = Form("(ped_60m[0]-pedt_60m[0])>-80&&(ped_60m[0]-pedt_60m[0])<20");
        cout << "5" << endl;
 
      
      TCut min60m_ch0_cut =  Form("min_60m[0]>-240.0 && min_60m[0]<-220.0");
      TCut ped60m_ch0_cut =  Form("ped_60m[0]>-201.0 && ped_60m[0]<-180.0");
      TCut pedt60m_ch0_cut =  Form("pedt_60m[0]>-195 && pedt_60m[0]<-185.0");
      TCut ped_diff60m_ch0_cut = Form("(ped_60m[0]-pedt_60m[0])>-20&&(ped_60m[0]-pedt_60m[0])<20");
       */
        
      TCut offset_cut;

      offset_cut = os60m_ch0_cut&&min60m_ch0_cut&&ped60m_ch0_cut&&pedt60m_ch0_cut&&ped_diff60m_ch0_cut;

        cout << "7" << endl;

      //*************

      char hist_fname[100];
      sprintf(hist_fname,"output_root/offset_cut_%d_2.root",file_date[n],file_date[n]);
      TFile *fwrite=new TFile(hist_fname,"recreate");

        cout << "7" << endl;

      //*************
      printf("Channel0 cuts\n");

        
    TH2F *q60m_ch0_hist = new TH2F("q60m_ch0_hist","",10000,0,50000000,3500,-2000,33000);
    tr->Project("q60m_ch0_hist","max_60m[0]:q_60m[0]",garan&&os60m_ch0_cut);
    q60m_ch0_hist->SetMarkerStyle(21);
    q60m_ch0_hist->SetMarkerSize(0.2);
    q60m_ch0_hist->Write();
    q60m_ch0_hist->Delete();
    
        cout << "8" << endl;

    TH2F *q60m_ch0_ran_hist = new TH2F("q60m_ch0_ran_hist","",10000,0,50000000,3500,-2000,33000);
    tr->Project("q60m_ch0_ran_hist","max_60m[0]:q_60m[0]",gran&&os60m_ch0_cut);
    q60m_ch0_ran_hist->SetMarkerStyle(21);
    q60m_ch0_ran_hist->SetMarkerSize(0.2);
    q60m_ch0_ran_hist->SetMarkerColor(2);
    q60m_ch0_ran_hist->Write();
    q60m_ch0_ran_hist->Delete();
    
    TH2F *q60m_ch0_af_hist = new TH2F("q60m_ch0_af_hist","",10000,0,50000000,3500,-2000,33000);
    tr->Project("q60m_ch0_af_hist","max_60m[0]:q_60m[0]",garan&&os60m_ch0_cut&&q60m_ch0_cut);
    q60m_ch0_af_hist->SetMarkerStyle(21);
    q60m_ch0_af_hist->SetMarkerSize(0.2);
    q60m_ch0_af_hist->SetMarkerColor(4);
    q60m_ch0_af_hist->Write();
    q60m_ch0_af_hist->Delete();
    
    TH2F *q60m_ch0_af_ran_hist = new TH2F("q60m_ch0_af_ran_hist","",10000,0,50000000,3500,-2000,33000);
    tr->Project("q60m_ch0_af_ran_hist","max_60m[0]:q_60m[0]",gran&&os60m_ch0_cut&&q60m_ch0_cut);
    q60m_ch0_af_ran_hist->SetMarkerStyle(21);
    q60m_ch0_af_ran_hist->SetMarkerSize(0.2);
    q60m_ch0_af_ran_hist->SetMarkerColor(2);
    q60m_ch0_af_ran_hist->Write();
    q60m_ch0_af_ran_hist->Delete();
    
    cout << "8" << endl;
        
        /*************/

      TH2F *min60m_ch0_hist = new TH2F("min60m_ch0_hist","",3500,-2000,33000,200,-1000,4000);
      tr->Project("min60m_ch0_hist","min_60m[0]:max_60m[0]",garan&&os60m_ch0_cut);
      min60m_ch0_hist->SetMarkerStyle(21);
      min60m_ch0_hist->SetMarkerSize(0.2);
      min60m_ch0_hist->Write();
      min60m_ch0_hist->Delete();

      TH2F *min60m_ch0_ran_hist = new TH2F("min60m_ch0_ran_hist","",3500,-2000,33000,200,-1000,4000);
      tr->Project("min60m_ch0_ran_hist","min_60m[0]:max_60m[0]",gran&&os60m_ch0_cut);
      min60m_ch0_ran_hist->SetMarkerStyle(21);
      min60m_ch0_ran_hist->SetMarkerSize(0.2);
      min60m_ch0_ran_hist->SetMarkerColor(2);
      min60m_ch0_ran_hist->Write();
      min60m_ch0_ran_hist->Delete();

      TH2F *min60m_ch0_af_hist = new TH2F("min60m_ch0_af_hist","",3500,-2000,33000,200,-1000,4000);
      tr->Project("min60m_ch0_af_hist","min_60m[0]:max_60m[0]",garan&&os60m_ch0_cut&&min60m_ch0_cut);
      min60m_ch0_af_hist->SetMarkerStyle(21);
      min60m_ch0_af_hist->SetMarkerSize(0.2);
      min60m_ch0_af_hist->SetMarkerColor(4);
      min60m_ch0_af_hist->Write();
      min60m_ch0_af_hist->Delete();

        TH2F *min60m_ch0_af_ran_hist = new TH2F("min60m_ch0_af_ran_hist","",3500,-2000,33000,200,-1000,4000);
        tr->Project("min60m_ch0_af_ran_hist","min_60m[0]:max_60m[0]",gran&&os60m_ch0_cut&&min60m_ch0_cut);
        min60m_ch0_af_ran_hist->SetMarkerStyle(21);
        min60m_ch0_af_ran_hist->SetMarkerSize(0.2);
        min60m_ch0_af_ran_hist->SetMarkerColor(2);
        min60m_ch0_af_ran_hist->Write();
        min60m_ch0_af_ran_hist->Delete();

        cout << "9" << endl;

      /*************/

      TH2F *ped60m_ch0_hist = new TH2F("ped60m_ch0_hist","",3500,-2000,33000,1000,-1000,4000);
      tr->Project("ped60m_ch0_hist","ped_60m[0]:max_60m[0]",garan&&os60m_ch0_cut&&min60m_ch0_cut);
      ped60m_ch0_hist->SetMarkerStyle(21);
      ped60m_ch0_hist->SetMarkerSize(0.2);
      ped60m_ch0_hist->Write();
      ped60m_ch0_hist->Delete();

      TH2F *ped60m_ch0_ran_hist = new TH2F("ped60m_ch0_ran_hist","",3500,-2000,33000,1000,-1000,4000);
      tr->Project("ped60m_ch0_ran_hist","ped_60m[0]:max_60m[0]",gran&&os60m_ch0_cut&&min60m_ch0_cut);
      ped60m_ch0_ran_hist->SetMarkerStyle(21);
      ped60m_ch0_ran_hist->SetMarkerSize(0.2);
      ped60m_ch0_ran_hist->SetMarkerColor(2);
      ped60m_ch0_ran_hist->Write();
      ped60m_ch0_ran_hist->Delete();

      TH2F *ped60m_ch0_af_hist = new TH2F("ped60m_ch0_af_hist","",3500,-2000,33000,1000,-1000,4000);
      tr->Project("ped60m_ch0_af_hist","ped_60m[0]:max_60m[0]",garan&&os60m_ch0_cut&&min60m_ch0_cut&&ped60m_ch0_cut);
      ped60m_ch0_af_hist->SetMarkerStyle(21);
      ped60m_ch0_af_hist->SetMarkerSize(0.2);
      ped60m_ch0_af_hist->SetMarkerColor(4);
      ped60m_ch0_af_hist->Write();
      ped60m_ch0_af_hist->Delete();

      TH2F *ped60m_ch0_af_ran_hist = new TH2F("ped60m_ch0_af_ran_hist","",3500,-2000,33000,1000,-1000,4000);
      tr->Project("ped60m_ch0_af_ran_hist","ped_60m[0]:max_60m[0]",gran&&os60m_ch0_cut&&min60m_ch0_cut&&ped60m_ch0_cut);
      ped60m_ch0_af_ran_hist->SetMarkerStyle(21);
      ped60m_ch0_af_ran_hist->SetMarkerSize(0.2);
      ped60m_ch0_af_ran_hist->SetMarkerColor(2);
      ped60m_ch0_af_ran_hist->Write();
      ped60m_ch0_af_ran_hist->Delete();

        cout << "11" << endl;

      /*************/

      TH2F *pedt60m_ch0_hist = new TH2F("pedt60m_ch0_hist","",3500,-2000,33000,1000,-1000,4000);
      tr->Project("pedt60m_ch0_hist","pedt_60m[0]:max_60m[0]",garan&&os60m_ch0_cut&&min60m_ch0_cut&&ped60m_ch0_cut);
      pedt60m_ch0_hist->SetMarkerStyle(21);
      pedt60m_ch0_hist->SetMarkerSize(0.2);
      pedt60m_ch0_hist->Write();
      pedt60m_ch0_hist->Delete();

      TH2F *pedt60m_ch0_ran_hist = new TH2F("pedt60m_ch0_ran_hist","",3500,-2000,33000,1000,-1000,4000);
      tr->Project("pedt60m_ch0_ran_hist","pedt_60m[0]:max_60m[0]",gran&&os60m_ch0_cut&&min60m_ch0_cut&&ped60m_ch0_cut);
      pedt60m_ch0_ran_hist->SetMarkerStyle(21);
      pedt60m_ch0_ran_hist->SetMarkerSize(0.2);
      pedt60m_ch0_ran_hist->SetMarkerColor(2);
      pedt60m_ch0_ran_hist->Write();
      pedt60m_ch0_ran_hist->Delete();

      TH2F *pedt60m_ch0_af_hist = new TH2F("pedt60m_ch0_af_hist","",3500,-2000,33000,1000,-1000,4000);
      tr->Project("pedt60m_ch0_af_hist","pedt_60m[0]:max_60m[0]",garan&&os60m_ch0_cut&&min60m_ch0_cut&&ped60m_ch0_cut&&pedt60m_ch0_cut);
      pedt60m_ch0_af_hist->SetMarkerStyle(21);
      pedt60m_ch0_af_hist->SetMarkerSize(0.2);
      pedt60m_ch0_af_hist->SetMarkerColor(4);
      pedt60m_ch0_af_hist->Write();
      pedt60m_ch0_af_hist->Delete();

        TH2F *pedt60m_ch0_af_ran_hist = new TH2F("pedt60m_ch0_af_ran_hist","",3500,-2000,33000,1000,-1000,4000);
        tr->Project("pedt60m_ch0_af_ran_hist","ped_60m[0]:max_60m[0]",gran&&os60m_ch0_cut&&min60m_ch0_cut&&ped60m_ch0_cut&&pedt60m_ch0_cut);
        pedt60m_ch0_af_ran_hist->SetMarkerStyle(21);
        pedt60m_ch0_af_ran_hist->SetMarkerSize(0.2);
        pedt60m_ch0_af_ran_hist->SetMarkerColor(2);
        pedt60m_ch0_af_ran_hist->Write();
        pedt60m_ch0_af_ran_hist->Delete();
        cout << "13" << endl;


      /*************/
      printf("ped-diff cuts\n");

      TH2F *diff_ped60m_ch0_pedt60m_ch0_hist = new TH2F("diff_ped60m_ch0_pedt60m_ch0_hist","",3500,-2000,33000,1000,-1000,1000);
      tr->Project("diff_ped60m_ch0_pedt60m_ch0_hist","ped_60m[0]-pedt_60m[0]:max_60m[0]",garan&&os60m_ch0_cut&&min60m_ch0_cut&&ped60m_ch0_cut&&pedt60m_ch0_cut);
      diff_ped60m_ch0_pedt60m_ch0_hist->SetMarkerStyle(21);
      diff_ped60m_ch0_pedt60m_ch0_hist->SetMarkerSize(0.2);
      diff_ped60m_ch0_pedt60m_ch0_hist->Write();
      diff_ped60m_ch0_pedt60m_ch0_hist->Delete();

      TH2F *diff_ped60m_ch0_pedt60m_ch0_ran_hist = new TH2F("diff_ped60m_ch0_pedt60m_ch0_ran_hist","",3500,-2000,33000,1000,-1000,1000);
      tr->Project("diff_ped60m_ch0_pedt60m_ch0_ran_hist","ped_60m[0]-pedt_60m[0]:max_60m[0]",gran&&os60m_ch0_cut&&min60m_ch0_cut&&ped60m_ch0_cut&&pedt60m_ch0_cut);
      diff_ped60m_ch0_pedt60m_ch0_ran_hist->SetMarkerStyle(21);
      diff_ped60m_ch0_pedt60m_ch0_ran_hist->SetMarkerSize(0.2);
      diff_ped60m_ch0_pedt60m_ch0_ran_hist->SetMarkerColor(2);
      diff_ped60m_ch0_pedt60m_ch0_ran_hist->Write();
      diff_ped60m_ch0_pedt60m_ch0_ran_hist->Delete();


      TH2F *diff_ped60m_ch0_pedt60m_ch0_af_hist = new TH2F("diff_ped60m_ch0_pedt60m_ch0_af_hist","",3500,-2000,33000,1000,-1000,1000);
      tr->Project("diff_ped60m_ch0_pedt60m_ch0_af_hist","ped_60m[0]-pedt_60m[0]:max_60m[0]",garan&&os60m_ch0_cut&&min60m_ch0_cut&&ped60m_ch0_cut&&pedt60m_ch0_cut&&ped_diff60m_ch0_cut);
      diff_ped60m_ch0_pedt60m_ch0_af_hist->SetMarkerStyle(21);
      diff_ped60m_ch0_pedt60m_ch0_af_hist->SetMarkerSize(0.2);
      diff_ped60m_ch0_pedt60m_ch0_af_hist->SetMarkerColor(4);
      diff_ped60m_ch0_pedt60m_ch0_af_hist->Write();
      diff_ped60m_ch0_pedt60m_ch0_af_hist->Delete();

        TH2F *diff_ped60m_ch0_pedt60m_ch0_af_ran_hist = new TH2F("diff_ped60m_ch0_pedt60m_ch0_af_ran_hist","",3500,-2000,33000,1000,-1000,1000);
        tr->Project("diff_ped60m_ch0_pedt60m_ch0_af_ran_hist","ped_60m[0]-pedt_60m[0]:max_60m[0]",gran&&os60m_ch0_cut&&min60m_ch0_cut&&ped60m_ch0_cut&&pedt60m_ch0_cut&&ped_diff60m_ch0_cut);
        diff_ped60m_ch0_pedt60m_ch0_af_ran_hist->SetMarkerStyle(21);
        diff_ped60m_ch0_pedt60m_ch0_af_ran_hist->SetMarkerSize(0.2);
        diff_ped60m_ch0_pedt60m_ch0_af_ran_hist->SetMarkerColor(4);
        diff_ped60m_ch0_pedt60m_ch0_af_ran_hist->Write();
        diff_ped60m_ch0_pedt60m_ch0_af_ran_hist->Delete();
        cout << "15" << endl;

      /*************/

      fwrite->Close();
    }

}
