#include "TMath.h"


void manual_fit()
{ 
 
  const int data_set = 1;
  char rootfile_name[data_set][100];
  int file_date[data_set];

  FILE *lhs;
  lhs = fopen("rootfile_names.txt","r");
  for(int jj= 0; jj<data_set; jj++)
    {
      fscanf(lhs,"%s  %d\n",&rootfile_name[jj],&file_date[jj]);
    }

    cout << "1: " << endl;
  for (Int_t n=0; n<data_set; n++) 
    { 
      char fname[100];
      sprintf(fname,"%s.root",rootfile_name[n]);
      printf("\nopen....%s date %d\n",rootfile_name[n],file_date[n]);
      TFile *fin=new TFile(fname);
    
      TTree *tr=(TTree*)fin->Get("tr");
        cout << "2: " << endl;

 
      TCut gran = Form("random_trig_on_off==1");  
      TCut garan = Form("random_trig_on_off==0"); 

      TCut os60m_ch0_cut = Form("max_60m[0]<32750");
        TCut min60m_ch0_cut =  Form("");
        TCut ped60m_ch0_cut =  Form("");
        TCut pedt60m_ch0_cut =  Form("");
        TCut ped_diff60m_ch0_cut = Form("");
        cout << "3: " << endl;

/*
        TCut min60m_ch0_cut =  Form("min_60m[0]>-255.0 && min_60m[0]<-245.0");
        TCut ped60m_ch0_cut =  Form("ped_60m[0]>0.0 && ped_60m[0]<-250.0");
        TCut pedt60m_ch0_cut =  Form("pedt_60m[0]>-220 && pedt_60m[0]<0.0");
        TCut ped_diff60m_ch0_cut = Form("(ped_60m[0]-pedt_60m[0])>-200&&(ped_60m[0]-pedt_60m[0])<0");
*/
/*
        TCut min60m_ch0_cut =  Form("min_60m[0]>-300.0 && min_60m[0]<-200.0");
        TCut ped60m_ch0_cut =  Form("ped_60m[0]>-240.0 && ped_60m[0]<-210.0");
        TCut pedt60m_ch0_cut =  Form("pedt_60m[0]>-220 && pedt_60m[0]<0.0");
        TCut ped_diff60m_ch0_cut = Form("(ped_60m[0]-pedt_60m[0])>-220&&(ped_60m[0]-pedt_60m[0])<0");
*/
        /*
        TCut min60m_ch0_cut =  Form("min_60m[0]>-500.0 && min_60m[0]<-300.0");
        TCut ped60m_ch0_cut =  Form("ped_60m[0]>-230.0 && ped_60m[0]<-160.0");
        TCut pedt60m_ch0_cut =  Form("pedt_60m[0]>-200 && pedt_60m[0]<0.0");
        TCut ped_diff60m_ch0_cut = Form("(ped_60m[0]-pedt_60m[0])>-200&&(ped_60m[0]-pedt_60m[0])<0");
*/
        
        /*
      TCut min60m_ch0_cut =  Form("min_60m[0]>-250.0 && min_60m[0]<-200.0");
      TCut ped60m_ch0_cut =  Form("ped_60m[0]>-220.0 && ped_60m[0]<-160.0");  
      TCut pedt60m_ch0_cut =  Form("pedt_60m[0]>-250 && pedt_60m[0]<100.0"); 
      TCut ped_diff60m_ch0_cut = Form("(ped_60m[0]-pedt_60m[0])>-200&&(ped_60m[0]-pedt_60m[0])<100");
*/
        cout << "4: " << endl;

      TCut offset_cut;

      offset_cut = os60m_ch0_cut&&min60m_ch0_cut&&ped60m_ch0_cut&&pedt60m_ch0_cut&&ped_diff60m_ch0_cut;



      TH1F *clbr_ped_he_ran = new TH1F("clbr_ped_he_ran","",75,0,40000);
      tr->Project("clbr_ped_he_ran","max_200m1[0]",garan);
      clbr_ped_he_ran->SetLineColor(6);
        cout << "5: " << endl;

      //TH1F *clbr_max_he = new TH1F("clbr_max_he","",500,0,10000000);
      TH2F *clbr_max_he = new TH2F("clbr_max_he","clbr_max_he",10,-0.25,0.25,100,0,30000);
      //tr->Project("clbr_max_he","max_60m[0]/q_60m[0]",offset_cut);//
      //tr->Project("clbr_max_he","fslope_200m[0]:max_60m[0]",offset_cut);
       // tr->Project("clbr_max_he","max_60m[0]:fslope_200m[0]",offset_cut);
        cout << "6: " << endl;

      clbr_max_he->SetLineColor(4);

      char fout_name[100];
      sprintf(fout_name,"for_fit/for_calib_%d_max_he_New.root",file_date[n]);

      TFile *fout=new TFile(fout_name,"recreate");

        cout << "7: " << endl;

      clbr_ped_he_ran->Write();
      clbr_max_he->Write();
      fout->Close();


    }

}
