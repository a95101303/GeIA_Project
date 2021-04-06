void spec0()
{
  TFile *file = new TFile("170701_1845_Ba_top12cm_thres7mV_FADC_RAW_Data_1_223.root");
  TTree *tr = (TTree*)file->Get("tr");

  TCut offset_cut = Form("min_200m1[1]>-30000 && ped_200m1[1]<-29680 && random_trig_on_off==0");
  
  
  TCut offset_cut_ch3 = Form("max_60m[3]<32750 &&  min_60m[3]>-300 && min_60m[3]<-100 && ped_60m[3]<0 && random_trig_on_off==0");
  
  TCut offset_cut_ch2 = Form("max_60m[2]<32750 &&  min_60m[2]>-400 && min_60m[2]<-100 && ped_60m[2]<0 && random_trig_on_off==0");
  
  TCut offset_cut_ch0 = Form("max_60m[0]<32750 &&  min_60m[0]>-600 && min_60m[0]<-200 && ped_60m[0]<0 && random_trig_on_off==0");
  TCut offset_cut_200m2ch1 = Form("min_200m2[1]>-14700 && min_200m2[1]<-14500 && ped_200m2[1]>-14500 && ped_200m2[1]<-14400 &&  max_time_bin_200m2[1]>4000&& max_time_bin_200m2[1]<6100&& pedt_200m2[1]<-14450");

  //QA
      Double_t  qa_200m2ch0_up_x1 = 32324.6;
      Double_t  qa_200m2ch0_up_y1 = -5.26606e+07;
      Double_t  qa_200m2ch0_up_x2 = -15182.7;
      Double_t  qa_200m2ch0_up_y2 = -1.73879e+08;
      Double_t  qa_200m2ch0_up_slope;
      Double_t  qa_200m2ch0_up_intercept;
      qa_200m2ch0_up_slope = ( qa_200m2ch0_up_y2- qa_200m2ch0_up_y1)/( qa_200m2ch0_up_x2- qa_200m2ch0_up_x1);
      qa_200m2ch0_up_intercept = (( qa_200m2ch0_up_y2- qa_200m2ch0_up_y1)/( qa_200m2ch0_up_x2- qa_200m2ch0_up_x1))*(-1* qa_200m2ch0_up_x1)+ qa_200m2ch0_up_y1;
      TCut  qa_200m2ch0_up = Form("q_200m2[1]<(%f*(max_200m2[1])+%f)", qa_200m2ch0_up_slope, qa_200m2ch0_up_intercept,qa_200m2ch0_up_x1,qa_200m2ch0_up_x1);
      



  TCanvas *plot1 = new TCanvas("plot1");
  // gPad->SetLogy(1);
 
  /*
  char E_def[200];
  sprintf(E_def,"-2.11117605e+01+ 1.50264e-02*(max_200m1[1]+31000)");
  TH1F *ba_hist = new TH1F("ba_hist","",8500,-10,750);
  tr->Project("ba_hist",E_def,offset_cut);
  ba_hist->Draw();
  */
  
  /*
  char E_def[200];
  sprintf(E_def,"7.70647823e-01+ 8.23776805e-03*(max_60m[3])");
  TH1F *ba_hist = new TH1F("ba_hist","",3100,-10,300);
  tr->Project("ba_hist",E_def,offset_cut_ch3 );
  ba_hist->Draw();
  */

   
  // char E_def[200];
  // sprintf(E_def,"8.93900870e-01+ 4.07626072e-03*(max_60m[2])");
  // TH1F *ba_hist = new TH1F("ba_hist","",350,-10,140);
  // tr->Project("ba_hist",E_def,offset_cut_ch2 );
  // ba_hist->Draw();
  

 
   // char E_def[200];
   // sprintf(E_def,"(1.70832448e-01+ 8.08395296e-04*(max_60m[0]))/2.0");
   // TH1F *ba_hist = new TH1F("ba_hist","",280,0,14);
   // tr->Project("ba_hist",E_def,offset_cut_ch0);
   // ba_hist->Draw();
 
   /*
   TFile *file = new TFile("170630_1950_Cs_thres4mV_FADC_RAW_Data_1_239.root");
   TTree *tr = (TTree*)file->Get("tr");
 
   char E_def1[200];
   sprintf(E_def1,"(1.70832448e-01+ 8.08395296e-04*(max_60m[0]))");
   TH1F *ba_hist1 = new TH1F("ba_hist1","",600,0,30);
   tr->Project("ba_hist1",E_def1,offset_cut_ch0);
   ba_hist1->SetLineColor(2);
   ba_hist1->Scale(1.0/31.0);
   ba_hist1->Draw("same");
   */

   TH1F *rs = new TH1F("rs","",675,-20000,34000);
   tr->Project("rs","max_200m2[1]",offset_cut_200m2ch1&&qa_200m2ch0_up);
   rs->Draw();
    
  // TH1F *rs = new TH1F("rs","",3400,-1000,33000);
  // tr->Project("rs","max_60m[2]");
  // rs->Draw();
    


  //TH1F *rs = new TH1F("rs","",350,-10,140);
  // TH1F *rs = new TH1F("rs","",2100,-10,370);
  //tr->Project("rs","(81.0/25310.0)*max_60m[2]");
  //rs->Draw();

  // pedt_200m1[0]-ped_200m1[0]

  //tr->Draw("pedt_200m2[1]-ped_200m2[1]:max_200m2[1]",offset_cut_200m2ch1&&qa_200m2ch0_up);
  }

