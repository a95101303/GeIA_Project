
void final_spec()
{ 

  char fname[100];
  sprintf(fname,"190627_1828_Ba_Gain20_FADC_RAW_Data_1_150.root");
  // printf("\nopen....%s date %d\n",rootfile_name[n],file_date[n]);
  TFile *fin=new TFile(fname);
    
  TTree *tr=(TTree*)fin->Get("tr");
 
 
  TCut gran = Form("random_trig_on_off==1");  
  TCut garan = Form("random_trig_on_off==0"); 

  TCut os60m_ch0_cut = Form("max_60m[0]<32750");

  TCut min60m_ch0_cut =  Form("min_60m[0]>-250.0 && min_60m[0]<-200.0");  
  TCut ped60m_ch0_cut =  Form("ped_60m[0]>-220.0 && ped_60m[0]<-160.0");  
  TCut pedt60m_ch0_cut =  Form("pedt_60m[0]>-250 && pedt_60m[0]<100.0"); 
  TCut ped_diff60m_ch0_cut = Form("(ped_60m[0]-pedt_60m[0])>-200&&(ped_60m[0]-pedt_60m[0])<100");

  TCut offset_cut;

  offset_cut = os60m_ch0_cut&&min60m_ch0_cut&&ped60m_ch0_cut&&pedt60m_ch0_cut&&ped_diff60m_ch0_cut;



  TH1F *clbr_max_he_optq = new TH1F("clbr_max_he_optq","",1800,0,1800);
  tr->Project("clbr_max_he_optq","7.01704143e-05*optq_60m[0]+1.90551059e+01",garan);//
  clbr_max_he_optq->SetLineColor(4);

  TH1F *clbr_max_he_before_pq = new TH1F("clbr_max_he_before_pq","",1800,0,1800);
  tr->Project("clbr_max_he_before_pq","1.28191572e-04*before_pq_60m[0]+1.48431412e+02",garan);//
  clbr_max_he_before_pq->SetLineColor(4);

  TH1F *clbr_max_he_after_pq = new TH1F("clbr_max_he_after_pq","",1800,0,1800);
  tr->Project("clbr_max_he_after_pq","1.18738431e-04*after_pq_60m[0]+7.82492821e+01",garan);//
  clbr_max_he_after_pq->SetLineColor(4);

  TH1F *clbr_max_he_max = new TH1F("clbr_max_he_max","",1800,0,1800);
  tr->Project("clbr_max_he_max","5.28340350e-02*max_60m[0]+8.83059656e+00",garan);//
  clbr_max_he_max->SetLineColor(4);

  TH1F *clbr_max_he_q = new TH1F("clbr_max_he_q","",1800,0,1800);
  tr->Project("clbr_max_he_q","6.75631888e-05*q_60m[0]+7.23595847e+01",garan);//
  clbr_max_he_q->SetLineColor(4);

  char fout_name[100];
  sprintf(fout_name,"for_fit/190627_small_run_final_calib_all.root");

  TFile *fout=new TFile(fout_name,"recreate");


    clbr_max_he_optq->Write();
    clbr_max_he_before_pq->Write();
    clbr_max_he_after_pq->Write();
    clbr_max_he_max->Write();
    clbr_max_he_q->Write();

  fout->Close();

}
