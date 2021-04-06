#include "TCanvas.h"
#include "txn_ntupe_class_L114.h"

void run_header::read_rheader(FILE *raw)
{
  fread(&program_start_time,sizeof(Double_t),1,raw);
  
  fread(&pos_triger_2g1,sizeof(Int_t),1,raw);
  fread(&pre_triger_2g1,sizeof(Int_t),1,raw);
  fread(&open_channel_2g1,sizeof(Int_t),1,raw);
  fread(&on_off_2g1,sizeof(Int_t),1,raw);

  fread(&pos_triger_200M1,sizeof(Int_t),1,raw);
  fread(&pre_triger_200M1,sizeof(Int_t),1,raw);
  fread(&open_channel_200M1,sizeof(Int_t),1,raw);
  fread(&on_off_200M1,sizeof(Int_t),1,raw);

  fread(&pos_triger_2g2,sizeof(Int_t),1,raw);
  fread(&pre_triger_2g2,sizeof(Int_t),1,raw);
  fread(&open_channel_2g2,sizeof(Int_t),1,raw);
  fread(&on_off_2g2,sizeof(Int_t),1,raw);

  fread(&pos_triger_200M2,sizeof(Int_t),1,raw);
  fread(&pre_triger_200M2,sizeof(Int_t),1,raw);
  fread(&open_channel_200M2,sizeof(Int_t),1,raw);
  fread(&on_off_200M2,sizeof(Int_t),1,raw);

  fread(&open_channel_60M,sizeof(Int_t),1,raw);
  fread(&pre_triger_60M,sizeof(Int_t),1,raw);
  fread(&pos_triger_60M,sizeof(Int_t),1,raw);
  fread(&on_off_60M,sizeof(Int_t),1,raw);
 
  fread(&open_channel_20M1,sizeof(Int_t),1,raw);
  fread(&pre_triger_20M1,sizeof(Int_t),1,raw);
  fread(&pos_triger_20M1,sizeof(Int_t),1,raw);
  fread(&on_off_20M1,sizeof(Int_t),1,raw);
  
  fread(&open_channel_20M2,sizeof(Int_t),1,raw);
  fread(&pre_triger_20M2,sizeof(Int_t),1,raw);  
  fread(&pos_triger_20M2,sizeof(Int_t),1,raw);
  fread(&on_off_20M2,sizeof(Int_t),1,raw);

  fread(&run_start_time,sizeof(Double_t),1,raw);
  fread(&FPGA_clock1,sizeof(ULong64_t),1,raw);
  fread(&FPGA_clock2,sizeof(ULong64_t),1,raw);
  fread(&FPGA_clock3,sizeof(ULong64_t),1,raw);

  length_2g1=pre_triger_2g1+pos_triger_2g1;
  length_2g2=pre_triger_2g2+pos_triger_2g2;
  length_60m=pre_triger_60M+pos_triger_60M;
  length_20m1=pre_triger_20M1+pos_triger_20M1;
  length_20m2=pre_triger_20M2+pos_triger_20M2;
  length_200m1=pre_triger_200M1+pos_triger_200M1;
  length_200m2=pre_triger_200M2+pos_triger_200M2;
  
  open_ch_2g1=ascii_to_binary(31,open_channel_2g1);      if(on_off_2g1==0)  { open_ch_2g1=0; }
  open_ch_2g2=ascii_to_binary(31,open_channel_2g2);      if(on_off_2g2==0)  { open_ch_2g2=0; }
  open_ch_60m=ascii_to_binary(31,open_channel_60M);      if(on_off_60M==0)  { open_ch_60m=0; }
  open_ch_20m1=ascii_to_binary(31,open_channel_20M1);    if(on_off_20M1==0) { open_ch_20m1=0; }
  open_ch_20m2=ascii_to_binary(31,open_channel_20M2);    if(on_off_20M2==0) { open_ch_20m2=0; }
  open_ch_200m1=ascii_to_binary(31,open_channel_200M1);  if(on_off_200M1==0){ open_ch_200m1=0; }
  open_ch_200m2=ascii_to_binary(31,open_channel_200M2);  if(on_off_200M2==0){ open_ch_200m2=0; }

  
  printf("2G1: "); if(on_off_2g1==1) {printf(" ON  ");} else {printf(" OFF  ");}
  printf("OPEN CHANNELS=%d  LENGTH=%d\n",open_ch_2g1,length_2g1);
  printf("2G2: "); if(on_off_2g2==1) {printf(" ON  ");} else {printf(" OFF  ");}
  printf("OPEN CHANNELS=%d  LENGTH=%d\n",open_ch_2g2,length_2g2);
  printf("60M: "); if(on_off_60M==1) {printf(" ON  ");} else {printf(" OFF  ");}
  printf(" OPEN CHANNELS=%d  LENGTH=%d\n",open_ch_60m,length_60m);
  printf("20M1: "); if(on_off_20M1==1) {printf("ON  ");} else {printf("OFF  ");}
  printf(" OPEN CHANNELS=%d  LENGTH=%d\n",open_ch_20m1,length_20m1);
  printf("20M2: "); if(on_off_20M2==1) {printf("ON  ");} else {printf("OFF  ");}
  printf("OPEN CHANNELS=%d  LENGTH=%d\n",open_ch_20m2,length_20m2);
  printf("200M1:"); if(on_off_200M1==1) {printf("ON  ");} else {printf("OFF  ");}
  printf(" OPEN CHANNELS=%d  LENGTH=%d\n",open_ch_200m1,length_200m1);
  printf("200M2:"); if(on_off_200M2==1) {printf("ON  ");} else {printf("OFF  ");}
  printf("OPEN CHANNELS=%d  LENGTH=%d\n",open_ch_200m2,length_200m2);
 

}

Int_t run_header::ascii_to_binary(Int_t bit_num,Int_t num)
{
  Int_t open_channels=0;

  for(Int_t bit=bit_num;bit>=0;bit--)
  {
    if(num&(1<<bit))
    { 
      open_channels++;
    }
  }
  
  return open_channels;
}

void event_header::read_eheader(FILE *raw)
{
    fread(&event_start_time1,sizeof(Long64_t),1,raw);
    fread(&event_start_time2,sizeof(Long64_t),1,raw);
    fread(&event_start_time3,sizeof(Long64_t),1,raw);
    fread(&veto_tdc,sizeof(Long64_t),32,raw);
    fread(&former_veto_tdc,sizeof(Long64_t),32,raw);
    fread(&late_former_veto_tdc,sizeof(Long64_t),32,raw);
    fread(&late_former_veto_tdc2,sizeof(Long64_t),32,raw);

    fread(&veto_tdc_or,sizeof(Long64_t),13,raw);
    fread(&former_veto_tdc_or,sizeof(Long64_t),13,raw);
    fread(&late_former_veto_tdc_or,sizeof(Long64_t),13,raw);
    fread(&late_former_veto_tdc2_or,sizeof(Long64_t),13,raw);
    
    fread(&signal_inhibit_time,sizeof(Long64_t),8,raw);  
    fread(&hitpattern,sizeof(Long64_t),1,raw);
    fread(&upper_veto_hp,sizeof(Long64_t),1,raw);
    fread(&other_veto_hp,sizeof(Long64_t),1,raw);
    fread(&other_inner_veto_hp,sizeof(Long64_t),1,raw);

    fread(&event_end_time,sizeof(Long64_t),1,raw);

    fread(&mins_OV,sizeof(Short_t),16,raw);
    fread(&mins_OV2,sizeof(Short_t),8,raw);
        
    
}

void event_body::read_ebody(FILE *raw)
{ 
  
  fread(&pulse_2g1,sizeof(Char_t),open_ch_2g1*length_2g1,raw);
  fread(&pulse_60m,sizeof(Short_t),open_ch_60m*length_60m,raw);
  fread(&pulse_20m1,sizeof(Short_t),open_ch_20m1*length_20m1,raw);
  fread(&pulse_20m2,sizeof(Short_t),open_ch_20m2*length_20m2,raw);
  fread(&pulse_200m1,sizeof(Short_t),open_ch_200m1*length_200m1,raw);
  fread(&pulse_200m2,sizeof(Short_t),open_ch_200m2*length_200m2,raw);
  fread(&pulse_2g2,sizeof(Char_t),open_ch_2g2*length_2g2,raw);
}

void run_tail::read_rtail(FILE *raw)
{
  fseek(raw,-(sizeof(UInt_t)*107+sizeof(Double_t)),SEEK_END);
  fread(&scalar,sizeof(UInt_t),107,raw);
  fread(&run_end_time,sizeof(Double_t),1,raw);
}

UInt_t eheader_size(void)
{
  return sizeof(Long64_t)*(196)+sizeof(Short_t)*24;
}

UInt_t ebody_size(void)
{
  total_ebody_size=0;

  total_ebody_size+=sizeof(Char_t)*open_ch_2g1*length_2g1;
  total_ebody_size+=sizeof(Char_t)*open_ch_2g2*length_2g2;
  total_ebody_size+=sizeof(Short_t)*open_ch_60m*length_60m;
  total_ebody_size+=sizeof(Short_t)*open_ch_20m1*length_20m1;
  total_ebody_size+=sizeof(Short_t)*open_ch_20m2*length_20m2;
  total_ebody_size+=sizeof(Short_t)*open_ch_200m1*length_200m1;
  total_ebody_size+=sizeof(Short_t)*open_ch_200m2*length_200m2;
  
  return total_ebody_size;//+8;
}



unsigned long event_header::filesize(FILE *raw)
{
  unsigned long a;

  fseek(raw,0,SEEK_END);
  a=ftell(raw);
  fseek(raw,0,SEEK_SET);
  
  return a;
}

unsigned long filesize(FILE *raw)
{
  unsigned long a;

  fseek(raw,0,SEEK_END);
  a=ftell(raw);
  fseek(raw,0,SEEK_SET);
  
  return a;
}




void get_name(char *name, char *filename)
{
    unsigned int len = strlen(filename);
    unsigned int start_count = 0;
    unsigned int end_count = 0;

    for(unsigned int i=0;i<len;i++)
    {
        if(filename[i]=='/') { start_count = i+1; }
    }
    
    end_count = len;
    strcpy(name,"NoName");

    if((end_count-start_count)>0)
    {
        for(unsigned int i=start_count;i<end_count;i++)
        { name[i-start_count]=filename[i];}
        name[end_count-start_count]='\0';
    }
}

void get_date(char *date, char *filename)
{
    unsigned int len = strlen(filename);
    unsigned int check = 0;
    unsigned int i = 0;

    while ((i<(len-6))&&(check!=7))
    {
        check = 0;
        date[0] = filename[i];
      
        if ((date[0]=='0')||(date[0]=='1')) { check = check + 1; }

        for(unsigned int j=1;j<6;j++)
        {
	  date[j]=filename[j+i];
	  if ((date[j]>='0')&&(date[j]<='9')) { check = check + 1; }
        }
	
        if((filename[i+6]<'0')||(filename[i+6]>'9')) { check = check + 1; }
        i++;
    }
    
    if (check!=7)
      { strcpy(date,"NoDate"); }
    else
    { date[6] = '\0'; }
}




