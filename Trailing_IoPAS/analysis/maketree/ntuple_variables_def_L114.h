#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TStyle.h"
#include "avmw_filter_lakhs.h"
#include "lakhs_sgolay_200m.h"
#include "tanh_fit.h"
// EVENT HEADER RELEVENT VARIABLES
void event_header::eh_variables(Int_t run,Int_t initial_run,Int_t block)
{
  for(Int_t sub_block=0;sub_block<100;sub_block++)
    {
      // TRIGGER HITPATTERN ON OFF
      if(hitpattern&(1<<0))
	{ pulser_on_off=1; }
      else
	{ pulser_on_off=0; }

      if(hitpattern&(1<<1))
	{ ge_on_off=1; }
      else
	{ ge_on_off=0; }

      if(hitpattern&(1<<2))
	{ random_trig_on_off=1; }
      else
	{ random_trig_on_off=0; }

      if(hitpattern&(1<<3))
	{ nai_1_on_off=1; }
      else
	{ nai_1_on_off=0; }

      if(hitpattern&(1<<4))
	{ nai_2_on_off=1; }
      else
	{ nai_2_on_off=0; }

      if(hitpattern&(1<<5))
	{ nai_3_on_off=1; }
      else
	{ nai_3_on_off=0; }

      if(hitpattern&(1<<6))
	{ nai_5_on_off=1; }
      else
	{ nai_5_on_off=0; }

      // VETO ON OFF
      for(Int_t bit=0;bit<32;bit++)
	{
	  if(upper_veto_hp&(1<<bit))
	    { upper_veto_on_off[bit]=1; }
	  else
	    { upper_veto_on_off[bit]=0; }

	  if(other_veto_hp&(1<<bit))
	    { other_veto_on_off[bit]=1; }
	  else
	    { other_veto_on_off[bit]=0; }

	  if(other_inner_veto_hp&(1<<bit))
	    { other_inner_on_off[bit]=1; }
	  else
	    { other_inner_on_off[bit]=0; }

	} 

    } // END of sub_block loop

}



// NT Variables for 60M
void event_body::nt_variables_60m(int random_bit)
{
  for(Int_t ch=0;ch<OPEN_CHANNEL_60M;ch++)
    {
      q_60m[ch]=0;
      max_60m[ch]=-100000; 
      omax_60m[ch]=-100000;
      tmax_60m[ch]=-100000;
      min_60m[ch]=100000;
      max_time_bin_60m[ch]=0xFFFF;
      omax_time_bin_60m[ch]=0xFFFF;
      tmax_time_bin_60m[ch]=0xFFFF;
      min_time_bin_60m[ch]=0xFFFF;
      ped_60m[ch]=0.0;
      pedt_60m[ch]=0.0;
      before_pq_60m[ch]=0.0;
      after_pq_60m[ch]=0.0;
      optq_60m[ch]=0;  
      q2_60m[ch]=0;
      t0_60m[ch]=0.0;
      T2_60m[ch]=0.0;
      
      for(unsigned int bin=0;bin<length_60m;bin++)
	{ 
	  q_60m[ch]+=(double)pulse_60m[ch][bin];
	
	  if(pulse_60m[ch][bin]>max_60m[ch])
	    { 
	      max_60m[ch]=pulse_60m[ch][bin];
	      max_time_bin_60m[ch]=bin;
	    }
	
	  if((pulse_60m[ch][bin]>omax_60m[ch])&&(bin>=800)&&(bin<3500))
	    { 
	      omax_60m[ch]=pulse_60m[ch][bin];
	      omax_time_bin_60m[ch]=bin;
	    } 
	  
	  if(pulse_60m[ch][bin]<min_60m[ch])
	    {
	      min_60m[ch]=pulse_60m[ch][bin];
	      min_time_bin_60m[ch]=bin;
	    }
	}

      for(unsigned int bin=0;bin<500;bin++)
	{ 
	  ped_60m[ch]+=pulse_60m[ch][bin]; 
	}
      ped_60m[ch]=ped_60m[ch]/500.0;
    
      for(unsigned int bin=length_60m-500;bin<length_60m;bin++)
	{
	  pedt_60m[ch]+=pulse_60m[ch][bin]; 
	}
      pedt_60m[ch]=pedt_60m[ch]/500.0;
  
 
      for(unsigned int bin=3500;bin<length_60m;bin++)
	{ 
	  if(pulse_60m[ch][bin]>tmax_60m[ch])
	    {
	      tmax_60m[ch]=pulse_60m[ch][bin];
	      tmax_time_bin_60m[ch]=bin;
	    }
	}
    
      for(int bin=0;bin<omax_time_bin_60m[ch];bin++)
	{ 
	  before_pq_60m[ch]+=(double)pulse_60m[ch][bin];
	}
    
      for(unsigned int bin=omax_time_bin_60m[ch];bin<length_60m;bin++)
	{ 
	  after_pq_60m[ch]+=(double)pulse_60m[ch][bin];
	} 
      
      if(random_bit==0)
	{
	  if((max_time_bin_60m[ch]>=800&&max_time_bin_60m[ch]<=3500))
	    {
	      for(int bin=(max_time_bin_60m[ch]-800);bin<(max_time_bin_60m[ch]+800);bin++)
		{
		  optq_60m[ch]+=(double)pulse_60m[ch][bin];
		}
	    }
	  else
	    {	
	      for(int bin=0;bin<1600;bin++) 
		{
		  optq_60m[ch]+=(double)pulse_60m[ch][bin];
		}
	    }
	}
      else
	{
	  for(int bin=0;bin<1600;bin++) 
	    {
	      optq_60m[ch]+=(double)pulse_60m[ch][bin];
	    }
	}

       for(unsigned int bin=0;bin<length_60m;bin++)
        {
          q2_60m[ch] = q2_60m[ch]
             + ((double)pulse_60m[ch][bin]-(double)ped_60m[ch])*((double)pulse_60m[ch][bin]-(double)ped_60m[ch]);
          t0_60m[ch] = t0_60m[ch]
             + (double)bin*((double)pulse_60m[ch][bin]-(double)ped_60m[ch])*((double)pulse_60m[ch][bin]-(double)ped_60m[ch]);
        }
      t0_60m[ch] = t0_60m[ch]/q2_60m[ch];

      for(unsigned int bin=0;bin<length_60m;bin++)
        {
          T2_60m[ch] = T2_60m[ch]
            + ((double)bin-t0_60m[ch])*((double)bin-t0_60m[ch])*((double)pulse_60m[ch][bin]-(double)ped_60m[ch])*((double)pulse_60m[ch][bin]-(double)ped_60m[ch]);
        }
      T2_60m[ch] = T2_60m[ch]/q2_60m[ch];


    }
}

// NT Variables for 20M1
void event_body::nt_variables_20m1(void)
{
  for(Int_t ch=0;ch<OPEN_CHANNEL_20M1;ch++)
    {
      q_20m1[ch]=0;
      ped_20m1[ch]=0;
      max_20m1[ch]=-100000;
      min_20m1[ch]=100000;
      max_timebin_t1_20m1[ch]=0xFFFF;
      max_timebin_t2_20m1[ch]=0xFFFF; 
      max_timebin_t3_20m1[ch]=0xFFFF;
      
      int th_20m=7000;
      unsigned int t1_temp=0;
      unsigned int  t2_temp=0;

      for(unsigned int bin=0;bin<length_20m1;bin++)
	{
	  q_20m1[ch]+=(double)pulse_20m1[ch][bin];

	  if(pulse_20m1[ch][bin]>max_20m1[ch])
	    {
	      max_20m1[ch]=pulse_20m1[ch][bin];
	      max_time_bin_20m1[ch]=bin;
	    }
      
	  if(pulse_20m1[ch][bin]<min_20m1[ch])
	    {
	      min_20m1[ch]=pulse_20m1[ch][bin];
	      min_time_bin_20m1[ch]=bin;
	    }
	}

      for(unsigned int bin=0;bin<200;bin++)
	{ ped_20m1[ch]+=pulse_20m1[ch][bin]; }
      ped_20m1[ch]=ped_20m1[ch]/200;

      for(unsigned int bin=length_20m1-500;bin<length_20m1;bin++)
	{ pedt_20m1[ch]+=pulse_20m1[ch][bin]; }
      pedt_20m1[ch]=pedt_20m1[ch]/500;
    

      
      for(unsigned int bin=360;bin>0;bin--)
	{
	  if(pulse_20m1[ch][bin]>th_20m&&pulse_20m1[ch][bin-1]<th_20m) 
	    { 
	      max_timebin_t1_20m1[ch]=bin;
	      t1_temp=bin;
	      break;
	    }
      }

      t1_temp = t1_temp-2;
      if(t1_temp>5)
	{
	  for(int bina = t1_temp; bina>0; bina--)
	    {
	      if(pulse_20m1[ch][bina]>th_20m&&pulse_20m1[ch][bina-1]<th_20m) 
		{ 
		  max_timebin_t2_20m1[ch]=bina;
		  t2_temp=bina;
		  break;
		}
	    }
	}
     
      t2_temp = t2_temp-2;
      if(t2_temp>10)
	{
	  for(int bina=t2_temp;bina>0;bina--)
	    {
	      if(pulse_20m1[ch][bina]>th_20m&&pulse_20m1[ch][bina-1]<th_20m) 
		{ 
		  max_timebin_t3_20m1[ch]=bina;
		  break;
		}
	    }
	}
  

    }
}


// NT Variables for 20M2
void event_body::nt_variables_20m2(void)
{
  for(Int_t ch=0;ch<OPEN_CHANNEL_20M2;ch++)
    {
      q_20m2[ch]=0;
      ped_20m2[ch]=0;
      max_20m2[ch]=-100000;
      min_20m2[ch]=100000;
      int th_20m=7000;
      max_timebin_t1_20m2[ch]=0xFFFF;
      max_timebin_t2_20m2[ch]=0xFFFF; 
      max_timebin_t3_20m2[ch]=0xFFFF;
      
      unsigned int t1_temp=0;
      unsigned int  t2_temp=0;

      for(unsigned int bin=0;bin<length_20m2;bin++)
	{
	  q_20m2[ch]+=(double)pulse_20m2[ch][bin];

	  if(pulse_20m2[ch][bin]>max_20m2[ch])
	    {
	      max_20m2[ch]=pulse_20m2[ch][bin];
	      max_time_bin_20m2[ch]=bin;
	    }
      
	  if(pulse_20m2[ch][bin]<min_20m2[ch])
	    {
	      min_20m2[ch]=pulse_20m2[ch][bin];
	      min_time_bin_20m2[ch]=bin;
	    }
	}

      for(unsigned int bin=0;bin<200;bin++)
	{ ped_20m2[ch]+=pulse_20m2[ch][bin]; }
      ped_20m2[ch]=ped_20m2[ch]/200;

      for(unsigned int bin=length_20m2-500;bin<length_20m2;bin++)
	{ pedt_20m2[ch]+=pulse_20m2[ch][bin]; }
      pedt_20m2[ch]=pedt_20m2[ch]/500;
    

      
      for(unsigned int bin=360;bin>0;bin--)
	{
	  if(pulse_20m2[ch][bin]>th_20m&&pulse_20m2[ch][bin-1]<th_20m) 
	    { 
	      max_timebin_t1_20m2[ch]=bin;
	      t1_temp=bin;
	      break;
	    }
      }

      t1_temp = t1_temp-2;
      if(t1_temp>5)
	{
	  for(int bina = t1_temp; bina>0; bina--)
	    {
	      if(pulse_20m2[ch][bina]>th_20m&&pulse_20m2[ch][bina-1]<th_20m) 
		{ 
		  max_timebin_t2_20m2[ch]=bina;
		  t2_temp=bina;
		  break;
		}
	    }
	}
     
      t2_temp = t2_temp-2;
      if(t2_temp>10)
	{
	  for(int bina=t2_temp;bina>0;bina--)
	    {
	      if(pulse_20m2[ch][bina]>th_20m&&pulse_20m2[ch][bina-1]<th_20m) 
		{ 
		  max_timebin_t3_20m2[ch]=bina;
		  break;
		}
	    }
	}
  

    }
}


// NT Variables for 200M1
void event_body::nt_variables_200m1(void)
{
  
  for(Int_t ch=0;ch<OPEN_CHANNEL_200M1;ch++)
    {
      q_200m1[ch]=0;
      ped_200m1[ch]=0;
      pedt_200m1[ch]=0;
      max_200m1[ch]=-100000;
      min_200m1[ch]=100000;
      q2_200m1[ch]=0;
      t0_200m1[ch]=0.0;
      T2_200m1[ch]=0.0;

      for(unsigned int bin=0;bin<length_200m1;bin++)
	{
	  q_200m1[ch]+=(double)pulse_200m1[ch][bin];
	  
	  if(pulse_200m1[ch][bin]>max_200m1[ch])
	    {
	      max_200m1[ch]=pulse_200m1[ch][bin];
	      max_time_bin_200m1[ch]=bin;
	    }
	  
	  if(pulse_200m1[ch][bin]<min_200m1[ch])
	    {
	      min_200m1[ch]=pulse_200m1[ch][bin];
	      min_time_bin_200m1[ch]=bin;
	    }
	}
      
      for(unsigned int bin=0;bin<500;bin++)
	{ped_200m1[ch]+=pulse_200m1[ch][bin]; }
      ped_200m1[ch]=ped_200m1[ch]/500.0;
      
      for(unsigned int bin=length_200m1-500;bin<length_200m1;bin++)
	{ pedt_200m1[ch]+=pulse_200m1[ch][bin]; }
      pedt_200m1[ch]=pedt_200m1[ch]/500.0;
 
      
      for(unsigned int bin=0;bin<length_200m1;bin++)
        {
          q2_200m1[ch] = q2_200m1[ch]
	    + ((double)pulse_200m1[ch][bin]-(double)ped_200m1[ch])*((double)pulse_200m1[ch][bin]-(double)ped_200m1[ch]);
          t0_200m1[ch] = t0_200m1[ch]
	    + (double)bin*((double)pulse_200m1[ch][bin]-(double)ped_200m1[ch])*((double)pulse_200m1[ch][bin]-(double)ped_200m1[ch]);
        }
      t0_200m1[ch] = t0_200m1[ch]/q2_200m1[ch];

      for(unsigned int bin=0;bin<length_200m1;bin++)
        {
          T2_200m1[ch] = T2_200m1[ch]
            + ((double)bin-t0_200m1[ch])*((double)bin-t0_200m1[ch])*((double)pulse_200m1[ch][bin]-(double)ped_200m1[ch])*((double)pulse_200m1[ch][bin]-(double)ped_200m1[ch]);
        }
      T2_200m1[ch] = T2_200m1[ch]/q2_200m1[ch];
      
    }
}


// NT Variables for 200M2
void event_body::nt_variables_200m2(void)
{
  
  for(Int_t ch=0;ch<OPEN_CHANNEL_200M2;ch++)
    {
      q_200m2[ch]=0;
      ped_200m2[ch]=0;
      pedt_200m2[ch]=0;
      max_200m2[ch]=-100000;
      min_200m2[ch]=100000;
      q2_200m2[ch]=0;
      t0_200m2[ch]=0.0;
      T2_200m2[ch]=0.0;

      for(unsigned int bin=0;bin<length_200m2;bin++)
	{
	  q_200m2[ch]+=(double)pulse_200m2[ch][bin];
	  
	  if(pulse_200m2[ch][bin]>max_200m2[ch])
	    {
	      max_200m2[ch]=pulse_200m2[ch][bin];
	      max_time_bin_200m2[ch]=bin;
	    }
	  
	  if(pulse_200m2[ch][bin]<min_200m2[ch])
	    {
	      min_200m2[ch]=pulse_200m2[ch][bin];
	      min_time_bin_200m2[ch]=bin;
	    }
	}
      
      for(unsigned int bin=0;bin<500;bin++)
	{ped_200m2[ch]+=pulse_200m2[ch][bin]; }
      ped_200m2[ch]=ped_200m2[ch]/500.0;
      
      for(unsigned int bin=length_200m2-500;bin<length_200m2;bin++)
	{ pedt_200m2[ch]+=pulse_200m2[ch][bin]; }
      pedt_200m2[ch]=pedt_200m2[ch]/500.0;
 
      
      for(unsigned int bin=0;bin<length_200m2;bin++)
        {
          q2_200m2[ch] = q2_200m2[ch]
	    + ((double)pulse_200m2[ch][bin]-(double)ped_200m2[ch])*((double)pulse_200m2[ch][bin]-(double)ped_200m2[ch]);
          t0_200m2[ch] = t0_200m2[ch]
	    + (double)bin*((double)pulse_200m2[ch][bin]-(double)ped_200m2[ch])*((double)pulse_200m2[ch][bin]-(double)ped_200m2[ch]);
        }
      t0_200m2[ch] = t0_200m2[ch]/q2_200m2[ch];

      for(unsigned int bin=0;bin<length_200m2;bin++)
        {
          T2_200m2[ch] = T2_200m2[ch]
            + ((double)bin-t0_200m2[ch])*((double)bin-t0_200m2[ch])*((double)pulse_200m2[ch][bin]-(double)ped_200m2[ch])*((double)pulse_200m2[ch][bin]-(double)ped_200m2[ch]);
        }
      T2_200m2[ch] = T2_200m2[ch]/q2_200m2[ch];
      
    }
}



void event_body::fit_variable_200m()
{
  // 200m
  chi2_200m =  0.0 ;
  famp_200m = -999.0 ;
  er_famp_200m =  0.0;
  fped_200m = -999.0 ;
  er_fped_200m =  0.0;
  fcross_200m  = -9999.0;
  er_fcross_200m =  0.0;
  fslope_200m  = -9.0;
  er_fslope_200m =  0.0;
  fmid_200m   = -9999.0;
      
      
  if(max_60m[0]>120&&max_60m[2]<15000 && min_60m[0]>-1000.0) 
    {
      int raw_pulse_200m[LENGTH_200M1];
      Int_t x_axis[LENGTH_200M1];
      float xaxis_lakh[LENGTH_200M1];
      
     
      //200M Fit Start 
      //smooth pulse and fit pulse
      int   i,j,m,nl,nr;
      for (i=0; i<NMAX_200M; i++) 
	{
	  ysave_200m[i]=pulse_200m1[0][i]; 
	  signal_200m[i]=pulse_200m1[0][i]; 
	}  
      
      nl=150; nr=150; m=4;                            
      
      index_pt_200m[1]=0;
      j=3;
      for (i=2; i<=nl+1; i++) {
	index_pt_200m[i]=i-j;
	j += 2;
      }
      j=2;
      for (i=nl+2; i<=nl+nr+1; i++) {
	index_pt_200m[i]=i-j;
	j += 2;
      }
      
      savgol_200m(c_200m,nl+nr+1,nl,nr,0,m);
      
      for (i=1; i<=NMAX_200M-nr; i++) 
	{  
	  signal_200m[i] =0.0;
	  
	  if(i>nr&&i<(NMAX_200M-nr))
	    {
	      for (j=1; j<=nl+nr+1; j++)
		{
		  if (i+index_pt_200m[j]>0) 
		    {
		      signal_200m[i] += c_200m[j]*ysave_200m[i+index_pt_200m[j]];
		      
		    }
		  
		}
	    }
	  else
	    {
	      signal_200m[i] = ysave_200m[i];
	    }
	}
      
      
      //Plot 200m pulses
      float fit_ped_200m =0.0;
      float fit_xed_200m =0.0;
      int mid_value_200m = 3000;
      
      
       const float a0_200m = 5.35419e+03; 
       const float a1_200m = 1.10697e-02; 
       const float a2_200m = 1.01078e+03;
      
      
      
      if(max_60m[0]>300)
	{
	  

	  mid_value_200m = (int)(a0_200m+(a1_200m*(max_60m[0]))-(a2_200m*exp(-6.0E-04*max_60m[0])));
	}
      
      if(max_60m[0]<=300)
	{
	  mid_value_200m = 4000;  
	}
      
	
	
	for(Int_t i=0;i<LENGTH_200M1;i++)
	  {
	    if(i<500) 
	      { fit_ped_200m +=signal_200m[i]; }
	    
	    if(i> LENGTH_200M1-600 && i< LENGTH_200M1-100) 
	      { fit_xed_200m +=signal_200m[i]; }
	  }
	fit_ped_200m = fit_ped_200m/(Float_t)500;
	fit_xed_200m = fit_xed_200m/(Float_t)500;
	
	
	for(int bin=0;bin<LENGTH_200M1;bin++)
	  { 
	    x_axis[bin]=bin; 
	    xaxis_lakh[bin]=(float)bin;
	}
      
	for(int i=0;i<LENGTH_200M1;i++)
	  { 
	    raw_pulse_200m[i] = pulse_200m1[0][i]; 
	  }
	
	
      TF1 *smooth_fitfun = new TF1("smooth_fitfun",tanh_fit,0,LENGTH_200M1,4);
      smooth_fitfun->FixParameter(0,(fit_xed_200m-fit_ped_200m));
      smooth_fitfun->FixParameter(1,(0.5*fit_xed_200m + 0.5*fit_ped_200m));
      smooth_fitfun->SetParameter(2,mid_value_200m);
      smooth_fitfun->SetParameter(3,0.0001);
     
     
      TGraph *graph_sfit = new TGraph(LENGTH_200M1,xaxis_lakh,signal_200m);
      graph_sfit->Fit("smooth_fitfun","Q");
     
      TF1 *fit_func_200m = new TF1("fit_func_200m",tanh_fit,0,LENGTH_200M1,4);
      fit_func_200m->SetParNames("Amp_200m","Ped_200m","Cross_200m","Slope_200m");
      
      fit_func_200m->FixParameter(0,smooth_fitfun->GetParameter(0));
      fit_func_200m->FixParameter(1,smooth_fitfun->GetParameter(1));
      fit_func_200m->SetParameter(2,smooth_fitfun->GetParameter(2));
      fit_func_200m->SetParameter(3,smooth_fitfun->GetParameter(3));
     
      TGraph *graph_rfit=new TGraph(LENGTH_200M1,x_axis,raw_pulse_200m);
      graph_rfit->Fit("fit_func_200m","Q");
      
      chi2_200m = (Float_t)fit_func_200m->GetChisquare();
      famp_200m = (Float_t)fit_func_200m->GetParameter(0);
      er_famp_200m =  (Float_t)fit_func_200m->GetParError(0);
      fped_200m = (Float_t)fit_func_200m->GetParameter(1);
      er_fped_200m =  (Float_t)fit_func_200m->GetParError(1);
      fcross_200m  = (Float_t)fit_func_200m->GetParameter(2);
      er_fcross_200m =  (Float_t)fit_func_200m->GetParError(2);
      fslope_200m  = (Float_t)fit_func_200m->GetParameter(3);
      er_fslope_200m = (Float_t)fit_func_200m->GetParError(3);
      fmid_200m   = (Float_t)mid_value_200m;
      
      graph_sfit->Delete();
      graph_rfit->Delete();
      fit_func_200m->Delete(); 
      smooth_fitfun->Delete(); 
    }
  else
    {
      // 200m
      chi2_200m =      -310000;
      famp_200m =      -310000 ;
      er_famp_200m =   -310000;
      fped_200m =      -310000;
      er_fped_200m =   -310000;
      fcross_200m  =   -310000;
      er_fcross_200m = -310000;
      fslope_200m  =   -310000;
      er_fslope_200m = -310000;
      fmid_200m   =    -310000;
    }
}



void event_body::trapez_filter_200m()
{
  trapez_q_200m = 0.0;
  trapez_pq_200m = 0.0;
  trapez_ped_200m = 0.0;
  trapez_pedt_200m = 0.0;
  trapez_max_200m = -99999.0;
  trapez_omax_200m = -99999.0;
  trapez_min_200m = 99999.0;
  trapez_max_time_bin_200m = 0xFFFF;
  trapez_min_time_bin_200m= 0xFFFF;
  trapez_omax_time_bin_200m = 0xFFFF;
  trapez_tprime_200m = 0.0;
  
  trapez_rms_200m = 0.0;
  trapez_local_max_200m = -99999.0;
  trapez_local_max_time_bin_200m = 0xFFFF;
  trapez_local_bfmax_200m = -99999.0;
  trapez_local_bfmax_time_bin_200m = 0xFFFF;
  trapez_local_afmax_200m = -99999.0;
  trapez_local_afmax_time_bin_200m = 0xFFFF;
  
  if(max_60m[0]>0 && min_60m[0]>-1000.0) 
    {
      float input_pulse[LENGTH_200M1] = {0.0};
      for(int bin=0;bin<LENGTH_200M1;bin++) 
	{ 
	  input_pulse[bin] = pulse_200m1[0][bin];
	}
      
      float first_amw_smooth_pulse_2g[LENGTH_200M1]; 
      float amw_smooth_pulse_2g[LENGTH_200M1]; 
      
      
      avmw_filter_lakhs(LENGTH_200M1,200,input_pulse,first_amw_smooth_pulse_2g) ;
      avmw_filter_lakhs(LENGTH_200M1,200,first_amw_smooth_pulse_2g,amw_smooth_pulse_2g) ;
      
      
      const int BIN_FAST_PULSE = 7500;
      //these number, L and G, are optimize for k8.1 fast pulse
      int L = 400;
      int G = 5;
      
      float pulse[BIN_FAST_PULSE] = {0.0};
      float filter[BIN_FAST_PULSE] = {0.0};
      
      for(int bin=0;bin<BIN_FAST_PULSE;bin++) 
	{ 
	  pulse[bin] = amw_smooth_pulse_2g[bin];
	}
      
      double pulse_n, pulse_p;
      for(int bin=0;bin<BIN_FAST_PULSE;bin++)
	{
	  pulse_n = 0;
	  for(int i=(bin-2*L-G+1);i<(bin-L-G+1);i++)
	    {
	      if((i>=0)&&(i<BIN_FAST_PULSE)) { pulse_n = pulse_n + pulse[i]; }
	      if(i<0) { pulse_n = pulse_n + pulse[0]; }
	      if(i>=BIN_FAST_PULSE) { pulse_n = pulse_n + pulse[BIN_FAST_PULSE-1]; }
	    }
	  
	  pulse_p = 0;
	  for(int i=(bin-L+1);i<(bin+1);i++)
	    {
	      if((i>=0)&&(i<BIN_FAST_PULSE)) { pulse_p = pulse_p + pulse[i]; } 
	      if(i<0) { pulse_p = pulse_p + pulse[0]; }
	      if(i>=BIN_FAST_PULSE) { pulse_p = pulse_p + pulse[BIN_FAST_PULSE-1]; }
	    }
	  
	  filter[bin] = (float)(pulse_p-pulse_n)/(float)L;
	}
      
      for(int bin=0;bin<BIN_FAST_PULSE;bin++)
	{
	  trapez_q_200m = trapez_q_200m + filter[bin];
	  if(filter[bin]>trapez_max_200m) 
	    { trapez_max_200m = filter[bin]; 
	      trapez_max_time_bin_200m =(int)bin; 
	    }
	  
	  if(filter[bin]>trapez_omax_200m&&(bin>=800)&&(bin<3000)) 
	    { trapez_omax_200m = filter[bin]; 
	      trapez_omax_time_bin_200m = (int)bin; 
	    }
	  if(filter[bin]<trapez_min_200m) 
	    { trapez_min_200m = filter[bin]; 
	      trapez_min_time_bin_200m = (int)bin; 
	    }
	}
      
      for(int bin=trapez_omax_time_bin_200m-800; bin<trapez_omax_time_bin_200m+800;bin++)
	{
	  trapez_pq_200m = trapez_pq_200m + filter[bin];
	}
            
      
      for(unsigned int bin=0;bin<500;bin++)
	{
	  trapez_ped_200m += filter[bin]; 
	}
      trapez_ped_200m = trapez_ped_200m/500.0;
      
      for(unsigned int bin=length_200m1-500;bin<length_200m1;bin++)
	{ 
	  trapez_pedt_200m += filter[bin]; 
	}
      trapez_pedt_200m = trapez_pedt_200m/500.0;
      
      for(int bin=0;bin<BIN_FAST_PULSE;bin++)
	{
	  trapez_tprime_200m = trapez_tprime_200m + (double)bin*filter[bin];
	}
      
      trapez_tprime_200m = trapez_tprime_200m/trapez_q_200m;
     
      double sq =0;
      
      for(int bin=0;bin<BIN_FAST_PULSE;bin++)
	{
	  sq = sq + filter[bin]*filter[bin];
	}
      if(((trapez_q_200m*trapez_q_200m)-sq)>0)
	{
	  trapez_rms_200m = (sqrt((trapez_q_200m*trapez_q_200m)-sq)/double(BIN_FAST_PULSE));
	}
      else
	{
	  trapez_rms_200m = -1.0*(sqrt(sq-(trapez_q_200m*trapez_q_200m))/double(BIN_FAST_PULSE));	  
	}

      
      for(int bin=0;bin<BIN_FAST_PULSE;bin++)
	{
	  if(filter[bin]>trapez_local_max_200m && filter[bin-2]<filter[bin] && filter[bin+2]<filter[bin]  ) 
	    { trapez_local_max_200m = filter[bin]; 
	      trapez_local_max_time_bin_200m =(int)bin; 
	    }
	}
      
      if(trapez_local_max_time_bin_200m>100 && trapez_local_max_time_bin_200m<7400)
	{
	  
	  for(int bin=trapez_local_max_time_bin_200m+90;bin<BIN_FAST_PULSE-10;bin++)
	    {
	      if(filter[bin]>trapez_local_afmax_200m &&filter[bin-2]<filter[bin] && filter[bin+2]<filter[bin]  ) 
		{ trapez_local_afmax_200m = filter[bin]; 
		  trapez_local_afmax_time_bin_200m =(int)bin; 
		}
	    }
	  
	  for(int bin = trapez_local_max_time_bin_200m-90;  bin > 3; bin--)
	    {
	      if(filter[bin]>trapez_local_bfmax_200m &&filter[bin-2]<filter[bin] && filter[bin+2]<filter[bin] ) 
		{ 
		  trapez_local_bfmax_200m = filter[bin]; 
		  trapez_local_bfmax_time_bin_200m =bin; 
		}
	    }
	}
      else
	{
	  for(int bin=trapez_omax_time_bin_200m;bin<BIN_FAST_PULSE-10;bin++)
	    {
	      if(filter[bin]>trapez_local_afmax_200m &&filter[bin-2]<filter[bin] && filter[bin+2]<filter[bin]  ) 
		{ trapez_local_afmax_200m = filter[bin]; 
		  trapez_local_afmax_time_bin_200m =(int)bin; 
		}
	    }
	  
	  for(int bin = trapez_omax_time_bin_200m;  bin > 3; bin--)
	    {
	      if(filter[bin]>trapez_local_bfmax_200m &&filter[bin-2]<filter[bin] && filter[bin+2]<filter[bin] ) 
		{ 
		  trapez_local_bfmax_200m = filter[bin]; 
		  trapez_local_bfmax_time_bin_200m =bin; 
		}
	    }
	}
      
      if(trapez_local_bfmax_200m == -99999)  
	{
	  trapez_local_bfmax_200m = filter[0];
	  trapez_local_bfmax_time_bin_200m = 0; 
	}  
	
      if(trapez_local_afmax_200m == -99999)  
	{
	  trapez_local_afmax_200m = filter[7450];
	  trapez_local_afmax_time_bin_200m = 7450; 
	}
    }
  else
    {
      trapez_q_200m =    -310000.0;
      trapez_pq_200m =   -310000.0;
      trapez_ped_200m =  -310000.0;
      trapez_pedt_200m = -310000.0;
      trapez_max_200m =  -310000.0;
      trapez_omax_200m = -310000.0;
      trapez_min_200m =  -310000.0 ;
      trapez_max_time_bin_200m = -310000.0;
      trapez_min_time_bin_200m= -310000.0;
      trapez_omax_time_bin_200m = -310000.0;
      trapez_tprime_200m =    -310000.0;
      trapez_rms_200m = -310000.0;
      trapez_local_max_200m = -310000.0;
      trapez_local_max_time_bin_200m = -310000.0;
      trapez_local_bfmax_200m = -310000.0;
      trapez_local_bfmax_time_bin_200m = -310000.0;
      trapez_local_afmax_200m = -310000.0;
      trapez_local_afmax_time_bin_200m = -310000.0;
      
    }
  
}

void event_body::fftw_variable_60m()
{

  double *fft_in;
  double *fft_result;
  fftw_complex *fft_out;
  fftw_complex *fft_cut;
  fftw_plan fft_forward;
  //fftw_plan fft_backward;

  fft_in = (double*)fftw_malloc(sizeof(double)*length_60m);
  fft_result  = (double*)fftw_malloc(sizeof(double)*length_60m);
  fft_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*length_60m);
  fft_cut  = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*length_60m);

  fft_W2_60m = 0.0;
  fft_q2_60m = 0.0;
  fft_max_60m = 0.0;
  fft_max_bin_60m = 0;
  fft1_W2_60m = 0.0;
  fft1_q2_60m = 0.0;
  fft1_max_60m = 0.0;
  fft1_max_bin_60m = 0;
  fft100_W2_60m = 0.0;
  fft100_q2_60m = 0.0;
  fft100_max_60m = 0.0;
  fft100_max_bin_60m = 0;
 
  for(unsigned int bin=0;bin<length_60m;bin++) { fft_in[bin] = (double)pulse_60m[0][bin]; }

  fft_forward = fftw_plan_dft_r2c_1d(length_60m,fft_in,fft_out,FFTW_ESTIMATE);
  fftw_execute(fft_forward);

  for(unsigned int bin=0;bin<length_60m;bin++)
    {
      if(fft_max_60m<(fft_out[bin][0]*fft_out[bin][0]+fft_out[bin][1]*fft_out[bin][1]))
        { 
          fft_max_60m = (fft_out[bin][0]*fft_out[bin][0]+fft_out[bin][1]*fft_out[bin][1]);
          fft_max_bin_60m = (double)bin;
        }

      if((bin!=0)&&(fft1_max_60m<(fft_out[bin][0]*fft_out[bin][0]+fft_out[bin][1]*fft_out[bin][1])))
        {
          fft1_max_60m = (fft_out[bin][0]*fft_out[bin][0]+fft_out[bin][1]*fft_out[bin][1]);
          fft1_max_bin_60m = (double)bin;
        }
        
      if((bin>=100)&&(fft100_max_60m<(fft_out[bin][0]*fft_out[bin][0]+fft_out[bin][1]*fft_out[bin][1])))
        {
          fft100_max_60m = (fft_out[bin][0]*fft_out[bin][0]+fft_out[bin][1]*fft_out[bin][1]);
          fft100_max_bin_60m = (double)bin;
        }

      fft_q2_60m = fft_q2_60m + (fft_out[bin][0]*fft_out[bin][0]+fft_out[bin][1]*fft_out[bin][1]);
      fft_W2_60m = fft_W2_60m + (double)(bin*bin)*(fft_out[bin][0]*fft_out[bin][0]+fft_out[bin][1]*fft_out[bin][1]);

      if(bin!=0)
        {
          fft1_q2_60m = fft1_q2_60m + (fft_out[bin][0]*fft_out[bin][0]+fft_out[bin][1]*fft_out[bin][1]);
          fft1_W2_60m = fft1_W2_60m + (double)(bin*bin)*(fft_out[bin][0]*fft_out[bin][0]+fft_out[bin][1]*fft_out[bin][1]);
        }

      if(bin>=100)
        {
          fft100_q2_60m = fft100_q2_60m + (fft_out[bin][0]*fft_out[bin][0]+fft_out[bin][1]*fft_out[bin][1]);
          fft100_W2_60m = fft100_W2_60m + (double)(bin*bin)*(fft_out[bin][0]*fft_out[bin][0]+fft_out[bin][1]*fft_out[bin][1]);
        }
    } 
      
  fft_W2_60m = fft_W2_60m/fft_q2_60m;
  fft1_W2_60m = fft1_W2_60m/fft1_q2_60m;
  fft100_W2_60m = fft100_W2_60m/fft100_q2_60m;

  fftw_destroy_plan(fft_forward);
    

      

  fftw_free(fft_in);
  fftw_free(fft_out);
  fftw_free(fft_result);
  fftw_free(fft_cut);
}


void event_body::fftw_variable_200m()
{

  double *fft_in;
  double *fft_result;
  fftw_complex *fft_out;
  fftw_complex *fft_cut;
  fftw_plan fft_forward;


  fft_in = (double*)fftw_malloc(sizeof(double)*LENGTH_200M1);
  fft_result  = (double*)fftw_malloc(sizeof(double)*LENGTH_200M1);
  fft_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*LENGTH_200M1);
  fft_cut  = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*LENGTH_200M1);

  
  fft_W2_200m = 0.0;
  fft_q2_200m = 0.0;
  fft_max_200m = 0.0;
  fft_max_bin_200m = 0;

  fft1_W2_200m = 0.0;
  fft1_q2_200m = 0.0;
  fft1_max_200m = 0.0;
  fft1_max_bin_200m = 0;

  fft100_W2_200m = 0.0;
  fft100_q2_200m = 0.0;
  fft100_max_200m = 0.0;
  fft100_max_bin_200m = 0;

  for(unsigned int bin=0;bin<LENGTH_200M1;bin++) { fft_in[bin] = (double)pulse_200m1[0][bin]; }

  fft_forward = fftw_plan_dft_r2c_1d(LENGTH_200M1,fft_in,fft_out,FFTW_ESTIMATE);
  fftw_execute(fft_forward);

  for(unsigned int bin=0;bin<LENGTH_200M1;bin++)
    {
      if(fft_max_200m<(fft_out[bin][0]*fft_out[bin][0]+fft_out[bin][1]*fft_out[bin][1]))
        { 
          fft_max_200m = (fft_out[bin][0]*fft_out[bin][0]+fft_out[bin][1]*fft_out[bin][1]);
          fft_max_bin_200m = (double)bin;
        }

      if((bin!=0)&&(fft1_max_200m<(fft_out[bin][0]*fft_out[bin][0]+fft_out[bin][1]*fft_out[bin][1])))
        {
          fft1_max_200m = (fft_out[bin][0]*fft_out[bin][0]+fft_out[bin][1]*fft_out[bin][1]);
          fft1_max_bin_200m = (double)bin;
        }

      if((bin>=100)&&(fft100_max_200m<(fft_out[bin][0]*fft_out[bin][0]+fft_out[bin][1]*fft_out[bin][1])))
        {
          fft100_max_200m = (fft_out[bin][0]*fft_out[bin][0]+fft_out[bin][1]*fft_out[bin][1]);
          fft100_max_bin_200m = (double)bin;
        }

      fft_q2_200m = fft_q2_200m + (fft_out[bin][0]*fft_out[bin][0]+fft_out[bin][1]*fft_out[bin][1]);
      fft_W2_200m = fft_W2_200m + (double)(bin*bin)*(fft_out[bin][0]*fft_out[bin][0]+fft_out[bin][1]*fft_out[bin][1]);

      if(bin!=0)
        {
          fft1_q2_200m = fft1_q2_200m + (fft_out[bin][0]*fft_out[bin][0]+fft_out[bin][1]*fft_out[bin][1]);
          fft1_W2_200m = fft1_W2_200m + (double)(bin*bin)*(fft_out[bin][0]*fft_out[bin][0]+fft_out[bin][1]*fft_out[bin][1]);
        }

      if(bin>=100)
        {
          fft100_q2_200m = fft100_q2_200m + (fft_out[bin][0]*fft_out[bin][0]+fft_out[bin][1]*fft_out[bin][1]);
          fft100_W2_200m = fft100_W2_200m + (double)(bin*bin)*(fft_out[bin][0]*fft_out[bin][0]+fft_out[bin][1]*fft_out[bin][1]);
        }
    }

  fft_W2_200m = fft_W2_200m/fft_q2_200m;
  fft1_W2_200m = fft1_W2_200m/fft1_q2_200m;
  fft100_W2_200m = fft100_W2_200m/fft100_q2_200m;


  fftw_destroy_plan(fft_forward);
      
      
  fftw_free(fft_in);
  fftw_free(fft_out);
  fftw_free(fft_result);
  fftw_free(fft_cut);
}
