#define OPEN_CHANNEL_2G1 1
#define OPEN_CHANNEL_2G2 1
#define OPEN_CHANNEL_200M1 1 
#define OPEN_CHANNEL_200M2 2
#define OPEN_CHANNEL_20M1 7
#define OPEN_CHANNEL_20M2 5
#define OPEN_CHANNEL_60M 4
#define LENGTH_2G1 31000
#define LENGTH_2G2 31000
#define LENGTH_200M1 7500
#define LENGTH_200M2 12000
#define LENGTH_20M1 1000
#define LENGTH_20M2 1000
#define LENGTH_60M 5760


UInt_t length_2g1,length_2g2,length_200m1,length_200m2;
UInt_t length_20m1,length_20m2,length_60m;
UInt_t open_ch_2g1,open_ch_2g2,open_ch_200m1,open_ch_200m2;
UInt_t open_ch_20m1,open_ch_20m2,open_ch_60m;
UInt_t total_eheader_size,total_ebody_size,block_size;

ULong64_t pointerposition;
ULong64_t last_start_time;
ULong64_t last_end_time;
ULong64_t FPGA_reference_time1;
ULong64_t FPGA_reference_time2;
ULong64_t FPGA_reference_time3;


class run_header
{
 public:

  Double_t program_start_time;

  Int_t pos_triger_2g1;
  Int_t pre_triger_2g1;
  Int_t open_channel_2g1;
  Int_t on_off_2g1;

  Int_t pos_triger_200M1;
  Int_t pre_triger_200M1;
  Int_t open_channel_200M1;
  Int_t on_off_200M1;

  Int_t pos_triger_2g2;
  Int_t pre_triger_2g2;
  Int_t open_channel_2g2;
  Int_t on_off_2g2;

  Int_t pos_triger_200M2;
  Int_t pre_triger_200M2;
  Int_t open_channel_200M2;
  Int_t on_off_200M2;

  Int_t pos_triger_60M;
  Int_t pre_triger_60M;
  Int_t open_channel_60M;
  Int_t on_off_60M;

  Int_t pos_triger_20M1;
  Int_t pre_triger_20M1;
  Int_t open_channel_20M1;
  Int_t on_off_20M1;

  Int_t pos_triger_20M2;
  Int_t pre_triger_20M2;
  Int_t open_channel_20M2;
  Int_t on_off_20M2;

  Double_t run_start_time;

  //  ULong64_t FPGA_clock;

  ULong64_t FPGA_clock1;
  ULong64_t FPGA_clock2;
  ULong64_t FPGA_clock3;

  void read_rheader(FILE *raw);
  Int_t ascii_to_binary(Int_t bit_num,Int_t num);

}rheader;

class event_header
{
 public:

  ULong64_t event_start_time1;
  ULong64_t event_start_time2;
  ULong64_t event_start_time3;
  ULong64_t idevt;

  ULong64_t veto_tdc[32];
  ULong64_t former_veto_tdc[32];
  ULong64_t late_former_veto_tdc[32];
  ULong64_t late_former_veto_tdc2[32];

  ULong64_t veto_tdc_or[13];
  ULong64_t former_veto_tdc_or[13];
  ULong64_t late_former_veto_tdc_or[13];
  ULong64_t late_former_veto_tdc2_or[13];

  ULong64_t signal_inhibit_time[8];
  ULong64_t hitpattern;
  ULong64_t upper_veto_hp;
  ULong64_t other_veto_hp;
  ULong64_t other_inner_veto_hp;

  ULong64_t event_end_time;

  Short_t mins_OV[16];
  Short_t mins_OV2[8];
  
  Int_t pulser_on_off;
  Int_t ge_on_off;
  Int_t random_trig_on_off;
  Int_t nai_1_on_off;
  Int_t nai_2_on_off;
  Int_t nai_3_on_off;
  Int_t nai_5_on_off;
  Int_t upper_veto_on_off[32];
  Int_t other_veto_on_off[32];
  Int_t other_inner_on_off[32];

  void read_eheader(FILE *raw);
  void eh_variables(Int_t run,Int_t initial_run,Int_t block);
  unsigned long filesize(FILE *raw);

}eheader;

class event_body
{
 public:

  Char_t pulse_2g1[OPEN_CHANNEL_2G1][LENGTH_2G1];
  Short_t pulse_60m[OPEN_CHANNEL_60M][LENGTH_60M];
  Short_t pulse_20m1[OPEN_CHANNEL_20M1][LENGTH_20M1];
  Short_t pulse_20m2[OPEN_CHANNEL_20M2][LENGTH_20M2];
  Short_t pulse_200m1[OPEN_CHANNEL_200M1][LENGTH_200M1];
  Short_t pulse_200m2[OPEN_CHANNEL_200M2][LENGTH_200M2];
  Char_t pulse_2g2[OPEN_CHANNEL_2G2][LENGTH_2G2];

  
    // 60M N-Tuple Variables
  double q_60m[OPEN_CHANNEL_60M]; 
  int max_60m[OPEN_CHANNEL_60M];
  int max_time_bin_60m[OPEN_CHANNEL_60M];
  int omax_60m[OPEN_CHANNEL_60M];
  int omax_time_bin_60m[OPEN_CHANNEL_60M];
  int tmax_60m[OPEN_CHANNEL_60M];
  int tmax_time_bin_60m[OPEN_CHANNEL_60M];
  int min_60m[OPEN_CHANNEL_60M];
  int min_time_bin_60m[OPEN_CHANNEL_60M];
  float ped_60m[OPEN_CHANNEL_60M];
  float pedt_60m[OPEN_CHANNEL_60M];
  double before_pq_60m[OPEN_CHANNEL_60M];
  double after_pq_60m[OPEN_CHANNEL_60M];
  double pq800_3500_60m[OPEN_CHANNEL_60M];
  double optq_60m[OPEN_CHANNEL_60M];
  double q2_60m[OPEN_CHANNEL_60M];
  double t0_60m[OPEN_CHANNEL_60M];
  double T2_60m[OPEN_CHANNEL_60M];

  // 200M 1 N-Tuple Variables
  double q_200m1[OPEN_CHANNEL_200M1]; 
  float ped_200m1[OPEN_CHANNEL_200M1];
  float pedt_200m1[OPEN_CHANNEL_200M1];
  int max_200m1[OPEN_CHANNEL_200M1];
  int min_200m1[OPEN_CHANNEL_200M1];
  int max_time_bin_200m1[OPEN_CHANNEL_200M1];
  int min_time_bin_200m1[OPEN_CHANNEL_200M1];
  double q2_200m1[OPEN_CHANNEL_200M1]; 
  double t0_200m1[OPEN_CHANNEL_200M1];
  double T2_200m1[OPEN_CHANNEL_200M1];

  // 200M 2 N-Tuple Variables
  double q_200m2[OPEN_CHANNEL_200M2]; 
  float ped_200m2[OPEN_CHANNEL_200M2];
  float pedt_200m2[OPEN_CHANNEL_200M2];
  int max_200m2[OPEN_CHANNEL_200M2];
  int min_200m2[OPEN_CHANNEL_200M2];
  int max_time_bin_200m2[OPEN_CHANNEL_200M2];
  int min_time_bin_200m2[OPEN_CHANNEL_200M2];
  double q2_200m2[OPEN_CHANNEL_200M2]; 
  double t0_200m2[OPEN_CHANNEL_200M2];
  double T2_200m2[OPEN_CHANNEL_200M2];

  // 20M 1 N-Tuple Variables
  double q_20m1[OPEN_CHANNEL_20M1]; 
  float ped_20m1[OPEN_CHANNEL_20M1];
  float pedt_20m1[OPEN_CHANNEL_20M1];
  int max_20m1[OPEN_CHANNEL_20M1];
  int min_20m1[OPEN_CHANNEL_20M1];
  int max_time_bin_20m1[OPEN_CHANNEL_20M1];
  int min_time_bin_20m1[OPEN_CHANNEL_20M1];
  int max_timebin_t1_20m1[OPEN_CHANNEL_20M1];
  int max_timebin_t2_20m1[OPEN_CHANNEL_20M1];
  int max_timebin_t3_20m1[OPEN_CHANNEL_20M1];
   
  // 20M 2 N-Tuple Variables
  double q_20m2[OPEN_CHANNEL_20M2]; 
  float ped_20m2[OPEN_CHANNEL_20M2];
  float pedt_20m2[OPEN_CHANNEL_20M2];
  int max_20m2[OPEN_CHANNEL_20M2];
  int min_20m2[OPEN_CHANNEL_20M2];
  int max_time_bin_20m2[OPEN_CHANNEL_20M2];
  int min_time_bin_20m2[OPEN_CHANNEL_20M2];
  int max_timebin_t1_20m2[OPEN_CHANNEL_20M2];
  int max_timebin_t2_20m2[OPEN_CHANNEL_20M2];
  int max_timebin_t3_20m2[OPEN_CHANNEL_20M2];
   
  // variables_fitting_200m
  Float_t chi2_200m;
  Float_t famp_200m;
  Float_t er_famp_200m;
  Float_t fped_200m; 
  Float_t er_fped_200m;
  Float_t fcross_200m;
  Float_t er_fcross_200m;
  Float_t fslope_200m;
  Float_t er_fslope_200m;
  Float_t fmid_200m;

  // variables_trapez_filter_200m
   double trapez_q_200m; 
   double trapez_pq_200m; 
   float trapez_ped_200m;
   float trapez_pedt_200m;
   float trapez_max_200m;
   float trapez_min_200m;
   float trapez_omax_200m;
   int trapez_max_time_bin_200m;
   int trapez_min_time_bin_200m;
   int trapez_omax_time_bin_200m;
   double trapez_tprime_200m; 
   double trapez_rms_200m;
   float trapez_local_max_200m;
   int trapez_local_max_time_bin_200m;
   float trapez_local_bfmax_200m;
   int trapez_local_bfmax_time_bin_200m;
   float trapez_local_afmax_200m;
   int trapez_local_afmax_time_bin_200m;
  
   //fftw
   Float_t fft_W2_60m;
   Float_t fft_q2_60m;
   Float_t fft_max_60m;
   Float_t fft_max_bin_60m;
   Float_t fft1_W2_60m;
   Float_t fft1_q2_60m;
   Float_t fft1_max_60m;
   Float_t fft1_max_bin_60m;
   Float_t fft100_W2_60m;
   Float_t fft100_q2_60m;
   Float_t fft100_max_60m;
   Float_t fft100_max_bin_60m;
   
   Float_t fft_W2_200m;
   Float_t fft_q2_200m;
   Float_t fft_max_200m;
   Float_t fft_max_bin_200m;
   Float_t fft1_W2_200m;
   Float_t fft1_q2_200m;
   Float_t fft1_max_200m;
   Float_t fft1_max_bin_200m;
   Float_t fft100_W2_200m;
   Float_t fft100_q2_200m;
   Float_t fft100_max_200m;
   Float_t fft100_max_bin_200m;



  void read_ebody(FILE *raw);
  void nt_variables_60m(int random_bit);
  void nt_variables_200m1(void);
  void nt_variables_200m2(void);
  void nt_variables_20m1(void);
  void nt_variables_20m2(void);
  void fit_variable_200m(void); 
  void trapez_filter_200m(void);
  void fftw_variable_60m(void);
  void fftw_variable_200m(void);
  void diff0_variable_200m(void);

}ebody;

class run_tail
{
 public:

  UInt_t scalar[107];
  Double_t run_end_time;

  void read_rtail(FILE *raw);

}rtail;


