
#include <iostream>
#include <fstream>

//#define  HARD_INPUT
#include <immintrin.h>
#include <malloc.h>
#include <math.h>
#include <ctime>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>

#ifndef NO_VERILATOR
#define DUMP_VCD
#include "svdpi.h"
//#include "Vtb_top_verilator__Dpi.h"
#include "Vldpc_wsh.h"
#ifdef DUMP_VCD
    #include <verilated_vcd_c.h>
#endif
#include "verilated.h"
#include <iomanip>
#include <exception>
#include <cstdio>
#include <cstdint>
#include <cerrno>

static vluint64_t main_time = 0;
#ifdef DUMP_VCD
VerilatedVcdC* tfp;
#endif
Vldpc_wsh     *top;
#endif

#define RSER         0.22  
#define RSCR         0.91   

#ifndef __WSH__
#define __WSH__

#include <stdio.h>

int ms2p_avx(char lratio[], char *dec_bit, int &cycle, int &num_dec) ;
void encoder(const char *g_name, unsigned long long *cw_bit, int mode_sel) ;

#endif  // __WSH__

extern int  var_N   ;
extern int  check_M ;
extern int  max_iter ;
extern char  **llrmap; // 2b or 3b
// ----------------
//
// ----------------
int   min_loop     = 1000 ;
int   min_error    = 300  ;
int   test_per_msg = 10000 ;

float sn_start     = 5.20;
float sn_end       = 14.95;
float sn_step      = 0.10;
int   binary_flip  ;
int   seflp        ;
int   scflp        ;
long  *idum        ;
long  longi        ;
double  var        ;

long  snr_runs     = 6000000000;
int    max_hist    = 500; 

int   *cnt_cor    ;
int   *cnt_uncor  ;  
int   *cnt_alias  ;
// ----------------
//
// ----------------
using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
double ran2(long *idum)
{
    int j;
    long k;
    static long idum2=123456789;
    static long iy=0;
    static long iv[NTAB];
    float temp;
    if (*idum <= 0) { 
        if (-(*idum) < 1) *idum=1; 
        else *idum = -(*idum);
        idum2 = (*idum);
        for (j=NTAB+7;j>=0;j--) { 
            k = (*idum)/IQ1;
            *idum = IA1*(*idum-k*IQ1)-k*IR1;
            if (*idum < 0) *idum += IM1;
            if (j < NTAB) iv[j] = *idum;
        }
        iy = iv[0];
    }
    k = (*idum)/IQ1;
    *idum = IA1*(*idum-k*IQ1)-k*IR1;
    if (*idum < 0) *idum += IM1; 
    
    k = idum2/IQ2;
    idum2 = IA2*(idum2-k*IQ2)-k*IR2;
    if (idum2 < 0) idum2 += IM2;
    
    j  = iy/NDIV; 
    iy = iv[j]-idum2; 
    iv[j] = *idum; 
    if (iy < 1) iy += IMM1;
    if ((temp=AM*iy) > RNMX) return RNMX; 
    else return temp;
}


void normal(long *idum, double var, double &n1, double &n2 )
{
    double x1, x2, s, t;
    do{
        x1 = ran2(idum);
        x2 = ran2(idum);
        x1 = 2 * x1 - 1;
        x2 = 2 * x2 - 1;
        s = x1 * x1 + x2 * x2;
    }while( s >= 1.0 );
    t = ( -2 * log(s) )/s;
    t = sqrt(t);
    n1 = var*x1*t;
    n2 = var*x2*t;
}

void write_report(ofstream & outFile, double sn, int errbit_cnt, int block_err, long i, int flag)
{
    int    j;
    long   all;
    double avg_hard_error;
    double block_err_rate;
    double out_BER;
    
    time_t now = time(0);
    tm *ltm = localtime(&now);
    outFile << "now is "   << 1900+ltm->tm_year<<"/"<< 1+ltm->tm_mon <<"/";
    outFile <<ltm->tm_mday <<" "<<ltm->tm_hour <<":"<<ltm->tm_min    <<":"<<ltm->tm_sec<< endl;      
    
    outFile << "NUM_VAR       = " << var_N          << " ;" << endl;
    outFile << "NUM_CHK       = " << check_M        << " ;" << endl;
    outFile << "max_iter      = " << max_iter       << " ;" << endl;
    outFile << "SNR           = " << sn             << " ;(dB)" << endl;
    outFile << "Error_Bit     = " << errbit_cnt     << " ;"<< endl;
    outFile << "Error_Block   = " << block_err      << " ;"<< endl;
    outFile << "Total_Block   = " << i              << " ;"<< endl;
    for( j = 0, all = 0; j < max_hist; j ++ ) 
    {
        if( cnt_cor[j] >= all )  
        {
            avg_hard_error = j+1 ;
            all = cnt_cor[j] ;
        }
    }
    outFile << "Average_Error = " << avg_hard_error << " ; // " << binary_flip << "/65535" ;
    outFile << endl ;
    
    block_err_rate = static_cast<double>(block_err)/static_cast<double>(i);
    out_BER        = static_cast<double>(errbit_cnt)/static_cast<double>(i*var_N);
    
    if(block_err > 0)   outFile << "BLER          = "<< (block_err_rate) << " ; " << 10*log10(block_err_rate) << "dB" <<endl;
    else                outFile << "BLER = " << (block_err_rate) << " ;" << endl;
    
    //outFile << "out_BER = " << out_BER << " ;" << endl;
    outFile << "  ----  Correctable initial Error Count ----- " << endl;
    for( j = 0; j < max_hist; j ++ ) 
    {
        if( cnt_cor[j] !=0 )       outFile << "Correctable_Error_" << j+1 <<" = "<< cnt_cor[j] <<" ;" << endl;
    }
    outFile << "  ----  Uncorrectable initial Error Count --- " << endl;
    for( j = 0; j < max_hist; j ++ ) 
    {
        if( cnt_uncor[j] !=0 )     outFile << "Uncorrectable_Error_" << j+1 <<" = "<< cnt_uncor[j] <<" ;" << endl;
    }    
    outFile << "------------------------------------------------------- " << endl;
    if(flag==1) 
    {
        cout << "SNR           = " << sn               << " dB" << endl;
        cout << "Average_Error = " << avg_hard_error   << " ;"<< endl;   
        cout << "BLER          = " << (block_err_rate) <<" = "<< block_err << "/" << i << ", "<< errbit_cnt << endl;    
        //cout << "out_BER       = " << out_BER << " ;" << endl;    
    }
}

void bpsk_awgn( __m256i m256z[], long *idum, double var, int *tx_code, char *lratio, int &num_awgn)
{
    #ifndef HARD_INPUT
    const char ch_idx = 1;   // default 
    #else
    const char ch_idx = 2;  
    #endif
    //
    //  if change ch_idx to 2 when hard-input decoding, 
    //      decoder need to change tuning parameter for extra correction gain
    //
    
    int i, sign;
    double n1, n2, a1, a2;
    num_awgn = 0;
    for( i = 0; i < var_N/2; i++ )
    {
        a1 = (tx_code[2*i  ]) ? -1 : 1;
        a2 = (tx_code[2*i+1]) ? -1 : 1;        
        normal( idum, var, n1, n2 );
        
        sign = ((a1+n1)>=0) ? 1 : -1;	
		lratio[2*i  ] = ( ((a1+n1)>0.5)||((a1+n1)<-0.5) )? sign*2 : sign*ch_idx;

        sign = ((a2+n2)>=0) ? 1 : -1;
		lratio[2*i+1] = ( ((a2+n2)>0.5)||((a2+n2)<-0.5) )? sign*2 : sign*ch_idx;
		
        if( (lratio  [2*i  ]<0)^(tx_code[2*i  ]==1) )  num_awgn++;
        if( (lratio  [2*i+1]<0)^(tx_code[2*i+1]==1) )  num_awgn++;

    }

    if(num_awgn>=max_hist)
    {
        num_awgn = max_hist-1;
        printf("Warning: num_awgn is larger than %d\n", max_hist);
    }  
}

unsigned long lfsr113(void) {
  // Generates random 32 bit numbers.
  // NOTE: the seed MUST satisfy
  // z1 > 1, z2 > 7, z3 > 15, and z4 > 127 
  static unsigned long z1 = 0x1234;
  static unsigned long z2 = 0x5678;
  static unsigned long z3 = 0x9abc;
  static unsigned long z4 = 0xdef0;
  unsigned long b;
  b  = (((z1 <<  6) ^ z1)  >> 13);
  z1 = (((z1 & 4294967294) << 18) ^ b);
  b  = (((z2 <<  2) ^ z2)  >> 27);
  z2 = (((z2 & 4294967288) <<  2) ^ b);
  b  = (((z3 << 13) ^ z3)  >> 21);
  z3 = (((z3 & 4294967280) <<  7) ^ b);
  b  = (((z4 <<  3) ^ z4)  >> 12);
  z4 = (((z4 & 4294967168) << 13) ^ b);
  return (z1 ^ z2 ^ z3 ^ z4);
}

void info_source(unsigned long long * cw_bit, int info_length) 
{
  int i, j, q;
  unsigned long long ran_num, qc;
  q  = info_length / 16;
  qc = 0 ;
  for (i = 0; i < q; i++) 
  {
    ran_num = lfsr113() % 0xffff;
    qc = ( ran_num << 48 )|( qc >> 16 );
    
    if( (i%4)==3 )
    {    
        cw_bit[ (int)(i/4)+(int)(check_M/64) ] = qc ; 
    }
  }   
  // parity is in front of info bits  
}

unsigned long long * cw_bit ; 
#ifndef NO_VERILATOR
void tick()
{
    top->clk = 0;
    top->eval();
    #ifdef DUMP_VCD
    tfp->dump(main_time);
    #endif
    main_time++;

    top->clk = 1;
    top->eval();
    #ifdef DUMP_VCD
    tfp->dump(main_time);
    #endif
    main_time++;
}

void top_init(int mode_select)
{
    // ===== Reset =====
    top->rst_n = 0;
    top->clear        = 0;
    top->enc0_h1s2    = 0;
    top->which_code   = mode_select;
    top->max_iter     = max_iter; 
    top->llr_sb       = 0xEFEF2121;   // only 1,2,3,F,E,D available
    top->hb_in        = 0;
    for(int ii=0; ii<4; ii++)   top->sb_in[ii]     = 0;
    top->info_q       = 0;
    for(int ii=0; ii<8; ii++)   top->sign4d_q[ii]  = 0;
    top->enc_vdin     = 0;
    top->enc_din      = 0;
    tick();
    tick();
    top->rst_n = 1;
    tick();  
}

void top_enc(unsigned long long *cw_bit, int mode_select)
{
    top->clear = 1;
    tick();  
    top->clear = 0;
    tick();  
    printf("num_var = %d\n", top->num_var);
    
    for(int ii=0, adr=0; ii<=var_N; ii+=64)
    {
        adr = (int)(((ii+check_M)%var_N)/64);
        top->enc_vdin = 1;
        top->enc_din  = (ii<(var_N-check_M)) ? (uint64_t) cw_bit[adr] : (uint64_t) 0;
        //printf("%d, %lx  \n", ii, top->enc_din);
        if((top->is_parity==1)&&(top->enc_vdout==1))
        {

            if(top->enc_dout!=cw_bit[adr-1])      printf("Mismatch Encoder !");
            //printf("%d, %lx vs gold %llx, \n", ii, top->enc_dout, cw_bit[adr-1]);
        }
        tick();
    }
    top->enc_vdin = 0;
}

void top_dec(char *lratio, int &rtl_cycle, int &num_rtl)
{
	static int init = 0;
    static uint64_t *hb_in;
    static uint32_t *sb_in;
    static uint64_t *info_b;
    static uint32_t *sign4d_b;	
	
	if(init==0)
	{
        hb_in       = (uint64_t *)malloc( sizeof(uint64_t)*(int)(var_N/64) );
        sb_in       = (uint32_t *)malloc( sizeof(uint32_t)*(int)(var_N/16) );
        info_b      = (uint64_t *)malloc( sizeof(uint64_t)*(int)(var_N/64) );
        sign4d_b    = (uint32_t *)malloc( sizeof(uint32_t)*(int)(var_N/8 ) );		
		init = 1;
	}
    for(int ii=0; ii<(int)(var_N/64); ii++)
    {
        int idx = (ii*64+check_M)%var_N ; 
        uint64_t hb_tmp   = 0;
        uint32_t sb_tmp[2]= {0 ,0};
        for(int jj=0; jj<64; jj++)
        {
            if( lratio[idx+jj]<0 )       hb_tmp |= (0x1ULL<<jj); 
            if((lratio[idx+jj]>1)||(lratio[idx+jj]<-1))
            {
                if(jj<32)   sb_tmp[0] |= (0x1ULL<<jj);
                else        sb_tmp[1] |= (0x1ULL<<(jj-32)); 
            }
        }
        hb_in[ii]  = hb_tmp ;
        for(int jj=0; jj<4; jj++)     sb_in[4*ii+jj] = (jj<2) ? sb_tmp[jj] : 0 ;
    }

    // ----------
    //
    // ----------
    uint64_t  hb_qin    = 0;
    uint32_t  sb_qin[4] = {0,0,0,0}; 
    uint64_t  info_qin  = 0;
    uint32_t  sign4d_qin[8] = {0,0,0,0,0,0,0,0};
    top->enc0_h1s2 = 2;   // hard-input or soft-input
    tick();  
    top->clear = 1;
    tick();  
    top->clear = 0;
    rtl_cycle  = 0;
    while(!top->out_done)
    {
        //  input
        top->hb_in = hb_qin;
        hb_qin     = hb_in[top->hb_addr];
        for(int ii = 0; ii<4; ii++)   top->sb_in[ii] = sb_qin[ii];
        for(int ii = 0; ii<4; ii++)   sb_qin[ii]     = sb_in[top->sb_addr*4+ii];
        //  output
        top->info_q = info_qin;
        if(top->info_cs)
        {
           info_qin = info_b[top->info_ar];
           if(top->info_we)
           {
               info_b[top->info_aw] = top->info_d;
           }
        } 
        //  sign4d
        for(int ii=0; ii<8; ii++)
        {
            top->sign4d_q[ii] = sign4d_qin[ii];
            if(top->sign4d_cs)
            {
               sign4d_qin[ii] = sign4d_b[top->sign4d_ar*8+ii];
               if(top->sign4d_we)
               {
                   sign4d_b[top->sign4d_aw*8+ii] = top->sign4d_d[ii];
               }
            } 
        }
        tick();
        rtl_cycle++;
    }
	tick();
    num_rtl = top->out_ecc ;
    tick();
}
#endif


typedef struct input_struct{
int      idx;
int     *encoded_bit;
char    *lratio;
char    *dec_bit;
int     *cnt_cor;
int     *cnt_uncor;
int     *cnt_alias;
int      block_i;
int      block_err;
int      errbit_cnt;
} pthrd_input;

void *job_thread(void *args)
{
    pthrd_input *input = static_cast<pthrd_input *>(args);
    int          idx         = input->idx;
    int         *encoded_bit = input->encoded_bit;
    char        *lratio      = input->lratio;
    char        *dec_bit     = input->dec_bit;
    int         *cnt_cor     = input->cnt_cor;
    int         *cnt_uncor   = input->cnt_uncor;
    int         *cnt_alias   = input->cnt_alias;    
    
    __m256i     *m256z       ; 
    int          cycle , rtl_cycle, num_rtl;
    int          ii, jj, this_err, ans_wrong=0, num_awgn, dec_ok, num_dec; 
    
    for( ii=0; ii<max_hist; ii++ )
    {
        cnt_cor[ii]   = 0;
        cnt_uncor[ii] = 0;
        cnt_alias[ii] = 0;
    }
    input->block_err  = 0;
    input->errbit_cnt = 0;

    for( ii=0; (ii<(test_per_msg*(1+1)))&&(input->block_err<min_error); ii++)
    {
         bpsk_awgn( m256z, idum, var, encoded_bit, lratio, num_awgn);  

         dec_ok = ms2p_avx( lratio, dec_bit, cycle, num_dec);
         
         #ifndef NO_VERILATOR
         if(ii>888)                      printf("sim end of  H%d_%d_QC64_DEG4\n", var_N, var_N-check_M);
         if(ii>888)                      exit(0);      // finite testing
    
         top_dec( lratio, rtl_cycle, num_rtl);
         
         // use ip run clock = 200M Hz
         printf("%3d th, ecc=%3d, throughput= %4d MBs (%4d cycles), ", ii, num_dec, (int)(var_N*200/8/cycle), cycle);
         if(cycle!=rtl_cycle)            printf("\n** Mismatch cycle  : c %d vs rtl %d \n", cycle, rtl_cycle);
         else if(num_dec!=num_rtl  )     printf("\n** Mismatch num_ecc: c %d vs rtl %d \n", num_dec, num_rtl);
		 else {
	        if(num_dec==65535)           printf("init_err= %3d", num_awgn);
			else                         printf("eq init");
		 }
		 printf("\n");
         #endif
         //  ---------------------
         //   statistics summary
         //  ---------------------      
         if((dec_ok==1)&&(num_dec!=num_awgn))
         {
            ans_wrong = 1 ;
            dec_ok    = 0 ;
         } 
         else ans_wrong = 0 ;
         
         if(ans_wrong==1)
         {
             for( jj=0, this_err=0; jj<var_N; jj++ )
             {
                 if( dec_bit[jj]!=encoded_bit[jj] ) 
                 {
                     this_err++;
                 }
             }    
             //if(this_err<(int)(var_N*0.02))
             {
                  for( jj=0; jj<var_N; jj++ )
                  {
                      if( dec_bit[jj]!=encoded_bit[jj] ) 
                      {
                          cout << jj << ",";
                      }
                  }    
                  cout <<"dist = "<< this_err <<", alias"<< endl;                   
             }
             if( this_err >= max_hist )     cnt_alias[max_hist-1]++;
             else if( this_err > 0    )     cnt_alias[this_err-1]++;            
         }
         
         if( dec_ok==0 )
         {
             if(longi>test_per_msg)     cout<<num_awgn<<", "<<input->block_err<<'/'<<ii<<endl;
             input->errbit_cnt += num_awgn;
             input->block_err ++;
             cnt_uncor[num_awgn-1]++;
         }
         else  cnt_cor[num_awgn-1]++;  
         
         if(ii==(test_per_msg*(idx+1)-1))   cout<<input->errbit_cnt<<':'<<input->block_err<<':'<<ii<<endl;      
    }
    input->block_i = ii ;
        
    return NULL;
}

int main(int argc, char** argv)
{
    #ifndef NO_VERILATOR
    Verilated::commandArgs(argc, argv);
    #ifdef DUMP_VCD
    Verilated::traceEverOn(true);
    #endif  
    #endif
    
    int *encoded_bit;
    int  errbit_cnt , block_err , error_limit ;
    int  mode_select, j , k ;
    double  all, sn , noise_var2;  
    char  report[100];
    ofstream  outFile , hstFile ;   
    pthrd_input  pd_in[1];
    pthread_t  threads[1];
    int   rc;
    void *result = NULL;
    //  --------------
    //    init
    //  --------------
    if(argc<2)     mode_select = 3 ;
    else
    {
        char* endp = nullptr;
        long v = std::strtol(argv[1], &endp, 10); 
        if (*argv[1] != '\0' && *endp == '\0') {
            mode_select = static_cast<int>(v);
        } else {
            mode_select = 7;
            std::cerr << "Invalid number: '" << argv[1] << "'. Using default " << mode_select << "\n";
        }
    }
    max_iter    = 48;

    // ------------------
    //     
    // ------------------
    encoder ( "LDPC_QLC_TLC_SHWu.txt", 0, mode_select);    
    idum   = (long *)malloc(   1*sizeof(long) ); 
    *idum  = -83168 ;  // rand_seed 
    sprintf( report, "report%d_%d_QC%d_DEG4.txt", var_N, var_N - check_M, 64); 
    outFile.open(report);
    if(!outFile)
    {
        cout << "Failed to open "<< report << endl;
        exit(0);   
    } 
    cw_bit      = (unsigned long long *) malloc(sizeof(unsigned long long) * (int)(var_N/64) );     
    encoded_bit = (int    *)malloc( var_N*sizeof(int)  );   
    cnt_cor     = (int    *)malloc( max_hist *sizeof(int ) );
    cnt_uncor   = (int    *)malloc( max_hist *sizeof(int ) );  
    cnt_alias   = (int    *)malloc( max_hist *sizeof(int ) );

    for( k=0; k<1; k++)
    {
        pd_in[k].idx         = k;
        pd_in[k].encoded_bit = encoded_bit ;
        pd_in[k].lratio      = (char     *)malloc( var_N*sizeof(char) );
        pd_in[k].dec_bit     = (char     *)malloc( var_N*sizeof(char) );
        pd_in[k].cnt_cor     = (int      *)malloc( max_hist *sizeof(int ) );
        pd_in[k].cnt_uncor   = (int      *)malloc( max_hist *sizeof(int ) );  
        pd_in[k].cnt_alias   = (int      *)malloc( max_hist *sizeof(int ) );
    }
    // ------------------
    //     
    // ------------------
    info_source( cw_bit, var_N - check_M);
    // ------------------
    //     
    // ------------------   
    // Notice : parity bits are put in the front of cw_bit.
    //       so data  bits start from cw_bit[ check_M/64 ] to cw_bit[ var_N/64 -1 ]
    //
    // testing encoding throughput
    //for(j=0; j<1000000; j++)//enc_testing
    encoder( "LDPC_QLC_TLC_SHWu.txt", cw_bit, mode_select); 
    //printf("testing encoding 1000000 cws\n");//enc_testing
    //exit(0);//enc_testing
    
    for (j = 0; j < var_N; j = j + 64) 
    {
        for (k = 0; k < 64; k++) 
        {
            encoded_bit[j+k] = (int)((cw_bit[(int)(j/64)]>>k)&1);
        }
    }   

    #ifndef NO_VERILATOR
    top = new Vldpc_wsh;
    #ifdef DUMP_VCD
    tfp = new VerilatedVcdC;
    top->trace(tfp, 99);
    tfp->open("ldpc_wsh.vcd");
    #endif
    top_init(mode_select);
    top_enc (cw_bit, mode_select);
    #endif

    for( sn=sn_start; sn<=sn_end; sn=sn+sn_step )
    {
        noise_var2  = 0.5*pow(10,(sn/-10))*(static_cast<double>(var_N)/static_cast<double>(var_N-check_M));
        var         = sqrt( noise_var2 );
        binary_flip = (int)(65536.0*0.5*erfc(1.0/var/sqrt(2)));    
        seflp       = (int)(binary_flip*( 1 - RSER)); 
        scflp       = (int)(binary_flip+(65535-binary_flip)*RSCR);  
        errbit_cnt  = 0;
        block_err   = 0;    
        longi       = 0;
        error_limit = min_error;             
        
        for( j=0; j<max_hist; j++ )
        {
            cnt_cor[j]    = 0;
            cnt_uncor[j]  = 0;
            cnt_alias[j]  = 0;
        }
        
        cout <<"max_iter= " <<max_iter << ", SNR= "<< sn <<" db"<<endl;
        //  ---------------
        //    LDPC main
        //  ---------------        
        while( (longi<min_loop)||(block_err<error_limit) )
        {
            rc = pthread_create(&threads[0], NULL, job_thread, (void *)&(pd_in[0]));
            if (rc) 
            {
                printf("ERORR; return code from pthread_create() is %d\n", rc);
                exit(0);
            }

            pthread_join(threads[0], &result);

            for( j=0; j<max_hist; j++ )
            {
                cnt_cor[j]   += pd_in[0].cnt_cor[j];
                cnt_uncor[j] += pd_in[0].cnt_uncor[j];
                cnt_alias[j] += pd_in[0].cnt_alias[j];
            }
            errbit_cnt += pd_in[0].errbit_cnt;
            block_err  += pd_in[0].block_err;
            longi      += pd_in[0].block_i;

            if(  longi >  snr_runs     )      break;                          
            if(  longi >= test_per_msg )      cout << errbit_cnt << ':' << block_err << ':' << longi << endl;
            if( (longi%(test_per_msg*(1+1)*5))==0 )
            {
                error_limit = 40 ; 
                
                hstFile.open ("./hstg.txt");
                write_report(hstFile, sn, errbit_cnt, block_err, longi, 0);
                hstFile.close(); 
            }
        } 
        write_report(outFile, sn, errbit_cnt, block_err, longi, 1);
        printf( "H_file  = H%d_%d_QC%d_DEG4\n", var_N, var_N - check_M, 64); 
    }
    outFile.close();
    delete idum;
    delete [] encoded_bit;
    delete [] cnt_cor;
    delete [] cnt_uncor;

    return 0;
}


