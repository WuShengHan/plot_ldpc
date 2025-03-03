//
// g++ -mavx2 -O3 -march=native -static -fvisibility=hidden -c -o wsh.o wsh.cpp; ar rcs libwsh.a wsh.o; strip --strip-unneeded libwsh.a; cd ..; g++ -o run_wsh -mavx2 -O3 -march=native main_wsh.cpp ./dec/libwsh.a -lpthread
//
#include <immintrin.h>
#include <malloc.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <ctime>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>

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
extern int  qcLen   ;
extern int  max_iter ;
extern char  **llrmap; // 3b0

// ----------------
//
// ----------------
int   min_loop     = 1000 ;
int   min_error    = 300  ;
int   test_per_msg = 10000 ;

float sn_start     = 2.10;
float sn_end       = 14.95;
float sn_step      = 0.10;
int   binary_flip  ;
int   seflp        ;
int   scflp        ;
long  *idum        ;
long  longi        ;
double  var        ;

long  snr_runs     = 6000000000;
int    max_hist    = 1000; 

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
	const char ch_idx = 4;  
	//
	//  change ch_idx to 7 when hard decoding, but need change decoder parameter for better quality
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
        if( ((a1+n1)>0.5)||((a1+n1)<-0.5) )   lratio[2*i  ] = sign*7;        // 7 is channel-bit-index, not llr for calculation
        else                                  lratio[2*i  ] = sign*ch_idx;   // 4 is channel-bit-index, not llr for calculation
        
        sign = ((a2+n2)>=0) ? 1 : -1;
        if( ((a2+n2)>0.5)||((a2+n2)<-0.5) )   lratio[2*i+1] = sign*7;        // 7 is channel-bit-index, not llr for calculation
        else                                  lratio[2*i+1] = sign*ch_idx;   // 4 is channel-bit-index, not llr for calculation

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
	int          cycle    ;
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

#if 1

int main(void)
{
    int *encoded_bit   ;
    int  errbit_cnt    ,  block_err    ,  error_limit ;
    int  j, k, num_gltn,  rnd256[8]    ;
    double  all, sn    ,  noise_var2   ;  
    char    str[94]    , tmp[100]      , report[100] ;
	unsigned long long * cw_bit   ; 
    ofstream  outFile , hstFile   ;   
	pthrd_input  pd_in[1];
	pthread_t  threads[1];
	int   rc;
	void *result = NULL;
    //  --------------
    //    init
    //  --------------
	int mode_select = 6 ;
	
	
	encoder ( "LDPC_QLC_TLC_SHWu.txt", 0, mode_select);    
    idum   = (long *)malloc(   1*sizeof(long) ); 
    *idum  = -83168 ;  // rand_seed	
    sprintf( report, "report%d_%d_QC%d_DEG4.txt", var_N, var_N - check_M, qcLen); 
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
		pd_in[k].lratio      = (char   *)malloc( var_N*sizeof(char) );
		pd_in[k].dec_bit     = (char   *)malloc( var_N*sizeof(char) );
        pd_in[k].cnt_cor     = (int    *)malloc( max_hist *sizeof(int ) );
        pd_in[k].cnt_uncor   = (int    *)malloc( max_hist *sizeof(int ) );  
        pd_in[k].cnt_alias   = (int    *)malloc( max_hist *sizeof(int ) );
	}

    // ------------------
	//    LLR setup for decoder
	// ------------------
	//max_iter     = 48;

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
	//for(j=0; j<10000000; j++)      // testing encoding throughput
	encoder( "LDPC_QLC_TLC_SHWu.txt", cw_bit, mode_select); 
    //printf("testing encoding 10000000 cws\n");
	//exit(0);
	
	for (j = 0; j < var_N; j = j + 64) 
	{
		for (k = 0; k < 64; k++) 
		{
			encoded_bit[j+k] = (int)((cw_bit[(int)(j/64)]>>k)&1);
		}
	}	

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
		num_gltn    = 0;
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

            if( (longi >  snr_runs)    )      break;                          
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
		printf( "H_file  = H%d_%d_QC%d_DEG4\n", var_N, var_N - check_M, qcLen); 
    }
    outFile.close();
    delete idum;
    delete [] encoded_bit;
    delete [] cnt_cor;
    delete [] cnt_uncor;

    return 0;
}


#else
	

int main(void)
{
    int *encoded_bit   ;
    int  errbit_cnt    ,  block_err    ,  error_limit ;
    int  j, k, num_gltn,  rnd256[8]    , ii;
    double  all, sn    ,  noise_var2   ;  
    char    str[94]    , tmp[100]      , report[100] ;
    unsigned long long * cw_bit ;
	char gname[128] = "LDPC_QLC_TLC_SHWu.txt";
	
    char        *lratio   ;   
    char        *dec_bit  ;
    int         *llr_pfct ;	
    __m256i     *m256z    ; 
	int          cycle    ;
    int   num_awgn, dec_ok, num_dec; 


    idum   = (long *)malloc(   1*sizeof(long) ); 
    *idum  = -83168 ;  // rand_seed	

    sn = 5.5 ;  
    noise_var2  = 0.5*pow(10,(sn/-10));
    var         = sqrt( noise_var2 );

    cw_bit      = (unsigned long long *) malloc(sizeof(unsigned long long) * (int)(9126/64) );	
    encoded_bit = (int    *)malloc( 9126*sizeof(int)  );	
	lratio      = (char   *)malloc( 9126*sizeof(char) );
    dec_bit     = (char   *)malloc( 9126*sizeof(char) );

    // init
	int mode_select = 0 ;


for(ii=0; ii<20; ii++)
{
    // -----------
	//
	// -----------
    mode_select = ii%7 ;
	encoder( gname, cw_bit, mode_select);	     // init if mode_select change	
	
	info_source(cw_bit, var_N - check_M);
    encoder( gname, cw_bit, mode_select);	
	for (j = 0; j < var_N; j = j + 64) 
	{
		for (k = 0; k < 64; k++) 
		{
			encoded_bit[j+k] = (cw_bit[(int)(j/64)]&(1ULL<<k)) ? 1 : 0 ;
		}
	}
    bpsk_awgn( m256z, idum, var, encoded_bit, lratio, num_awgn); 

    dec_ok = ms2p_avx( lratio, dec_bit, cycle, num_dec);	

    if (dec_ok ==1) 
	{
        for (j = 0, k = 0; j < var_N; j++)   
			if(dec_bit[k]!=encoded_bit[k])    k++ ; 
			
		if(k!=0)    			printf("miscorrection! %d\n", k);		
		if(num_dec!=num_awgn)   printf("wrong number %d, %d!\n", num_dec, num_awgn);
    }
    printf("%d cw, num_awgn = %d, num_dec = %d, using %d cycle\n", ii+1, num_awgn, num_dec, cycle);    
}

}


#endif


