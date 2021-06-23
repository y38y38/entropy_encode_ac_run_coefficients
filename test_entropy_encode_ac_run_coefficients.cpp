
#include <iostream>
#include <verilated.h>
#include "Ventropy_encode_ac_run_coefficients.h"
#include <verilated_fst_c.h> 



int16_t vals[2016];

#define MAX_COEFFICIENT_NUM_PER_BLOCK (64)

int g_print_flag = 0;
void setBit(uint32_t buf, uint32_t size_of_bit)
{
	if (g_print_flag == 1)
		printf("s %x %d\n", buf, size_of_bit);
}


static void golomb_rice_code(int32_t k, uint32_t val)
{
    int32_t q  = val >> k;

    if (k ==0) {
        if (q != 0) {
            setBit(1,q+1);
 //		printf("s %x %d\n", 1,q+1);
       } else {
//		printf("s %x %d\n", 1,1);
        setBit( 1,1);
		}
    } else {
        uint32_t tmp = (k==0) ? 1 : (2<<(k-1));
        uint32_t r = val & (tmp -1 );

        uint32_t codeword = (1 << k) | r;
        setBit(codeword, q + 1 + k );
//		printf("s %x %d\n", codeword, 1+1+k);
    }
    return;
}
int set_bit = 0;
int set_bit_size = 0;
//1 2
static void exp_golomb_code(int32_t k, uint32_t val) {

	//LOG
    int32_t q = floor(log2(val + ((k==0) ? 1 : (2<<(k-1))))) - k;
//printf("q=%d k=%d %d %f %f\n",q,k, val + ((k==0) ? 1 : (2<<(k-1))),
// log2(2),floor(1));
    uint32_t sum = val + ((k==0) ? 1 : (2<<(k-1)));

    int32_t codeword_length = (2 * q) + k + 1;

  setBit(sum, codeword_length+set_bit_size);
  
// 		printf("s %x %d\n", sum, codeword_length+set_bit_size);
     return;
}
static void rice_exp_combo_code(int32_t last_rice_q, int32_t k_rice, int32_t k_exp, uint32_t val)
{
    uint32_t value = (last_rice_q + 1) << k_rice;

    if (val < value) {
        golomb_rice_code(k_rice, val);
    } else {
//		set_bit = 1;
		set_bit_size = last_rice_q + 1;
//        setBit(bitstream, 0,last_rice_q + 1);
        exp_golomb_code(k_exp, val - value);
//		set_bit = 0;
		set_bit_size = 0;
    }
    return;
}
static void encode_vlc_codeword_ac_run(int32_t previousRun, int32_t val)
{
    if ((previousRun== 0)||(previousRun== 1)) {
        rice_exp_combo_code(2,0,1, val);
    } else if ((previousRun== 2)||(previousRun== 3)) {
        rice_exp_combo_code(1,0,1, val);
    } else if (previousRun== 4) {
        exp_golomb_code(0, val);
    } else if ((previousRun>= 5) && (previousRun <= 8))  {
        rice_exp_combo_code(1,1,2, val);
    } else if ((previousRun>= 9) && (previousRun <= 14))  {
        exp_golomb_code(1, val);
    } else {
        exp_golomb_code(2, val);
    }
    return;

}
static void encode_vlc_codeword_ac_level(int32_t previousLevel, int32_t val)
{
    if (previousLevel== 0) {
        rice_exp_combo_code(2,0,2, val);
    } else if (previousLevel== 1) {
        rice_exp_combo_code(1,0,1, val);
    } else if (previousLevel== 2) {
        rice_exp_combo_code(2,0,1, val);
    } else if (previousLevel == 3)  {
        exp_golomb_code(0, val);
    } else if ((previousLevel>= 4) && (previousLevel<= 7))  {
        exp_golomb_code(1, val);
    } else {
        exp_golomb_code(2, val);
    }
    return;

}
static void entropy_encode_dc_coefficient(bool first, int32_t abs_previousDCDiff , int val)
{
    if (first) {
        exp_golomb_code(5, val);
    } else if (abs_previousDCDiff == 0) {
        exp_golomb_code(0, val);
    } else if (abs_previousDCDiff == 1) {
        exp_golomb_code(1, val);
    } else if (abs_previousDCDiff == 2) {
        rice_exp_combo_code(1,2,3, val);
    } else {
        exp_golomb_code(3, val);
    }
    return;

}

static int32_t GetAbs(int32_t val)
{
    if (val < 0) {
        //printf("m\n");
        return val * -1;
    } else {
        //printf("p\n");
        return val;
    }
}



static int32_t Signedintegertosymbolmapping(int32_t val)
{
    uint32_t sn;
    if (val >=0 ) {
        sn = GetAbs(val) << 1;

    } else {
        sn = (GetAbs(val) << 1) - 1;
    }
    return sn;
}
void entropy_encode_dc_coefficients(int16_t*coefficients, int32_t numBlocks)
{
    int32_t DcCoeff;
    int32_t val;
    int32_t previousDCCoeff;
    int32_t previousDCDiff;
    int32_t n;
    int32_t dc_coeff_difference;
    int32_t abs_previousDCDiff;

    DcCoeff = (coefficients[0]) ;
    val = Signedintegertosymbolmapping(DcCoeff);
    entropy_encode_dc_coefficient(true, 0, val);

    
    previousDCCoeff= DcCoeff;
    previousDCDiff = 3;
    n = 1;
    while( n <numBlocks) {
        DcCoeff = (coefficients[n++ * 64]); 
        dc_coeff_difference = DcCoeff - previousDCCoeff;
        if (previousDCDiff < 0) {
            dc_coeff_difference *= -1;
        }
        val = Signedintegertosymbolmapping(dc_coeff_difference);
        abs_previousDCDiff = GetAbs(previousDCDiff );
        entropy_encode_dc_coefficient(false, abs_previousDCDiff, val);
        previousDCDiff = DcCoeff - previousDCCoeff;
        previousDCCoeff= DcCoeff;

    }
    return ;
}

//from figure 4
static const uint8_t block_pattern_scan_table[64] = {
     0,  1,  4,  5, 16, 17, 21, 22,
     2,  3,  6,  7, 18, 20, 23, 28,
     8,  9, 12, 13, 19, 24, 27, 29,
    10, 11, 14, 15, 25, 26, 30, 31,
    32, 33, 37, 38, 45, 46, 53, 54,
    34, 36, 39, 44, 47, 52, 55, 60,
    35, 40, 43, 48, 51, 56, 59, 61,
    41, 42, 49, 50, 57, 58, 62, 63,
};
static uint8_t block_pattern_scan_read_order_table[64];


uint32_t entropy_encode_ac_coefficients(int16_t*coefficients, int32_t numBlocks)
{
    int32_t block;
    int32_t conefficient;
    int32_t run;
    int32_t level;
    int32_t abs_level_minus_1;
    int32_t previousRun = 4;
    int32_t previousLevelSymbol = 1;
    int32_t position;

    run = 0;

int i=0;
    //start is 1 because 0 equal dc position
    for (conefficient = 1; conefficient< MAX_COEFFICIENT_NUM_PER_BLOCK; conefficient++) {
        position = block_pattern_scan_read_order_table[conefficient];
        for (block=0; block < numBlocks; block++) {
//            level = coefficients[(block * MAX_COEFFICIENT_NUM_PER_BLOCK) + position] ;
            level = coefficients[i] ;
//			printf("%d\n", level);
			i++;

            if (level != 0) {
				g_print_flag = 1;
				printf("previousRun %d %d\n", previousRun, i);
				printf("run %d\n", run);
                encode_vlc_codeword_ac_run(previousRun, run);
				g_print_flag = 0;

                abs_level_minus_1 = GetAbs(level) - 1;
                encode_vlc_codeword_ac_level( previousLevelSymbol, abs_level_minus_1);
                if (level >=0) {
                    //setBit(bitstream, 0,1);
                } else {
                    //setBit(bitstream, 1,1);
                }

                previousRun = run;
                previousLevelSymbol = abs_level_minus_1;
                run    = 0;

            } else {
                run++;
            }
        }
    }
    return 0;
}


int time_counter = 0;

int main(int argc, char** argv) {

//	Verilated::commandArgs(argc, argv);
	FILE *in = fopen(argv[1], "r");
	if (in==NULL) {
		printf("err\n");
	}
	int i;
	for (i=0;i<2016;i++) {
		fscanf(in, "%d", &vals[i]);
	}
	fclose(in);

//	entropy_encode_ac_coefficients(vals, 32);
//	return 0;
//entropy_encode_dc_coefficients(vals);
	// Instantiate DUT
	Ventropy_encode_ac_run_coefficients *dut = new Ventropy_encode_ac_run_coefficients();
	// Trace DUMP ON
	Verilated::traceEverOn(true);
	VerilatedFstC* tfp = new VerilatedFstC;

	dut->trace(tfp, 100);  // Trace 100 levels of hierarchy
	tfp->open("simx.fst");

	// Format
	dut->reset_n = 0;
	dut->clk = 0;

	// Reset Time
	while (time_counter < 10) {
		dut->clk = !dut->clk; // Toggle clock
		dut->eval();
		tfp->dump(time_counter);  // 波形ダンプ用の記述を追加
		time_counter++;
	}
	// Release reset
	dut->reset_n = 1;

	int state = 0;
//	int val[12] = {5,-6,-8,15,15,-9,-3,-5,0,0,0,0};
//	int i =0;
//	int k,val,codeword_length,sum;
//			printf("previousDCDiff %d\n", dut->previousDCDiff);
//			printf("abs_previousDCDiff %d\n", dut->abs_previousDCDiff);
//printf("\n\n");
/*
			printf("%x %d\n", dut->sum, dut->LENGTH);
			printf("abs_previousDCDiff %d\n", dut->abs_previousDCDiff);
			printf("abs_previousDCDiff_next %d\n", dut->abs_previousDCDiff_next);
			printf("previousDCCoeff %d\n", dut->previousDCCoeff);
			printf("previousDCDiff %d\n", dut->previousDCDiff);
			printf("dc_coeff_difference %d\n", dut->dc_coeff_difference);
			printf("val %d\n", dut->val);
			printf("val_n %d\n", dut->val_n);
			printf("DcCoeff %d\n", dut->DcCoeff);
*/
//	while (time_counter < 68 && !Verilated::gotFinish()) {
	i=0;
//	while (time_counter < 82 && !Verilated::gotFinish()) {
	while (!Verilated::gotFinish()) {
		dut->clk = !dut->clk; // Toggle clock
		if (dut->clk) {
//			fscanf(in, "%d,%d,%d,%d",&k, &val, &sum, &codeword_length);
//			if (state == 0) {
//				dut->k = k;
//				if (i>=12) {
//					printf("end\n");
//					break;
//				}
				dut->Coeff = vals[i];
				i++;
				if (i==2016) {
					break;
				}
//				state = 1;
//			}
#if 0

printf("Coeff= %d %d\n", dut->Coeff, i);
printf("previousRun= %d\n", dut->previousRun);
printf("run =%d\n",dut->run);
printf("run_n =%d\n",dut->run_n);
printf("is_expo_golomb_code =%d\n",dut->is_expo_golomb_code);
printf("is_add_setbit =%d\n",dut->is_add_setbit);
printf("k =%d\n",dut->k);
printf("q =%d\n",dut->q);
#endif
		}

		// Evaluate DUT
		dut->eval();
		if (dut->clk) {
			
//			uint32_t sum;
//			int32_t codword_length;
//			exp_golomb_code(dut->k, dut->input_data, &sum, &codword_length);
//			if ((dut->sum != sum) || (dut->CODEWORD_LENGTH != codeword_length)) {
//				printf("q=%d\n",dut->Q);
//				printf("k=%d,input_data=%d sum=%d len=%d %x, sum=%d len=%d\n", 
//				dut->k, dut->input_data, dut->sum, dut->CODEWORD_LENGTH, dut->output_enable,  sum, codeword_length);

//			}
			if (dut->codeword_length != 0)
			 {
//			printf("s %x %d %d\n", dut->sum, dut->codeword_length, i);
			printf("s %x %d\n", dut->sum, dut->codeword_length);
//			printf("s %x %d %x\n", dut->sum, dut->codeword_length,dut->output_enable );
//printf("q =%d\n",dut->q);
//printf("is_expo_golomb_code =%d\n",dut->is_expo_golomb_code);

			}
		}
		tfp->dump(time_counter);  // 波形ダンプ用の記述を追加
		time_counter++;
//		break;
	}

	dut->final();
	tfp->close(); 
}