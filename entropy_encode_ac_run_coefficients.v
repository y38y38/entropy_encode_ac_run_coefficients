module entropy_encode_ac_run_coefficients(
	input clk,
	input reset_n,

	//本当は19bitで足りるが、本関数の処理上桁溢れする可能性があるので、
	//1bit多く用意しておく。
	input [19:0] Coeff,
	output reg [23:0] output_enable,//mask
	output reg [23:0] sum,
	output reg [31:0] codeword_length,

	output reg [19:0] previousRun,
output reg signed [19:0] run,
output reg signed [19:0] run_n,
output reg [1:0] is_expo_golomb_code,
output reg [1:0] is_add_setbit,
output reg [2:0] k,
output reg [31:0] q



);

function [19:0] getabs;
	input [19:0] value;
	begin
		if (value[19] != 1'b1) begin
			getabs = value;
		end else begin
			getabs =  (~(value - 1));
		end
	end
endfunction

function [31:0] getfloorclog2;
	input [19:0] val;
	begin
		reg [19:0] in_val;
		in_val = val;
		for (getfloorclog2=0; in_val>0; getfloorclog2=getfloorclog2+1) begin
			in_val = in_val>>1;
		end
		getfloorclog2 = getfloorclog2 - 1;
	end
endfunction



function [23:0] bitmask;
	input [5:0] val;
	reg [5:0] index = 6'h0;
	begin
		bitmask = 24'h1;
		for(index=1;index<val;index=index+1) begin
			bitmask = (bitmask<<1) | 1;
		end
	end
endfunction



always @(posedge clk, negedge reset_n) begin
	if (!reset_n) begin
		previousRun <= 20'h4;
		run <= 20'h0;
		is_expo_golomb_code <= 2'h2;

	end else begin
		if (Coeff != 0) begin
			if ((previousRun == 0) || (previousRun == 1)) begin
				if (run < 3) begin
					is_expo_golomb_code <= 2'b0;
					is_add_setbit<=2'h0;
					k <= 0;
//					run_n <= 100;
					run_n <= run;
					q = run;

				end else begin
					is_expo_golomb_code <= 2'b1;
					is_add_setbit<=2'h3;
					k <= 1;
//					run_n <= 100;
					run_n <= run - 3;
//run_n <= -4;
				end
			end else if ((previousRun == 2) || (previousRun == 3)) begin
				if (run < 2) begin
					is_expo_golomb_code <= 2'b0;
					is_add_setbit<=2'h0;
					k <= 0;
					run_n <= run;
					q = run;

				end else begin
					is_expo_golomb_code <= 2'b1;
					is_add_setbit<=2'h2;
					k <= 1;
					run_n <= run - 2;
				end
			end else if ((previousRun == 4)) begin
				is_expo_golomb_code <= 2'b1;
				is_add_setbit <= 2'h0;
				k <= 0;
				run_n <= run;
			end else if ((previousRun >= 5) && (previousRun <= 8)) begin
				if (run < 4) begin
					is_expo_golomb_code <= 2'b0;
					is_add_setbit<=2'h0;
					k <= 1;
					run_n <= run;
					q = run>>1;

				end else begin
					is_expo_golomb_code <= 2'b1;
					is_add_setbit<=2'h2;
					k <= 2;
					run_n <= run - 4;
				end
			end else if ((previousRun >= 9) && (previousRun <= 14)) begin
				is_expo_golomb_code <= 2'b1;
				is_add_setbit <= 2'h0;
				k <= 1;
				run_n <= run;
			end else begin
				is_expo_golomb_code <= 2'b1;
				is_add_setbit <= 2'h0;
				k <= 2;
				run_n <= run;
			end
			previousRun <= run;
			run <= 0;
		end else begin
			run <= run + 1;
			is_expo_golomb_code <= 2'b10;

		end

	end
end





//assign LENGTH = codeword_length;


//reg [31:0] codeword_length = 32'h0;

//exp_golomb_code
always @(posedge clk, negedge reset_n) begin
	if (!reset_n) begin
		output_enable = 24'h0;
		sum = 24'h0;
	end else begin
		if (is_expo_golomb_code == 2'b1) begin
			q = getfloorclog2((run_n + (1<<(k)))) - k;
			//q =  input_data + 16'h1;
			sum[19:0] = run_n + (1<<k);
//			if (is_add_setbit == 1'b1) begin
				//dd
				codeword_length = (2 * q) + k + 1 + is_add_setbit;
//			end else begin
//				codeword_length = (2 * q) + k + 1;
//			end
			output_enable = bitmask(codeword_length);
		end
	end
end


//golomb_rice_code
always @(posedge clk, negedge reset_n) begin
	if (!reset_n) begin
		output_enable = 24'h0;
		sum = 24'h0;
		codeword_length = 32'h0;
	end else begin
		if (is_expo_golomb_code == 2'b0) begin
			//q=と、if(qのタイミング)
//			q = run_n >> k;
			if (k==0) begin
				if(q!=0) begin
					sum = 1;
					codeword_length = q+1;
					output_enable = bitmask(codeword_length);
				end else begin
					sum = 1;
					codeword_length = 1;
					output_enable = 1;
				end
			end else begin
				// 0x4 | 1 & 0x3
				// 0x5 
				sum = (1<<k) | (run_n & ((1<<k) - 1));
//				sum = 20'h111 ;//(1<<k) | (run_n & ((1<<k) - 1));
				codeword_length = q + 1 + k;
				output_enable = bitmask( codeword_length);	
			end
		end
	end
end

always @(posedge clk, negedge reset_n) begin
	if (!reset_n) begin
	end else begin
		if (is_expo_golomb_code == 2'h2) begin
			sum = 0;
			codeword_length = 0;
			output_enable = 0;
		end
	end
end



endmodule;
