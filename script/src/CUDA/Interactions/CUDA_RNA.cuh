/* System constants */
__constant__ int MD_N[1];
__constant__ int MD_n_forces[1];

__constant__ float MD_hb_multi[1];
__constant__ float MD_F1_A[2];
__constant__ float MD_F1_RC[2];
__constant__ float MD_F1_R0[2];
__constant__ float MD_F1_BLOW[2];
__constant__ float MD_F1_BHIGH[2];
__constant__ float MD_F1_RLOW[2];
__constant__ float MD_F1_RHIGH[2];
__constant__ float MD_F1_RCLOW[2];
__constant__ float MD_F1_RCHIGH[2];
// 50 = 2 * 5 * 5
__constant__ float MD_F1_EPS[50];
__constant__ float MD_F1_SHIFT[50];

__constant__ float MD_F2_K[2];
__constant__ float MD_F2_RC[2];
__constant__ float MD_F2_R0[2];
__constant__ float MD_F2_BLOW[2];
__constant__ float MD_F2_RLOW[2];
__constant__ float MD_F2_RCLOW[2];
__constant__ float MD_F2_BHIGH[2];
__constant__ float MD_F2_RCHIGH[2];
__constant__ float MD_F2_RHIGH[2];

__constant__ float MD_F5_PHI_A[4];
__constant__ float MD_F5_PHI_B[4];
__constant__ float MD_F5_PHI_XC[4];
__constant__ float MD_F5_PHI_XS[4];

__constant__ float MD_dh_RC[1];
__constant__ float MD_dh_RHIGH[1];
__constant__ float MD_dh_prefactor[1];
__constant__ float MD_dh_B[1];
__constant__ float MD_dh_minus_kappa[1];
__constant__ bool MD_dh_half_charged_ends[1];

//__constant__ struct Model RNA_MODEL;

struct CUDAModel {
	float RNA_POS_BACK;
	float RNA_POS_STACK;
	float RNA_POS_BASE;
	// RNA_POS_STACK - RNA_POS_BACK
	float RNA_GAMMA;

	//this is stacking site with your 3-neighbor
	float RNA_POS_STACK_3_a1;
	float RNA_POS_STACK_3_a2;
	//this is stacking site for interaction with your 5-neighbor
	float RNA_POS_STACK_5_a1;
	float RNA_POS_STACK_5_a2;
	/*
	 * FENE
	 */
	float RNA_FENE_EPS;
	float RNA_FENE_R0;
	float RNA_FENE_DELTA;
	float RNA_FENE_DELTA2;
	/*
	 * EXCLUDED VOLUME
	 */
	float RNA_EXCL_EPS;
	float RNA_EXCL_S1;
	float RNA_EXCL_S2;
	float RNA_EXCL_S3;
	float RNA_EXCL_S4;
	float RNA_EXCL_R1;
	float RNA_EXCL_R2;
	float RNA_EXCL_R3;
	float RNA_EXCL_R4;
	float RNA_EXCL_B1;
	float RNA_EXCL_B2;
	float RNA_EXCL_B3;
	float RNA_EXCL_B4;
	float RNA_EXCL_RC1;
	float RNA_EXCL_RC2;
	float RNA_EXCL_RC3;
	float RNA_EXCL_RC4;

	/*
	 * HYDROGEN BONDING
	 */
	// radial part
	float RNA_HYDR_EPS;
	float RNA_HYDR_A;
	float RNA_HYDR_RC;
	float RNA_HYDR_R0;
	float RNA_HYDR_BLOW;
	float RNA_HYDR_BHIGH;
	float RNA_HYDR_RLOW;
	float RNA_HYDR_RHIGH;
	float RNA_HYDR_RCLOW;
	float RNA_HYDR_RCHIGH;

	// angular part

	float RNA_HYDR_THETA1_A;
	float RNA_HYDR_THETA1_B;
	float RNA_HYDR_THETA1_T0;
	float RNA_HYDR_THETA1_TS;
	float RNA_HYDR_THETA1_TC;

	float RNA_HYDR_THETA2_A;
	float RNA_HYDR_THETA2_B;
	float RNA_HYDR_THETA2_T0;
	float RNA_HYDR_THETA2_TS;
	float RNA_HYDR_THETA2_TC;

	float RNA_HYDR_THETA3_A;
	float RNA_HYDR_THETA3_B;
	float RNA_HYDR_THETA3_T0;
	float RNA_HYDR_THETA3_TS;
	float RNA_HYDR_THETA3_TC;

	float RNA_HYDR_THETA4_A;
	float RNA_HYDR_THETA4_B;
	float RNA_HYDR_THETA4_T0;
	float RNA_HYDR_THETA4_TS;
	float RNA_HYDR_THETA4_TC;

	float RNA_HYDR_THETA7_A;
	float RNA_HYDR_THETA7_B;
	float RNA_HYDR_THETA7_T0;
	float RNA_HYDR_THETA7_TS;
	float RNA_HYDR_THETA7_TC;

	float RNA_HYDR_THETA8_A;
	float RNA_HYDR_THETA8_B;
	float RNA_HYDR_THETA8_T0;
	float RNA_HYDR_THETA8_TS;
	float RNA_HYDR_THETA8_TC;

	/*
	 * STACKING
	 */
	// radial part
	float RNA_STCK_BASE_EPS;
	float RNA_STCK_FACT_EPS;
	float RNA_STCK_A;
	float RNA_STCK_RC;
	float RNA_STCK_R0;
	float RNA_STCK_BLOW;
	float RNA_STCK_BHIGH;
	float RNA_STCK_RLOW;
	float RNA_STCK_RHIGH;
	float RNA_STCK_RCLOW;
	float RNA_STCK_RCHIGH;

	// angular part, theta
	float RNA_STCK_THETA4_A;
	float RNA_STCK_THETA4_B;
	float RNA_STCK_THETA4_T0;
	float RNA_STCK_THETA4_TS;
	float RNA_STCK_THETA4_TC;
	float RNA_STCK_THETA5_A;
	float RNA_STCK_THETA5_B;
	float RNA_STCK_THETA5_T0;
	float RNA_STCK_THETA5_TS;
	float RNA_STCK_THETA5_TC;
	float RNA_STCK_THETA6_A;
	float RNA_STCK_THETA6_B;
	float RNA_STCK_THETA6_T0;
	float RNA_STCK_THETA6_TS;
	float RNA_STCK_THETA6_TC;

	float STCK_THETAB1_A;
	float STCK_THETAB1_B;
	float STCK_THETAB1_T0;
	float STCK_THETAB1_TS;
	float STCK_THETAB1_TC;

	float STCK_THETAB2_A;
	float STCK_THETAB2_B;
	float STCK_THETAB2_T0;
	float STCK_THETAB2_TS;
	float STCK_THETAB2_TC;

	// angular part, phi
	float RNA_STCK_PHI1_A;
	float RNA_STCK_PHI1_B;
	float RNA_STCK_PHI1_XC;
	float RNA_STCK_PHI1_XS;
	float RNA_STCK_PHI2_A;
	float RNA_STCK_PHI2_B;
	float RNA_STCK_PHI2_XC;
	float RNA_STCK_PHI2_XS;

	/*
	 *  CROSS STACKING
	 */
	// radial part
	float RNA_CRST_R0;
	float RNA_CRST_RC;
	float RNA_CRST_K;
	float RNA_CRST_BLOW;
	float RNA_CRST_RLOW;
	float RNA_CRST_RCLOW;
	float RNA_CRST_BHIGH;
	float RNA_CRST_RHIGH;
	float RNA_CRST_RCHIGH;
	// angular part;
	float RNA_CRST_THETA1_A;
	float RNA_CRST_THETA1_B;
	float RNA_CRST_THETA1_T0;
	float RNA_CRST_THETA1_TS;
	float RNA_CRST_THETA1_TC;
	float RNA_CRST_THETA2_A;
	float RNA_CRST_THETA2_B;
	float RNA_CRST_THETA2_T0;
	float RNA_CRST_THETA2_TS;
	float RNA_CRST_THETA2_TC;
	float RNA_CRST_THETA3_A;
	float RNA_CRST_THETA3_B;
	float RNA_CRST_THETA3_T0;
	float RNA_CRST_THETA3_TS;
	float RNA_CRST_THETA3_TC;
	float RNA_CRST_THETA4_A;
	float RNA_CRST_THETA4_B;
	float RNA_CRST_THETA4_T0;
	float RNA_CRST_THETA4_TS;
	float RNA_CRST_THETA4_TC;
	float RNA_CRST_THETA7_A;
	float RNA_CRST_THETA7_B;
	float RNA_CRST_THETA7_T0;
	float RNA_CRST_THETA7_TS;
	float RNA_CRST_THETA7_TC;
	float RNA_CRST_THETA8_A;
	float RNA_CRST_THETA8_B;
	float RNA_CRST_THETA8_T0;
	float RNA_CRST_THETA8_TS;
	float RNA_CRST_THETA8_TC;

	/*
	 *  COAXIAL STACKING
	 */
	// radial part
	float RNA_CXST_R0;
	float RNA_CXST_RC;
	float RNA_CXST_K;
	float RNA_CXST_BLOW;
	float RNA_CXST_RLOW;
	float RNA_CXST_RCLOW;
	float RNA_CXST_BHIGH;
	float RNA_CXST_RHIGH;
	float RNA_CXST_RCHIGH;
	// angular part;

	float RNA_CXST_THETA1_A;
	float RNA_CXST_THETA1_B;
	float RNA_CXST_THETA1_T0;
	float RNA_CXST_THETA1_TS;
	float RNA_CXST_THETA1_TC;
	float RNA_CXST_THETA4_A;
	float RNA_CXST_THETA4_B;
	float RNA_CXST_THETA4_T0;
	float RNA_CXST_THETA4_TS;
	float RNA_CXST_THETA4_TC;
	float RNA_CXST_THETA5_A;
	float RNA_CXST_THETA5_B;
	float RNA_CXST_THETA5_T0;
	float RNA_CXST_THETA5_TS;
	float RNA_CXST_THETA5_TC;
	float RNA_CXST_THETA6_A;
	float RNA_CXST_THETA6_B;
	float RNA_CXST_THETA6_T0;
	float RNA_CXST_THETA6_TS;
	float RNA_CXST_THETA6_TC;
	// angular part, phi3 and phi4

	float RNA_CXST_PHI3_A;
	float RNA_CXST_PHI3_B;
	float RNA_CXST_PHI3_XC;
	float RNA_CXST_PHI3_XS;
	float RNA_CXST_PHI4_A;
	float RNA_CXST_PHI4_B;
	float RNA_CXST_PHI4_XC;
	float RNA_CXST_PHI4_XS;

//new stuff
	//float INCLINATION;
	//float BASE_SIZE_R;
	//float ROT_SUGAR;
	float p3_x, p3_y, p3_z, p5_x, p5_y, p5_z;
	float RNA_POS_BACK_a1, RNA_POS_BACK_a2, RNA_POS_BACK_a3;
};

__constant__ CUDAModel rnamodel;

#include "../cuda_utils/CUDA_lr_common.cuh"

__forceinline__ __device__ void _excluded_volume(const c_number4 &r, c_number4 &F, c_number sigma, c_number rstar, c_number b, c_number rc) {
	c_number rsqr = CUDA_DOT(r, r);

	F.x = F.y = F.z = F.w = (c_number) 0.f;
	if(rsqr < SQR(rc)) {
		if(rsqr > SQR(rstar)) {
			c_number rmod = sqrt(rsqr);
			c_number rrc = rmod - rc;
			c_number fmod = 2.f * rnamodel.RNA_EXCL_EPS * b * rrc / rmod;
			F.x = r.x * fmod;
			F.y = r.y * fmod;
			F.z = r.z * fmod;
			F.w = rnamodel.RNA_EXCL_EPS * b * SQR(rrc);
		}
		else {
			c_number lj_part = CUB(SQR(sigma)/rsqr);
			c_number fmod = 24.f * rnamodel.RNA_EXCL_EPS * (lj_part - 2.f * SQR(lj_part)) / rsqr;
			F.x = r.x * fmod;
			F.y = r.y * fmod;
			F.z = r.z * fmod;
			F.w = 4.f * rnamodel.RNA_EXCL_EPS * (SQR(lj_part) - lj_part);
		}
	}
}

__forceinline__ __device__ c_number _f1(c_number r, int type, int n3, int n5) {
	c_number val = (c_number) 0.f;
	if(r < MD_F1_RCHIGH[type]) {
		int eps_index = 25 * type + n3 * 5 + n5;
		if(r > MD_F1_RHIGH[type]) {
			val = MD_F1_EPS[eps_index] * MD_F1_BHIGH[type] * SQR(r - MD_F1_RCHIGH[type]);
		}
		else if(r > MD_F1_RLOW[type]) {
			c_number tmp = 1.f - expf(-(r - MD_F1_R0[type]) * MD_F1_A[type]);
			val = MD_F1_EPS[eps_index] * SQR(tmp) - MD_F1_SHIFT[eps_index];
		}
		else if(r > MD_F1_RCLOW[type]) {
			val = MD_F1_EPS[eps_index] * MD_F1_BLOW[type] * SQR(r - MD_F1_RCLOW[type]);
		}
	}

	return val;
}

__forceinline__ __device__ c_number _f1D(c_number r, int type, int n3, int n5) {
	c_number val = (c_number) 0.f;
	int eps_index = 0;
	if(r < MD_F1_RCHIGH[type]) {
		eps_index = 25 * type + n3 * 5 + n5;
		if(r > MD_F1_RHIGH[type]) {
			val = 2.f * MD_F1_BHIGH[type] * (r - MD_F1_RCHIGH[type]);
		}
		else if(r > MD_F1_RLOW[type]) {
			c_number tmp = expf(-(r - MD_F1_R0[type]) * MD_F1_A[type]);
			val = 2.f * (1.f - tmp) * tmp * MD_F1_A[type];
		}
		else if(r > MD_F1_RCLOW[type]) {
			val = 2.f * MD_F1_BLOW[type] * (r - MD_F1_RCLOW[type]);
		}
	}

	return MD_F1_EPS[eps_index] * val;
}

__forceinline__ __device__ c_number _fX(c_number r, int type, int n3, int n5) {
	c_number val = (c_number) 0.f;
	if(r < MD_F1_RCHIGH[type]) {
		int eps_index = 25 * type + n3 * 5 + n5;
		if(r > MD_F1_RHIGH[type]) {
			val = MD_F1_EPS[eps_index] * MD_F1_BHIGH[type] * SQR(r - MD_F1_RCHIGH[type]);
		}
		else if(r > MD_F1_R0[type]) {
			c_number tmp = 1.f - expf(-(r - MD_F1_R0[type]) * MD_F1_A[type]);
			val = MD_F1_EPS[eps_index] * SQR(tmp) - MD_F1_SHIFT[eps_index];
		}
		else {
			val = -MD_F1_SHIFT[eps_index];
			//val = MD_F1_EPS[eps_index] * MD_F1_BLOW[type] * SQR(r - MD_F1_RCLOW[type]);
		}
	}

	return val;
}

__forceinline__ __device__ c_number _fXD(c_number r, int type, int n3, int n5) {
	c_number val = (c_number) 0.f;
	int eps_index = 0;
	if(r < MD_F1_RCHIGH[type]) {
		eps_index = 25 * type + n3 * 5 + n5;
		if(r > MD_F1_RHIGH[type]) {
			val = 2.f * MD_F1_BHIGH[type] * (r - MD_F1_RCHIGH[type]);
		}
		else if(r > MD_F1_R0[type]) {
			c_number tmp = expf(-(r - MD_F1_R0[type]) * MD_F1_A[type]);
			val = 2.f * (1.f - tmp) * tmp * MD_F1_A[type];
		}

	}

	return MD_F1_EPS[eps_index] * val;
}

__forceinline__ __device__ c_number _f2(c_number r, int type) {
	c_number val = (c_number) 0.f;
	if(r < MD_F2_RCHIGH[type]) {
		if(r > MD_F2_RHIGH[type]) {
			val = MD_F2_K[type] * MD_F2_BHIGH[type] * SQR(r - MD_F2_RCHIGH[type]);
		}
		else if(r > MD_F2_RLOW[type]) {
			val = (MD_F2_K[type] * 0.5f) * (SQR(r - MD_F2_R0[type]) - SQR(MD_F2_RC[type] - MD_F2_R0[type]));
		}
		else if(r > MD_F2_RCLOW[type]) {
			val = MD_F2_K[type] * MD_F2_BLOW[type] * SQR(r - MD_F2_RCLOW[type]);
		}
	}
	return val;
}

__forceinline__ __device__ c_number _f2D(c_number r, int type) {
	c_number val = (c_number) 0.f;
	if(r < MD_F2_RCHIGH[type]) {
		if(r > MD_F2_RHIGH[type]) {
			val = 2.f * MD_F2_K[type] * MD_F2_BHIGH[type] * (r - MD_F2_RCHIGH[type]);
		}
		else if(r > MD_F2_RLOW[type]) {
			val = MD_F2_K[type] * (r - MD_F2_R0[type]);
		}
		else if(r > MD_F2_RCLOW[type]) {
			val = 2.f * MD_F2_K[type] * MD_F2_BLOW[type] * (r - MD_F2_RCLOW[type]);
		}
	}
	/*
	 if(type == 1)
	 printf("DEVICE _f2D: %f,  Running with %f %f %f %f %f %f %f\n",val,MD_F2_RCHIGH[type],MD_F2_RHIGH[type],MD_F2_BHIGH[type],MD_F2_BHIGH[type],MD_F2_RCLOW[type],MD_F2_BLOW[type], MD_F2_K[type]);
	 */
	return val;
}

__forceinline__ __device__ c_number _f4(c_number t, float t0, float ts, float tc, float a, float b) {
	c_number val = (c_number) 0.f;
	t -= t0;
	if(t < 0) t = -t;

	if(t < tc) {
		if(t > ts) {
			// smoothing
			val = b * SQR(tc - t);
		}
		else val = (c_number) 1.f - a * SQR(t);
	}

	return val;
}

__forceinline__ __device__ c_number _f4Dsin(c_number t, float t0, float ts, float tc, float a, float b) {
	c_number val = (c_number) 0.f;
	c_number tt0 = t - t0;
	// this function is a parabola centered in t0. If tt0 < 0 then the value of the function
	// is the same but the value of its derivative has the opposite sign, so m = -1
	c_number m = copysignf((c_number) 1.f, tt0);
	tt0 = copysignf(tt0, (c_number) 1.f);

	if(tt0 < tc) {
		c_number sint = sinf(t);
		if(tt0 > ts) {
			// smoothing
			val = b * (tt0 - tc) / sint;
		}
		else {
			if(SQR(sint) > 1e-12f) val = -a * tt0 / sint;
			else val = -a;
		}
	}

	return 2.f * m * val;
}

__forceinline__ __device__ c_number _f5(c_number f, int type) {
	c_number val = (c_number) 0.f;

	if(f > MD_F5_PHI_XC[type]) {
		if(f < MD_F5_PHI_XS[type]) {
			val = MD_F5_PHI_B[type] * SQR(MD_F5_PHI_XC[type] - f);
		}
		else if(f < 0.f) {
			val = (c_number) 1.f - MD_F5_PHI_A[type] * SQR(f);
		}
		else val = 1.f;
	}

	return val;
}

__forceinline__ __device__ c_number _f5D(c_number f, int type) {
	c_number val = (c_number) 0.f;

	if(f > MD_F5_PHI_XC[type]) {
		if(f < MD_F5_PHI_XS[type]) {
			val = 2.f * MD_F5_PHI_B[type] * (f - MD_F5_PHI_XC[type]);
		}
		else if(f < 0.f) {
			val = (c_number) -2.f * MD_F5_PHI_A[type] * f;
		}
	}

	return val;
}

template<bool qIsN3>
__device__ void _bonded_excluded_volume(c_number4 &r, c_number4 &n3pos_base, c_number4 &n3pos_back, c_number4 &n5pos_base, c_number4 &n5pos_back, c_number4 &F, c_number4 &T) {
	c_number4 Ftmp;
	// BASE-BASE
	c_number4 rcenter = r + n3pos_base - n5pos_base;
	_excluded_volume(rcenter, Ftmp, rnamodel.RNA_EXCL_S2, rnamodel.RNA_EXCL_R2, rnamodel.RNA_EXCL_B2, rnamodel.RNA_EXCL_RC2);
	c_number4 torquep1 = (qIsN3) ? _cross(n5pos_base, Ftmp) : _cross(n3pos_base, Ftmp);
	F += Ftmp;

	// n5-BASE vs. n3-BACK
	rcenter = r + n3pos_back - n5pos_base;
	_excluded_volume(rcenter, Ftmp, rnamodel.RNA_EXCL_S3, rnamodel.RNA_EXCL_R3, rnamodel.RNA_EXCL_B3, rnamodel.RNA_EXCL_RC3);
	c_number4 torquep2 = (qIsN3) ? _cross(n5pos_base, Ftmp) : _cross(n3pos_back, Ftmp);
	F += Ftmp;

	// n5-BACK vs. n3-BASE
	rcenter = r + n3pos_base - n5pos_back;
	_excluded_volume(rcenter, Ftmp, rnamodel.RNA_EXCL_S4, rnamodel.RNA_EXCL_R4, rnamodel.RNA_EXCL_B4, rnamodel.RNA_EXCL_RC4);
	c_number4 torquep3 = (qIsN3) ? _cross(n5pos_back, Ftmp) : _cross(n3pos_base, Ftmp);
	F += Ftmp;

	T += torquep1 + torquep2 + torquep3;
}

template<bool qIsN3>
__device__ void _bonded_part(c_number4 &n5pos, c_number4 &n5x, c_number4 &n5y, c_number4 &n5z, c_number4 &n3pos, c_number4 &n3x, c_number4 &n3y, c_number4 &n3z, c_number4 &F, c_number4 &T, bool average, bool use_mbf, c_number mbf_xmax, c_number mbf_finf) {
	//printf("Hello from bonded part function \n");
	int n3type = get_particle_type(n3pos);
	int n5type = get_particle_type(n5pos);

	c_number4 r = make_c_number4(n3pos.x - n5pos.x, n3pos.y - n5pos.y, n3pos.z - n5pos.z, (c_number) 0);

	c_number4 n5pos_back;
	c_number4 n3pos_back;
	//if(grooving) n5pos_back = n5x * POS_MM_BACK1 + n5y * POS_MM_BACK2;
	//else n5pos_back = n5x * rnamodel.RNA_POS_BACK;

	n3pos_back = n3x * rnamodel.RNA_POS_BACK_a1 + n3y * rnamodel.RNA_POS_BACK_a2 + n3z * rnamodel.RNA_POS_BACK_a3;
	n5pos_back = n5x * rnamodel.RNA_POS_BACK_a1 + n5y * rnamodel.RNA_POS_BACK_a2 + n5z * rnamodel.RNA_POS_BACK_a3;

	c_number4 n5pos_base = n5x * rnamodel.RNA_POS_BASE;

	//if(grooving) n3pos_back = n3x * POS_MM_BACK1 + n3y * POS_MM_BACK2;
	//else n3pos_back = n3x * rnamodel.RNA_POS_BACK;
	c_number4 n3pos_base = n3x * rnamodel.RNA_POS_BASE;
	c_number4 n3pos_stack_5 = n3x * rnamodel.RNA_POS_STACK_5_a1 + n3y * rnamodel.RNA_POS_STACK_5_a2;  //n3x * rnamodel.RNA_POS_STACK;
	//c_number4 n3pos_stack_3 = n3x * rnamodel.RNA_POS_STACK_3_a1   + n3y * rnamodel.RNA_POS_STACK_3_a2;  //n3x * rnamodel.RNA_POS_STACK;

	//c_number4 n5pos_stack_5 = n5x * rnamodel.RNA_POS_STACK_5_a1   + n5y * rnamodel.RNA_POS_STACK_5_a2;  //n3x * rnamodel.RNA_POS_STACK;
	c_number4 n5pos_stack_3 = n5x * rnamodel.RNA_POS_STACK_3_a1 + n5y * rnamodel.RNA_POS_STACK_3_a2;  //n3x * rnamodel.RNA_POS_STACK;

	c_number4 rback = r + n3pos_back - n5pos_back;

	c_number rbackmod = _module(rback);

	c_number4 rbackdir = make_c_number4(rback.x / rbackmod, rback.y / rbackmod, rback.z / rbackmod, 0);

	c_number rbackr0 = rbackmod - rnamodel.RNA_FENE_R0;

	c_number4 Ftmp;

	if(use_mbf == true && fabsf(rbackr0) > mbf_xmax) {
		// this is the "relax" potential, i.e. the standard FENE up to xmax and then something like A + B log(r) for r>xmax
		c_number fene_xmax = -(rnamodel.RNA_FENE_EPS / 2.f) * logf(1.f - mbf_xmax * mbf_xmax / rnamodel.RNA_FENE_DELTA2);
		c_number mbf_fmax = (rnamodel.RNA_FENE_EPS * mbf_xmax / (rnamodel.RNA_FENE_DELTA2 - SQR(mbf_xmax)));
		c_number long_xmax = (mbf_fmax - mbf_finf) * mbf_xmax * logf(mbf_xmax) + mbf_finf * mbf_xmax;
		Ftmp = rback * (copysignf(1.f, rbackr0) * ((mbf_fmax - mbf_finf) * mbf_xmax / fabsf(rbackr0) + mbf_finf) / rbackmod);
		Ftmp.w = (mbf_fmax - mbf_finf) * mbf_xmax * logf(fabsf(rbackr0)) + mbf_finf * fabsf(rbackr0) - long_xmax + fene_xmax;
	}
	else {
		Ftmp = rback * ((rnamodel.RNA_FENE_EPS * rbackr0 / (rnamodel.RNA_FENE_DELTA2 - SQR(rbackr0))) / rbackmod);
		Ftmp.w = -rnamodel.RNA_FENE_EPS * ((c_number) 0.5f) * logf(1 - SQR(rbackr0) / rnamodel.RNA_FENE_DELTA2);
	}
	c_number4 Ttmp = (qIsN3) ? _cross(n5pos_back, Ftmp) : _cross(n3pos_back, Ftmp);
	// EXCLUDED VOLUME
	_bonded_excluded_volume<qIsN3>(r, n3pos_base, n3pos_back, n5pos_base, n5pos_back, Ftmp, Ttmp);

	if(qIsN3) {
		F += Ftmp;
		T += Ttmp;
	}
	else {
		F -= Ftmp;
		T -= Ttmp;
	}

	// STACKING
	c_number4 rstack = r + n3pos_stack_5 - n5pos_stack_3;

	c_number rstackmod = _module(rstack);
	c_number4 rstackdir = make_c_number4(rstack.x / rstackmod, rstack.y / rstackmod, rstack.z / rstackmod, 0);
	// This is the position the backbone would have with major-minor grooves the same width.
	// We need to do this to implement different major-minor groove widths because rback is
	// used as a reference point for things that have nothing to do with the actual backbone
	// position (in this case, the stacking interaction).
	//c_number4 rbackref = r + n3x * rnamodel.RNA_POS_BACK - n5x * rnamodel.RNA_POS_BACK;
	//c_number rbackrefmod = _module(rbackref);

	c_number t4 = CUDA_LRACOS(CUDA_DOT(n3z, n5z));

	c_number t5 = CUDA_LRACOS(CUDA_DOT(n5z, rstackdir));
	c_number t6 = CUDA_LRACOS(-CUDA_DOT(n3z, rstackdir));

	c_number cosphi1 = CUDA_DOT(n5y, rbackdir); // / rbackrefmod;
	c_number cosphi2 = CUDA_DOT(n3y, rbackdir); // / rbackrefmod;

	//c_number costB1 = -rback_back_dir * p->int_centers[RNANucleotide::BBVECTOR_3];
	//c_number costB2 = -rback_back_dir * q->int_centers[RNANucleotide::BBVECTOR_5];

	c_number4 n3bbvector_5 = (n3x * rnamodel.p5_x + n3y * rnamodel.p5_y + n3z * rnamodel.p5_z);
	c_number4 n5bbvector_3 = (n5x * rnamodel.p3_x + n5y * rnamodel.p3_y + n5z * rnamodel.p3_z);

	c_number costB1 = -CUDA_DOT(rbackdir, n5bbvector_3);
	c_number costB2 = -CUDA_DOT(rbackdir, n3bbvector_5);
	c_number tB1 = CUDA_LRACOS(costB1);
	c_number tB2 = CUDA_LRACOS(costB2);

	// functions
	c_number f1 = _f1(rstackmod, STCK_F1, n3type, n5type);
	//c_number f4t4 = _f4(t4, rnamodel.RNA_STCK_THETA4_T0, rnamodel.RNA_STCK_THETA4_TS, rnamodel.RNA_STCK_THETA4_TC, rnamodel.RNA_STCK_THETA4_A, rnamodel.RNA_STCK_THETA4_B);
	c_number f4t5 = _f4(PI - t5, rnamodel.RNA_STCK_THETA5_T0, rnamodel.RNA_STCK_THETA5_TS, rnamodel.RNA_STCK_THETA5_TC, rnamodel.RNA_STCK_THETA5_A, rnamodel.RNA_STCK_THETA5_B);
	c_number f4t6 = _f4(t6, rnamodel.RNA_STCK_THETA6_T0, rnamodel.RNA_STCK_THETA6_TS, rnamodel.RNA_STCK_THETA6_TC, rnamodel.RNA_STCK_THETA6_A, rnamodel.RNA_STCK_THETA6_B);

	c_number f4tB1 = _f4(tB1, rnamodel.STCK_THETAB1_T0, rnamodel.STCK_THETAB1_TS, rnamodel.STCK_THETAB1_TC, rnamodel.STCK_THETAB1_A, rnamodel.STCK_THETAB1_B);
	c_number f4tB2 = _f4(tB2, rnamodel.STCK_THETAB2_T0, rnamodel.STCK_THETAB2_TS, rnamodel.STCK_THETAB2_TC, rnamodel.STCK_THETAB2_A, rnamodel.STCK_THETAB2_B);

	c_number f5phi1 = _f5(cosphi1, STCK_F5_PHI1);
	c_number f5phi2 = _f5(cosphi2, STCK_F5_PHI2);

	c_number energy = f1 * f4tB1 * f4tB2 * f4t5 * f4t6 * f5phi1 * f5phi2;

	//c_number4 pre_Ttmp = Ttmp;
	//c_number4 pre_Ftmp = Ftmp;

	if(energy != (c_number) 0) {
		// and their derivatives
		c_number f1D = _f1D(rstackmod, STCK_F1, n3type, n5type);
		//c_number f4t4Dsin = _f4Dsin(t4, rnamodel.RNA_STCK_THETA4_T0, rnamodel.RNA_STCK_THETA4_TS, rnamodel.RNA_STCK_THETA4_TC, rnamodel.RNA_STCK_THETA4_A, rnamodel.RNA_STCK_THETA4_B);
		c_number f4t5Dsin = _f4Dsin(PI - t5, rnamodel.RNA_STCK_THETA5_T0, rnamodel.RNA_STCK_THETA5_TS, rnamodel.RNA_STCK_THETA5_TC, rnamodel.RNA_STCK_THETA5_A, rnamodel.RNA_STCK_THETA5_B);
		c_number f4t6Dsin = _f4Dsin(t6, rnamodel.RNA_STCK_THETA6_T0, rnamodel.RNA_STCK_THETA6_TS, rnamodel.RNA_STCK_THETA6_TC, rnamodel.RNA_STCK_THETA6_A, rnamodel.RNA_STCK_THETA6_B);

		c_number f5phi1D = _f5D(cosphi1, STCK_F5_PHI1);
		c_number f5phi2D = _f5D(cosphi2, STCK_F5_PHI2);

		c_number f4tB1Dsin = _f4Dsin(tB1, rnamodel.STCK_THETAB1_T0, rnamodel.STCK_THETAB1_TS, rnamodel.STCK_THETAB1_TC, rnamodel.STCK_THETAB1_A, rnamodel.STCK_THETAB1_B);
		c_number f4tB2Dsin = _f4Dsin(tB2, rnamodel.STCK_THETAB2_T0, rnamodel.STCK_THETAB2_TS, rnamodel.STCK_THETAB2_TC, rnamodel.STCK_THETAB2_A, rnamodel.STCK_THETAB2_B);

		// RADIAL
		Ftmp = -rstackdir * (energy * f1D / f1);

		//if(grooving) printf("Amydevice: (qisN3 = %d) calculated Ftmp=(%f,%f,%f) \n",qIsN3,Ftmp.x,Ftmp.y,Ftmp.z);
		// THETA 5
		Ftmp -= (n5z - cosf(t5) * rstackdir) * (f1 * f4t5Dsin * f4t6 * f5phi1 * f5phi2 * f4tB1 * f4tB2 / rstackmod); //(n5z - cosf(t5) * rstackdir) * (energy * f4t5Dsin / (f4t5 * rstackmod));
		// THETA 6
		Ftmp -= (n3z + cosf(t6) * rstackdir) * (energy * f4t6Dsin / (f4t6 * rstackmod));
		//if(grooving) printf("mydevice: (qisN3 = %d) calculated Ftmp=(%f,%f,%f) \n ",qIsN3,Ftmp.x,Ftmp.y,Ftmp.z);

		c_number4 Fposback = -(n5bbvector_3 + rbackdir * costB1) * (f1 * f4t5 * f4t6 * f5phi1 * f5phi2 * f4tB1Dsin * f4tB2 / rbackmod);
		Fposback -= (n3bbvector_5 + rbackdir * costB2) * (f1 * f4t5 * f4t6 * f5phi1 * f5phi2 * f4tB1 * f4tB2Dsin / rbackmod);

		//force acting on the backbone site

		// COS PHI 1
		// here particle p is referred to using the a while particle q is referred with the b
		/*
		 c_number ra2 = CUDA_DOT(rstackdir, n5y);
		 c_number ra1 = CUDA_DOT(rstackdir, n5x);
		 c_number rb1 = CUDA_DOT(rstackdir, n3x);
		 c_number a2b1 = CUDA_DOT(n5y, n3x);
		 c_number dcosphi1dr = (SQR(rstackmod)*ra2 - ra2*SQR(rbackrefmod) - rstackmod*(a2b1 + ra2*(-ra1 + rb1))*rnamodel.RNA_GAMMA + a2b1*(-ra1 + rb1)*SQR(rnamodel.RNA_GAMMA))/(SQR(rbackrefmod)*rbackrefmod);
		 c_number dcosphi1dra1 = rstackmod*rnamodel.RNA_GAMMA*(rstackmod*ra2 - a2b1*rnamodel.RNA_GAMMA)/(SQR(rbackrefmod)*rbackrefmod);
		 c_number dcosphi1dra2 = -rstackmod / rbackrefmod;
		 c_number dcosphi1drb1 = -(rstackmod*rnamodel.RNA_GAMMA*(rstackmod*ra2 - a2b1*rnamodel.RNA_GAMMA))/(SQR(rbackrefmod)*rbackrefmod);
		 c_number dcosphi1da1b1 = SQR(rnamodel.RNA_GAMMA)*(-rstackmod*ra2 + a2b1*rnamodel.RNA_GAMMA)/(SQR(rbackrefmod)*rbackrefmod);
		 c_number dcosphi1da2b1 = rnamodel.RNA_GAMMA / rbackrefmod;
		 */
		//c_number force_part_phi1 = energy * f5phi1D / f5phi1;
		c_number force_part_phi1 = -f1 * f4t5 * f4t6 * f5phi1D * f5phi2 * f4tB1 * f4tB2;
		Fposback += (n5y - rback * (cosphi1 / rbackmod)) * (force_part_phi1 / rbackmod);

		/*
		 Ftmp -= (rstackdir * dcosphi1dr +
		 ((n5y - ra2*rstackdir) * dcosphi1dra2 +
		 (n5x - ra1*rstackdir) * dcosphi1dra1 +
		 (n3x - rb1*rstackdir) * dcosphi1drb1) / rstackmod) * force_part_phi1;

		 // COS PHI 2
		 // here particle p -> b, particle q -> a
		 ra2 = CUDA_DOT(rstackdir, n3y);
		 ra1 = rb1;
		 rb1 = CUDA_DOT(rstackdir, n5x);
		 a2b1 = CUDA_DOT(n3y, n5x);
		 c_number dcosphi2dr = ((rstackmod*ra2 + a2b1*rnamodel.RNA_GAMMA)*(rstackmod + (rb1 - ra1)*rnamodel.RNA_GAMMA) - ra2*SQR(rbackrefmod))/(SQR(rbackrefmod)*rbackrefmod);
		 c_number dcosphi2dra1 = -rstackmod*rnamodel.RNA_GAMMA*(rstackmod*ra2 + a2b1*rnamodel.RNA_GAMMA)/(SQR(rbackrefmod)*rbackrefmod);
		 c_number dcosphi2dra2 = -rstackmod / rbackrefmod;
		 c_number dcosphi2drb1 = (rstackmod*rnamodel.RNA_GAMMA*(rstackmod*ra2 + a2b1*rnamodel.RNA_GAMMA))/(SQR(rbackrefmod)*rbackrefmod);
		 c_number dcosphi2da1b1 = -SQR(rnamodel.RNA_GAMMA)*(rstackmod*ra2 + a2b1*rnamodel.RNA_GAMMA)/(SQR(rbackrefmod)*rbackrefmod);
		 c_number dcosphi2da2b1 = -rnamodel.RNA_GAMMA / rbackrefmod;

		 c_number force_part_phi2 = energy * f5phi2D / f5phi2;

		 Ftmp -= (rstackdir * dcosphi2dr +
		 ((n3y - rstackdir * ra2) * dcosphi2dra2 +
		 (n3x - rstackdir * ra1) * dcosphi2dra1 +
		 (n5x - rstackdir * rb1) * dcosphi2drb1) / rstackmod) * force_part_phi2;

		 */

		c_number force_part_phi2 = -f1 * f4t5 * f4t6 * f5phi1 * f5phi2D;
		Fposback += (n3y - rback * (cosphi2 / rbackmod)) * (force_part_phi2 / rbackmod);
		//if(grooving) printf("mydevice: (qisN3 = %d) calculated Fposback=(%f,%f,%f) \n ",qIsN3,Fposback.x,Fposback.y,Fposback.z);

		// TORQUE
		if(qIsN3) {
			Ttmp = -_cross(n5pos_stack_3, Ftmp);
			Ttmp -= _cross(n5pos_back, Fposback);

			// THETA 5
			Ttmp += _cross(rstackdir, n5z) * (f1 * f4t5Dsin * f4t6 * f5phi1 * f5phi2 * f4tB1 * f4tB2);

			//thetaB1
			Ttmp += _cross(rbackdir, n5bbvector_3) * (f1 * f4t5 * f4t6 * f5phi1 * f5phi2 * f4tB1Dsin * f4tB2);

			Ttmp += _cross(n5y, rbackdir) * force_part_phi1;

		}
		else {

			Ttmp = _cross(n3pos_stack_5, Ftmp);
			Ttmp += _cross(n3pos_back, Fposback);
			Ttmp += _cross(n3y, rbackdir) * force_part_phi2;

			// THETA 6
			//printf("Entered function with qIsN3 = %d\n",(int)qIsN3);
			Ttmp += _cross(rstackdir, n3z) * f1 * f4t5 * f4t6Dsin * f5phi1 * f5phi2 * f4tB1 * f4tB2;

			//theta B2
			Ttmp += _cross(rbackdir, n3bbvector_5) * (f1 * f4t5 * f4t6 * f5phi1 * f5phi2 * f4tB1 * f4tB2Dsin);

		}

		Ftmp.w = energy;

		//printf("Approaching decision loop with qIsN3 = %d\n",(int)qIsN3);
		if(qIsN3) {

			Ftmp += Fposback;

			T += Ttmp;
			F -= Ftmp;

		}
		else {
			Ftmp += Fposback;

			T += Ttmp;
			F += Ftmp;
		}

	}

	//Ftmp = F;
	//if(grooving)
	//  printf("This is bonded function that calculated F=(%f,%f,%f) and T=(%f,%f,%f)  and qIsN3 is %d \n",Ftmp.x,Ftmp.y,Ftmp.z,Ttmp.x,Ttmp.y,Ttmp.z,(int)qIsN3);

	//printf("Goodbye from bonded part function \n");
}

__device__
void _particle_particle_interaction(c_number4 ppos, c_number4 a1, c_number4 a2, c_number4 a3, c_number4 qpos, c_number4 b1, c_number4 b2, c_number4 b3, c_number4 &F, c_number4 &T, bool average, bool use_debye_huckel, bool mismatch_repulsion, LR_bonds pbonds, LR_bonds qbonds, CUDABox *box) {
	int ptype = get_particle_type(ppos);
	int qtype = get_particle_type(qpos);
	int pbtype = get_particle_btype(ppos);
	int qbtype = get_particle_btype(qpos);
	int int_type = pbtype + qbtype;

	c_number4 r = box->minimum_image(ppos, qpos);

	c_number4 ppos_back = a1 * rnamodel.RNA_POS_BACK_a1 + a2 * rnamodel.RNA_POS_BACK_a2 + a3 * rnamodel.RNA_POS_BACK_a3;
	//if(grooving) ppos_back = POS_MM_BACK1 * a1 + POS_MM_BACK2 * a2;
	//else ppos_back = rnamodel.RNA_POS_BACK * a1;

	c_number4 ppos_base = rnamodel.RNA_POS_BASE * a1;
	c_number4 ppos_stack = rnamodel.RNA_POS_STACK * a1;

	c_number4 qpos_back = b1 * rnamodel.RNA_POS_BACK_a1 + b2 * rnamodel.RNA_POS_BACK_a2 + b3 * rnamodel.RNA_POS_BACK_a3;
	//if(grooving) qpos_back = POS_MM_BACK1 * b1 + POS_MM_BACK2 * b2;
	//else qpos_back = rnamodel.RNA_POS_BACK * b1;

	c_number4 qpos_base = rnamodel.RNA_POS_BASE * b1;
	c_number4 qpos_stack = rnamodel.RNA_POS_STACK * b1;

	c_number old_Tw = T.w;

	// excluded volume
	// BACK-BACK
	c_number4 Ftmp = make_c_number4(0, 0, 0, 0);
	c_number4 rbackbone = r + qpos_back - ppos_back;
	_excluded_volume(rbackbone, Ftmp, rnamodel.RNA_EXCL_S1, rnamodel.RNA_EXCL_R1, rnamodel.RNA_EXCL_B1, rnamodel.RNA_EXCL_RC1);
	c_number4 Ttmp = _cross(ppos_back, Ftmp);
	_bonded_excluded_volume<true>(r, qpos_base, qpos_back, ppos_base, ppos_back, Ftmp, Ttmp);

	F += Ftmp;

	// HYDROGEN BONDING
	int is_pair = int(int_type == 3);
	if(!average) //allow for wobble bp
	{
		if(int_type == 4 && ((qtype == N_T && ptype == N_G) || (qtype == N_G && ptype == N_T))) is_pair = 1;
	}

	c_number hb_energy = (c_number) 0;
	c_number4 rhydro = r + qpos_base - ppos_base;
	c_number rhydromodsqr = CUDA_DOT(rhydro, rhydro);
	if((is_pair || mismatch_repulsion) && SQR(rnamodel.RNA_HYDR_RCLOW) < rhydromodsqr && rhydromodsqr < SQR(rnamodel.RNA_HYDR_RCHIGH)) {
		c_number hb_multi = (abs(qbtype) >= 300 && abs(pbtype) >= 300) ? MD_hb_multi[0] : 1.f;
		// versor and magnitude of the base-base separation
		c_number rhydromod = sqrtf(rhydromodsqr);
		c_number4 rhydrodir = rhydro / rhydromod;

		// angles involved in the HB interaction
		c_number t1 = CUDA_LRACOS(-CUDA_DOT(a1, b1));
		c_number t2 = CUDA_LRACOS(-CUDA_DOT(b1, rhydrodir));
		c_number t3 = CUDA_LRACOS(CUDA_DOT(a1, rhydrodir));
		c_number t4 = CUDA_LRACOS(CUDA_DOT(a3, b3));
		c_number t7 = CUDA_LRACOS(-CUDA_DOT(rhydrodir, b3));
		c_number t8 = CUDA_LRACOS(CUDA_DOT(rhydrodir, a3));

		// functions called at their relevant arguments
		c_number f1 = hb_multi * _f1(rhydromod, RNA_HYDR_F1, ptype, qtype);
		if(mismatch_repulsion && !is_pair) {
			f1 = _fX(rhydromod, RNA_HYDR_F1, 0, 0);
		}
		c_number f4t1 = _f4(t1, rnamodel.RNA_HYDR_THETA1_T0, rnamodel.RNA_HYDR_THETA1_TS, rnamodel.RNA_HYDR_THETA1_TC, rnamodel.RNA_HYDR_THETA1_A, rnamodel.RNA_HYDR_THETA1_B);
		c_number f4t2 = _f4(t2, rnamodel.RNA_HYDR_THETA2_T0, rnamodel.RNA_HYDR_THETA2_TS, rnamodel.RNA_HYDR_THETA2_TC, rnamodel.RNA_HYDR_THETA2_A, rnamodel.RNA_HYDR_THETA2_B);
		c_number f4t3 = _f4(t3, rnamodel.RNA_HYDR_THETA3_T0, rnamodel.RNA_HYDR_THETA3_TS, rnamodel.RNA_HYDR_THETA3_TC, rnamodel.RNA_HYDR_THETA3_A, rnamodel.RNA_HYDR_THETA3_B);
		c_number f4t4 = _f4(t4, rnamodel.RNA_HYDR_THETA4_T0, rnamodel.RNA_HYDR_THETA4_TS, rnamodel.RNA_HYDR_THETA4_TC, rnamodel.RNA_HYDR_THETA4_A, rnamodel.RNA_HYDR_THETA4_B);
		c_number f4t7 = _f4(t7, rnamodel.RNA_HYDR_THETA7_T0, rnamodel.RNA_HYDR_THETA7_TS, rnamodel.RNA_HYDR_THETA7_TC, rnamodel.RNA_HYDR_THETA7_A, rnamodel.RNA_HYDR_THETA7_B);
		c_number f4t8 = _f4(t8, rnamodel.RNA_HYDR_THETA8_T0, rnamodel.RNA_HYDR_THETA8_TS, rnamodel.RNA_HYDR_THETA8_TC, rnamodel.RNA_HYDR_THETA8_A, rnamodel.RNA_HYDR_THETA8_B);

		hb_energy = f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;

		if(hb_energy < (c_number) 0 || (hb_energy > 0 && mismatch_repulsion && !is_pair)) {
			// derivatives called at the relevant arguments
			c_number f1D = hb_multi * _f1D(rhydromod, HYDR_F1, ptype, qtype);
			if(mismatch_repulsion && !is_pair) {
				f1D = _fXD(rhydromod, RNA_HYDR_F1, 0, 0);
			}
			c_number f4t1Dsin = -_f4Dsin(t1, rnamodel.RNA_HYDR_THETA1_T0, rnamodel.RNA_HYDR_THETA1_TS, rnamodel.RNA_HYDR_THETA1_TC, rnamodel.RNA_HYDR_THETA1_A, rnamodel.RNA_HYDR_THETA1_B);
			c_number f4t2Dsin = -_f4Dsin(t2, rnamodel.RNA_HYDR_THETA2_T0, rnamodel.RNA_HYDR_THETA2_TS, rnamodel.RNA_HYDR_THETA2_TC, rnamodel.RNA_HYDR_THETA2_A, rnamodel.RNA_HYDR_THETA2_B);
			c_number f4t3Dsin = _f4Dsin(t3, rnamodel.RNA_HYDR_THETA3_T0, rnamodel.RNA_HYDR_THETA3_TS, rnamodel.RNA_HYDR_THETA3_TC, rnamodel.RNA_HYDR_THETA3_A, rnamodel.RNA_HYDR_THETA3_B);
			c_number f4t4Dsin = _f4Dsin(t4, rnamodel.RNA_HYDR_THETA4_T0, rnamodel.RNA_HYDR_THETA4_TS, rnamodel.RNA_HYDR_THETA4_TC, rnamodel.RNA_HYDR_THETA4_A, rnamodel.RNA_HYDR_THETA4_B);
			c_number f4t7Dsin = -_f4Dsin(t7, rnamodel.RNA_HYDR_THETA7_T0, rnamodel.RNA_HYDR_THETA7_TS, rnamodel.RNA_HYDR_THETA7_TC, rnamodel.RNA_HYDR_THETA7_A, rnamodel.RNA_HYDR_THETA7_B);
			c_number f4t8Dsin = _f4Dsin(t8, rnamodel.RNA_HYDR_THETA8_T0, rnamodel.RNA_HYDR_THETA8_TS, rnamodel.RNA_HYDR_THETA8_TC, rnamodel.RNA_HYDR_THETA8_A, rnamodel.RNA_HYDR_THETA8_B);

			// RADIAL PART
			Ftmp = rhydrodir * hb_energy * f1D / f1;

			// TETA4; t4 = LRACOS (a3 * b3);
			Ttmp -= _cross(a3, b3) * (-hb_energy * f4t4Dsin / f4t4);

			// TETA1; t1 = LRACOS (-a1 * b1);
			Ttmp -= _cross(a1, b1) * (-hb_energy * f4t1Dsin / f4t1);

			// TETA2; t2 = LRACOS (-b1 * rhydrodir);
			Ftmp -= (b1 + rhydrodir * cosf(t2)) * (hb_energy * f4t2Dsin / (f4t2 * rhydromod));

			// TETA3; t3 = LRACOS (a1 * rhydrodir);
			c_number part = -hb_energy * f4t3Dsin / f4t3;
			Ftmp -= (a1 - rhydrodir * cosf(t3)) * (-part / rhydromod);
			Ttmp += _cross(rhydrodir, a1) * part;

			// THETA7; t7 = LRACOS (-rhydrodir * b3);
			Ftmp -= (b3 + rhydrodir * cosf(t7)) * (hb_energy * f4t7Dsin / (f4t7 * rhydromod));

			// THETA 8; t8 = LRACOS (rhydrodir * a3);
			part = -hb_energy * f4t8Dsin / f4t8;
			Ftmp -= (a3 - rhydrodir * cosf(t8)) * (-part / rhydromod);
			Ttmp += _cross(rhydrodir, a3) * part;

			Ttmp += _cross(ppos_base, Ftmp);

			Ftmp.w = hb_energy;
			F += Ftmp;
		}
	}
	// END HYDROGEN BONDING

	// CROSS STACKING
	c_number4 rcstack = rhydro;
	c_number rcstackmodsqr = rhydromodsqr;
	if(SQR(rnamodel.RNA_CRST_RCLOW) < rcstackmodsqr && rcstackmodsqr < SQR(rnamodel.RNA_CRST_RCHIGH)) {
		c_number rcstackmod = sqrtf(rcstackmodsqr);
		c_number4 rcstackdir = rcstack / rcstackmod;

		// angles involved in the CSTCK interaction
		c_number t1 = CUDA_LRACOS(-CUDA_DOT(a1, b1));

		c_number t2 = CUDA_LRACOS(-CUDA_DOT(b1, rcstackdir));
		c_number t3 = CUDA_LRACOS(CUDA_DOT(a1, rcstackdir));

		//c_number t4 = CUDA_LRACOS ( CUDA_DOT(a3, b3));

		c_number t7 = CUDA_LRACOS(-CUDA_DOT(rcstackdir, b3));
		c_number t8 = CUDA_LRACOS(CUDA_DOT(rcstackdir, a3));

		// functions called at their relevant arguments
		c_number f2 = _f2(rcstackmod, RNA_CRST_F2);
		c_number f4t1 = _f4(t1, rnamodel.RNA_CRST_THETA1_T0, rnamodel.RNA_CRST_THETA1_TS, rnamodel.RNA_CRST_THETA1_TC, rnamodel.RNA_CRST_THETA1_A, rnamodel.RNA_CRST_THETA1_B);
		c_number f4t2 = _f4(t2, rnamodel.RNA_CRST_THETA2_T0, rnamodel.RNA_CRST_THETA2_TS, rnamodel.RNA_CRST_THETA2_TC, rnamodel.RNA_CRST_THETA2_A, rnamodel.RNA_CRST_THETA2_B);
		c_number f4t3 = _f4(t3, rnamodel.RNA_CRST_THETA3_T0, rnamodel.RNA_CRST_THETA3_TS, rnamodel.RNA_CRST_THETA3_TC, rnamodel.RNA_CRST_THETA3_A, rnamodel.RNA_CRST_THETA3_B);
		//	c_number f4t4 = _f4(t4, rnamodel.RNA_CRST_THETA4_T0, rnamodel.RNA_CRST_THETA4_TS, rnamodel.RNA_CRST_THETA4_TC, rnamodel.RNA_CRST_THETA4_A, rnamodel.RNA_CRST_THETA4_B) +
		//		_f4(PI - t4, rnamodel.RNA_CRST_THETA4_T0, rnamodel.RNA_CRST_THETA4_TS, rnamodel.RNA_CRST_THETA4_TC, rnamodel.RNA_CRST_THETA4_A, rnamodel.RNA_CRST_THETA4_B);
		c_number f4t7 = _f4(t7, rnamodel.RNA_CRST_THETA7_T0, rnamodel.RNA_CRST_THETA7_TS, rnamodel.RNA_CRST_THETA7_TC, rnamodel.RNA_CRST_THETA7_A, rnamodel.RNA_CRST_THETA7_B) + _f4(PI - t7, rnamodel.RNA_CRST_THETA7_T0, rnamodel.RNA_CRST_THETA7_TS, rnamodel.RNA_CRST_THETA7_TC, rnamodel.RNA_CRST_THETA7_A, rnamodel.RNA_CRST_THETA7_B);
		c_number f4t8 = _f4(t8, rnamodel.RNA_CRST_THETA8_T0, rnamodel.RNA_CRST_THETA8_TS, rnamodel.RNA_CRST_THETA8_TC, rnamodel.RNA_CRST_THETA8_A, rnamodel.RNA_CRST_THETA8_B) + _f4(PI - t8, rnamodel.RNA_CRST_THETA8_T0, rnamodel.RNA_CRST_THETA8_TS, rnamodel.RNA_CRST_THETA8_TC, rnamodel.RNA_CRST_THETA8_A, rnamodel.RNA_CRST_THETA8_B);

		c_number cstk_energy = f2 * f4t1 * f4t2 * f4t3 * f4t7 * f4t8;

		if(cstk_energy < (c_number) 0) {
			// derivatives called at the relevant arguments
			c_number f2D = _f2D(rcstackmod, CRST_F2);
			c_number f4t1Dsin = -_f4Dsin(t1, rnamodel.RNA_CRST_THETA1_T0, rnamodel.RNA_CRST_THETA1_TS, rnamodel.RNA_CRST_THETA1_TC, rnamodel.RNA_CRST_THETA1_A, rnamodel.RNA_CRST_THETA1_B);
			c_number f4t2Dsin = -_f4Dsin(t2, rnamodel.RNA_CRST_THETA2_T0, rnamodel.RNA_CRST_THETA2_TS, rnamodel.RNA_CRST_THETA2_TC, rnamodel.RNA_CRST_THETA2_A, rnamodel.RNA_CRST_THETA2_B);
			c_number f4t3Dsin = _f4Dsin(t3, rnamodel.RNA_CRST_THETA3_T0, rnamodel.RNA_CRST_THETA3_TS, rnamodel.RNA_CRST_THETA3_TC, rnamodel.RNA_CRST_THETA3_A, rnamodel.RNA_CRST_THETA3_B);
			//c_number f4t4Dsin = _f4Dsin(t4, rnamodel.RNA_CRST_THETA4_T0, rnamodel.RNA_CRST_THETA4_TS, rnamodel.RNA_CRST_THETA4_TC, rnamodel.RNA_CRST_THETA4_A, rnamodel.RNA_CRST_THETA4_B) -
			//		_f4Dsin(PI - t4, rnamodel.RNA_CRST_THETA4_T0, rnamodel.RNA_CRST_THETA4_TS, rnamodel.RNA_CRST_THETA4_TC, rnamodel.RNA_CRST_THETA4_A, rnamodel.RNA_CRST_THETA4_B);
			c_number f4t7Dsin = -_f4Dsin(t7, rnamodel.RNA_CRST_THETA7_T0, rnamodel.RNA_CRST_THETA7_TS, rnamodel.RNA_CRST_THETA7_TC, rnamodel.RNA_CRST_THETA7_A, rnamodel.RNA_CRST_THETA7_B) + _f4Dsin(PI - t7, rnamodel.RNA_CRST_THETA7_T0, rnamodel.RNA_CRST_THETA7_TS, rnamodel.RNA_CRST_THETA7_TC, rnamodel.RNA_CRST_THETA7_A, rnamodel.RNA_CRST_THETA7_B);
			c_number f4t8Dsin = _f4Dsin(t8, rnamodel.RNA_CRST_THETA8_T0, rnamodel.RNA_CRST_THETA8_TS, rnamodel.RNA_CRST_THETA8_TC, rnamodel.RNA_CRST_THETA8_A, rnamodel.RNA_CRST_THETA8_B) - _f4Dsin(PI - t8, rnamodel.RNA_CRST_THETA8_T0, rnamodel.RNA_CRST_THETA8_TS, rnamodel.RNA_CRST_THETA8_TC, rnamodel.RNA_CRST_THETA8_A, rnamodel.RNA_CRST_THETA8_B);

			// RADIAL PART
			Ftmp = rcstackdir * (cstk_energy * f2D / f2);

			// THETA1; t1 = LRACOS (-a1 * b1);
			Ttmp -= _cross(a1, b1) * (-cstk_energy * f4t1Dsin / f4t1);

			// TETA2; t2 = LRACOS (-b1 * rhydrodir);
			Ftmp -= (b1 + rcstackdir * cosf(t2)) * (cstk_energy * f4t2Dsin / (f4t2 * rcstackmod));

			// TETA3; t3 = LRACOS (a1 * rhydrodir);
			c_number part = -cstk_energy * f4t3Dsin / f4t3;
			Ftmp -= (a1 - rcstackdir * cosf(t3)) * (-part / rcstackmod);
			Ttmp += _cross(rcstackdir, a1) * part;

			// TETA4; t4 = LRACOS (a3 * b3);
			//Ttmp -= _cross(a3, b3) * (-cstk_energy * f4t4Dsin / f4t4);

			// THETA7; t7 = LRACOS (-rcsrackir * b3);
			Ftmp -= (b3 + rcstackdir * cosf(t7)) * (cstk_energy * f4t7Dsin / (f4t7 * rcstackmod));

			// THETA 8; t8 = LRACOS (rhydrodir * a3);
			part = -cstk_energy * f4t8Dsin / f4t8;
			Ftmp -= (a3 - rcstackdir * cosf(t8)) * (-part / rcstackmod);
			Ttmp += _cross(rcstackdir, a3) * part;

			Ttmp += _cross(ppos_base, Ftmp);

			Ftmp.w = cstk_energy;
			F += Ftmp;
		}
	}

	// COAXIAL STACKING
	c_number4 rstack = r + qpos_stack - ppos_stack;
	c_number rstackmodsqr = CUDA_DOT(rstack, rstack);
	if(SQR(rnamodel.RNA_CXST_RCLOW) < rstackmodsqr && rstackmodsqr < SQR(rnamodel.RNA_CXST_RCHIGH)) {
		c_number rstackmod = sqrtf(rstackmodsqr);
		c_number4 rstackdir = rstack / rstackmod;

		// angles involved in the CXST interaction
		c_number t1 = CUDA_LRACOS(-CUDA_DOT(a1, b1));
		c_number t4 = CUDA_LRACOS(CUDA_DOT(a3, b3));
		c_number t5 = CUDA_LRACOS(CUDA_DOT(a3, rstackdir));
		c_number t6 = CUDA_LRACOS(-CUDA_DOT(b3, rstackdir));

		// This is the position the backbone would have with major-minor grooves the same width.
		// We need to do this to implement different major-minor groove widths because rback is
		// used as a reference point for things that have nothing to do with the actual backbone
		// position (in this case, the coaxial stacking interaction).
		//c_number4 rbackboneref = r + rnamodel.RNA_POS_BACK * b1 - rnamodel.RNA_POS_BACK * a1;
		c_number rbackmod = _module(rbackbone);
		c_number4 rbackbonedir = rbackbone / rbackmod;

		c_number cosphi3 = CUDA_DOT(rstackdir, (_cross(rbackbonedir, a1)));
		c_number cosphi4 = CUDA_DOT(rstackdir, (_cross(rbackbonedir, b1)));

		// functions called at their relevant arguments
		c_number f2 = _f2(rstackmod, CXST_F2);
		c_number f4t1 = _f4(t1, rnamodel.RNA_CXST_THETA1_T0, rnamodel.RNA_CXST_THETA1_TS, rnamodel.RNA_CXST_THETA1_TC, rnamodel.RNA_CXST_THETA1_A, rnamodel.RNA_CXST_THETA1_B) + _f4(2 * PI - t1, rnamodel.RNA_CXST_THETA1_T0, rnamodel.RNA_CXST_THETA1_TS, rnamodel.RNA_CXST_THETA1_TC, rnamodel.RNA_CXST_THETA1_A, rnamodel.RNA_CXST_THETA1_B);
		c_number f4t4 = _f4(t4, rnamodel.RNA_CXST_THETA4_T0, rnamodel.RNA_CXST_THETA4_TS, rnamodel.RNA_CXST_THETA4_TC, rnamodel.RNA_CXST_THETA4_A, rnamodel.RNA_CXST_THETA4_B);
		c_number f4t5 = _f4(t5, rnamodel.RNA_CXST_THETA5_T0, rnamodel.RNA_CXST_THETA5_TS, rnamodel.RNA_CXST_THETA5_TC, rnamodel.RNA_CXST_THETA5_A, rnamodel.RNA_CXST_THETA5_B) + _f4(PI - t5, rnamodel.RNA_CXST_THETA5_T0, rnamodel.RNA_CXST_THETA5_TS, rnamodel.RNA_CXST_THETA5_TC, rnamodel.RNA_CXST_THETA5_A, rnamodel.RNA_CXST_THETA5_B);
		c_number f4t6 = _f4(t6, rnamodel.RNA_CXST_THETA6_T0, rnamodel.RNA_CXST_THETA6_TS, rnamodel.RNA_CXST_THETA6_TC, rnamodel.RNA_CXST_THETA6_A, rnamodel.RNA_CXST_THETA6_B) + _f4(PI - t6, rnamodel.RNA_CXST_THETA6_T0, rnamodel.RNA_CXST_THETA6_TS, rnamodel.RNA_CXST_THETA6_TC, rnamodel.RNA_CXST_THETA6_A, rnamodel.RNA_CXST_THETA6_B);
		c_number f5cosphi3 = _f5(cosphi3, RNA_CXST_F5_PHI3);

		c_number f5cosphi4 = _f5(cosphi4, RNA_CXST_F5_PHI4);

		c_number cxst_energy = f2 * f4t1 * f4t4 * f4t5 * f4t6 * f5cosphi3 * f5cosphi4;

		//if(grooving)  printf("AADEVICE: Energy =  %f,  f2 = %f, f4t1=%f, f4t4=%f, f4t5=%f, f4t6=%f, f5cosphi3=%f, f5cosphi4=%f \n",cxst_energy,f2 , f4t1 , f4t4 ,f4t5 , f4t6 ,f5cosphi3 , f5cosphi4);

		//
		if(cxst_energy < (c_number) 0) {
			// derivatives called at the relevant arguments
			c_number f2D = _f2D(rstackmod, RNA_CXST_F2);
			c_number f4t1Dsin = -_f4Dsin(t1, rnamodel.RNA_CXST_THETA1_T0, rnamodel.RNA_CXST_THETA1_TS, rnamodel.RNA_CXST_THETA1_TC, rnamodel.RNA_CXST_THETA1_A, rnamodel.RNA_CXST_THETA1_B) + _f4Dsin(2 * PI - t1, rnamodel.RNA_CXST_THETA1_T0, rnamodel.RNA_CXST_THETA1_TS, rnamodel.RNA_CXST_THETA1_TC, rnamodel.RNA_CXST_THETA1_A, rnamodel.RNA_CXST_THETA1_B);
			c_number f4t4Dsin = _f4Dsin(t4, rnamodel.RNA_CXST_THETA4_T0, rnamodel.RNA_CXST_THETA4_TS, rnamodel.RNA_CXST_THETA4_TC, rnamodel.RNA_CXST_THETA4_A, rnamodel.RNA_CXST_THETA4_B);
			c_number f4t5Dsin = _f4Dsin(t5, rnamodel.RNA_CXST_THETA5_T0, rnamodel.RNA_CXST_THETA5_TS, rnamodel.RNA_CXST_THETA5_TC, rnamodel.RNA_CXST_THETA5_A, rnamodel.RNA_CXST_THETA5_B) - _f4Dsin(PI - t5, rnamodel.RNA_CXST_THETA5_T0, rnamodel.RNA_CXST_THETA5_TS, rnamodel.RNA_CXST_THETA5_TC, rnamodel.RNA_CXST_THETA5_A, rnamodel.RNA_CXST_THETA5_B);
			c_number f4t6Dsin = -_f4Dsin(t6, rnamodel.RNA_CXST_THETA6_T0, rnamodel.RNA_CXST_THETA6_TS, rnamodel.RNA_CXST_THETA6_TC, rnamodel.RNA_CXST_THETA6_A, rnamodel.RNA_CXST_THETA6_B) + _f4Dsin(PI - t6, rnamodel.RNA_CXST_THETA6_T0, rnamodel.RNA_CXST_THETA6_TS, rnamodel.RNA_CXST_THETA6_TC, rnamodel.RNA_CXST_THETA6_A, rnamodel.RNA_CXST_THETA6_B);
			c_number f5Dcosphi3 = _f5D(cosphi3, RNA_CXST_F5_PHI3);
			c_number f5Dcosphi4 = _f5D(cosphi4, RNA_CXST_F5_PHI4);

			// RADIAL PART
			Ftmp = -rstackdir * (cxst_energy * f2D / f2);

			//if(grooving)  printf("DEVICE radial force: Ftmp = (%f,%f,%f) with stackdir = (%f,%f,%f), a3= (%f,%f,%f) \n",Ftmp.x,Ftmp.y,Ftmp.z,rstackdir.x,rstackdir.y,rstackdir.z,a3.x,a3.y,a3.z);
			//if(grooving)  printf("AADEVICE:, with rstackmod=%g, Force terms f2 = %g, f4t1=%g, f4t4=%g, f4t5=%g, f4t6=%g, f5cosphi3=%g, f5cosphi4=%g \n",rstackmod,f2D , f4t1Dsin , f4t4Dsin ,f4t5Dsin , f4t6Dsin ,f5Dcosphi3 , f5Dcosphi4);

			// THETA1; t1 = LRACOS (-a1 * b1);
			Ttmp -= _cross(a1, b1) * (-cxst_energy * f4t1Dsin / f4t1);

			// TETA4; t4 = LRACOS (a3 * b3);
			Ttmp -= _cross(a3, b3) * (-cxst_energy * f4t4Dsin / f4t4);

			// THETA5; t5 = LRACOS ( a3 * rstackdir);
			c_number part = cxst_energy * f4t5Dsin / f4t5;
			Ftmp += (a3 - rstackdir * cosf(t5)) / rstackmod * part;
			Ttmp -= _cross(rstackdir, a3) * part;
			//if(grooving)  printf("DEVICE theta5 force: Ftmp = (%f,%f,%f) with stackdir = (%f,%f,%f), a3= (%f,%f,%f) \n",Ftmp.x,Ftmp.y,Ftmp.z,rstackdir.x,rstackdir.y,rstackdir.z,a3.x,a3.y,a3.z);

			// THETA6; t6 = LRACOS (-b3 * rstackdir);
			Ftmp += (b3 + rstackdir * cosf(t6)) * (cxst_energy * f4t6Dsin / (f4t6 * rstackmod));

			//if(grooving)  printf("Before: Ttmp = (%f,%f,%f)\n",Ttmp.x,Ttmp.y,Ttmp.z);

			Ttmp -= _cross(ppos_stack, Ftmp);
			//if(grooving)  printf("AADEVICE: Ttmp = (%f,%f,%f)\n",Ttmp.x,Ttmp.y,Ttmp.z);

			c_number4 rba1 = _cross(rbackbonedir, a1);
			c_number4 a1rs = _cross(a1, rstackdir);
			c_number rb_dot_a1rs = CUDA_DOT(rbackbonedir, a1rs);
			c_number rs_dot_rba1 = CUDA_DOT(rstackdir, rba1);

			c_number force_c = f2 * f4t1 * f4t4 * f4t5 * f4t6 * f5cosphi4 * f5Dcosphi3;

			c_number4 forcestack = -(rba1 - rstackdir * rs_dot_rba1) * (force_c / rstackmod);
			c_number4 forceback = -(a1rs - rbackbonedir * rb_dot_a1rs) * (force_c / rbackmod);

			c_number4 myforce = forcestack;
			myforce += forceback;

			// for p
			c_number4 mytorque1 = -_cross(ppos_stack, forcestack); //a1 * (p->int_centers[RNANucleotide::STACK] * rbackbonedir) - rbackbonedir * (p->int_centers[RNANucleotide::STACK] * a1 )
			mytorque1 -= _cross(ppos_back, forceback);  // p->int_centers[RNANucleotide::BACK].cross(forceback);

			// for q
			c_number4 mytorque2 = _cross(qpos_stack, forcestack); //  q->int_centers[RNANucleotide::STACK].cross( forcestack ); //a1 * (p->int_centers[RNANucleotide::STACK] * rbackbonedir) - rbackbonedir * (p->int_centers[RNANucleotide::STACK] * a1 )
			mytorque2 += _cross(qpos_back, forceback); // q->int_centers[RNANucleotide::BACK].cross(forceback);

			mytorque1 -= force_c * _cross(a1, _cross(rstackdir, rbackbonedir)); // a1.cross(rstackdir.cross(rbackbonedir));  // rstackdir * (rbackbonedir * a1) - rbackbonedir * (rstackdir * a1);

			force_c = f2 * f4t1 * f4t4 * f4t5 * f4t6 * f5cosphi3 * f5Dcosphi4;
			rba1 = _cross(rbackbonedir, b1); // rbackbonedir.cross(b1);
			a1rs = _cross(b1, rstackdir);   //    b1.cross(rstackdir);
			rb_dot_a1rs = CUDA_DOT(rbackbonedir, a1rs);
			rs_dot_rba1 = CUDA_DOT(rstackdir, rba1);

			forcestack = -(rba1 - rstackdir * rs_dot_rba1) * (force_c / rstackmod);
			forceback = -(a1rs - rbackbonedir * rb_dot_a1rs) * (force_c / rbackmod);

			myforce += forcestack;
			myforce += forceback;

			// for p
			mytorque2 += _cross(qpos_stack, forcestack); // q->int_centers[RNANucleotide::STACK].cross( forcestack ); //a1 * (p->int_centers[RNANucleotide::STACK] * rbackbonedir) - rbackbonedir * (p->int_centers[RNANucleotide::STACK] * a1 )
			mytorque2 += _cross(qpos_back, forceback); //  q->int_centers[RNANucleotide::BACK].cross(forceback);

			// for q
			mytorque1 -= _cross(ppos_stack, forcestack); // p->int_centers[RNANucleotide::STACK].cross( forcestack ); //a1 * (p->int_centers[RNANucleotide::STACK] * rbackbonedir) - rbackbonedir * (p->int_centers[RNANucleotide::STACK] * a1 )
			mytorque1 -= _cross(ppos_back, forceback); //p->int_centers[RNANucleotide::BACK].cross(forceback);

			mytorque2 -= force_c * _cross(b1, rbackbonedir); // b1.cross(rstackdir.cross(rbackbonedir));  // rstackdir * (rbackbonedir * a1) - rbackbonedir * (rstackdir * a1);

			Ftmp -= myforce;
			Ttmp += mytorque1;

			//if(grooving)  printf("XXX DEVICE: Ftmp=( %f %f %f ), Ttmp=( %f %f %f ), Energy =  %f,  f2 = %f, f4t1=%f, f4t4=%f, f4t5=%f, f4t6=%f, f5cosphi3=%f, f5cosphi4=%f \n",Ftmp.x,Ftmp.y,Ftmp.z, Ttmp.x, Ttmp.y, Ttmp.z,cxst_energy,f2 , f4t1 , f4t4 ,f4t5 , f4t6 ,f5cosphi3 , f5cosphi4);

			//if(grooving)  printf("XXX DEVICE: Ftmp=( %f %f %f ), Ttmp=( %f %f %f )\n",Ftmp.x,Ftmp.y,Ftmp.z, Ttmp.x, Ttmp.y, Ttmp.z);

			//if(grooving && (f5Dcosphi4 != 0 || f5Dcosphi3 != 0))  printf("DEVICE: f5Dcosphi4: %f f5Dcosphi3: %f  \n",f5Dcosphi4,f5Dcosphi3);

			//if(grooving)  printf("Energy =  %f,  f2 = %f, f4t1=%f, f4t4=%f, f4t5=%f, f4t6=%f, f5cosphi3=%f, f5cosphi4=%f \n",Ftmp.x,Ftmp.y,Ftmp.z, Ttmp.x, Ttmp.y, Ttmp.z,cxst_energy,f2 , f4t1 , f4t4 ,f4t5 , f4t6 ,f5cosphi3 , f5cosphi4);

//
//
//			// COSPHI3
//			c_number rbackrefmodcub = rbackrefmod * rbackrefmod * rbackrefmod;
//
//			//c_number a1b1 = a1 * b1;
//			c_number a2b1 = CUDA_DOT(a2, b1);
//			c_number a3b1 = CUDA_DOT(a3, b1);
//			c_number ra1 = CUDA_DOT(rstackdir, a1);
//			c_number ra2 = CUDA_DOT(rstackdir, a2);
//			c_number ra3 = CUDA_DOT(rstackdir, a3);
//			c_number rb1 = CUDA_DOT(rstackdir, b1);
//
//			c_number parentesi = (ra3 * a2b1 - ra2 * a3b1);
//			c_number dcdr    = -rnamodel.RNA_GAMMA * parentesi * (rnamodel.RNA_GAMMA * (ra1 - rb1) + rstackmod) / rbackrefmodcub;
//			c_number dcda1b1 =  rnamodel.RNA_GAMMA * SQR(rnamodel.RNA_GAMMA) * parentesi / rbackrefmodcub;
//			c_number dcda2b1 =  rnamodel.RNA_GAMMA * ra3 / rbackrefmod;
//			c_number dcda3b1 = -rnamodel.RNA_GAMMA * ra2 / rbackrefmod;
//			c_number dcdra1  = -SQR(rnamodel.RNA_GAMMA) * parentesi * rstackmod / rbackrefmodcub;
//			c_number dcdra2  = -rnamodel.RNA_GAMMA * a3b1 / rbackrefmod;
//			c_number dcdra3  =  rnamodel.RNA_GAMMA * a2b1 / rbackrefmod;
//			c_number dcdrb1  = -dcdra1;
//
//			part = cxst_energy * 2 * f5cosphi3D / f5cosphi3;
//
//			Ftmp -= part * (rstackdir * dcdr +
//						    ((a1 - rstackdir * ra1) * dcdra1 +
//							(a2 - rstackdir * ra2) * dcdra2 +
//							(a3 - rstackdir * ra3) * dcdra3 +
//							(b1 - rstackdir * rb1) * dcdrb1) / rstackmod);
//
//			Ttmp += part * (_cross(rstackdir, a1) * dcdra1 +
//						    _cross(rstackdir, a2) * dcdra2 +
//						    _cross(rstackdir ,a3) * dcdra3);
//
//			Ttmp -= part * (_cross(a1, b1) * dcda1b1 +
//						    _cross(a2, b1) * dcda2b1 +
//						    _cross(a3, b1) * dcda3b1);
//

			Ftmp.w = cxst_energy;
			F -= Ftmp;

		}
	}

	//DEBYE-HUCKEL
	if(use_debye_huckel) {
		c_number rbackmod = _module(rbackbone);
		if(rbackmod < MD_dh_RC[0]) {
			c_number4 rbackdir = rbackbone / rbackmod;
			if(rbackmod < MD_dh_RHIGH[0]) {
				Ftmp = rbackdir * (-MD_dh_prefactor[0] * expf(MD_dh_minus_kappa[0] * rbackmod) * (MD_dh_minus_kappa[0] / rbackmod - 1.0f / SQR(rbackmod)));
				Ftmp.w = expf(rbackmod * MD_dh_minus_kappa[0]) * (MD_dh_prefactor[0] / rbackmod);
			}
			else {
				Ftmp = rbackdir * (-2.0f * MD_dh_B[0] * (rbackmod - MD_dh_RC[0]));
				Ftmp.w = MD_dh_B[0] * SQR(rbackmod - MD_dh_RC[0]);
			}

			// check for half-charge strand ends
			if(MD_dh_half_charged_ends[0] && (pbonds.n3 == P_INVALID || pbonds.n5 == P_INVALID)) {
				Ftmp *= 0.5f;
			}
			if(MD_dh_half_charged_ends[0] && (qbonds.n3 == P_INVALID || qbonds.n5 == P_INVALID)) {
				Ftmp *= 0.5f;
			}

			Ttmp -= _cross(ppos_back, Ftmp);
			F -= Ftmp;
		}
	}

	T += Ttmp;

	//if(grooving)  printf("XXX DEVICE all pairs: Ftmp=( %f %f %f ), Ttmp=( %f %f %f ), Energy =  %f \n",Ftmp.x,Ftmp.y,Ftmp.z, Ttmp.x, Ttmp.y, Ttmp.z,Ftmp.w);

	// this component stores the energy due to hydrogen bonding
	T.w = old_Tw + hb_energy;
}

// forces + second step without lists
__global__ void rna_forces(c_number4 *poss, GPU_quat *orientations, c_number4 *forces, c_number4 *torques, LR_bonds *bonds, bool average, bool use_debye_huckel, bool mismatch_repulsion, bool use_mbf, c_number mbf_xmax, c_number mbf_finf, CUDABox *box) {
	if(IND >= MD_N[0]) return;

	c_number4 F = forces[IND];
	c_number4 T = make_c_number4(0, 0, 0, 0);
	LR_bonds bs = bonds[IND];
	c_number4 ppos = poss[IND];
	c_number4 a1, a2, a3;
	get_vectors_from_quat(orientations[IND], a1, a2, a3); //Returns vectors a1,a2 and a3 as they would be in the GPU matrix. These are necessary even in pure quaternion dynamics

	if(bs.n3 != P_INVALID) {
		c_number4 qpos = poss[bs.n3];
		c_number4 b1, b2, b3;
		get_vectors_from_quat(orientations[bs.n3], b1, b2, b3);
		_bonded_part<true>(ppos, a1, a2, a3, qpos, b1, b2, b3, F, T, average, use_mbf, mbf_xmax, mbf_finf);

	}

	if(bs.n5 != P_INVALID) {
		c_number4 qpos = poss[bs.n5];

		c_number4 b1, b2, b3;
		get_vectors_from_quat(orientations[bs.n5], b1, b2, b3);
		_bonded_part<false>(qpos, b1, b2, b3, ppos, a1, a2, a3, F, T, average, use_mbf, mbf_xmax, mbf_finf);

	}

	const int type = get_particle_type(ppos);
	T.w = (c_number) 0;
	for(int j = 0; j < MD_N[0]; j++) {
		if(j != IND && bs.n3 != j && bs.n5 != j) {
			const c_number4 qpos = poss[j];
			c_number4 b1, b2, b3;
			get_vectors_from_quat(orientations[j], b1, b2, b3);
			LR_bonds qbonds = bonds[j];

			_particle_particle_interaction(ppos, a1, a2, a3, qpos, b1, b2, b3, F, T, average, use_debye_huckel, mismatch_repulsion, bs, qbonds, box);
		}
	}

	T = _vectors_transpose_c_number4_product(a1, a2, a3, T);

	forces[IND] = F;
	torques[IND] = T;
}

__global__ void rna_forces_edge_nonbonded(c_number4 *poss, GPU_quat *orientations, c_number4 *forces, c_number4 *torques, edge_bond *edge_list, int n_edges, LR_bonds *bonds, bool average, bool use_debye_huckel, bool mismatch_repulsion, CUDABox *box) {
	if(IND >= n_edges) return;

	c_number4 dF = make_c_number4(0, 0, 0, 0);
	c_number4 dT = make_c_number4(0, 0, 0, 0);

	edge_bond b = edge_list[IND];

	// get info for particle 1
	c_number4 ppos = poss[b.from];
	c_number4 a1, a2, a3;
	get_vectors_from_quat(orientations[b.from], a1, a2, a3);

	// get info for particle 2
	c_number4 qpos = poss[b.to];
	c_number4 b1, b2, b3;
	get_vectors_from_quat(orientations[b.to], b1, b2, b3);

	LR_bonds pbonds = bonds[b.from];
	LR_bonds qbonds = bonds[b.to];
	_particle_particle_interaction(ppos, a1, a2, a3, qpos, b1, b2, b3, dF, dT, average, use_debye_huckel, mismatch_repulsion, pbonds, qbonds, box);

	int from_index = MD_N[0] * (IND % MD_n_forces[0]) + b.from;
	//int from_index = MD_N[0]*(b.n_from % MD_n_forces[0]) + b.from;
	if((dF.x * dF.x + dF.y * dF.y + dF.z * dF.z + dF.w * dF.w) > (c_number) 0.f) LR_atomicAddXYZ(&(forces[from_index]), dF);
	if((dT.x * dT.x + dT.y * dT.y + dT.z * dT.z + dT.w * dT.w) > (c_number) 0.f) LR_atomicAddXYZ(&(torques[from_index]), dT);

	// Allen Eq. 6 pag 3:
	c_number4 dr = box->minimum_image(ppos, qpos); // returns qpos-ppos
	c_number4 crx = _cross(dr, dF);
	dT.x = -dT.x + crx.x;
	dT.y = -dT.y + crx.y;
	dT.z = -dT.z + crx.z;

	dF.x = -dF.x;
	dF.y = -dF.y;
	dF.z = -dF.z;

	int to_index = MD_N[0] * (IND % MD_n_forces[0]) + b.to;
	//int to_index = MD_N[0]*(b.n_to % MD_n_forces[0]) + b.to;
	if((dF.x * dF.x + dF.y * dF.y + dF.z * dF.z + dF.w * dF.w) > (c_number) 0.f) LR_atomicAddXYZ(&(forces[to_index]), dF);
	if((dT.x * dT.x + dT.y * dT.y + dT.z * dT.z + dT.w * dT.w) > (c_number) 0.f) LR_atomicAddXYZ(&(torques[to_index]), dT);
}

// bonded interactions for edge-based approach
__global__ void rna_forces_edge_bonded(c_number4 *poss, GPU_quat *orientations, c_number4 *forces, c_number4 *torques, LR_bonds *bonds, bool average, bool use_mbf, c_number mbf_xmax, c_number mbf_finf) {
	if(IND >= MD_N[0]) return;

	c_number4 F0, T0;

	F0.x = forces[IND].x;
	F0.y = forces[IND].y;
	F0.z = forces[IND].z;
	F0.w = forces[IND].w;
	T0.x = torques[IND].x;
	T0.y = torques[IND].y;
	T0.z = torques[IND].z;
	T0.w = torques[IND].w;

	c_number4 dF = make_c_number4(0, 0, 0, 0);
	c_number4 dT = make_c_number4(0, 0, 0, 0);
	c_number4 ppos = poss[IND];
	LR_bonds bs = bonds[IND];
	c_number4 a1, a2, a3;
	get_vectors_from_quat(orientations[IND], a1, a2, a3);

	if(bs.n3 != P_INVALID) {
		c_number4 qpos = poss[bs.n3];
		c_number4 b1, b2, b3;
		get_vectors_from_quat(orientations[bs.n3], b1, b2, b3);

		//grooving =false;
		_bonded_part<true>(ppos, a1, a2, a3, qpos, b1, b2, b3, dF, dT, average, use_mbf, mbf_xmax, mbf_finf);
	}
	if(bs.n5 != P_INVALID) {
		c_number4 qpos = poss[bs.n5];
		c_number4 b1, b2, b3;
		get_vectors_from_quat(orientations[bs.n5], b1, b2, b3);

		_bonded_part<false>(qpos, b1, b2, b3, ppos, a1, a2, a3, dF, dT, average, use_mbf, mbf_xmax, mbf_finf);

	}

	forces[IND] = (dF + F0);
	torques[IND] = (dT + T0);

	torques[IND] = _vectors_transpose_c_number4_product(a1, a2, a3, torques[IND]);
}

// forces + second step with verlet lists
__global__ void rna_forces(c_number4 *poss, GPU_quat *orientations, c_number4 *forces, c_number4 *torques, int *matrix_neighs, int *c_number_neighs, LR_bonds *bonds, bool average, bool use_debye_huckel, bool mismatch_repulsion, bool use_mbf, c_number mbf_xmax, c_number mbf_finf, CUDABox *box) {
	if(IND >= MD_N[0]) return;

	//c_number4 F = make_c_number4(0, 0, 0, 0);
	c_number4 F = forces[IND];
	c_number4 T = make_c_number4(0, 0, 0, 0);
	c_number4 ppos = poss[IND];
	LR_bonds bs = bonds[IND];
	// particle axes according to Allen's paper
	c_number4 a1, a2, a3;
	get_vectors_from_quat(orientations[IND], a1, a2, a3);

	if(bs.n3 != P_INVALID) {
		c_number4 qpos = poss[bs.n3];
		c_number4 b1, b2, b3;
		get_vectors_from_quat(orientations[bs.n3], b1, b2, b3);

		_bonded_part<true>(ppos, a1, a2, a3, qpos, b1, b2, b3, F, T, average, use_mbf, mbf_xmax, mbf_finf);

	}
	if(bs.n5 != P_INVALID) {
		c_number4 qpos = poss[bs.n5];
		c_number4 b1, b2, b3;
		get_vectors_from_quat(orientations[bs.n5], b1, b2, b3);

		_bonded_part<false>(qpos, b1, b2, b3, ppos, a1, a2, a3, F, T, average, use_mbf, mbf_xmax, mbf_finf);
	}

	const int type = get_particle_type(ppos);
	const int num_neighs = c_number_neighs[IND];

	T.w = (c_number) 0;
	for(int j = 0; j < num_neighs; j++) {
		const int k_index = matrix_neighs[j * MD_N[0] + IND];

		const c_number4 qpos = poss[k_index];
		c_number4 b1, b2, b3;
		get_vectors_from_quat(orientations[k_index], b1, b2, b3);
		LR_bonds pbonds = bonds[IND];
		LR_bonds qbonds = bonds[k_index];
		_particle_particle_interaction(ppos, a1, a2, a3, qpos, b1, b2, b3, F, T, average, use_debye_huckel, mismatch_repulsion, pbonds, qbonds, box);
	}

	T = _vectors_transpose_c_number4_product(a1, a2, a3, T);

	forces[IND] = F;
	torques[IND] = T;
}

//FFS order parameter pre-calculations

// check whether a particular pair of particles have hydrogen bonding energy lower than a given threshold hb_threshold (which may vary)
__global__ void rna_hb_op_precalc(c_number4 *poss, GPU_quat *orientations, int *op_pairs1, int *op_pairs2, float *hb_energies, int n_threads, bool *region_is_nearhb, CUDABox *box) {
	if(IND >= n_threads) return;

	int pind = op_pairs1[IND];
	int qind = op_pairs2[IND];
	// get distance between this nucleotide pair's "com"s
	c_number4 ppos = poss[pind];
	c_number4 qpos = poss[qind];
	c_number4 r = box->minimum_image(ppos, qpos);

	// check whether hb energy is below a certain threshold for this nucleotide pair
	int ptype = get_particle_type(ppos);
	int qtype = get_particle_type(qpos);
	int pbtype = get_particle_btype(ppos);
	int qbtype = get_particle_btype(qpos);
	int int_type = pbtype + qbtype;

	GPU_quat po = orientations[pind];
	GPU_quat qo = orientations[qind];

	//This gets an extra two vectors that are not needed, but the function doesn't seem to be called at all, so should make no difference. 
	c_number4 a1, a2, a3, b1, b2, b3;
	get_vectors_from_quat(po, a1, a2, a3);

	get_vectors_from_quat(qo, b1, b2, b3);

	c_number4 ppos_base = rnamodel.RNA_POS_BASE * a1;
	c_number4 qpos_base = rnamodel.RNA_POS_BASE * b1;

	// HYDROGEN BONDING

	/*
	 if(!average) //allow for wobble bp
	 {
	 if( int_type == 4 && ( (qtype == N_T && ptype == N_G ) || (qtype == N_G && ptype == N_T)  ))
	 is_pair = 1;
	 }*/

	int is_pair = int(int_type == 3);
	c_number hb_energy = (c_number) 0;
	c_number4 rhydro = r + qpos_base - ppos_base;
	c_number rhydromodsqr = CUDA_DOT(rhydro, rhydro);
	if(is_pair && SQR(rnamodel.RNA_HYDR_RCLOW) < rhydromodsqr && rhydromodsqr < SQR(rnamodel.RNA_HYDR_RCHIGH)) {
		c_number hb_multi = (abs(qbtype) >= 300 && abs(pbtype) >= 300) ? MD_hb_multi[0] : 1.f;
		// versor and magnitude of the base-base separation
		c_number rhydromod = sqrtf(rhydromodsqr);
		c_number4 rhydrodir = rhydro / rhydromod;

		// angles involved in the HB interaction
		c_number t1 = CUDA_LRACOS(-CUDA_DOT(a1, b1));
		c_number t2 = CUDA_LRACOS(-CUDA_DOT(b1, rhydrodir));
		c_number t3 = CUDA_LRACOS(CUDA_DOT(a1, rhydrodir));
		c_number t4 = CUDA_LRACOS(CUDA_DOT(a3, b3));
		c_number t7 = CUDA_LRACOS(-CUDA_DOT(rhydrodir, b3));
		c_number t8 = CUDA_LRACOS(CUDA_DOT(rhydrodir, a3));

		// functions called at their relevant arguments
		c_number f1 = hb_multi * _f1(rhydromod, RNA_HYDR_F1, ptype, qtype);
		c_number f4t1 = _f4(t1, rnamodel.RNA_HYDR_THETA1_T0, rnamodel.RNA_HYDR_THETA1_TS, rnamodel.RNA_HYDR_THETA1_TC, rnamodel.RNA_HYDR_THETA1_A, rnamodel.RNA_HYDR_THETA1_B);
		c_number f4t2 = _f4(t2, rnamodel.RNA_HYDR_THETA2_T0, rnamodel.RNA_HYDR_THETA2_TS, rnamodel.RNA_HYDR_THETA2_TC, rnamodel.RNA_HYDR_THETA2_A, rnamodel.RNA_HYDR_THETA2_B);
		c_number f4t3 = _f4(t3, rnamodel.RNA_HYDR_THETA3_T0, rnamodel.RNA_HYDR_THETA3_TS, rnamodel.RNA_HYDR_THETA3_TC, rnamodel.RNA_HYDR_THETA3_A, rnamodel.RNA_HYDR_THETA3_B);
		c_number f4t4 = _f4(t4, rnamodel.RNA_HYDR_THETA4_T0, rnamodel.RNA_HYDR_THETA4_TS, rnamodel.RNA_HYDR_THETA4_TC, rnamodel.RNA_HYDR_THETA4_A, rnamodel.RNA_HYDR_THETA4_B);
		c_number f4t7 = _f4(t7, rnamodel.RNA_HYDR_THETA7_T0, rnamodel.RNA_HYDR_THETA7_TS, rnamodel.RNA_HYDR_THETA7_TC, rnamodel.RNA_HYDR_THETA7_A, rnamodel.RNA_HYDR_THETA7_B);
		c_number f4t8 = _f4(t8, rnamodel.RNA_HYDR_THETA8_T0, rnamodel.RNA_HYDR_THETA8_TS, rnamodel.RNA_HYDR_THETA8_TC, rnamodel.RNA_HYDR_THETA8_A, rnamodel.RNA_HYDR_THETA8_B);

		hb_energy = f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;
	}
	//END HYDROGEN BONDING

	hb_energies[IND] = hb_energy;
}

__global__ void rna_near_hb_op_precalc(c_number4 *poss, GPU_quat *orientations, int *op_pairs1, int *op_pairs2, bool *nearly_bonded_array, int n_threads, bool *region_is_nearhb, CUDABox *box) {
	if(IND >= n_threads) return;

	int pind = op_pairs1[IND];
	int qind = op_pairs2[IND];
	// get distance between this nucleotide pair's "com"s
	c_number4 ppos = poss[pind];
	c_number4 qpos = poss[qind];
	c_number4 r = box->minimum_image(ppos, qpos);

	// check whether hb energy is below a certain threshold for this nucleotide pair
	int ptype = get_particle_type(ppos);
	int qtype = get_particle_type(qpos);
	int pbtype = get_particle_btype(ppos);
	int qbtype = get_particle_btype(qpos);
	int int_type = pbtype + qbtype;

	GPU_quat po = orientations[pind];
	GPU_quat qo = orientations[qind];

	//This gets extra a2 and b2 vectors that aren't needed. get_vectors_from_quat could easily be modified to only return the relevant vectors, but it will make computationally very little difference, since most of the same c_numbers need to be calculated anyway. Perhaps worth changing for memory considerations, however
	c_number4 a1, a2, a3, b1, b2, b3;
	get_vectors_from_quat(po, a1, a2, a3);
	get_vectors_from_quat(qo, b1, b2, b3);

	c_number4 ppos_base = rnamodel.RNA_POS_BASE * a1;
	c_number4 qpos_base = rnamodel.RNA_POS_BASE * b1;

	// HYDROGEN BONDING

	int is_pair = int(int_type == 3);
	/*if(!average) //allow for wobble bp
	 {
	 if( int_type == 4 && ( (qtype == N_T && ptype == N_G ) || (qtype == N_G && ptype == N_T)  ))
	 is_pair = 1;
	 }*/

	c_number4 rhydro = r + qpos_base - ppos_base;
	c_number rhydromodsqr = CUDA_DOT(rhydro, rhydro);

	int total_nonzero;
	bool nearly_bonded = false;

	if(is_pair && SQR(rnamodel.RNA_HYDR_RCLOW) < rhydromodsqr && rhydromodsqr < SQR(rnamodel.RNA_HYDR_RCHIGH)) {
		c_number hb_multi = (abs(qbtype) >= 300 && abs(pbtype) >= 300) ? MD_hb_multi[0] : 1.f;
		// versor and magnitude of the base-base separation
		c_number rhydromod = sqrtf(rhydromodsqr);
		c_number4 rhydrodir = rhydro / rhydromod;

		// angles involved in the HB interaction
		c_number t1 = CUDA_LRACOS(-CUDA_DOT(a1, b1));
		c_number t2 = CUDA_LRACOS(-CUDA_DOT(b1, rhydrodir));
		c_number t3 = CUDA_LRACOS(CUDA_DOT(a1, rhydrodir));
		c_number t4 = CUDA_LRACOS(CUDA_DOT(a3, b3));
		c_number t7 = CUDA_LRACOS(-CUDA_DOT(rhydrodir, b3));
		c_number t8 = CUDA_LRACOS(CUDA_DOT(rhydrodir, a3));

		// functions called at their relevant arguments
		c_number f1 = hb_multi * _f1(rhydromod, RNA_HYDR_F1, ptype, qtype);
		c_number f4t1 = _f4(t1, rnamodel.RNA_HYDR_THETA1_T0, rnamodel.RNA_HYDR_THETA1_TS, rnamodel.RNA_HYDR_THETA1_TC, rnamodel.RNA_HYDR_THETA1_A, rnamodel.RNA_HYDR_THETA1_B);
		c_number f4t2 = _f4(t2, rnamodel.RNA_HYDR_THETA2_T0, rnamodel.RNA_HYDR_THETA2_TS, rnamodel.RNA_HYDR_THETA2_TC, rnamodel.RNA_HYDR_THETA2_A, rnamodel.RNA_HYDR_THETA2_B);
		c_number f4t3 = _f4(t3, rnamodel.RNA_HYDR_THETA3_T0, rnamodel.RNA_HYDR_THETA3_TS, rnamodel.RNA_HYDR_THETA3_TC, rnamodel.RNA_HYDR_THETA3_A, rnamodel.RNA_HYDR_THETA3_B);
		c_number f4t4 = _f4(t4, rnamodel.RNA_HYDR_THETA4_T0, rnamodel.RNA_HYDR_THETA4_TS, rnamodel.RNA_HYDR_THETA4_TC, rnamodel.RNA_HYDR_THETA4_A, rnamodel.RNA_HYDR_THETA4_B);
		c_number f4t7 = _f4(t7, rnamodel.RNA_HYDR_THETA7_T0, rnamodel.RNA_HYDR_THETA7_TS, rnamodel.RNA_HYDR_THETA7_TC, rnamodel.RNA_HYDR_THETA7_A, rnamodel.RNA_HYDR_THETA7_B);
		c_number f4t8 = _f4(t8, rnamodel.RNA_HYDR_THETA8_T0, rnamodel.RNA_HYDR_THETA8_TS, rnamodel.RNA_HYDR_THETA8_TC, rnamodel.RNA_HYDR_THETA8_A, rnamodel.RNA_HYDR_THETA8_B);

		// the nearly bonded order parameter requires either all or all but one of these factors to be non-zero
		bool f1not0 = f1 < 0 ? 1 : 0;
		bool f4t1not0 = f4t1 > 0 ? 1 : 0;
		bool f4t2not0 = f4t2 > 0 ? 1 : 0;
		bool f4t3not0 = f4t3 > 0 ? 1 : 0;
		bool f4t4not0 = f4t4 > 0 ? 1 : 0;
		bool f4t7not0 = f4t7 > 0 ? 1 : 0;
		bool f4t8not0 = f4t8 > 0 ? 1 : 0;
		total_nonzero = f1not0 + f4t1not0 + f4t2not0 + f4t3not0 + f4t4not0 + f4t7not0 + f4t8not0;
		if(total_nonzero >= 6) nearly_bonded = true;
	}

	// END HYDROGEN BONDING

	nearly_bonded_array[IND] = nearly_bonded;
}

// compute the distance between a pair of particles
__global__ void rna_dist_op_precalc(c_number4 *poss, GPU_quat *orientations, int *op_pairs1, int *op_pairs2, c_number *op_dists, int n_threads, CUDABox *box) {
	if(IND >= n_threads) return;

	int pind = op_pairs1[IND];
	int qind = op_pairs2[IND];

	// get distance between this nucleotide pair's "com"s
	c_number4 ppos = poss[pind];
	c_number4 qpos = poss[qind];
	c_number4 r = box->minimum_image(ppos, qpos);

	GPU_quat po = orientations[pind];
	GPU_quat qo = orientations[qind];

	//This gets extra a2 and b2 vectors that aren't needed. get_vectors_from_quat could easily be modified to only return the relevant vectors, but it will make computationally very little difference, since most of the same c_numbers need to be calculated anyway. Perhaps worth changing for memory considerations, however 
	c_number4 a1, a2, a3, b1, b2, b3;
	get_vectors_from_quat(po, a1, a2, a3);
	get_vectors_from_quat(qo, b1, b2, b3);

	c_number4 ppos_base = rnamodel.RNA_POS_BASE * a1;
	c_number4 qpos_base = rnamodel.RNA_POS_BASE * b1;

	c_number4 rbase = r + qpos_base - ppos_base;
	op_dists[IND] = _module(rbase);
}
