#ifndef __MKL_APPEND_H__
#define __MKL_APPEND_H__
#include <mkl.h>
#include <iostream>

#ifndef MKL_ALIGN
#define MKL_ALIGN 64
#endif

extern const lapack_complex_double COMPLEX_ONE;
extern const lapack_complex_double COMPLEX_MINUS_ONE;
extern const lapack_complex_double COMPLEX_I;
extern const lapack_complex_double COMPLEX_MINUS_I;
extern const lapack_complex_double COMPLEX_ZERO;

inline const lapack_complex_double conj(const lapack_complex_double a);
inline const lapack_complex_double operator+(const lapack_complex_double a, const lapack_complex_double b);
inline void operator+=(lapack_complex_double& a, const lapack_complex_double b);
inline const lapack_complex_double operator-(const lapack_complex_double a, const lapack_complex_double b);
inline const lapack_complex_double operator*(const lapack_complex_double a, const lapack_complex_double b);

inline std::ostream& operator <<(std::ostream& os, const lapack_complex_double z);

//construct the identity matrix
//in *packed* storage
class Zgem_I {
private:
	lapack_complex_double* p;
	MKL_INT N;
public:
	Zgem_I();
	bool Prepare_I(MKL_INT _N);	//若内存分配失败，返回false
	const lapack_complex_double* const Get_I() const;
	MKL_INT Get_N() const;
	~Zgem_I();
};

//output the matrix(vector) G(n*m) with ostream
void zgem_out(const lapack_complex_double* G, int m, int n, int ld);

//output the Hermitian matrix in packed storage, with ostream
void zhpm_out(const lapack_complex_double* H, int n);

//R(n*m)=conjugate transpose of G(m*n)
void zgem_ct(const lapack_complex_double* G, int m, int n, lapack_complex_double* R);

//Hermitan matrix, packed storage, add
void zhpma(const lapack_complex_double* H1, const lapack_complex_double* H2, int n, lapack_complex_double* Hr);

//DirectProduct for general matrixs(vectors)
void zgemk(const lapack_complex_double* G1, int m, int n, const lapack_complex_double* G2, int p, int q, lapack_complex_double* R);

//DirectProduct for vectors
void zvek(const lapack_complex_double* v1, int m, int inc_v1, const lapack_complex_double* v2, int n, int inc_v2, lapack_complex_double* vr, int inc_vr);
void zvek_append(const lapack_complex_double* v1, int m, int inc_v1, const lapack_complex_double* v2, int n, int inc_v2, lapack_complex_double* vr, int inc_vr);

//Compute Hr(mn*mn):=H1(m*m)(DirectProduct)H2(n*n)
//full storage, L, row major,
//H1 should be stored in lower triangle, row major, Hermitian
void zhelrmk(const lapack_complex_double* H1, int m, int ld1, const lapack_complex_double* H2, int n, int ld2, lapack_complex_double* Hr, int ldr);

//Compute Hr(mn*mn):=H1(m*m)(DirectProduct)H2(n*n)
//full storage, L, row major,
//H1 should be stored in lower triangle, row major, Hermitian
//appending mode, meaning result will be added to Hr
void zhelrmk_append(const lapack_complex_double* H1, int m, int ld1, const lapack_complex_double* H2, int n, int ld2, lapack_complex_double* Hr, int ldr);

//Compute Hr(mn*mn):=H1(m*m)(DirectProduct)H2(n*n)
//*packed* storage, L, row major,
void zhplrmk(const lapack_complex_double* H1, int m, const lapack_complex_double* H2, int n, lapack_complex_double* Hr);

//Compute Hr(mn*mn):=H1(m*m)(DirectProduct)H2(n*n)
//*packed* storage, L, row major,
//appending mode, meaning result will be added to Hr
void zhplrmk_append(const lapack_complex_double* H1, int m, const lapack_complex_double* H2, int n, lapack_complex_double* Hr);


//将该矢量归一化，返回归一化前的模
double zvenml(int n, lapack_complex_double* v, int incv);

#endif