#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <emscripten/emscripten.h>

#define THRESHOLD 64  // Ngưỡng chuyển sang nhân thường

// ======================= HÀM TẠO MA TRẬN NGẪU NHIÊN ==========================
void taoMaTranNgauNhien(int n, int *A) {
	for (int i = 0; i < n * n; i++)
		A[i] = rand() % 10;
}

// ======================= NHÂN MA TRẬN THƯỜNG ==============================
void nhanMaTranThuong(int n, int *A, int *B, int *C) {
	for (int i = 0; i<n; i++)
		for (int j = 0; j<n; j++){
			int sum = 0;
			for (int k = 0; k<n; k++)
				sum += A[i*n + k] * B[k*n + j];
			C[i*n + j] = sum;
		}
}

// ======================= CỘNG / TRỪ MA TRẬN ================================
void congMaTran(int n, int *A, int *B, int *C){
	for (int i = 0; i<n*n; i++)
		C[i] = A[i] + B[i];
}

void truMaTran(int n, int *A, int *B, int *C){
	for (int i = 0; i<n*n; i++)
		C[i] = A[i] - B[i];
}

// ======================= STRASSEN TỐI ƯU ==================================
void strassen_opt(int n, int *A, int *B, int *C, int *temp) {
	if (n <= THRESHOLD){
		nhanMaTranThuong(n, A, B, C);
		return;
	}

	int k = n / 2;
	int size = k*k;

	int *A11 = A;
	int *A12 = A + k;
	int *A21 = A + k*n;
	int *A22 = A + k*n + k;

	int *B11 = B;
	int *B12 = B + k;
	int *B21 = B + k*n;
	int *B22 = B + k*n + k;

	int *C11 = C;
	int *C12 = C + k;
	int *C21 = C + k*n;
	int *C22 = C + k*n + k;

	int *M1 = temp;
	int *M2 = temp + size;
	int *M3 = temp + 2 * size;
	int *M4 = temp + 3 * size;
	int *M5 = temp + 4 * size;
	int *M6 = temp + 5 * size;
	int *M7 = temp + 6 * size;
	int *T1 = temp + 7 * size;
	int *T2 = temp + 8 * size;

	// M1 = (A11+A22)*(B11+B22)
	congMaTran(k, A11, A22, T1);
	congMaTran(k, B11, B22, T2);
	strassen_opt(k, T1, T2, M1, temp + 9 * size);

	// M2 = (A21+A22)*B11
	congMaTran(k, A21, A22, T1);
	strassen_opt(k, T1, B11, M2, temp + 9 * size);

	// M3 = A11*(B12-B22)
	truMaTran(k, B12, B22, T2);
	strassen_opt(k, A11, T2, M3, temp + 9 * size);

	// M4 = A22*(B21-B11)
	truMaTran(k, B21, B11, T2);
	strassen_opt(k, A22, T2, M4, temp + 9 * size);

	// M5 = (A11+A12)*B22
	congMaTran(k, A11, A12, T1);
	strassen_opt(k, T1, B22, M5, temp + 9 * size);

	// M6 = (A21-A11)*(B11+B12)
	truMaTran(k, A21, A11, T1);
	congMaTran(k, B11, B12, T2);
	strassen_opt(k, T1, T2, M6, temp + 9 * size);

	// M7 = (A12-A22)*(B21+B22)
	truMaTran(k, A12, A22, T1);
	congMaTran(k, B21, B22, T2);
	strassen_opt(k, T1, T2, M7, temp + 9 * size);

	// Tính C11..C22
	for (int i = 0; i<k; i++)
		for (int j = 0; j<k; j++){
			C11[i*n + j] = M1[i*k + j] + M4[i*k + j] - M5[i*k + j] + M7[i*k + j];
			C12[i*n + j] = M3[i*k + j] + M5[i*k + j];
			C21[i*n + j] = M2[i*k + j] + M4[i*k + j];
			C22[i*n + j] = M1[i*k + j] - M2[i*k + j] + M3[i*k + j] + M6[i*k + j];
		}
}


// ======================= STRASSEN WINOGRAD TỐI ƯU ========================
void strassenWinograd_opt(int n, int *A, int *B, int *C, int *temp){
		if (n <= THRESHOLD){
			nhanMaTranThuong(n, A, B, C);
			return;
		}

		int k = n / 2;
		int size = k*k;

		int *A11 = A;
		int *A12 = A + k;
		int *A21 = A + k*n;
		int *A22 = A + k*n + k;

		int *B11 = B;
		int *B12 = B + k;
		int *B21 = B + k*n;
		int *B22 = B + k*n + k;

		int *C11 = C;
		int *C12 = C + k;
		int *C21 = C + k*n;
		int *C22 = C + k*n + k;

		int *S1 = temp;          // A21 + A22
		int *S2 = temp + size;   // S1 - A11
		int *S3 = temp + 2 * size; // B12 - B11
		int *S4 = temp + 3 * size; // B22 - S3

		int *M1 = temp + 4 * size;
		int *M2 = temp + 5 * size;
		int *M3 = temp + 6 * size;
		int *M4 = temp + 7 * size;
		int *M5 = temp + 8 * size;
		int *M6 = temp + 9 * size;
		int *M7 = temp + 10 * size;

		int *T1 = temp + 11 * size;
		int *T2 = temp + 12 * size;

		// S1 = A21 + A22
		congMaTran(k, A21, A22, S1);

		// S2 = S1 - A11
		truMaTran(k, S1, A11, S2);

		// S3 = B12 - B11
		truMaTran(k, B12, B11, S3);

		// S4 = B22 - S3
		truMaTran(k, B22, S3, S4);

		// M1 = S2 * S4
		strassen_opt(k, S2, S4, M1, temp + 13 * size);

		// M2 = A11 * B11
		strassen_opt(k, A11, B11, M2, temp + 13 * size);

		// M3 = A12 * B21
		strassen_opt(k, A12, B21, M3, temp + 13 * size);

		// M4 = (A11 - A21) * (B22 - B12)
		truMaTran(k, A11, A21, T1);
		truMaTran(k, B22, B12, T2);
		strassen_opt(k, T1, T2, M4, temp + 14 * size);

		// M5 = S1 * S3
		strassen_opt(k, S1, S3, M5, temp + 14 * size);

		// M6 = (A12 - S2) * B22
		truMaTran(k, A12, S2, T1);
		strassen_opt(k, T1, B22, M6, temp + 14 * size);

		// M7 = A22 * (S4 - B21)
		truMaTran(k, S4, B21, T1);
		strassen_opt(k, A22, T1, M7, temp + 14 * size);

		// T1 = M1 + M2
		congMaTran(k, M1, M2, T1);

		// T2 = T1 + M4
		congMaTran(k, T1, M4, T2);

		// C11 = M2 + M3
		congMaTran(k, M2, M3, C11);

		// C12 = T1 + M5
		congMaTran(k, T2, M5, C12);

		// C21 = T2 - M7
		truMaTran(k, T2, M7, C21);

		// C22 = T2 + M5
		congMaTran(k, T2, M5, C22);
	}
// ======================= MAIN =============================================
EMSCRIPTEN_KEEPALIVE
void runSimulation(int k){
	int n = 1 << k;

	int *A = (int*)malloc(n*n*sizeof(int));
	int *B = (int*)malloc(n*n*sizeof(int));
	int *C1 = (int*)malloc(n*n*sizeof(int));
	int *C2 = (int*)malloc(n*n*sizeof(int));
	int *C3 = (int*)malloc(n*n*sizeof(int));
	int *buffer = (int*)malloc(n*n * 16 * sizeof(int));

	srand(time(NULL));
	taoMaTranNgauNhien(n, A);
	taoMaTranNgauNhien(n, B);

	clock_t start, end;
	double t1, t2, t3;

	start = clock(); nhanMaTranThuong(n, A, B, C1); end = clock(); t1 = (double)(end - start) / CLOCKS_PER_SEC;
	start = clock(); strassen_opt(n, A, B, C2, buffer); end = clock(); t2 = (double)(end - start) / CLOCKS_PER_SEC;
	start = clock(); strassenWinograd_opt(n, A, B, C3, buffer); end = clock(); t3 = (double)(end - start) / CLOCKS_PER_SEC;

	printf("===== KET QUA THOI GIAN TINH =====\n");
	printf("Nhan thuong       : %.6f giay\n", t1);
	printf("Strassen          : %.6f giay\n", t2);
	printf("Strassen Winograd : %.6f giay\n", t3);

	free(A); free(B); free(C1); free(C2); free(C3); free(buffer);
}