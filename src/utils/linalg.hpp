//////////////////////////////////////////////////////////////////////////
//  linalg.hpp -- column-major dense matrix + thin BLAS/LAPACK wrappers
//
//
// It provides a small C++ matrix class and thin wrappers around BLAS/LAPACK routines for:
// dense matrix storage,
// matrix–matrix multiplication,
// solving dense linear systems,
// LU factorization,
// repeatedly solving systems with the same matrix.
//////////////////////////////////////////////////////////////////////////

#pragma once
#include <stdexcept>
#include <vector>

extern "C" {
void dgetrf_(const int* m, const int* n, double* a, const int* lda, int* ipiv, int* info); // solve using LU
void dgetrs_(const char* tr, const int* n, const int* nrhs, const double* a, const int* lda,// compute LU factorization
             const int* ipiv, double* b, const int* ldb, int* info);
void dgesv_(const int* n, const int* nrhs, double* a, const int* lda, int* ipiv, double* b, //solve Ax=V
            const int* ldb, int* info);
void dgemm_(const char* ta, const char* tb, const int* m, const int* n, const int* k, //matrix multiplication
            const double* alpha, const double* a, const int* lda, const double* b,
            const int* ldb, const double* beta, double* c, const int* ldc);
}

namespace bem {

// Column-major dense matrix: element (i, j) at a[i + j*n]
struct Mat {
    int n = 0, m = 0;
    std::vector<double> a;
    Mat() = default;
    Mat(int n_, int m_) : n(n_), m(m_), a((size_t)n_ * m_, 0.0) {}
    double& operator()(int i, int j) { return a[i + (size_t)j * n]; }
    double operator()(int i, int j) const { return a[i + (size_t)j * n]; }
};


// C = alpha*A*B + beta*C, matrix multiplication
inline void gemm(double alpha, const Mat& A, const Mat& B, double beta, Mat& C) {
    if (A.m != B.n || C.n != A.n || C.m != B.m) throw std::runtime_error("gemm: dimension mismatch");// dimension check: A: n × m, B: m × p, C: n × p
    if (A.n == 0 || B.m == 0 || A.m == 0) return;
    const char N = 'N';
    dgemm_(&N, &N, &A.n, &B.m, &A.m, &alpha, A.a.data(), &A.n, B.a.data(), &B.n, &beta,
           C.a.data(), &C.n);
}



// B <- A^{-1} B   (A is copied, so it is left unchanged)
inline void solve_dense(Mat A, Mat& B) {
    std::vector<int> ip(A.n);
    int info;
    dgesv_(&A.n, &B.m, A.a.data(), &A.n, ip.data(), B.a.data(), &B.n, &info);
    if (info) throw std::runtime_error("dgesv failed");
}


// LU factorisation kept for repeated solves with different right-hand sides
struct LUFactor {
    Mat LU;
    std::vector<int> piv;
    void factor(Mat A) {
        LU = std::move(A);
        piv.resize(LU.n);
        int info;
        dgetrf_(&LU.n, &LU.n, LU.a.data(), &LU.n, piv.data(), &info);
        if (info) throw std::runtime_error("dgetrf failed");
    }
    void solve(Mat& B) const {  // B <- A^{-1} B
        if (B.m == 0) return;
        int info;
        const char N = 'N';
        dgetrs_(&N, &LU.n, &B.m, LU.a.data(), &LU.n, piv.data(), B.a.data(), &B.n, &info);
    }
};

}  // namespace bem
//
//                  mesh.hpp
//                     │
//                     │ creates
//                     ▼
//           triangular boundary mesh
//                     │
//                     │ geometry()
//                     ▼
//        ┌────────────┼─────────────┐
//        │            │             │
//    centroid       normal        area
//        │            │             │
//        └────────────┼─────────────┘
//                     │
//                     ▼
//              BEM discretization
//                     │
//           evaluate influence
//              coefficients
//                     │
//                     ▼
//              dense matrix A
//                     │
//                     ▼
//                linalg.hpp
//                     │
//              ┌──────┴──────┐
//              │             │
//           dgemm          dgesv/
//                          dgetrf+
//                          dgetrs
//              │             │
//              ▼             ▼
//         matrix ops      solve Ax=b