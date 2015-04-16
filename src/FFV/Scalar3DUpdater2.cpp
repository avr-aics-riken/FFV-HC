#include "Scalar3DUpdater2.h"

#ifdef BCMT_NAMESPACE
namespace BCMT_NAMESPACE {
#endif

template <>
int Scalar3DUpdater2<int>::interpolateF2C(const int* fData, const Index3DS& fIndex, int I, int J, int K) {
					int i = 2 * I;
					int j = 2 * J;
					int k = 2 * K;
					if( fData[fIndex(i  ,j  ,k  )] == 1 ||
							fData[fIndex(i+1,j  ,k  )] == 1 ||
							fData[fIndex(i  ,j+1,k  )] == 1 ||
							fData[fIndex(i  ,j  ,k+1)] == 1 ||
							fData[fIndex(i  ,j+1,k+1)] == 1 ||
							fData[fIndex(i+1,j  ,k+1)] == 1 ||
							fData[fIndex(i+1,j+1,k  )] == 1 ||
							fData[fIndex(i+1,j+1,k+1)] == 1 ) {
						return 1;
					} 
					return -1;
				}


template <>
int Scalar3DUpdater2<int>::interpolateC2F(const int* cData, const Index3DS& cIndex, int i, int j, int k) {
					int I, J, K;
					double r, s, t;
					linearInterpolate(i, nx, I, r);
					linearInterpolate(j, ny, J, s);
					linearInterpolate(k, nz, K, t);

					if( cData[cIndex(I  ,J  ,K  )] == 1 ) {
						return 1;
					}
					return -1;
				}

#ifdef BCMT_NAMESPACE
} // namespace BCMT_NAMESPACE
#endif

