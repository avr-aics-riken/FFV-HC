#include "Scalar3DUpdater2.h"

#ifdef BCMT_NAMESPACE
namespace BCMT_NAMESPACE {
#endif

/*
template <>
int Scalar3DUpdater2<int>::interpolateF2C(const int* fData, const Index3DS& fIndex, int I, int J, int K) {
std::cout << "ok" << std::endl;
					int i = 2 * I;
					int j = 2 * J;
					int k = 2 * K;
					if( (double)(fData[fIndex(i  ,j  ,k  )]) > 0.5 &&
							(double)(fData[fIndex(i+1,j  ,k  )]) > 0.5 &&
							(double)(fData[fIndex(i  ,j+1,k  )]) > 0.5 &&
							(double)(fData[fIndex(i  ,j  ,k+1)]) > 0.5 &&
							(double)(fData[fIndex(i  ,j+1,k+1)]) > 0.5 &&
							(double)(fData[fIndex(i+1,j  ,k+1)]) > 0.5 &&
							(double)(fData[fIndex(i+1,j+1,k  )]) > 0.5 &&
							(double)(fData[fIndex(i+1,j+1,k+1)]) > 0.5 ) {
						return 1;
					} 
					return 0;

					return 0.125 * (fData[fIndex(i  ,j  ,k  )] + fData[fIndex(i+1,j  ,k  )]
							+ fData[fIndex(i  ,j+1,k  )] + fData[fIndex(i+1,j+1,k  )]
							+ fData[fIndex(i  ,j  ,k+1)] + fData[fIndex(i+1,j  ,k+1)]
							+ fData[fIndex(i  ,j+1,k+1)] + fData[fIndex(i+1,j+1,k+1)]);
				}


template <>
int Scalar3DUpdater2<int>::interpolateC2F(const int* cData, const Index3DS& cIndex, int i, int j, int k) {
std::cout << "ok2" << std::endl;
					int I, J, K;
					double r, s, t;
					linearInterpolate(i, nx, I, r);
					linearInterpolate(j, ny, J, s);
					linearInterpolate(k, nz, K, t);

					if( (double)(cData[cIndex(I  ,J  ,K  )]) > 0.5 ) {
						return 1;
					}
					return 0;

					return (1.0-t)*( 
							(1.0-s)*( (1.0-r)*cData[cIndex(I  ,J  ,K  )] + r*cData[cIndex(I+1,J  ,K  )] )
							+ s*( (1.0-r)*cData[cIndex(I  ,J+1,K  )] + r*cData[cIndex(I+1,J+1,K  )] )
							)
						+t*(
								(1.0-s)*( (1.0-r)*cData[cIndex(I  ,J  ,K+1)] + r*cData[cIndex(I+1,J  ,K+1)] )
								+ s*( (1.0-r)*cData[cIndex(I  ,J+1,K+1)] + r*cData[cIndex(I+1,J+1,K+1)] )
							 );
				}
*/

#ifdef BCMT_NAMESPACE
} // namespace BCMT_NAMESPACE
#endif

