///
/// @file Scalar3DUpdater4.h
/// @brief スカラデータクラス仮想セルアップデータ
///

#ifndef SCALAR_3D_UPDATER4_H
#define SCALAR_3D_UPDATER4_H

#include "BCMTools.h"
#include "VCUpdater.h"
#include "Scalar3D.h"
///#include "real.h"
#include "omp.h"

#ifdef BCMT_NAMESPACE
namespace BCMT_NAMESPACE {
#endif


/// スカラデータクラス仮想セルアップデータ.
///
///  @note 通信と補間の順序は，簡単のためL→L+1もL+1→Lも，
///        送信元で補間を行なってから通信．
///
///  @todo 補間計算部分をFortranで実装
///
///
template <typename T>
class Scalar3DUpdater4 : public VCUpdater {

public:

private:

  Scalar3D<T>* dataClass;   ///< 仮想セル同期対象データクラス

  T* sendBuffer[NUM_FACE][NUM_SUBFACE];  ///< 送信データバッファテーブル
  T* recvBuffer[NUM_FACE][NUM_SUBFACE];  ///< 受信データバッファテーブル

  Scalar3D<T>* neighborDataClass[NUM_FACE][NUM_SUBFACE];  ///< 隣接データクラステーブル

  int nx, ny, nz, vc;
public:
  /// コンストラクタ.
  ///
  ///  @param[in] neighborInfo 隣接情報配列
  ///  @param[in] comm MPIコミュニケータ(ディフォルトMPI::COMM_WORLD)
  ///
  Scalar3DUpdater4(const NeighborInfo* neighborInfo,
//                  const MPI_Comm comm = MPI_COMM_WORLD)
                  const MPI::Comm& comm = MPI::COMM_WORLD)
    : VCUpdater(neighborInfo, comm) {
    clearCommBufferPointer();
    clearNeighbor();

  }

  /// デストラクタ.
  ~Scalar3DUpdater4() {}

  /// 仮想セル同期対象データクラスを登録.
  void setDataClass(DataClass* dc) {
    dataClass = dynamic_cast<Scalar3D<T>*>(dc);
    nx = dataClass->getSizeX();
    ny = dataClass->getSizeY();
    nz = dataClass->getSizeZ();
    vc = dataClass->getVCSize();
  }

  /// 仮想セル同期データ送信に必要なバッファサイズを取得(同レベル間).
  size_t getSendBufferByteSize(Face face) const {
    return sizeof(T) * getCommBufferSize(face);
  }

  /// 仮想セル同期データ送信に必要なバッファサイズを取得(レベルL+1→L).
  size_t getSendBufferByteSizeF2C(Face face, Subface subface) const {
    return sizeof(T) * getCommBufferSize(face) / 4;
  }

  /// 仮想セル同期データ送信に必要なバッファサイズを取得(レベルL→L+1).
  size_t getSendBufferByteSizeC2F(Face face, Subface subface) const {
    return sizeof(T) * getCommBufferSize(face);
  }

  /// 仮想セル同期データ受信に必要なバッファサイズを取得(同レベル間).
  size_t getRecvBufferByteSize(Face face) const {
    return sizeof(T) * getCommBufferSize(face);
  }

  /// 仮想セル同期データ受信に必要なバッファサイズを取得(レベルL+1→L).
  size_t getRecvBufferByteSizeF2C(Face face, Subface subface) const {
    return sizeof(T) * getCommBufferSize(face) / 4;
  }

  /// 仮想セル同期データ受信に必要なバッファサイズを取得(レベルL→L+1).
  size_t getRecvBufferByteSizeC2F(Face face, Subface subface) const {
    return sizeof(T) * getCommBufferSize(face);
  }

  /// 仮想セル同期データ送信バッファ用PointerSetterオブジェクトを取得.
  PointerSetterBase* getSendBufferPointerSetter(Face face, Subface subface) {
    return new PointerSetter<T>(&sendBuffer[face][subface]);
  }

  /// 仮想セル同期データ受信バッファ用PointerSetterオブジェクトを取得.
  PointerSetterBase* getRecvBufferPointerSetter(Face face, Subface subface) {
    return new PointerSetter<T>(&recvBuffer[face][subface]);
  }


public:

  /// 同並列計算ノード内の隣接データクラスを登録.
  void setNeighbor(Face face, Subface subface, DataClass* dataClass) {
    neighborDataClass[face][subface] = dynamic_cast<Scalar3D<T>*>(dataClass);
  }

  /// 隣接データクラスの登録解除.
  void clearNeighbor(Face face, Subface subface) {
    neighborDataClass[face][subface] = 0;
  }

  /// 隣接データクラスの登録解除.
  void clearNeighbor() {
    for (int i = 0; i < NUM_FACE; ++i) {
      for (int j = 0; j < NUM_SUBFACE; ++j) {
        clearNeighbor(Face(i), Subface(j));
      }
    }
  }

  /// 通信バッファテーブルのエントリをクリア.
  void clearCommBufferPointer(Face face, Subface subface) {
    sendBuffer[face][subface] = recvBuffer[face][subface] = 0;
  }

  /// 通信バッファテーブルをクリア.
  void clearCommBufferPointer() {
    for (int i = 0; i < NUM_FACE; ++i) {
      for (int j = 0; j < NUM_SUBFACE; ++j) {
        clearCommBufferPointer(Face(i), Subface(j));
      }
    }
  }

private:

  /// 通信バッファサイズを計算.
  size_t getCommBufferSize(Face face) const {
    switch (face) {
      case X_M:
      case X_P:
        return ny * nz * vc;
      case Y_M:
      case Y_P:
        return nz * nx * vc;
      case Z_M:
      case Z_P:
        return nx * ny * vc;
      default:
        Exit(EX_FAILURE);
    }
    /* NOTREACHED */
  }

/*
  T interpolateC2F(const T* cData, const Index3DS& cIndex, int i, int j, int k) {
    int I, J, K;
    double r, s, t;
    linearInterpolate(i, nx, I, r);
    linearInterpolate(j, ny, J, s);
    linearInterpolate(k, nz, K, t);

    return (1.0-t)*(
         (1.0-s)*( (1.0-r)*cData[cIndex(I  ,J  ,K  )] + r*cData[cIndex(I+1,J  ,K  )] )
             + s*( (1.0-r)*cData[cIndex(I  ,J+1,K  )] + r*cData[cIndex(I+1,J+1,K  )] )
        )
       +t*(
         (1.0-s)*( (1.0-r)*cData[cIndex(I  ,J  ,K+1)] + r*cData[cIndex(I+1,J  ,K+1)] )
             + s*( (1.0-r)*cData[cIndex(I  ,J+1,K+1)] + r*cData[cIndex(I+1,J+1,K+1)] )
        );
  }

  void linearInterpolate(int i, int n, int& I, double& r) {
#if 1
    I = std::min(std::max(i/2 - 1 + i%2, 0), n - 2);
    r = -0.25 + 0.5 * i - double(I);
#else
    if (i == 0) {
      I = 0;
      r = -0.25;
    }
    else if (i == 2*n-1) {
      I = n - 2;
      r = 1.25;
    }
    else if (i%2 == 0) {
      I = i/2 - 1;
      r = 0.75;
    }
    else {
      I = i/2;
      r = 0.25;
    }
#endif
  }
*/

//FEAST.s
  /// レベルL+1→Lの線形補間 (2x2)
  T interpolateF2C_2x2(const T* fData, const Index3DS& fIndex, int I, int J, int K) {

    int is=2*I;
    int js=2*J;
    int ks=2*K; 
    int ie=is+1;
    int je=js+1;
    int ke=ks+1;
 
    T XT,YT,ZT;
    XT=(T)(is+ie)*0.5;
    YT=(T)(js+je)*0.5;
    ZT=(T)(ks+ke)*0.5;

    return interpolate_NxN(fData, fIndex, is, ie, js, je, ks, ke, XT, YT, ZT);

  }

  /// レベルL+1→Lの線形補間 (4x4)
  T interpolateF2C_4x4(const T* fData, const Index3DS& fIndex, int I, int J, int K) {

    int is=2*I-1;
    int js=2*J-1;
    int ks=2*K-1;
    int ie=is+3;
    int je=js+3;
    int ke=ks+3;

    T XT,YT,ZT;
    XT=(T)(is+ie)*0.5;
    YT=(T)(js+je)*0.5;
    ZT=(T)(ks+ke)*0.5;

    return interpolate_NxN(fData, fIndex, is, ie, js, je, ks, ke, XT, YT, ZT);

  }

  /// 高次の線形補間 (NxN)
  T interpolate_NxN(const T* Data, const Index3DS& Index, int is, int ie, int js, int je, int ks, int ke, T XT, T YT, T ZT) {

    if( is < 0   ) is=0;
    if( ie >= nx ) ie=nx-1;
    if( js < 0   ) js=0;
    if( je >= ny ) je=ny-1;
    if( ks < 0   ) ks=0;
    if( ke >= nz ) ke=nz-1;

    //距離による重み係数の計算
    T wx[4],wy[4],wz[4];
    //x方向
    for(int i=is; i<=ie; i++) {
      int ii=i-is;
      wx[ii]=1.0;
      for(int n=is; n<=ie; n++) { 
        if( i != n ) wx[ii]=wx[ii]*( (XT-(T)n)/((T)i-(T)n) );
      }
    }

    //y方向
    for(int j=js; j<=je; j++) {
      int jj=j-js;
      wy[jj]=1.0;
      for(int n=js; n<=je; n++) {
        if( j != n ) wy[jj]=wy[jj]*( (YT-(T)n)/((T)j-(T)n) );
      }
    } 

    //z方向
    for(int k=ks; k<=ke; k++) {
      int kk=k-ks;
      wz[kk]=1.0;
      for(int n=ks; n<=ke; n++) {
        if( k != n ) wz[kk]=wz[kk]*( (ZT-(T)n)/((T)k-(T)n) );
      }
    } 

    //QTの計算
    T QT=0.0;
    int i,j,k;
    for(int kk=ks; kk<=ke; kk++) {
    k=kk-ks;
    for(int jj=js; jj<=je; jj++) {
    j=jj-js;
    for(int ii=is; ii<=ie; ii++) {
      i=ii-is;
      QT=QT+( Data[Index(ii,jj,kk)]*wx[i]*wy[j]*wz[k] );
    }}}

    return QT;

  }

  /// レベルL+1→Lの線形補間 (2x2)
  T interpolateF2C_2x2(const Scalar3D<T>& f, int I, int J, int K) {

    int is=2*I;
    int js=2*J;
    int ks=2*K; 
    int ie=is+1;
    int je=js+1;
    int ke=ks+1;
 
    T XT,YT,ZT;
    XT=(T)(is+ie)*0.5;
    YT=(T)(js+je)*0.5;
    ZT=(T)(ks+ke)*0.5;

    return interpolate_NxN(f, is, ie, js, je, ks, ke, XT, YT, ZT);

  }

  /// レベルL+1→Lの線形補間 (4x4)
  T interpolateF2C_4x4(const Scalar3D<T>& f, const int I, int J, int K) {

    int is=2*I-1;
    int js=2*J-1;
    int ks=2*K-1;
    int ie=is+3;
    int je=js+3;
    int ke=ks+3;

    T XT,YT,ZT;
    XT=(T)(is+ie)*0.5;
    YT=(T)(js+je)*0.5;
    ZT=(T)(ks+ke)*0.5;

    return interpolate_NxN(f, is, ie, js, je, ks, ke, XT, YT, ZT);

  }

  /// 高次の線形補間 (NxN)
  T interpolate_NxN(const Scalar3D<T>& f, int is, int ie, int js, int je, int ks, int ke, T XT, T YT, T ZT) {

    if( is < 0   ) is=0;
    if( ie >= nx ) ie=nx-1;
    if( js < 0   ) js=0;
    if( je >= ny ) je=ny-1;
    if( ks < 0   ) ks=0;
    if( ke >= nz ) ke=nz-1;

    //距離による重み係数の計算
    T wx[4],wy[4],wz[4];
    //x方向
    for(int i=is; i<=ie; i++) {
      int ii=i-is;
      wx[ii]=1.0;
      for(int n=is; n<=ie; n++) { 
        if( i != n ) wx[ii]=wx[ii]*( (XT-(T)n)/((T)i-(T)n) );
      }
    }

    //y方向
    for(int j=js; j<=je; j++) {
      int jj=j-js;
      wy[jj]=1.0;
      for(int n=js; n<=je; n++) {
        if( j != n ) wy[jj]=wy[jj]*( (YT-(T)n)/((T)j-(T)n) );
      }
    } 

    //z方向
    for(int k=ks; k<=ke; k++) {
      int kk=k-ks;
      wz[kk]=1.0;
      for(int n=ks; n<=ke; n++) {
        if( k != n ) wz[kk]=wz[kk]*( (ZT-(T)n)/((T)k-(T)n) );
      }
    } 

    //QTの計算
    T QT=0.0;
    int i,j,k;
    for(int kk=ks; kk<=ke; kk++) {
    k=kk-ks;
    for(int jj=js; jj<=je; jj++) {
    j=jj-js;
    for(int ii=is; ii<=ie; ii++) {
      i=ii-is;
      QT=QT+( f(ii,jj,kk)*wx[i]*wy[j]*wz[k] );
    }}}

    return QT;

  }

//FEAST.e

//FEAST.s
  /// X C2F補間
  T interpolateC2F_X(const T* cData, const Index3DS& cIndex, int i, int j, int k) 
  {

    T XT,YT,ZT;
    int is,js,ks,ie,je,ke;
    targetXYZ(i,j,k,XT,YT,ZT,is,js,ks);

    //スタートエンドの補正
    if( j >= ny ) js=js+1; 
    else          js=js-1;
    if( k >= nz ) ks=ks+1;
    else          ks=ks-1;
    ie=is+2;
    je=js+2;
    ke=ks+2;

    return interpolate_NxN(cData, cIndex, is, ie, js, je, ks, ke, XT, YT, ZT);
  }

  /// X C2F補間
  T interpolateC2F_X(const Scalar3D<T>& c, int i, int j, int k) 
  {

    T XT,YT,ZT;
    int is,js,ks,ie,je,ke;
    targetXYZ(i,j,k,XT,YT,ZT,is,js,ks);

    //スタートエンドの補正
    if( j >= ny ) js=js+1; 
    else          js=js-1;
    if( k >= nz ) ks=ks+1;
    else          ks=ks-1;
    ie=is+2;
    je=js+2;
    ke=ks+2;

    return interpolate_NxN(c, is, ie, js, je, ks, ke, XT, YT, ZT);
  }

  /// Y C2F補間
  T interpolateC2F_Y(const T* cData, const Index3DS& cIndex, int i, int j, int k) 
  {

    T XT,YT,ZT;
    int is,js,ks,ie,je,ke;
    targetXYZ(i,j,k,XT,YT,ZT,is,js,ks);

    //スタートエンドの補正
    if( i >= nx ) is=is+1;
    else          is=is-1;
    if( k >= nz ) ks=ks+1;
    else          ks=ks-1;
    ie=is+2;
    je=js+2;
    ke=ks+2;

    return interpolate_NxN(cData, cIndex, is, ie, js, je, ks, ke, XT, YT, ZT);
  }

  T interpolateC2F_Y(const Scalar3D<T>& c, int i, int j, int k) 
  {

    T XT,YT,ZT;
    int is,js,ks,ie,je,ke;
    targetXYZ(i,j,k,XT,YT,ZT,is,js,ks);

    //スタートエンドの補正
    if( i >= nx ) is=is+1;
    else          is=is-1;
    if( k >= nz ) ks=ks+1;
    else          ks=ks-1;
    ie=is+2;
    je=js+2;
    ke=ks+2;

    return interpolate_NxN(c, is, ie, js, je, ks, ke, XT, YT, ZT);
  }

  /// Z C2F補間
  T interpolateC2F_Z(const T* cData, const Index3DS& cIndex, int i, int j, int k) 
  {

    T XT,YT,ZT;
    int is,js,ks,ie,je,ke;
    targetXYZ(i,j,k,XT,YT,ZT,is,js,ks);

    //スタートエンドの補正
    if( i >= nx ) is=is+1;
    else          is=is-1;
    if( j >= ny ) js=js+1;
    else          js=js-1;
    ie=is+2;
    je=js+2;
    ke=ks+2;

    return interpolate_NxN(cData, cIndex, is, ie, js, je, ks, ke, XT, YT, ZT);
  }

  T interpolateC2F_Z(const Scalar3D<T>& c, int i, int j, int k) 
  {

    T XT,YT,ZT;
    int is,js,ks,ie,je,ke;
    targetXYZ(i,j,k,XT,YT,ZT,is,js,ks);

    //スタートエンドの補正
    if( i >= nx ) is=is+1;
    else          is=is-1;
    if( j >= ny ) js=js+1;
    else          js=js-1;
    ie=is+2;
    je=js+2;
    ke=ks+2;

    return interpolate_NxN(c, is, ie, js, je, ks, ke, XT, YT, ZT);
  }

  /// target座標値を求める
  void targetXYZ(int i, int j, int k, T &XT, T &YT, T &ZT, int &is, int &js, int &ks)
  {

    XT = T(i/2);
    if( i%2 == 0 ) XT = XT-0.25;
    else           XT = XT+0.25;
    YT = T(j/2);
    if( j%2 == 0 ) YT = YT-0.25;
    else           YT = YT+0.25;
    ZT = T(k/2);
    if( k%2 == 0 ) ZT = ZT-0.25;
    else           ZT = ZT+0.25;

    if( i>=nx ) is=i/2-2;
    else        is=i/2;
    if( j>=ny ) js=j/2-2;
    else        js=j/2;
    if( k>=nz ) ks=k/2-2;
    else        ks=k/2;

  } 

//FEAST.e

/// 隣接データクラスから仮想セルデータをコピー(同レベル間).
void copyFromNeighbor(Face face)
{

  Scalar3D<T>* dc = neighborDataClass[face][0];
  if (!dc) return;
  switch (face) {
    case X_M:
      dataClass->copyFromDataClass(-vc, 0, 0,  dc->getSizeX()-vc, 0, 0, vc, ny, nz,  dc);
      break;
    case X_P:
      dataClass->copyFromDataClass(nx, 0, 0,  0, 0, 0,  vc, ny, nz,  dc);
      break;
    case Y_M:
      dataClass->copyFromDataClass(0, -vc, 0,  0, dc->getSizeY()-vc, 0, nx, vc, nz,  dc);
      break;
    case Y_P:
      dataClass->copyFromDataClass(0, ny, 0,  0, 0, 0,  nx, vc, nz,  dc);
      break;
    case Z_M:
      dataClass->copyFromDataClass(0, 0, -vc,  0, 0, dc->getSizeZ()-vc, nx, ny, vc,  dc);
      break;
    case Z_P:
      dataClass->copyFromDataClass(0, 0, nz,  0, 0, 0,  nx, ny, vc,  dc);
      break;
    default:
      break;
  }
}


/// 隣接データクラスから仮想セルデータをコピー(レベルL+1→L).
void copyFromNeighborF2C(Face face, Subface subface)
{

  T* cData = dataClass->getData();
  Index3DS cIndex = dataClass->getIndex();
  Scalar3D<T>* f = neighborDataClass[face][subface];
  T* fData = f->getData();
  Index3DS fIndex = f->getIndex();

  copyFromNeighborF2C_0(nx, ny, nz, vc, face, subface, fData, fIndex, cData, cIndex);
}


/// 隣接データクラスから仮想セルデータをコピー(レベルL→L+1).
void copyFromNeighborC2F(Face face, Subface subface)
{

  T* fData = dataClass->getData();
  Index3DS fIndex = dataClass->getIndex();
  Scalar3D<T>* c = neighborDataClass[face][0];
  T* cData = c->getData();
  Index3DS cIndex = c->getIndex();

  copyFromNeighborC2F_0(nx, ny, nz, vc, face, subface, cData, cIndex, fData, fIndex);
}


/// 送信バッファに仮想セルデータをコピー(同レベル間).
void copyToCommBuffer(Face face)
{

  T* buffer = sendBuffer[face][0];
  if (!buffer) return;
  switch (face) {
    case X_M:
      dataClass->copyToBuffer(0, 0, 0,  vc, ny, nz,  buffer);
      break;
    case X_P:
      dataClass->copyToBuffer(nx-vc, 0, 0,  vc, ny, nz,  buffer);
      break;
    case Y_M:
      dataClass->copyToBuffer(0, 0, 0,  nx, vc, nz,  buffer);
      break;
    case Y_P:
      dataClass->copyToBuffer(0, ny-vc, 0,  nx, vc, nz,  buffer);
      break;
    case Z_M:
      dataClass->copyToBuffer(0, 0, 0,  nx, ny, vc,  buffer);
      break;
    case Z_P:
      dataClass->copyToBuffer(0, 0, nz-vc,  nx, ny, vc,  buffer);
      break;
    default:
      break;
  }
}


/// 送信バッファに仮想セルデータをコピー(レベルL+1→L).
void copyToCommBufferF2C(Face face, Subface subface)
{
  T* buffer = sendBuffer[face][0];
  T* fData = dataClass->getData();
  Index3DS fIndex = dataClass->getIndex();

  copyToCommBufferF2C_0(nx, ny, nz, vc, face, subface, fData, fIndex, buffer);
}


/// 送信バッファに仮想セルデータをコピー(レベルL→L+1).
void copyToCommBufferC2F(Face face, Subface subface)
{
  T* cData = dataClass->getData();
  Index3DS cIndex = dataClass->getIndex();
  T* buffer = sendBuffer[face][subface];

  copyToCommBufferC2F_0(nx, ny, nz, vc, face, subface, cData, cIndex, buffer);
}


/// 受信バッファから仮想セルデータをコピー(同レベル間).
void copyFromCommBuffer(Face face)
{

  T* buffer = recvBuffer[face][0];
  if (!buffer) return;
  switch (face) {
    case X_M:
      dataClass->copyFromBuffer(-vc, 0, 0,  vc, ny, nz,  buffer);
      break;
    case X_P:
      dataClass->copyFromBuffer(nx, 0, 0,  vc, ny, nz,  buffer);
      break;
    case Y_M:
      dataClass->copyFromBuffer(0, -vc, 0,  nx, vc, nz,  buffer);
      break;
    case Y_P:
      dataClass->copyFromBuffer(0, ny, 0,  nx, vc, nz,  buffer);
      break;
    case Z_M:
      dataClass->copyFromBuffer(0, 0, -vc,  nx, ny, vc,  buffer);
      break;
    case Z_P:
      dataClass->copyFromBuffer(0, 0, nz,  nx, ny, vc,  buffer);
      break;
    default:
      break;
  }
}

/// 受信バッファから仮想セルデータをコピー(レベルL+1→L).
void copyFromCommBufferF2C(Face face, Subface subface)
{

  copyFromCommBufferF2C_0(face, subface);

/*
  T* buffer = recvBuffer[face][subface];
  switch (face) {
    case X_M:
    {
      int j0 = (ny/2) * subfaceOrigin0(subface);
      int k0 = (nz/2) * subfaceOrigin1(subface);
      dataClass->copyFromBuffer(-vc, j0, k0, vc, ny/2, nz/2, buffer);
      break;
    }
    case X_P:
    {
      int j0 = (ny/2) * subfaceOrigin0(subface);
      int k0 = (nz/2) * subfaceOrigin1(subface);
      dataClass->copyFromBuffer(nx, j0, k0, vc, ny/2, nz/2, buffer);
      break;
    }
    case Y_M:
    {
      int k0 = (nz/2) * subfaceOrigin0(subface);
      int i0 = (nx/2) * subfaceOrigin1(subface);
      dataClass->copyFromBuffer(i0, -vc, k0,  nx/2, vc, nz/2,  buffer);
      break;
    }
    case Y_P:
    {
      int k0 = (nz/2) * subfaceOrigin0(subface);
      int i0 = (nx/2) * subfaceOrigin1(subface);
      dataClass->copyFromBuffer(i0, ny, k0,  nx/2, vc, nz/2,  buffer);
      break;
    }
    case Z_M:
    {
      int i0 = (nx/2) * subfaceOrigin0(subface);
      int j0 = (ny/2) * subfaceOrigin1(subface);
      dataClass->copyFromBuffer(i0, j0, -vc,  nx/2, ny/2, vc,  buffer);
      break;
    }
    case Z_P:
    {
      int i0 = (nx/2) * subfaceOrigin0(subface);
      int j0 = (ny/2) * subfaceOrigin1(subface);
      dataClass->copyFromBuffer(i0, j0, nz,  nx/2, ny/2, vc,  buffer);
      break;
    }
    default:
      break;
  }
*/

}

/// 受信バッファから仮想セルデータをコピー(レベルL+1→L).
void copyFromCommBufferF2C_0(Face face, Subface subface)
{

  T* cData = dataClass->getData();
  Index3DS cIndex = dataClass->getIndex();

  T* buffer = recvBuffer[face][subface];
  if(!buffer) return;

  switch (face) {
    case X_M:
    {
      int j0 = (ny/2) * subfaceOrigin0(subface);
      int k0 = (nz/2) * subfaceOrigin1(subface);
      //dataClass->copyFromBuffer(-vc, j0, k0, vc, ny/2, nz/2, buffer);
#pragma omp parallel for if(nz >= 16)
      for (int k = k0; k < k0 + nz/2; ++k) {
      for (int j = j0; j < j0 + ny/2; ++j) {
      for (int i = 0; i < vc; ++i) {
        cData[cIndex(-1-i,j,k)] = buffer[i + vc*(j-j0) + (vc*ny/2)*(k-k0)];
      }}}
      break;
    }
    case X_P:
    {
      int j0 = (ny/2) * subfaceOrigin0(subface);
      int k0 = (nz/2) * subfaceOrigin1(subface);
      //dataClass->copyFromBuffer(nx, j0, k0, vc, ny/2, nz/2, buffer);
#pragma omp parallel for if(nz >= 16)
      for (int k = k0; k < k0 + nz/2; ++k) {
      for (int j = j0; j < j0 + ny/2; ++j) {
      for (int i = nx; i < vc + nx;   ++i) {
        cData[cIndex(i,j,k)] = buffer[i-nx + vc*(j-j0) + (vc*ny/2)*(k-k0)];
      }}}
      break;
    }
    case Y_M:
    {
      int k0 = (nz/2) * subfaceOrigin0(subface);
      int i0 = (nx/2) * subfaceOrigin1(subface);
      //dataClass->copyFromBuffer(i0, -vc, k0,  nx/2, vc, nz/2,  buffer);
#pragma omp parallel for if(nz >= 16)
      for (int k = k0; k < k0 + nz/2; ++k) {
      for (int j = 0;  j < vc;        ++j) {
      for (int i = i0; i < i0 + nx/2; ++i) {
        cData[cIndex(i,-1-j,k)] = buffer[i-i0 + nx/2*j + (nx/2*vc)*(k-k0)];
      }}}
      break;
    }
    case Y_P:
    {
      int k0 = (nz/2) * subfaceOrigin0(subface);
      int i0 = (nx/2) * subfaceOrigin1(subface);
      //dataClass->copyFromBuffer(i0, ny, k0,  nx/2, vc, nz/2,  buffer);
#pragma omp parallel for if(nz >= 16)
      for (int k = k0; k < k0 + nz/2; ++k) {
      for (int j = ny; j < vc + ny;   ++j) {
      for (int i = i0; i < i0 + nx/2; ++i) {
        cData[cIndex(i,j,k)] = buffer[i-i0 + nx/2*(j-ny) + (nx/2*vc)*(k-k0)];
      }}}
      break;
    }
    case Z_M:
    {
      int i0 = (nx/2) * subfaceOrigin0(subface);
      int j0 = (ny/2) * subfaceOrigin1(subface);
      //dataClass->copyFromBuffer(i0, j0, -vc,  nx/2, ny/2, vc,  buffer);
#pragma omp parallel for if(nz >= 16)
      for (int k = 0;  k < vc;        ++k) {
      for (int j = j0; j < j0 + ny/2; ++j) {
      for (int i = i0; i < i0 + nx/2; ++i) {
        cData[cIndex(i,j,-1-k)] = buffer[i-i0 + nx/2*(j-j0) + (nx/2*ny/2)*k];
      }}}
      break;
    }
    case Z_P:
    {
      int i0 = (nx/2) * subfaceOrigin0(subface);
      int j0 = (ny/2) * subfaceOrigin1(subface);
      //dataClass->copyFromBuffer(i0, j0, nz,  nx/2, ny/2, vc,  buffer);
#pragma omp parallel for if(nz >= 16)
      for (int k = nz; k < vc + nz;   ++k) {
      for (int j = j0; j < j0 + ny/2; ++j) {
      for (int i = i0; i < i0 + nx/2; ++i) {
        cData[cIndex(i,j,k)] = buffer[i-i0 + nx/2*(j-j0) + (nx/2*ny/2)*(k-nz)];
      }}}
      break;
    }
    default:
      break;
  }
}


/// 受信バッファから仮想セルデータをコピー(レベルL→L+1).
void copyFromCommBufferC2F(Face face, Subface subface)
{

  //copyFromCommBuffer(face);
  copyFromCommBufferC2F_0(face);

}

void copyFromCommBufferC2F_0(Face face)
{
  T* fData = dataClass->getData();
  Index3DS fIndex = dataClass->getIndex();

  T* buffer = recvBuffer[face][0];
  if (!buffer) return;
  switch (face) {
    case X_M:
#pragma omp parallel for if(nz >= 16)
      for (int k = 0; k < nz; ++k) {
      for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < vc; ++i) {
        fData[fIndex(-1-i,j,k)] = buffer[i + vc*j + (vc*ny)*k];
      }}}
      break;

    case X_P:
#pragma omp parallel for if(nz >= 16)
      for (int k = 0;  k < nz; ++k) {
      for (int j = 0;  j < ny; ++j) {
      for (int i = nx; i < nx+vc; ++i) {
        fData[fIndex(i,j,k)] = buffer[i-nx + vc*j + (vc*ny)*k];
      }}}
      break;

    case Y_M:
#pragma omp parallel for if(nz >= 16)
      for (int k = 0; k < nz; ++k) {
      for (int j = 0; j < vc; ++j) {
      for (int i = 0; i < nx; ++i) {
        fData[fIndex(i,-1-j,k)] = buffer[i + nx*j + (nx*vc)*k];
      }}}
      break;

    case Y_P:
#pragma omp parallel for if(nz >= 16)
      for (int k = 0;  k < nz; ++k) {
      for (int j = ny; j < ny+vc; ++j) {
      for (int i = 0;  i < nx; ++i) {
        fData[fIndex(i,j,k)] = buffer[i + nx*(j-ny) + (nx*vc)*k];
      }}}
      break;

    case Z_M:
#pragma omp parallel for if(nz >= 16)
      for (int k = 0; k < vc; ++k) {
      for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        fData[fIndex(i,j,-1-k)] = buffer[i + nx*j + (nx*ny)*k];
      }}}
      break;

    case Z_P:
#pragma omp parallel for if(nz >= 16)
      for (int k = nz; k < nz+vc; ++k) {
      for (int j = 0;  j < ny; ++j) {
      for (int i = 0;  i < nx; ++i) {
        fData[fIndex(i,j,k)] = buffer[i + nx*j + (nx*ny)*(k-nz)];
      }}}
      break;
    default:
      break;
  }

}

 
void copyFromNeighborF2C_0(int nx, int ny, int nz, int vc,
                                               Face face, Subface subface,
                                               const T* fData, Index3DS fIndex,
                                               T* cData, Index3DS cIndex)
{
  switch (face) {
//X
    case X_M:
    {
      int j0 = (ny/2) * subfaceOrigin0(subface);
      int k0 = (nz/2) * subfaceOrigin1(subface);
//#ifdef _BLOCK_IS_LARGE_
//#pragma omp parallel for collapse(3)
//#endif
      for (int k = 0; k < nz/2; k++) {
      for (int j = 0; j < ny/2; j++) {
        for (int i = 0; i < min(1,vc); i++ ) {
          cData[cIndex(-1-i, j+j0, k+k0)] = interpolateF2C_2x2(fData, fIndex, nx/2-1-i, j, k);
        }
        for (int i = 1; i < vc; i++) {
          cData[cIndex(-1-i, j+j0, k+k0)] = interpolateF2C_4x4(fData, fIndex, nx/2-1-i, j, k);
        }
      }}
      break;
    }
    case X_P:
    {
      int j0 = (ny/2) * subfaceOrigin0(subface);
      int k0 = (nz/2) * subfaceOrigin1(subface);
//#pragma omp parallel for collapse(3)
      for (int k = 0; k < nz/2; k++) {
      for (int j = 0; j < ny/2; j++) {
        for (int i = 0; i < min(1,vc); i++) {
          cData[cIndex(i+nx, j+j0, k+k0)] = interpolateF2C_2x2(fData, fIndex, i, j, k);
        }
        for (int i = 1; i < vc; i++) {
          cData[cIndex(i+nx, j+j0, k+k0)] = interpolateF2C_4x4(fData, fIndex, i, j, k);
        }
      }}
      break;
    }
//Y
    case Y_M:
    {
      int k0 = (nz/2) * subfaceOrigin0(subface);
      int i0 = (nx/2) * subfaceOrigin1(subface);
//#pragma omp parallel for collapse(3)
      for (int k = 0; k < nz/2; k++) {
        for (int j = 0; j < min(1,vc); j++) {
        for (int i = 0; i < nx/2; i++) {
          cData[cIndex(i+i0, -1-j, k+k0)] = interpolateF2C_2x2(fData, fIndex, i, ny/2-1-j, k);
        }}
        for (int j = 1; j < vc; j++) {
        for (int i = 0; i < nx/2; i++) {
          cData[cIndex(i+i0, -1-j, k+k0)] = interpolateF2C_4x4(fData, fIndex, i, ny/2-1-j, k);
        }}
      }
      break;
    }
    case Y_P:
    {
      int k0 = (nz/2) * subfaceOrigin0(subface);
      int i0 = (nx/2) * subfaceOrigin1(subface);
//#pragma omp parallel for collapse(3)
      for (int k = 0; k < nz/2; k++) {
        for (int j = 0; j < min(1,vc); j++) {
        for (int i = 0; i < nx/2; i++) {
          cData[cIndex(i+i0, j+ny, k+k0)] = interpolateF2C_2x2(fData, fIndex, i, j, k);
        }}
        for (int j = 1; j < vc; j++) {
        for (int i = 0; i < nx/2; i++) {
          cData[cIndex(i+i0, j+ny, k+k0)] = interpolateF2C_4x4(fData, fIndex, i, j, k);
        }}
      }
      break;
    }
//Z
    case Z_M:
    {
      int i0 = (nx/2) * subfaceOrigin0(subface);
      int j0 = (ny/2) * subfaceOrigin1(subface);
//#pragma omp parallel for collapse(3)
      for (int k = 0; k < min(1,vc); k++) {
      for (int j = 0; j < ny/2; j++) {
      for (int i = 0; i < nx/2; i++) {
        cData[cIndex(i+i0, j+j0, -1-k)] = interpolateF2C_2x2(fData, fIndex, i, j, nz/2-1-k);
      }}}
      for (int k = 1; k < vc; k++) {
      for (int j = 0; j < ny/2; j++) {
      for (int i = 0; i < nx/2; i++) {
        cData[cIndex(i+i0, j+j0, -1-k)] = interpolateF2C_4x4(fData, fIndex, i, j, nz/2-1-k);
      }}}
      break;
    }
    case Z_P:
    {
      int i0 = (nx/2) * subfaceOrigin0(subface);
      int j0 = (ny/2) * subfaceOrigin1(subface);
//#pragma omp parallel for collapse(3)
      for (int k = 0; k < min(1,vc); k++) {
      for (int j = 0; j < ny/2; j++) {
      for (int i = 0; i < nx/2; i++) {
        cData[cIndex(i+i0, j+j0, k+nz)] = interpolateF2C_2x2(fData, fIndex, i, j, k);
      }}}
      for (int k = 1; k < vc; k++) {
      for (int j = 0; j < ny/2; j++) {
      for (int i = 0; i < nx/2; i++) {
        cData[cIndex(i+i0, j+j0, k+nz)] = interpolateF2C_4x4(fData, fIndex, i, j, k);
      }}}
      break;
    }
    default:
      break;
  }
}


void copyFromNeighborC2F_0(int nx, int ny, int nz, int vc,
                                               Face face, Subface subface,
                                               const T* cData, Index3DS cIndex,
                                               T* fData, Index3DS fIndex)
{

  switch (face) {
    case X_M:
    {
      int J0 = ny * subfaceOrigin0(subface);
      int K0 = nz * subfaceOrigin1(subface);
//#pragma omp parallel for collapse(3)
      for (int K = 0; K < nz; K++) {
      for (int J = 0; J < ny; J++) {
      for (int I = 0; I < vc; I++) {
        fData[fIndex(-1-I, J, K)] = interpolateC2F_X(cData, cIndex, 2*nx-1-I, J+J0, K+K0);
      }}}
      break;
    }
    case X_P:
    {
      int J0 = ny * subfaceOrigin0(subface);
      int K0 = nz * subfaceOrigin1(subface);
//#pragma omp parallel for collapse(3)
      for (int K = 0; K < nz; K++) {
      for (int J = 0; J < ny; J++) {
      for (int I = 0; I < vc; I++) {
        fData[fIndex(I+nx, J, K)] = interpolateC2F_X(cData, cIndex, I, J+J0, K+K0);
      }}}
      break;
    }
    case Y_M:
    {
      int K0 = nz * subfaceOrigin0(subface);
      int I0 = nx * subfaceOrigin1(subface);
//#pragma omp parallel for collapse(3)
      for (int K = 0; K < nz; K++) {
      for (int J = 0; J < vc; J++) {
      for (int I = 0; I < nx; I++) {
        fData[fIndex(I, -1-J, K)] = interpolateC2F_Y(cData, cIndex, I+I0, 2*ny-1-J, K+K0);
      }}}
      break;
    }
    case Y_P:
    {
      int K0 = nz * subfaceOrigin0(subface);
      int I0 = nx * subfaceOrigin1(subface);
//#pragma omp parallel for collapse(3)
      for (int K = 0; K < nz; K++) {
      for (int J = 0; J < vc; J++) {
      for (int I = 0; I < nx; I++) {
        fData[fIndex(I, J+ny, K)] = interpolateC2F_Y(cData, cIndex, I+I0, J, K+K0);
      }}}
      break;
    }
    case Z_M:
    {
      int I0 = nx * subfaceOrigin0(subface);
      int J0 = ny * subfaceOrigin1(subface);
//#pragma omp parallel for collapse(3)
      for (int K = 0; K < vc; K++) {
      for (int J = 0; J < ny; J++) {
      for (int I = 0; I < nx; I++) {
        fData[fIndex(I, J, -1-K)] = interpolateC2F_Z(cData, cIndex, I+I0, J+J0, 2*nz-1-K);
      }}}
      break;
    }
    case Z_P:
    {
      int I0 = nx * subfaceOrigin0(subface);
      int J0 = ny * subfaceOrigin1(subface);
//#pragma omp parallel for collapse(3)
      for (int K = 0; K < vc; K++) {
      for (int J = 0; J < ny; J++) {
      for (int I = 0; I < nx; I++) {
        fData[fIndex(I, J, K+nz)] = interpolateC2F_Z(cData, cIndex, I+I0, J+J0, K);
      }}}
      break;
    }
    default:
      break;
  }
}


void copyToCommBufferC2F_0(int nx, int ny, int nz, int vc,
                                               Face face, Subface subface,
                                               const T* cData, Index3DS cIndex,
                                               T* buffer)
{
  int ii = 0;
  switch (face) {
    case X_M:
    {
      int J0 = ny * subfaceOrigin0(subface);
      int K0 = nz * subfaceOrigin1(subface);
//#pragma omp parallel for collapse(3)
      for (int K = 0; K < nz; K++) {
      for (int J = 0; J < ny; J++) {
      for (int I = 0; I < vc; I++) {
        int m = I + vc*(J + ny*K);
        buffer[m] = interpolateC2F_X(cData, cIndex, I, J+J0, K+K0);
      }}}
      break;
    }
    case X_P:
    {
      int J0 = ny * subfaceOrigin0(subface);
      int K0 = nz * subfaceOrigin1(subface);
//#pragma omp parallel for collapse(3)
      for (int K = 0; K < nz; K++) {
      for (int J = 0; J < ny; J++) {
      for (int I = 0; I < vc; I++) {
        int m = I + vc*(J + ny*K);
        buffer[m] = interpolateC2F_X(cData, cIndex, 2*nx-1-I, J+J0, K+K0);
      }}}
      break;
    }
    case Y_M:
    {
      int K0 = nz * subfaceOrigin0(subface);
      int I0 = nx * subfaceOrigin1(subface);
//#pragma omp parallel for collapse(3)
      for (int K = 0; K < nz; K++) {
      for (int J = 0; J < vc; J++) {
      for (int I = 0; I < nx; I++) {
        int m = I + nx*(J + vc*K);
        buffer[m] = interpolateC2F_Y(cData, cIndex, I+I0, J, K+K0);
      }}}
      break;
    }
    case Y_P:
    {
      int K0 = nz * subfaceOrigin0(subface);
      int I0 = nx * subfaceOrigin1(subface);
//#pragma omp parallel for collapse(3)
      for (int K = 0; K < nz; K++) {
      for (int J = 0; J < vc; J++) {
      for (int I = 0; I < nx; I++) {
        int m = I + nx*(J + vc*K);
        buffer[m] = interpolateC2F_Y(cData, cIndex, I+I0, 2*ny-1-J, K+K0);
      }}}
      break;
    }
    case Z_M:
    {
      int I0 = nx * subfaceOrigin0(subface);
      int J0 = ny * subfaceOrigin1(subface);
//#pragma omp parallel for collapse(3)
      for (int K = 0; K < vc; K++) {
      for (int J = 0; J < ny; J++) {
      for (int I = 0; I < nx; I++) {
        int m = I + nx*(J + ny*K);
        buffer[m] = interpolateC2F_Z(cData, cIndex, I+I0, J+J0, K);
      }}}
      break;
    }
    case Z_P:
    {
      int I0 = nx * subfaceOrigin0(subface);
      int J0 = ny * subfaceOrigin1(subface);
//#pragma omp parallel for collapse(3)
      for (int K = 0; K < vc; K++) {
      for (int J = 0; J < ny; J++) {
      for (int I = 0; I < nx; I++) {
        int m = I + nx*(J + ny*K);
        buffer[m] = interpolateC2F_Z(cData, cIndex, I+I0, J+J0, 2*nz-1-K);
      }}}
      break;
    }
    default:
      break;
  }
}

void copyToCommBufferF2C_0(int nx, int ny, int nz, int vc,
                                               Face face, Subface subface,
                                               const T* fData, Index3DS fIndex,
                                               T* buffer)
{
  int ii = 0;
  switch (face) {
    case X_M:
    {
//#pragma omp parallel for collapse(3)
      for (int k = 0; k < nz/2; k++) {
      for (int j = 0; j < ny/2; j++) {
        for (int i = 0; i < min(1,vc); i++) {
          int m = i + vc*(j + ny/2*k);
          buffer[m] = interpolateF2C_2x2(fData, fIndex, i, j, k);
        }
        for (int i = 1; i < vc; i++) {
          int m = i + vc*(j + ny/2*k);
          buffer[m] = interpolateF2C_4x4(fData, fIndex, i, j, k);
        }
      }}
      break;
    }
    case X_P:
    {
//#pragma omp parallel for collapse(3)
      for (int k = 0; k < nz/2; k++) {
      for (int j = 0; j < ny/2; j++) {
        for (int i = 0; i < min(1,vc); i++) {
          int m = i + vc*(j + ny/2*k);
          buffer[m] = interpolateF2C_2x2(fData, fIndex, nx/2-1-i, j, k);
        }
        for (int i = 1; i < vc; i++) {
          int m = i + vc*(j + ny/2*k);
          buffer[m] = interpolateF2C_4x4(fData, fIndex, nx/2-1-i, j, k);
        }
      }}
      break;
    }
    case Y_M:
    {
//#pragma omp parallel for collapse(3)
      for (int k = 0; k < nz/2; k++) {
        for (int j = 0; j < min(1,vc); j++) {
        for (int i = 0; i < nx/2; i++) {
          int m = i + nx/2*(j + vc*k);
          buffer[m] = interpolateF2C_2x2(fData, fIndex, i, j, k);
        }}
        for (int j = 1; j < vc; j++) {
        for (int i = 0; i < nx/2; i++) {
          int m = i + nx/2*(j + vc*k);
          buffer[m] = interpolateF2C_4x4(fData, fIndex, i, j, k);
        }}
      }
      break;
    }
    case Y_P:
    {
//#pragma omp parallel for collapse(3)
      for (int k = 0; k < nz/2; k++) {
        for (int j = 0; j < min(1,vc); j++) {
        for (int i = 0; i < nx/2; i++) {
          int m = i + nx/2*(j + vc*k);
          buffer[m] =  interpolateF2C_2x2(fData, fIndex, i, ny/2-1-j, k);
        }}
        for (int j = 1; j < vc; j++) {
        for (int i = 0; i < nx/2; i++) {
          int m = i + nx/2*(j + vc*k);
          buffer[m] =  interpolateF2C_4x4(fData, fIndex, i, ny/2-1-j, k);
        }}
      }
      break;
    }
    case Z_M:
    {
//#pragma omp parallel for collapse(3)
      for (int k = 0; k < min(1,vc); k++) {
      for (int j = 0; j < ny/2; j++) {
      for (int i = 0; i < nx/2; i++) {
        int m = i + nx/2*(j + ny/2*k);
        buffer[m] = interpolateF2C_2x2(fData, fIndex, i, j, k);
      }}}
      for (int k = 1; k < vc; k++) {
      for (int j = 0; j < ny/2; j++) {
      for (int i = 0; i < nx/2; i++) {
        int m = i + nx/2*(j + ny/2*k);
        buffer[m] = interpolateF2C_4x4(fData, fIndex, i, j, k);
      }}}
      break;
    }
    case Z_P:
    {
//#pragma omp parallel for collapse(3)
      for (int k = 0; k < min(1,vc); k++) {
      for (int j = 0; j < ny/2; j++) {
      for (int i = 0; i < nx/2; i++) {
        int m = i + nx/2*(j + ny/2*k);
        buffer[m] = interpolateF2C_2x2(fData, fIndex, i, j, nz/2-1-k);
      }}}
      for (int k = 1; k < vc; k++) {
      for (int j = 0; j < ny/2; j++) {
      for (int i = 0; i < nx/2; i++) {
        int m = i + nx/2*(j + ny/2*k);
        buffer[m] = interpolateF2C_4x4(fData, fIndex, i, j, nz/2-1-k);
      }}}
      break;
    }
    default:
      break;
  }
}

};

#ifdef BCMT_NAMESPACE
} // namespace BCMT_NAMESPACE
#endif

#endif // SCALAR_3D_UPDATER4_H
