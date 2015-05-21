///
/// @file VtkWriter.h
/// @brief VTKフォーマットによるBCMデータ出力クラス
/// 
/// @note 
///  
///

#ifndef VTK_WRITER_H
#define VTK_WRITER_H

#include <string>
#include <sstream>
#include <fstream>
#include <cstring>
#include <sys/stat.h>
#include <mpi.h>
#include <zlib.h>

#include "BCMTools.h"
#include "BlockManager.h"
#include "Scalar3D.h"
#include "Vector3D.h"
#include "BCMOctree.h"
#include "Partition.h"

using namespace std;

#ifdef BCMT_NAMESPACE
namespace BCMT_NAMESPACE {
#endif


	class VtkWriter {
		BlockManager& blockManager;
		const MPI::Intracomm& comm;

		public:
		VtkWriter()
			: blockManager(BlockManager::getInstance()),
				comm(blockManager.getCommunicator()) {
		}

		/// デストラクタ.
		~VtkWriter() {
		}

		template <typename T>
			void printVTIC(
					T* sData, 
					const char* path,
					const char* prefix,
					const char* label,
					int step, int rank, int block,
					int NX, int NY, int NZ,
					int vc,
					double ox, double oy, double oz,
					double dx) {
				ostringstream ossFileName;
				ossFileName << path;
				ossFileName << "/";
				ossFileName.width(10);
				ossFileName.setf(ios::fixed);
				ossFileName.fill('0');
				ossFileName << step;
				ossFileName << "/";
				ossFileName << prefix;
				ossFileName << label;
/*
				ossFileName << "-";
				ossFileName.width(5);
				ossFileName.setf(ios::fixed);
				ossFileName.fill('0');
				ossFileName << rank;
*/
				ossFileName << "-";
				ossFileName.width(5);
				ossFileName.setf(ios::fixed);
				ossFileName.fill('0');
				ossFileName << block;
				ossFileName << "-";
				ossFileName.width(10);
				ossFileName.setf(ios::fixed);
				ossFileName.fill('0');
				ossFileName << step;
				ossFileName << ".vti";

				int iNX1 = 0;
				int iNY1 = 0;
				int iNZ1 = 0;
				int iNXN = NX;
				int iNYN = NY;
				int iNZN = NZ;

				unsigned int nSizeX = NX;
				unsigned int nSizeY = NY;
				unsigned int nSizeZ = NZ;
				unsigned int nSize = nSizeX*nSizeY*nSizeZ;

				ofstream ofs;
				ofs.open(ossFileName.str().c_str(), ios::out);
				ofs << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"";
#ifdef __FUJITSU
				ofs << "BigEndian";
#else
				ofs << "LittleEndian";
#endif
				ofs << "\">" << endl;
				ofs << "<ImageData WholeExtent=\"";
				ofs << iNX1 << " ";
				ofs << iNXN << " ";
				ofs << iNY1 << " ";
				ofs << iNYN << " ";
				ofs << iNZ1 << " ";
				ofs << iNZN << "\" ";
				ofs << "Origin=\"";
				ofs.setf(std::ios::scientific, std::ios::floatfield);
				ofs.precision(16);
				ofs << ox << " ";
				ofs << oy << " ";
				ofs << oz << "\" ";
				ofs << "Spacing=\"";
				ofs << dx << " ";
				ofs << dx << " ";
				ofs << dx << "\">" << std::endl;
				ofs << "<Piece Extent=\"";
				ofs << iNX1 << " ";
				ofs << iNXN << " ";
				ofs << iNY1 << " ";
				ofs << iNYN << " ";
				ofs << iNZ1 << " ";
				ofs << iNZN << "\">" << std::endl;
				ofs << "<PointData>" << endl;
				ofs << "</PointData>" << endl;
				ofs << "<CellData>" << endl;

				ofs << "<DataArray type=\"Float32\" Name=\"";
				ofs << label;
				ofs << "\" format=\"appended\" offset=\"";
				ofs << (sizeof(int) + sizeof(float)*nSize)*0;
				ofs << "\"/>" << endl;

				ofs << "</CellData>" << endl;
				ofs << "<Coordinates>" << endl;
				ofs << "</Coordinates>" << endl;
				ofs << "</Piece>" << endl;
				ofs << "</ImageData>" << endl;
				ofs << "<AppendedData encoding=\"raw\">" << endl;
				ofs << "_";
				ofs.close();

				ofs.open(ossFileName.str().c_str(), ios::out | ios::app | ios::binary);

				float* pScalar = new float[nSize];
				int nBytes = sizeof(float)*nSize;

				for(int k=0; k<NZ; k++) {
					for(int j=0; j<NY; j++) {
						for(int i=0; i<NX; i++) {
							int i0 = i + vc;
							int j0 = j + vc;
							int k0 = k + vc;
							int m0 = i0 + (NX+2*vc)*( j0 + (NY+2*vc)*k0 );
							float data = sData[m0];
							int m = i + NX*( j + NY*k );
							pScalar[m] = data;
						}
					}
				}
				ofs.write((char*)&nBytes, sizeof(int));
				ofs.write((char*)pScalar, sizeof(float)*(nSize));

				ofs.close();

				ofs.open(ossFileName.str().c_str(), ios::out | ios::app);
				ofs << endl;
				ofs << "</AppendedData>" << endl;
				ofs << "</VTKFile>" << endl;
				ofs.close();

				delete [] pScalar;
			}

		template <typename T>
			void printVTIC(
					T* sDataP, 
					T* sDataUX, 
					T* sDataUY, 
					T* sDataUZ, 
					const char* path,
					const char* prefix,
					const char* label,
					int step, int rank, int block,
					int NX, int NY, int NZ,
					int vc,
					double ox, double oy, double oz,
					double dx) {
				ostringstream ossFileName;
				ossFileName << path;
				ossFileName << "/";
				ossFileName.width(10);
				ossFileName.setf(ios::fixed);
				ossFileName.fill('0');
				ossFileName << step;
				ossFileName << "/";
				ossFileName << prefix;
				ossFileName << label;
/*
				ossFileName << "-";
				ossFileName.width(5);
				ossFileName.setf(ios::fixed);
				ossFileName.fill('0');
				ossFileName << rank;
*/
				ossFileName << "-";
				ossFileName.width(5);
				ossFileName.setf(ios::fixed);
				ossFileName.fill('0');
				ossFileName << block;
				ossFileName << "-";
				ossFileName.width(10);
				ossFileName.setf(ios::fixed);
				ossFileName.fill('0');
				ossFileName << step;
				ossFileName << ".vti";

				int iNX1 = 0;
				int iNY1 = 0;
				int iNZ1 = 0;
				int iNXN = NX;
				int iNYN = NY;
				int iNZN = NZ;

				unsigned int nSizeX = NX;
				unsigned int nSizeY = NY;
				unsigned int nSizeZ = NZ;
				unsigned int nSize = nSizeX*nSizeY*nSizeZ;

				ofstream ofs;
				ofs.open(ossFileName.str().c_str(), ios::out);
				ofs << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"";
#ifdef __FUJITSU
				ofs << "BigEndian";
#else
				ofs << "LittleEndian";
#endif
				ofs << "\">" << endl;
				ofs << "<ImageData WholeExtent=\"";
				ofs << iNX1 << " ";
				ofs << iNXN << " ";
				ofs << iNY1 << " ";
				ofs << iNYN << " ";
				ofs << iNZ1 << " ";
				ofs << iNZN << "\" ";
				ofs << "Origin=\"";
				ofs.setf(std::ios::scientific, std::ios::floatfield);
				ofs.precision(16);
				ofs << ox << " ";
				ofs << oy << " ";
				ofs << oz << "\" ";
				ofs << "Spacing=\"";
				ofs << dx << " ";
				ofs << dx << " ";
				ofs << dx << "\">" << std::endl;
				ofs << "<Piece Extent=\"";
				ofs << iNX1 << " ";
				ofs << iNXN << " ";
				ofs << iNY1 << " ";
				ofs << iNYN << " ";
				ofs << iNZ1 << " ";
				ofs << iNZN << "\">" << std::endl;
				ofs << "<PointData>" << endl;
				ofs << "</PointData>" << endl;
				ofs << "<CellData>" << endl;

				ofs << "<DataArray type=\"Float32\" Name=\"";
				ofs << "p";
				ofs << "\" format=\"appended\" offset=\"";
				ofs << (sizeof(int) + sizeof(float)*nSize)*0;
				ofs << "\"/>" << endl;

				ofs << "<DataArray type=\"Float32\" Name=\"";
				ofs << "u";
				ofs << "\" NumberOfComponents=\"3\" format=\"appended\" offset=\"";
				ofs << (sizeof(int) + sizeof(float)*nSize)*1;
				ofs << "\"/>" << endl;

				ofs << "</CellData>" << endl;
				ofs << "<Coordinates>" << endl;
				ofs << "</Coordinates>" << endl;
				ofs << "</Piece>" << endl;
				ofs << "</ImageData>" << endl;
				ofs << "<AppendedData encoding=\"raw\">" << endl;
				ofs << "_";
				ofs.close();

				ofs.open(ossFileName.str().c_str(), ios::out | ios::app | ios::binary);

				float* pScalar = new float[nSize];
				float* pVector = new float[nSize*3];
				int nBytesS = sizeof(float)*nSize;
				int nBytesV = sizeof(float)*nSize*3;

				for(int k=0; k<NZ; k++) {
					for(int j=0; j<NY; j++) {
						for(int i=0; i<NX; i++) {
							int i0 = i + vc;
							int j0 = j + vc;
							int k0 = k + vc;
							int m0 = i0 + (NX+2*vc)*( j0 + (NY+2*vc)*k0 );
							float data = sDataP[m0];
							int m = i + NX*( j + NY*k );
							pScalar[m] = data;
						}
					}
				}
				ofs.write((char*)&nBytesS, sizeof(int));
				ofs.write((char*)pScalar, sizeof(float)*(nSize));

				for(int k=0; k<NZ; k++) {
					for(int j=0; j<NY; j++) {
						for(int i=0; i<NX; i++) {
							int i0 = i + vc;
							int j0 = j + vc;
							int k0 = k + vc;
							int m0 = i0 + (NX+2*vc)*( j0 + (NY+2*vc)*k0 );
							float data0 = sDataUX[m0];
							float data1 = sDataUY[m0];
							float data2 = sDataUZ[m0];
							int m = i + NX*( j + NY*k );
							pVector[3*m + 0] = data0;
							pVector[3*m + 1] = data1;
							pVector[3*m + 2] = data2;
						}
					}
				}
				ofs.write((char*)&nBytesV, sizeof(int));
				ofs.write((char*)pVector, sizeof(float)*(nSize*3));

				ofs.close();

				ofs.open(ossFileName.str().c_str(), ios::out | ios::app);
				ofs << endl;
				ofs << "</AppendedData>" << endl;
				ofs << "</VTKFile>" << endl;
				ofs.close();

				delete [] pScalar;
				delete [] pVector;
			}

		template <typename T>
			void printVTIC(
					T* sDataP, 
					T* sDataUX, 
					T* sDataUY, 
					T* sDataUZ, 
					T* sDataT, 
					const char* path,
					const char* prefix,
					const char* label,
					int step, int rank, int block,
					int NX, int NY, int NZ,
					int vc,
					double ox, double oy, double oz,
					double dx) {
				ostringstream ossFileName;
				ossFileName << path;
				ossFileName << "/";
				ossFileName.width(10);
				ossFileName.setf(ios::fixed);
				ossFileName.fill('0');
				ossFileName << step;
				ossFileName << "/";
				ossFileName << prefix;
				ossFileName << label;
/*
				ossFileName << "-";
				ossFileName.width(5);
				ossFileName.setf(ios::fixed);
				ossFileName.fill('0');
				ossFileName << rank;
*/
				ossFileName << "-";
				ossFileName.width(5);
				ossFileName.setf(ios::fixed);
				ossFileName.fill('0');
				ossFileName << block;
				ossFileName << "-";
				ossFileName.width(10);
				ossFileName.setf(ios::fixed);
				ossFileName.fill('0');
				ossFileName << step;
				ossFileName << ".vti";

				int iNX1 = 0;
				int iNY1 = 0;
				int iNZ1 = 0;
				int iNXN = NX;
				int iNYN = NY;
				int iNZN = NZ;

				unsigned int nSizeX = NX;
				unsigned int nSizeY = NY;
				unsigned int nSizeZ = NZ;
				unsigned int nSize = nSizeX*nSizeY*nSizeZ;

				ofstream ofs;
				ofs.open(ossFileName.str().c_str(), ios::out);
				ofs << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"";
#ifdef __FUJITSU
				ofs << "BigEndian";
#else
				ofs << "LittleEndian";
#endif
				ofs << "\">" << endl;
				ofs << "<ImageData WholeExtent=\"";
				ofs << iNX1 << " ";
				ofs << iNXN << " ";
				ofs << iNY1 << " ";
				ofs << iNYN << " ";
				ofs << iNZ1 << " ";
				ofs << iNZN << "\" ";
				ofs << "Origin=\"";
				ofs.setf(std::ios::scientific, std::ios::floatfield);
				ofs.precision(16);
				ofs << ox << " ";
				ofs << oy << " ";
				ofs << oz << "\" ";
				ofs << "Spacing=\"";
				ofs << dx << " ";
				ofs << dx << " ";
				ofs << dx << "\">" << std::endl;
				ofs << "<Piece Extent=\"";
				ofs << iNX1 << " ";
				ofs << iNXN << " ";
				ofs << iNY1 << " ";
				ofs << iNYN << " ";
				ofs << iNZ1 << " ";
				ofs << iNZN << "\">" << std::endl;
				ofs << "<PointData>" << endl;
				ofs << "</PointData>" << endl;
				ofs << "<CellData>" << endl;

				ofs << "<DataArray type=\"Float32\" Name=\"";
				ofs << "p";
				ofs << "\" format=\"appended\" offset=\"";
				ofs << sizeof(int)*0 + sizeof(float)*nSize*0;
				ofs << "\"/>" << endl;

				ofs << "<DataArray type=\"Float32\" Name=\"";
				ofs << "u";
				ofs << "\" NumberOfComponents=\"3\" format=\"appended\" offset=\"";
				ofs << sizeof(int)*1 + sizeof(float)*nSize*1;
				ofs << "\"/>" << endl;

				ofs << "<DataArray type=\"Float32\" Name=\"";
				ofs << "t";
				ofs << "\" format=\"appended\" offset=\"";
				ofs << sizeof(int)*2 + sizeof(float)*nSize*4;
				ofs << "\"/>" << endl;

				ofs << "</CellData>" << endl;
				ofs << "<Coordinates>" << endl;
				ofs << "</Coordinates>" << endl;
				ofs << "</Piece>" << endl;
				ofs << "</ImageData>" << endl;
				ofs << "<AppendedData encoding=\"raw\">" << endl;
				ofs << "_";
				ofs.close();

				ofs.open(ossFileName.str().c_str(), ios::out | ios::app | ios::binary);

				float* pScalar = new float[nSize];
				float* pVector = new float[nSize*3];
				int nBytesS = sizeof(float)*nSize;
				int nBytesV = sizeof(float)*nSize*3;

				for(int k=0; k<NZ; k++) {
					for(int j=0; j<NY; j++) {
						for(int i=0; i<NX; i++) {
							int i0 = i + vc;
							int j0 = j + vc;
							int k0 = k + vc;
							int m0 = i0 + (NX+2*vc)*( j0 + (NY+2*vc)*k0 );
							float data = sDataP[m0];
							int m = i + NX*( j + NY*k );
							pScalar[m] = data;
						}
					}
				}
				ofs.write((char*)&nBytesS, sizeof(int));
				ofs.write((char*)pScalar, sizeof(float)*(nSize));

				for(int k=0; k<NZ; k++) {
					for(int j=0; j<NY; j++) {
						for(int i=0; i<NX; i++) {
							int i0 = i + vc;
							int j0 = j + vc;
							int k0 = k + vc;
							int m0 = i0 + (NX+2*vc)*( j0 + (NY+2*vc)*k0 );
							float data0 = sDataUX[m0];
							float data1 = sDataUY[m0];
							float data2 = sDataUZ[m0];
							int m = i + NX*( j + NY*k );
							pVector[3*m + 0] = data0;
							pVector[3*m + 1] = data1;
							pVector[3*m + 2] = data2;
						}
					}
				}
				ofs.write((char*)&nBytesV, sizeof(int));
				ofs.write((char*)pVector, sizeof(float)*(nSize*3));

				for(int k=0; k<NZ; k++) {
					for(int j=0; j<NY; j++) {
						for(int i=0; i<NX; i++) {
							int i0 = i + vc;
							int j0 = j + vc;
							int k0 = k + vc;
							int m0 = i0 + (NX+2*vc)*( j0 + (NY+2*vc)*k0 );
							float data = sDataT[m0];
							int m = i + NX*( j + NY*k );
							pScalar[m] = data;
						}
					}
				}
				ofs.write((char*)&nBytesS, sizeof(int));
				ofs.write((char*)pScalar, sizeof(float)*(nSize));

				ofs.close();

				ofs.open(ossFileName.str().c_str(), ios::out | ios::app);
				ofs << endl;
				ofs << "</AppendedData>" << endl;
				ofs << "</VTKFile>" << endl;
				ofs.close();

				delete [] pScalar;
				delete [] pVector;
			}

		unsigned long GetMaxCompressionSpace(unsigned long size) {
			return size + (size + 999)/1000 + 12;
		}

		template <typename T>
			void printVTIC_Z(
					T* sDataP, 
					T* sDataUX, 
					T* sDataUY, 
					T* sDataUZ, 
					T* sDataT, 
					const char* path,
					const char* prefix,
					const char* label,
					int step, int rank, int block,
					int NX, int NY, int NZ,
					int vc,
					double ox, double oy, double oz,
					double dx) {
				ostringstream ossFileName;
				ossFileName << path;
				ossFileName << "/";
				ossFileName.width(10);
				ossFileName.setf(ios::fixed);
				ossFileName.fill('0');
				ossFileName << step;
				ossFileName << "/";
				ossFileName << prefix;
				ossFileName << label;
/*
				ossFileName << "-";
				ossFileName.width(5);
				ossFileName.setf(ios::fixed);
				ossFileName.fill('0');
				ossFileName << rank;
*/
				ossFileName << "-";
				ossFileName.width(5);
				ossFileName.setf(ios::fixed);
				ossFileName.fill('0');
				ossFileName << block;
				ossFileName << "-";
				ossFileName.width(10);
				ossFileName.setf(ios::fixed);
				ossFileName.fill('0');
				ossFileName << step;
				ossFileName << ".vti";

				int iNX1 = 0;
				int iNY1 = 0;
				int iNZ1 = 0;
				int iNXN = NX;
				int iNYN = NY;
				int iNZN = NZ;

				unsigned int nSizeX = NX;
				unsigned int nSizeY = NY;
				unsigned int nSizeZ = NZ;
				unsigned int nSize = nSizeX*nSizeY*nSizeZ;

				float* pP = new float[nSize];
				int nBytesP = sizeof(float)*nSize;
				for(int k=0; k<NZ; k++) {
					for(int j=0; j<NY; j++) {
						for(int i=0; i<NX; i++) {
							int i0 = i + vc;
							int j0 = j + vc;
							int k0 = k + vc;
							int m0 = i0 + (NX+2*vc)*( j0 + (NY+2*vc)*k0 );
							int m = i + NX*( j + NY*k );
							pP[m] = (float)sDataP[m0];
						}
					}
				}
				const unsigned char* pP_u = (const unsigned char*)pP;
				unsigned long        nP_u = nSize*sizeof(float);
				unsigned long        nP_c = GetMaxCompressionSpace(nP_u);
				unsigned char*       pP_c = new unsigned char [nP_c];
				if( compress2(pP_c, &nP_c, pP_u, nP_u, Z_DEFAULT_COMPRESSION) != Z_OK ) {
					std::cout << "Error: zlib" << std::endl;
				}
				delete [] pP;

				float* pU = new float[nSize*3];
				int nBytesU = sizeof(float)*nSize*3;
				for(int k=0; k<NZ; k++) {
					for(int j=0; j<NY; j++) {
						for(int i=0; i<NX; i++) {
							int i0 = i + vc;
							int j0 = j + vc;
							int k0 = k + vc;
							int m0 = i0 + (NX+2*vc)*( j0 + (NY+2*vc)*k0 );
							int m = i + NX*( j + NY*k );
							pU[3*m + 0] = (float)sDataUX[m0];
							pU[3*m + 1] = (float)sDataUY[m0];
							pU[3*m + 2] = (float)sDataUZ[m0];
						}
					}
				}
				const unsigned char* pU_u = (const unsigned char*)pU;
				unsigned long        nU_u = 3*nSize*sizeof(float);
				unsigned long        nU_c = GetMaxCompressionSpace(nU_u);
				unsigned char*       pU_c = new unsigned char [nU_c];
				if( compress2(pU_c, &nU_c, pU_u, nU_u, Z_DEFAULT_COMPRESSION) != Z_OK ) {
					std::cout << "Error: zlib" << std::endl;
				}
				delete [] pU;

				float* pT = new float[nSize];
				int nBytesT = sizeof(float)*nSize;
				for(int k=0; k<NZ; k++) {
					for(int j=0; j<NY; j++) {
						for(int i=0; i<NX; i++) {
							int i0 = i + vc;
							int j0 = j + vc;
							int k0 = k + vc;
							int m0 = i0 + (NX+2*vc)*( j0 + (NY+2*vc)*k0 );
							int m = i + NX*( j + NY*k );
							pT[m] = (float)sDataT[m0];
						}
					}
				}
				const unsigned char* pT_u = (const unsigned char*)pT;
				unsigned long        nT_u = nSize*sizeof(float);
				unsigned long        nT_c = GetMaxCompressionSpace(nT_u);
				unsigned char*       pT_c = new unsigned char [nT_c];
				if( compress2(pT_c, &nT_c, pT_u, nT_u, Z_DEFAULT_COMPRESSION) != Z_OK ) {
					std::cout << "Error: zlib" << std::endl;
				}
				delete [] pT;

				ofstream ofs;
				ofs.open(ossFileName.str().c_str(), ios::out);
				ofs << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"";
#ifdef __FUJITSU
				ofs << "BigEndian";
#else
				ofs << "LittleEndian";
#endif
				ofs << "\" compressor=\"vtkZLibDataCompressor\">" << endl;
				//		ofs << "\">" << endl;
				ofs << "<ImageData WholeExtent=\"";
				ofs << iNX1 << " ";
				ofs << iNXN << " ";
				ofs << iNY1 << " ";
				ofs << iNYN << " ";
				ofs << iNZ1 << " ";
				ofs << iNZN << "\" ";
				ofs << "Origin=\"";
				ofs.setf(std::ios::scientific, std::ios::floatfield);
				ofs.precision(16);
				ofs << ox << " ";
				ofs << oy << " ";
				ofs << oz << "\" ";
				ofs << "Spacing=\"";
				ofs << dx << " ";
				ofs << dx << " ";
				ofs << dx << "\">" << std::endl;
				ofs << "<Piece Extent=\"";
				ofs << iNX1 << " ";
				ofs << iNXN << " ";
				ofs << iNY1 << " ";
				ofs << iNYN << " ";
				ofs << iNZ1 << " ";
				ofs << iNZN << "\">" << std::endl;
				ofs << "<PointData>" << endl;
				ofs << "</PointData>" << endl;
				ofs << "<CellData>" << endl;

				ofs << "<DataArray type=\"Float32\" Name=\"";
				ofs << "p";
				ofs << "\" format=\"appended\" offset=\"";
				ofs << sizeof(int)*0 + 0;
				ofs << "\"/>" << endl;

				ofs << "<DataArray type=\"Float32\" Name=\"";
				ofs << "u";
				ofs << "\" NumberOfComponents=\"3\" format=\"appended\" offset=\"";
				ofs << sizeof(int)*4 + nP_c;
				ofs << "\"/>" << endl;

				ofs << "<DataArray type=\"Float32\" Name=\"";
				ofs << "t";
				ofs << "\" format=\"appended\" offset=\"";
				ofs << sizeof(int)*8 + nP_c + nU_c;
				ofs << "\"/>" << endl;

				ofs << "</CellData>" << endl;
				ofs << "<Coordinates>" << endl;
				ofs << "</Coordinates>" << endl;
				ofs << "</Piece>" << endl;
				ofs << "</ImageData>" << endl;
				ofs << "<AppendedData encoding=\"raw\">" << endl;
				ofs << "_";
				ofs.close();

				int nB = 1;
				int psize = 0;
				ofs.open(ossFileName.str().c_str(), ios::out | ios::app | ios::binary);

				ofs.write((const char*)&nB, sizeof(int));
				ofs.write((const char*)&nP_u, sizeof(int));
				ofs.write((const char*)&psize, sizeof(int));
				ofs.write((const char*)&nP_c, sizeof(int));
				ofs.write((const char*)pP_c, nP_c);

				ofs.write((const char*)&nB, sizeof(int));
				ofs.write((const char*)&nU_u, sizeof(int));
				ofs.write((const char*)&psize, sizeof(int));
				ofs.write((const char*)&nU_c, sizeof(int));
				ofs.write((const char*)pU_c, nU_c);

				ofs.write((const char*)&nB, sizeof(int));
				ofs.write((const char*)&nT_u, sizeof(int));
				ofs.write((const char*)&psize, sizeof(int));
				ofs.write((const char*)&nT_c, sizeof(int));
				ofs.write((const char*)pT_c, nT_c);

				ofs.close();

				ofs.open(ossFileName.str().c_str(), ios::out | ios::app);
				ofs << endl;
				ofs << "</AppendedData>" << endl;
				ofs << "</VTKFile>" << endl;
				ofs.close();

				delete [] pP_c;
				delete [] pU_c;
				delete [] pT_c;
			}

		template <typename T>
			void writeScalar_OverlappingAMR(
					int dataClassID_P,
					int vc,
					const std::string path,
					const std::string prefix,
					const std::string name,
					int step,
					int maxLevel,
					int minLevel,
					RootGrid* rootGrid,
					BCMOctree* tree,
					Partition* partition,
					Vec3r rootOrigin,
					double rootLength) {
				ostringstream ossFileNameTime;
				ossFileNameTime << path;
				mkdir(ossFileNameTime.str().c_str(), 0755);
				ossFileNameTime << "/";
				ossFileNameTime.width(10);
				ossFileNameTime.setf(ios::fixed);
				ossFileNameTime.fill('0');
				ossFileNameTime << step;
				mkdir(ossFileNameTime.str().c_str(), 0755);

				const Vec3i& size = blockManager.getSize();
				int myrank = comm.Get_rank();

				float* dataP  = new float[(size.x) * (size.y) * (size.z)];

				for (int id = 0; id < blockManager.getNumBlock(); ++id) {
					BlockBase* block = blockManager.getBlock(id);
					Vec3i size = block->getSize();
					Vec3r origin = block->getOrigin();
					Vec3r blockSize = block->getBlockSize();
					Vec3r cellSize = block->getCellSize();
					int level = block->getLevel();
					int blockid = id + blockManager.getStartID();

					Scalar3D<T>* sp = dynamic_cast<Scalar3D<T>*>(block->getDataClass(dataClassID_P));
					T* sDataP = sp->getData();

					//			printVTIC(sDataP, path.c_str(), prefix.c_str(), name.c_str(), step, myrank, blockid, size[0], size[1], size[2], vc, origin[0], origin[1], origin[2], cellSize[0]);
				}

				delete[] dataP;

				ostringstream ossFileName;
				ossFileName << path;
				ossFileName << "/";
				ossFileName << prefix;
				ossFileName << name.c_str();
				ossFileName << "-";
				ossFileName.width(10);
				ossFileName.setf(ios::fixed);
				ossFileName.fill('0');
				ossFileName << step;
				ossFileName << ".vth";

				if( myrank == 0 ) {
					ofstream ofs;
					ofs.open(ossFileName.str().c_str(), ios::out);
					ofs << "<VTKFile type=\"vtkOverlappingAMR\" version=\"1.1\">" << endl;
					ofs << "<vtkOverlappingAMR origin=\"0 0 0\" grid_description=\"XYZ\">" << endl;

					std::vector<Node*>& leafNodeArray = tree->getLeafNodeArray();
					//			int lmax = difflevel+1;
					int lmax = maxLevel+1;
					for(int n=0; n<lmax; n++) {
						double dx = rootLength/(1 << n);
						ofs << "\t<Block level=\"";
						ofs << n;
						ofs << "\" spacing=\"";
						ofs << dx;
						ofs << " ";
						ofs << dx;
						ofs << " ";
						ofs << dx;
						ofs << "\">";
						ofs << endl;
						for (int iRank = 0; iRank < comm.Get_size(); iRank++) {
							for (int id = partition->getStart(iRank); id < partition->getEnd(iRank); id++) {
								Node* node = leafNodeArray[id];
								Vec3r origin = tree->getOrigin(node) * rootLength;
								Vec3r blockSize = node->getBlockSize() * rootLength;
								Vec3r cellSize;
								cellSize.x = blockSize.x / size.x;
								cellSize.y = blockSize.y / size.y;
								cellSize.z = blockSize.z / size.z;
								int level = node->getLevel();

								ostringstream ossFileName2;
								ossFileName2 << "./";
								ossFileName2.width(10);
								ossFileName2.setf(ios::fixed);
								ossFileName2.fill('0');
								ossFileName2 << step;
								ossFileName2 << "/";
								ossFileName2 << prefix;
								ossFileName2 << name.c_str();
/*
								ossFileName2 << "-";
								ossFileName2.width(5);
								ossFileName2.setf(ios::fixed);
								ossFileName2.fill('0');
								ossFileName2 << iRank;
*/
								ossFileName2 << "-";
								ossFileName2.width(5);
								ossFileName2.setf(ios::fixed);
								ossFileName2.fill('0');
								ossFileName2 << id;
								ossFileName2 << "-";
								ossFileName2.width(10);
								ossFileName2.setf(ios::fixed);
								ossFileName2.fill('0');
								ossFileName2 << step;
								ossFileName2 << ".vti";

								int lx = size.x*(1 << level);
								int ly = size.y*(1 << level);
								int lz = size.z*(1 << level);

								Vec3r origin2 = tree->getOrigin(node);

								double nx0 = (origin2.x)*lx;
								double ny0 = (origin2.y)*ly;
								double nz0 = (origin2.z)*lz;

								double nx1 = nx0 + size.x - 1;
								double ny1 = ny0 + size.y - 1;
								double nz1 = nz0 + size.z - 1;

								if( level == n ) {
									ofs << "\t\t<DataSet index=\"";
									ofs << id;
									ofs << "\" amr_box=\"";
									ofs << nx0;
									ofs << " ";
									ofs << nx1;
									ofs << " ";
									ofs << ny0;
									ofs << " ";
									ofs << ny1;
									ofs << " ";
									ofs << nz0;
									ofs << " ";
									ofs << nz1;
									ofs << "\" file=\"";
									ofs << ossFileName2.str().c_str();
									ofs << "\">";
									ofs << endl;
									ofs << "\t\t</DataSet>" << endl;	
								}
							}
						}
						ofs << "\t</Block>";
						ofs << endl;
						ofs << endl;
					}

					ofs << "</vtkOverlappingAMR>" << endl;
					ofs << "</VTKFile>" << endl;
					ofs.close();
				}
			}

		template <typename T>
			void writeScalar(
					int dataClassID_P,
					int vc,
					const std::string path,
					const std::string prefix,
					const std::string name,
					int step,
					int difflevel,
					RootGrid* rootGrid,
					BCMOctree* tree,
					Partition* partition,
					Vec3r rootOrigin,
					double rootLength) {
				ostringstream ossFileNameTime;
				ossFileNameTime << path;
				mkdir(ossFileNameTime.str().c_str(), 0755);
				ossFileNameTime << "/";
				ossFileNameTime.width(10);
				ossFileNameTime.setf(ios::fixed);
				ossFileNameTime.fill('0');
				ossFileNameTime << step;
				mkdir(ossFileNameTime.str().c_str(), 0755);

				const Vec3i& size = blockManager.getSize();
				int myrank = comm.Get_rank();

				float* dataP  = new float[(size.x) * (size.y) * (size.z)];

				for (int id = 0; id < blockManager.getNumBlock(); ++id) {
					BlockBase* block = blockManager.getBlock(id);
					Vec3i size = block->getSize();
					Vec3r origin = block->getOrigin();
					Vec3r blockSize = block->getBlockSize();
					Vec3r cellSize = block->getCellSize();
					int level = block->getLevel();
					int blockid = id + blockManager.getStartID();

					Scalar3D<T>* sp = dynamic_cast<Scalar3D<T>*>(block->getDataClass(dataClassID_P));
					T* sDataP = sp->getData();

					printVTIC(sDataP, path.c_str(), prefix.c_str(), name.c_str(), step, myrank, blockid, size[0], size[1], size[2], vc, origin[0], origin[1], origin[2], cellSize[0]);
				}

				delete[] dataP;

				ostringstream ossFileName;
				ossFileName << path;
				ossFileName << "/";
				ossFileName << prefix;
				ossFileName << name.c_str();
				ossFileName << "-";
				ossFileName.width(10);
				ossFileName.setf(ios::fixed);
				ossFileName.fill('0');
				ossFileName << step;
				ossFileName << ".vthb";

				if( myrank == 0 ) {
					ofstream ofs;
					ofs.open(ossFileName.str().c_str(), ios::out);
					ofs << "<VTKFile type=\"vtkHierarchicalBoxDataSet\" version=\"0.1\">" << endl;
					ofs << "<vtkHierarchicalBoxDataSet>" << endl;

					for(int n=0; n<difflevel+1; n++) {
						ofs << "<RefinementRatio level=\"";
						ofs << n;
						ofs << "\" refinement=\"2\"/>";
						ofs << endl;
					}
					ofs << endl;

					std::vector<Node*>& leafNodeArray = tree->getLeafNodeArray();
					for (int iRank = 0; iRank < comm.Get_size(); iRank++) {
						for (int id = partition->getStart(iRank); id < partition->getEnd(iRank); id++) {
							Node* node = leafNodeArray[id];
							Vec3r origin = tree->getOrigin(node) * rootLength;
							Vec3r blockSize = node->getBlockSize() * rootLength;
							Vec3r cellSize;
							cellSize.x = blockSize.x / size.x;
							cellSize.y = blockSize.y / size.y;
							cellSize.z = blockSize.z / size.z;
							int level = node->getLevel();

							ostringstream ossFileName2;
							ossFileName2 << "./";
							ossFileName2.width(10);
							ossFileName2.setf(ios::fixed);
							ossFileName2.fill('0');
							ossFileName2 << step;
							ossFileName2 << "/";
							ossFileName2 << prefix;
							ossFileName2 << name.c_str();
/*
							ossFileName2 << "-";
							ossFileName2.width(5);
							ossFileName2.setf(ios::fixed);
							ossFileName2.fill('0');
							ossFileName2 << iRank;
*/
							ossFileName2 << "-";
							ossFileName2.width(5);
							ossFileName2.setf(ios::fixed);
							ossFileName2.fill('0');
							ossFileName2 << id;
							ossFileName2 << "-";
							ossFileName2.width(10);
							ossFileName2.setf(ios::fixed);
							ossFileName2.fill('0');
							ossFileName2 << step;
							ossFileName2 << ".vti";

							int lx = size.x*(1 << level);
							int ly = size.y*(1 << level);
							int lz = size.z*(1 << level);

							Vec3r origin2 = tree->getOrigin(node);

							double nx0 = (origin2.x)*lx;
							double ny0 = (origin2.y)*ly;
							double nz0 = (origin2.z)*lz;

							double nx1 = nx0 + size.x - 1;
							double ny1 = ny0 + size.y - 1;
							double nz1 = nz0 + size.z - 1;

							ofs << "<DataSet group=\"";
							ofs << level;
							ofs << "\" dataset=\"";
							ofs << id;
							ofs << "\" amr_box=\"";
							ofs << nx0;
							ofs << " ";
							ofs << nx1;
							ofs << " ";
							ofs << ny0;
							ofs << " ";
							ofs << ny1;
							ofs << " ";
							ofs << nz0;
							ofs << " ";
							ofs << nz1;
							ofs << "\" file=\"";
							ofs << ossFileName2.str().c_str();
							ofs << "\"/>";
							ofs << endl;
						}
					}

					ofs << "</vtkHierarchicalBoxDataSet>" << endl;
					ofs << "</VTKFile>" << endl;
					ofs.close();
				}
			}

		template <typename T>
			void writePU(
					int dataClassID_P,
					int dataClassID_UX, int dataClassID_UY, int dataClassID_UZ,
					int vc,
					const std::string path,
					const std::string prefix,
					const std::string name,
					int step,
					int difflevel,
					RootGrid* rootGrid,
					BCMOctree* tree,
					Partition* partition,
					Vec3r rootOrigin,
					double rootLength) {
				ostringstream ossFileNameTime;
				ossFileNameTime << path;
				mkdir(ossFileNameTime.str().c_str(), 0755);
				ossFileNameTime << "/";
				ossFileNameTime.width(10);
				ossFileNameTime.setf(ios::fixed);
				ossFileNameTime.fill('0');
				ossFileNameTime << step;
				mkdir(ossFileNameTime.str().c_str(), 0755);

				const Vec3i& size = blockManager.getSize();
				int myrank = comm.Get_rank();

				float* dataP  = new float[(size.x) * (size.y) * (size.z)];
				float* dataUX = new float[(size.x) * (size.y) * (size.z)];
				float* dataUY = new float[(size.x) * (size.y) * (size.z)];
				float* dataUZ = new float[(size.x) * (size.y) * (size.z)];

				for (int id = 0; id < blockManager.getNumBlock(); ++id) {
					BlockBase* block = blockManager.getBlock(id);
					Vec3i size = block->getSize();
					Vec3r origin = block->getOrigin();
					Vec3r blockSize = block->getBlockSize();
					Vec3r cellSize = block->getCellSize();
					int level = block->getLevel();
					int blockid = id + blockManager.getStartID();

					Scalar3D<T>* sp = dynamic_cast<Scalar3D<T>*>(block->getDataClass(dataClassID_P));
					T* sDataP = sp->getData();
					Scalar3D<T>* sux = dynamic_cast<Scalar3D<T>*>(block->getDataClass(dataClassID_UX));
					T* sDataUX = sux->getData();
					Scalar3D<T>* suy = dynamic_cast<Scalar3D<T>*>(block->getDataClass(dataClassID_UY));
					T* sDataUY = suy->getData();
					Scalar3D<T>* suz = dynamic_cast<Scalar3D<T>*>(block->getDataClass(dataClassID_UZ));
					T* sDataUZ = suz->getData();

					printVTIC(sDataP, sDataUX, sDataUY, sDataUZ, path.c_str(), prefix.c_str(), name.c_str(), step, myrank, blockid, size[0], size[1], size[2], vc, origin[0], origin[1], origin[2], cellSize[0]);
				}

				delete[] dataP;
				delete[] dataUX;
				delete[] dataUY;
				delete[] dataUZ;

				ostringstream ossFileName;
				ossFileName << path;
				ossFileName << "/";
				ossFileName << prefix;
				ossFileName << name.c_str();
				ossFileName << "-";
				ossFileName.width(10);
				ossFileName.setf(ios::fixed);
				ossFileName.fill('0');
				ossFileName << step;
				ossFileName << ".vthb";

				if( myrank == 0 ) {
					ofstream ofs;
					ofs.open(ossFileName.str().c_str(), ios::out);
					ofs << "<VTKFile type=\"vtkHierarchicalBoxDataSet\" version=\"0.1\">" << endl;
					ofs << "<vtkHierarchicalBoxDataSet>" << endl;

					for(int n=0; n<difflevel+1; n++) {
						ofs << "<RefinementRatio level=\"";
						ofs << n;
						ofs << "\" refinement=\"2\"/>";
						ofs << endl;
					}
					ofs << endl;

					std::vector<Node*>& leafNodeArray = tree->getLeafNodeArray();
					for (int iRank = 0; iRank < comm.Get_size(); iRank++) {
						for (int id = partition->getStart(iRank); id < partition->getEnd(iRank); id++) {
							Node* node = leafNodeArray[id];
							Vec3r origin = tree->getOrigin(node) * rootLength;
							Vec3r blockSize = node->getBlockSize() * rootLength;
							Vec3r cellSize;
							cellSize.x = blockSize.x / size.x;
							cellSize.y = blockSize.y / size.y;
							cellSize.z = blockSize.z / size.z;
							int level = node->getLevel();

							ostringstream ossFileName2;
							ossFileName2 << "./";
							ossFileName2.width(10);
							ossFileName2.setf(ios::fixed);
							ossFileName2.fill('0');
							ossFileName2 << step;
							ossFileName2 << "/";
							ossFileName2 << prefix;
							ossFileName2 << name.c_str();
/*
							ossFileName2 << "-";
							ossFileName2.width(5);
							ossFileName2.setf(ios::fixed);
							ossFileName2.fill('0');
							ossFileName2 << iRank;
*/
							ossFileName2 << "-";
							ossFileName2.width(5);
							ossFileName2.setf(ios::fixed);
							ossFileName2.fill('0');
							ossFileName2 << id;
							ossFileName2 << "-";
							ossFileName2.width(10);
							ossFileName2.setf(ios::fixed);
							ossFileName2.fill('0');
							ossFileName2 << step;
							ossFileName2 << ".vti";

							int lx = size.x*(1 << level);
							int ly = size.y*(1 << level);
							int lz = size.z*(1 << level);

							Vec3r origin2 = tree->getOrigin(node);

							double nx0 = (origin2.x)*lx;
							double ny0 = (origin2.y)*ly;
							double nz0 = (origin2.z)*lz;

							double nx1 = nx0 + size.x - 1;
							double ny1 = ny0 + size.y - 1;
							double nz1 = nz0 + size.z - 1;

							ofs << "<DataSet group=\"";
							ofs << level;
							ofs << "\" dataset=\"";
							ofs << id;
							ofs << "\" amr_box=\"";
							ofs << nx0;
							ofs << " ";
							ofs << nx1;
							ofs << " ";
							ofs << ny0;
							ofs << " ";
							ofs << ny1;
							ofs << " ";
							ofs << nz0;
							ofs << " ";
							ofs << nz1;
							ofs << "\" file=\"";
							ofs << ossFileName2.str().c_str();
							ofs << "\"/>";
							ofs << endl;
						}
					}

					ofs << "</vtkHierarchicalBoxDataSet>" << endl;
					ofs << "</VTKFile>" << endl;
					ofs.close();
				}
			}

		template <typename T>
			void writePUT(
					int dataClassID_P,
					int dataClassID_UX, int dataClassID_UY, int dataClassID_UZ,
					int dataClassID_T,
					int vc,
					const std::string path,
					const std::string prefix,
					const std::string name,
					int step,
					int difflevel,
					RootGrid* rootGrid,
					BCMOctree* tree,
					Partition* partition,
					Vec3r rootOrigin,
					double rootLength) {
				ostringstream ossFileNameTime;
				ossFileNameTime << path;
				mkdir(ossFileNameTime.str().c_str(), 0755);
				ossFileNameTime << "/";
				ossFileNameTime.width(10);
				ossFileNameTime.setf(ios::fixed);
				ossFileNameTime.fill('0');
				ossFileNameTime << step;
				mkdir(ossFileNameTime.str().c_str(), 0755);

				const Vec3i& size = blockManager.getSize();
				int myrank = comm.Get_rank();

				float* dataP  = new float[(size.x) * (size.y) * (size.z)];
				float* dataUX = new float[(size.x) * (size.y) * (size.z)];
				float* dataUY = new float[(size.x) * (size.y) * (size.z)];
				float* dataUZ = new float[(size.x) * (size.y) * (size.z)];
				float* dataT  = new float[(size.x) * (size.y) * (size.z)];

				for (int id = 0; id < blockManager.getNumBlock(); ++id) {
					BlockBase* block = blockManager.getBlock(id);
					Vec3i size = block->getSize();
					Vec3r origin = block->getOrigin();
					Vec3r blockSize = block->getBlockSize();
					Vec3r cellSize = block->getCellSize();
					int level = block->getLevel();
					int blockid = id + blockManager.getStartID();

					Scalar3D<T>* sp = dynamic_cast<Scalar3D<T>*>(block->getDataClass(dataClassID_P));
					T* sDataP = sp->getData();
					Scalar3D<T>* sux = dynamic_cast<Scalar3D<T>*>(block->getDataClass(dataClassID_UX));
					T* sDataUX = sux->getData();
					Scalar3D<T>* suy = dynamic_cast<Scalar3D<T>*>(block->getDataClass(dataClassID_UY));
					T* sDataUY = suy->getData();
					Scalar3D<T>* suz = dynamic_cast<Scalar3D<T>*>(block->getDataClass(dataClassID_UZ));
					T* sDataUZ = suz->getData();
					Scalar3D<T>* st = dynamic_cast<Scalar3D<T>*>(block->getDataClass(dataClassID_T));
					T* sDataT = st->getData();

					printVTIC(sDataP, sDataUX, sDataUY, sDataUZ, sDataT, path.c_str(), prefix.c_str(), name.c_str(), step, myrank, blockid, size[0], size[1], size[2], vc, origin[0], origin[1], origin[2], cellSize[0]);
				}

				delete[] dataP;
				delete[] dataUX;
				delete[] dataUY;
				delete[] dataUZ;
				delete[] dataT;

				ostringstream ossFileName;
				ossFileName << path;
				ossFileName << "/";
				ossFileName << prefix;
				ossFileName << name.c_str();
				ossFileName << "-";
				ossFileName.width(10);
				ossFileName.setf(ios::fixed);
				ossFileName.fill('0');
				ossFileName << step;
				ossFileName << ".vthb";

				if( myrank == 0 ) {
					ofstream ofs;
					ofs.open(ossFileName.str().c_str(), ios::out);
					ofs << "<VTKFile type=\"vtkHierarchicalBoxDataSet\" version=\"0.1\">" << endl;
					ofs << "<vtkHierarchicalBoxDataSet>" << endl;

					for(int n=0; n<difflevel+1; n++) {
						ofs << "<RefinementRatio level=\"";
						ofs << n;
						ofs << "\" refinement=\"2\"/>";
						ofs << endl;
					}
					ofs << endl;

					std::vector<Node*>& leafNodeArray = tree->getLeafNodeArray();
					for (int iRank = 0; iRank < comm.Get_size(); iRank++) {
						for (int id = partition->getStart(iRank); id < partition->getEnd(iRank); id++) {
							Node* node = leafNodeArray[id];
							Vec3r origin = tree->getOrigin(node) * rootLength;
							Vec3r blockSize = node->getBlockSize() * rootLength;
							Vec3r cellSize;
							cellSize.x = blockSize.x / size.x;
							cellSize.y = blockSize.y / size.y;
							cellSize.z = blockSize.z / size.z;
							int level = node->getLevel();

							ostringstream ossFileName2;
							ossFileName2 << "./";
							ossFileName2.width(10);
							ossFileName2.setf(ios::fixed);
							ossFileName2.fill('0');
							ossFileName2 << step;
							ossFileName2 << "/";
							ossFileName2 << prefix;
							ossFileName2 << name.c_str();
/*
							ossFileName2 << "-";
							ossFileName2.width(5);
							ossFileName2.setf(ios::fixed);
							ossFileName2.fill('0');
							ossFileName2 << iRank;
*/
							ossFileName2 << "-";
							ossFileName2.width(5);
							ossFileName2.setf(ios::fixed);
							ossFileName2.fill('0');
							ossFileName2 << id;
							ossFileName2 << "-";
							ossFileName2.width(10);
							ossFileName2.setf(ios::fixed);
							ossFileName2.fill('0');
							ossFileName2 << step;
							ossFileName2 << ".vti";

							int lx = size.x*(1 << level);
							int ly = size.y*(1 << level);
							int lz = size.z*(1 << level);

							Vec3r origin2 = tree->getOrigin(node);

							double nx0 = (origin2.x)*lx;
							double ny0 = (origin2.y)*ly;
							double nz0 = (origin2.z)*lz;

							double nx1 = nx0 + size.x - 1;
							double ny1 = ny0 + size.y - 1;
							double nz1 = nz0 + size.z - 1;

							ofs << "<DataSet group=\"";
							ofs << level;
							ofs << "\" dataset=\"";
							ofs << id;
							ofs << "\" amr_box=\"";
							ofs << nx0;
							ofs << " ";
							ofs << nx1;
							ofs << " ";
							ofs << ny0;
							ofs << " ";
							ofs << ny1;
							ofs << " ";
							ofs << nz0;
							ofs << " ";
							ofs << nz1;
							ofs << "\" file=\"";
							ofs << ossFileName2.str().c_str();
							ofs << "\"/>";
							ofs << endl;
						}
					}

					ofs << "</vtkHierarchicalBoxDataSet>" << endl;
					ofs << "</VTKFile>" << endl;
					ofs.close();
				}
			}

		template <typename T>
			void writeVtkOverlappingAMR_PUT(
					int dataClassID_P,
					int dataClassID_UX, int dataClassID_UY, int dataClassID_UZ,
					int dataClassID_T,
					int vc,
					const std::string path,
					const std::string prefix,
					const std::string name,
					int step,
					int maxLevel,
					int minLevel,
					RootGrid* rootGrid,
					BCMOctree* tree,
					Partition* partition,
					Vec3r rootOrigin,
					double rootLength) {
				ostringstream ossFileNameTime;
				ossFileNameTime << path;
				mkdir(ossFileNameTime.str().c_str(), 0755);
				ossFileNameTime << "/";
				ossFileNameTime.width(10);
				ossFileNameTime.setf(ios::fixed);
				ossFileNameTime.fill('0');
				ossFileNameTime << step;
				mkdir(ossFileNameTime.str().c_str(), 0755);

				const Vec3i& size = blockManager.getSize();
				int myrank = comm.Get_rank();

				float* dataP  = new float[(size.x) * (size.y) * (size.z)];
				float* dataUX = new float[(size.x) * (size.y) * (size.z)];
				float* dataUY = new float[(size.x) * (size.y) * (size.z)];
				float* dataUZ = new float[(size.x) * (size.y) * (size.z)];
				float* dataT  = new float[(size.x) * (size.y) * (size.z)];

				for (int id = 0; id < blockManager.getNumBlock(); ++id) {
					BlockBase* block = blockManager.getBlock(id);
					Vec3i size = block->getSize();
					Vec3r origin = block->getOrigin();
					Vec3r blockSize = block->getBlockSize();
					Vec3r cellSize = block->getCellSize();
					int level = block->getLevel();
					int blockid = id + blockManager.getStartID();

					Scalar3D<T>* sp = dynamic_cast<Scalar3D<T>*>(block->getDataClass(dataClassID_P));
					T* sDataP = sp->getData();
					Scalar3D<T>* sux = dynamic_cast<Scalar3D<T>*>(block->getDataClass(dataClassID_UX));
					T* sDataUX = sux->getData();
					Scalar3D<T>* suy = dynamic_cast<Scalar3D<T>*>(block->getDataClass(dataClassID_UY));
					T* sDataUY = suy->getData();
					Scalar3D<T>* suz = dynamic_cast<Scalar3D<T>*>(block->getDataClass(dataClassID_UZ));
					T* sDataUZ = suz->getData();
					Scalar3D<T>* st = dynamic_cast<Scalar3D<T>*>(block->getDataClass(dataClassID_T));
					T* sDataT = st->getData();

					//			printVTIC(sDataP, sDataUX, sDataUY, sDataUZ, sDataT, name.c_str(), step, myrank, blockid, size[0], size[1], size[2], vc, origin[0], origin[1], origin[2], cellSize[0]);
					//			printVTIP(sData, name.c_str(), step, myrank, blockid, size[0], size[1], size[2], vc, origin[0], origin[1], origin[2], cellSize[0]);
				}

				delete[] dataP;
				delete[] dataUX;
				delete[] dataUY;
				delete[] dataUZ;
				delete[] dataT;

				ostringstream ossFileName;
				ossFileName << path;
				ossFileName << "/";
				ossFileName << prefix;
				ossFileName << name.c_str();
				ossFileName << "-";
				ossFileName.width(10);
				ossFileName.setf(ios::fixed);
				ossFileName.fill('0');
				ossFileName << step;
				ossFileName << ".vth";

				if( myrank == 0 ) {
					ofstream ofs;
					ofs.open(ossFileName.str().c_str(), ios::out);
					ofs << "<VTKFile type=\"vtkOverlappingAMR\" version=\"1.1\">" << endl;
					ofs << "<vtkOverlappingAMR origin=\"0 0 0\" grid_description=\"XYZ\">" << endl;

					std::vector<Node*>& leafNodeArray = tree->getLeafNodeArray();
					//			int lmax = difflevel+1;
					int lmax = maxLevel+1;
					for(int n=0; n<lmax; n++) {
						double dx = rootLength/(1 << n);
						ofs << "\t<Block level=\"";
						ofs << n;
						ofs << "\" spacing=\"";
						ofs << dx;
						ofs << " ";
						ofs << dx;
						ofs << " ";
						ofs << dx;
						ofs << "\">";
						ofs << endl;
						for (int iRank = 0; iRank < comm.Get_size(); iRank++) {
							for (int id = partition->getStart(iRank); id < partition->getEnd(iRank); id++) {
								Node* node = leafNodeArray[id];
								Vec3r origin = tree->getOrigin(node) * rootLength;
								Vec3r blockSize = node->getBlockSize() * rootLength;
								Vec3r cellSize;
								cellSize.x = blockSize.x / size.x;
								cellSize.y = blockSize.y / size.y;
								cellSize.z = blockSize.z / size.z;
								int level = node->getLevel();

								ostringstream ossFileName2;
								ossFileName2 << "./";
								ossFileName2.width(10);
								ossFileName2.setf(ios::fixed);
								ossFileName2.fill('0');
								ossFileName2 << step;
								ossFileName2 << "/";
								ossFileName2 << prefix;
								ossFileName2 << name.c_str();
/*
								ossFileName2 << "-";
								ossFileName2.width(5);
								ossFileName2.setf(ios::fixed);
								ossFileName2.fill('0');
								ossFileName2 << iRank;
*/
								ossFileName2 << "-";
								ossFileName2.width(5);
								ossFileName2.setf(ios::fixed);
								ossFileName2.fill('0');
								ossFileName2 << id;
								ossFileName2 << "-";
								ossFileName2.width(10);
								ossFileName2.setf(ios::fixed);
								ossFileName2.fill('0');
								ossFileName2 << step;
								ossFileName2 << ".vti";

								int lx = size.x*(1 << level);
								int ly = size.y*(1 << level);
								int lz = size.z*(1 << level);

								Vec3r origin2 = tree->getOrigin(node);

								double nx0 = (origin2.x)*lx;
								double ny0 = (origin2.y)*ly;
								double nz0 = (origin2.z)*lz;

								double nx1 = nx0 + size.x - 1;
								double ny1 = ny0 + size.y - 1;
								double nz1 = nz0 + size.z - 1;

								if( level == n ) {
									ofs << "\t\t<DataSet index=\"";
									ofs << id;
									ofs << "\" amr_box=\"";
									ofs << nx0;
									ofs << " ";
									ofs << nx1;
									ofs << " ";
									ofs << ny0;
									ofs << " ";
									ofs << ny1;
									ofs << " ";
									ofs << nz0;
									ofs << " ";
									ofs << nz1;
									ofs << "\" file=\"";
									ofs << ossFileName2.str().c_str();
									ofs << "\">";
									ofs << endl;
									ofs << "\t\t</DataSet>" << endl;	
								}
							}
						}
						ofs << "\t</Block>";
						ofs << endl;
						ofs << endl;
					}

					ofs << "</vtkOverlappingAMR>" << endl;
					ofs << "</VTKFile>" << endl;
					ofs.close();
				}
			}
	};

#ifdef BCMT_NAMESPACE
} // namespace BCMT_NAMESPACE
#endif

#endif // VTK_WRITER_H

