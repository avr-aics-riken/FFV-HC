#ifndef FFVMC3D_H
#define FFVMC3D_H

#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>

#include "BCMTools.h"
#include "BlockManager.h"
#include "Scalar3D.h"
#include "BCMOctree.h"
#include "Partition.h"

typedef struct _Vertex {
	double x;
	double y;
	double z;
	double value;
}Vertex;

extern const int edgeTable[256];
extern const int triTable[256][16];
extern Vertex GetIntersection(Vertex v1, Vertex v2, double p1, double p2, double threshold);

class FFVMC {
public:
	FFVMC()
		: blockManager(BlockManager::getInstance()),
			comm(blockManager.getCommunicator()) {
		int myrank = comm.Get_rank();
		for(int id=0; id<blockManager.getNumBlock(); ++id) {
			BlockBase* block = blockManager.getBlock(id);
		}

		if( myrank == 0 ) {
		} else {
		}

		PNX = 0;
		PNY = 0;
		PNZ = 0;
		pPointData = 0;
		vVertexList.clear();
	}
	virtual ~FFVMC() {
	}

private:
	BlockManager& blockManager;
	const MPI::Intracomm& comm;

public:
template <typename T>
	void writeContour(
				int dataClassID,
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
				double rootLength,
				double threshold) {
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

		for (int id = 0; id < blockManager.getNumBlock(); ++id) {
			BlockBase* block = blockManager.getBlock(id);
			Vec3i size = block->getSize();
			Vec3r origin = block->getOrigin();
			Vec3r blockSize = block->getBlockSize();
			Vec3r cellSize = block->getCellSize();
			int level = block->getLevel();

			Scalar3D<T>* s = dynamic_cast<Scalar3D<T>*>(block->getDataClass(dataClassID));
			T* sData = s->getData();

			int nx = size.x;
			int ny = size.y;
			int nz = size.z;
			int nv = vc;
			double ox = origin.x;
			double oy = origin.y;
			double oz = origin.z;
			double dx = cellSize.x;
			writeContourLocal(
					sData,
					threshold,
					path,
					prefix,
					name,
					step, myrank, id,
					nx, ny, nz,
					nv,
					ox, oy, oz,
					dx);
		}

		if( myrank != 0 ) {
			return;
		}

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
		ossFileName << ".pvtp";

		std::ofstream ofs;
		ofs.open(ossFileName.str().c_str(), std::ios::out);
		ofs << "<VTKFile type=\"PPolyData\" version=\"0.1\">" << std::endl;
		ofs << "<PPolyData GhostLevel=\"0\">" << std::endl;

		ofs << "<PPointData>" << std::endl;
		ofs << "<PDataArray type=\"Float32\" Name=\"";
		ofs << name;
		ofs << "\" format=\"ascii\"/>" << std::endl;
		ofs << "</PPointData>" << std::endl;

		ofs << "<PCellData>" << std::endl;
		ofs << "</PCellData>" << std::endl;

		ofs << "<PPoints>" << std::endl;
		ofs << "<PDataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\"/>" << std::endl;
		ofs << "</PPoints>" << std::endl;

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
				ossFileName2 << "-";
				ossFileName2.width(5);
				ossFileName2.setf(ios::fixed);
				ossFileName2.fill('0');
				ossFileName2 << iRank;
				ossFileName2 << "-";
				ossFileName2.width(5);
				ossFileName2.setf(ios::fixed);
				ossFileName2.fill('0');
				ossFileName2 << id - partition->getStart(iRank);
				ossFileName2 << "-";
				ossFileName2.width(10);
				ossFileName2.setf(ios::fixed);
				ossFileName2.fill('0');
				ossFileName2 << step;
				ossFileName2 << ".vtp";

				ofs << "<Piece Source=\"";
				ofs << ossFileName2.str();
				ofs << "\"/>" << std::endl;

			}
		}

		ofs << "</PPolyData>" << std::endl;
		ofs << "</VTKFile>" << std::endl;
		ofs.close();
	}

template <typename T>
	void writeContourLocal(
			T* pData,
			double threshold,
			const std::string path,
			const std::string prefix,
			const std::string name,
			int step, int rank, int block,
			int NX, int NY, int NZ,
			int NV,
			double ox, double oy, double oz,
			double dx)
	{
		InitPointData(pData, NX, NY, NZ, NV);
		ClearTriangles();
		DetectTriangles(threshold, ox, oy, oz, dx);
		PrintVTP(path, prefix, name, step, rank, block, threshold);
	}

template <typename T>
	void InitPointData(
			T* pData,
			int NX, int NY, int NZ,
			int NV)
	{
		if( NX + 1 == PNX &&
				NY + 1 == PNY &&
				NZ + 1 == PNZ ) {
		} else {
			PNX = NX + 1;
			PNY = NY + 1;
			PNZ = NZ + 1;

			if( pPointData ) {
				delete [] pPointData;
				pPointData = 0;
			}
			pPointData = new float [PNX*PNY*PNZ];
		}

		int CX = NX + 2*NV;
		int CY = NY + 2*NV;
		int CZ = NZ + 2*NV;

#pragma omp parallel for
		for(int k=0; k<PNZ; k++) {
			for(int j=0; j<PNY; j++) {
				for(int i=0; i<PNX; i++) {
					int i0 = NV + i;
					int j0 = NV + j;
					int k0 = NV + k;
					int m0 = i0 + CX*(j0 + CY*k0);
					int m1 = i0-1 + CX*(j0 + CY*k0);
					int m2 = i0-1 + CX*((j0-1) + CY*k0);
					int m3 = i0 + CX*((j0-1) + CY*k0);
					int m4 = i0 + CX*(j0 + CY*(k0-1));
					int m5 = i0-1 + CX*(j0 + CY*(k0-1));
					int m6 = i0-1 + CX*((j0-1) + CY*(k0-1));
					int m7 = i0 + CX*((j0-1) + CY*(k0-1));
					T phi0 = pData[m0];
					T phi1 = pData[m1];
					T phi2 = pData[m2];
					T phi3 = pData[m3];
					T phi4 = pData[m4];
					T phi5 = pData[m5];
					T phi6 = pData[m6];
					T phi7 = pData[m7];

					int mp = i + PNX*(j + PNY*k);
					pPointData[mp] = 0.125*(phi0 + phi1 + phi2 + phi3 + phi4 + phi5 + phi6 + phi7);
				}
			}
		}
	}

	void ClearTriangles() {
		vVertexList.clear();
	}

	void DetectTriangles(
			double threshold,
			double ox, double oy, double oz,
			double dx)
	{
#pragma omp parallel for
		for(int k=0; k<PNZ-1; k++) {
			for(int j=0; j<PNY-1; j++) {
				for(int i=0; i<PNX-1; i++) {
					Vertex v[8];
					v[3].x = ox + dx*i;
					v[3].y = oy + dx*j;
					v[3].z = oz + dx*k;
					v[2].x = ox + dx*(i+1);
					v[2].y = oy + dx*j;
					v[2].z = oz + dx*k;
					v[1].x = ox + dx*(i+1);
					v[1].y = oy + dx*(j+1);
					v[1].z = oz + dx*k;
					v[0].x = ox + dx*i;
					v[0].y = oy + dx*(j+1);
					v[0].z = oz + dx*k;
					v[7].x = ox + dx*i;
					v[7].y = oy + dx*j;
					v[7].z = oz + dx*(k+1);
					v[6].x = ox + dx*(i+1);
					v[6].y = oy + dx*j;
					v[6].z = oz + dx*(k+1);
					v[5].x = ox + dx*(i+1);
					v[5].y = oy + dx*(j+1);
					v[5].z = oz + dx*(k+1);
					v[4].x = ox + dx*i;
					v[4].y = oy + dx*(j+1);
					v[4].z = oz + dx*(k+1);

					double p[8];
					p[3] = pPointData[i   + PNX*(j   + PNY*k)];
					p[2] = pPointData[i+1 + PNX*(j   + PNY*k)];
					p[1] = pPointData[i+1 + PNX*(j+1 + PNY*k)];
					p[0] = pPointData[i   + PNX*(j+1 + PNY*k)];
					p[7] = pPointData[i   + PNX*(j   + PNY*(k+1))];
					p[6] = pPointData[i+1 + PNX*(j   + PNY*(k+1))];
					p[5] = pPointData[i+1 + PNX*(j+1 + PNY*(k+1))];
					p[4] = pPointData[i   + PNX*(j+1 + PNY*(k+1))];

					int cubeindex = 0;
					if( p[0] < threshold ) {
						cubeindex |= 1;
					}
					if( p[1] < threshold ) {
						cubeindex |= 2;
					}
					if( p[2] < threshold ) {
						cubeindex |= 4;
					}
					if( p[3] < threshold ) {
						cubeindex |= 8;
					}
					if( p[4] < threshold ) {
						cubeindex |= 16;
					}
					if( p[5] < threshold ) {
						cubeindex |= 32;
					}
					if( p[6] < threshold ) {
						cubeindex |= 64;
					}
					if( p[7] < threshold ) {
						cubeindex |= 128;
					}

					Vertex u[12];
					if( edgeTable[cubeindex] & 1 ) {
						u[0] = GetIntersection(v[0], v[1], p[0], p[1], threshold);
					}
					if( edgeTable[cubeindex] & 2 ) {
						u[1] = GetIntersection(v[1], v[2], p[1], p[2], threshold);
					}
					if( edgeTable[cubeindex] & 4 ) {
						u[2] = GetIntersection(v[2], v[3], p[2], p[3], threshold);
					}
					if( edgeTable[cubeindex] & 8 ) {
						u[3] = GetIntersection(v[3], v[0], p[3], p[0], threshold);
					}
					if( edgeTable[cubeindex] & 16 ) {
						u[4] = GetIntersection(v[4], v[5], p[4], p[5], threshold);
					}
					if( edgeTable[cubeindex] & 32 ) {
						u[5] = GetIntersection(v[5], v[6], p[5], p[6], threshold);
					}
					if( edgeTable[cubeindex] & 64 ) {
						u[6] = GetIntersection(v[6], v[7], p[6], p[7], threshold);
					}
					if( edgeTable[cubeindex] & 128 ) {
						u[7] = GetIntersection(v[7], v[4], p[7], p[4], threshold);
					}
					if( edgeTable[cubeindex] & 256 ) {
						u[8] = GetIntersection(v[0], v[4], p[0], p[4], threshold);
					}
					if( edgeTable[cubeindex] & 512 ) {
						u[9] = GetIntersection(v[1], v[5], p[1], p[5], threshold);
					}
					if( edgeTable[cubeindex] & 1024 ) {
						u[10] = GetIntersection(v[2], v[6], p[2], p[6], threshold);
					}
					if( edgeTable[cubeindex] & 2048 ) {
						u[11] = GetIntersection(v[3], v[7], p[3], p[7], threshold);
					}

					for(int n=0; n<12; n++) {
						u[n].value = threshold;
					}

#pragma omp critical
{
					int nTriangle = 0;
					for(int i=0; triTable[cubeindex][i] != -1; i+=3) {
						vVertexList.push_back( u[triTable[cubeindex][i]] );
						vVertexList.push_back( u[triTable[cubeindex][i+1]] );
						vVertexList.push_back( u[triTable[cubeindex][i+2]] );
						nTriangle++;
					}
}
//				std::cout << nTriangle << std::endl;
				}
			}
		}
	}

	void PrintVTP(
			const std::string path,
			const std::string prefix,
			const std::string name,
			int step, int rank, int block,
			double threshold)
	{
		std::ostringstream ossFileName2;
		ossFileName2 << path;
		ossFileName2 << "/";
		ossFileName2.width(10);
		ossFileName2.setf(std::ios::fixed);
		ossFileName2.fill('0');
		ossFileName2 << step;
		ossFileName2 << "/";
		ossFileName2 << prefix;
		ossFileName2 << name.c_str();
		ossFileName2 << "-";
		ossFileName2.width(5);
		ossFileName2.setf(std::ios::fixed);
		ossFileName2.fill('0');
		ossFileName2 << rank;
		ossFileName2 << "-";
		ossFileName2.width(5);
		ossFileName2.setf(std::ios::fixed);
		ossFileName2.fill('0');
		ossFileName2 << block;
		ossFileName2 << "-";
		ossFileName2.width(10);
		ossFileName2.setf(std::ios::fixed);
		ossFileName2.fill('0');
		ossFileName2 << step;
		ossFileName2 << ".vtp";

		int nPoints = vVertexList.size();
		int nPolys = vVertexList.size()/3;
		int nLines = 0;

		std::ofstream ofs;
		ofs.open(ossFileName2.str().c_str(), std::ios::out);
		ofs << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"";
#ifdef __FUJITSU
		ofs << "BigEndian";
#else
		ofs << "LittleEndian";
#endif
		ofs << "\">" << std::endl;
		ofs << "<PolyData>" << std::endl;

		ofs << "<Piece NumberOfPoints=\"";
		ofs << nPoints;
		ofs << "\" NumberOfVerts=\"0\" NumberOfLines=\"";
		ofs << nLines;
		ofs << "\" NumberOfStrips=\"0\" NumberOfPolys=\"";
		ofs << nPolys;
		ofs << "\">" << std::endl;

		ofs << "<Points>" << std::endl;
		ofs << "<DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
		for(int i=0; i<vVertexList.size(); i++) {
			ofs << vVertexList[i].x << " " << vVertexList[i].y << " " << vVertexList[i].z << std::endl;
		}
		ofs << "</DataArray>" << std::endl;
		ofs << "</Points>" << std::endl;

		ofs << "<Polys>" << std::endl;
		ofs << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << std::endl;
		for(int i=0; i<vVertexList.size()/3; i++) {
			ofs << 3*i << " " << 3*i+1 << " " << 3*i+2 << std::endl;
		}
		ofs << "</DataArray>" << std::endl;
		ofs << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << std::endl;
		for(int i=0; i<vVertexList.size()/3; i++) {
			ofs << 3*(i+1) << std::endl;
		}
		ofs << "</DataArray>" << std::endl;
		ofs << "</Polys>" << std::endl;

		ofs << "<Lines>" << std::endl;
		ofs << "</Lines>" << std::endl;

		ofs << "<Verts>" << std::endl;
		ofs << "</Verts>" << std::endl;

		ofs << "<Strips>" << std::endl;
		ofs << "</Strips>" << std::endl;

		ofs << "<PointData>" << std::endl;
		ofs << "<DataArray type=\"Float32\" Name=\"";
		ofs << name;
		ofs << "\" format=\"ascii\">" << std::endl;
		for(int i=0; i<vVertexList.size(); i++) {
			ofs << vVertexList[i].value << std::endl;
		}
		ofs << "</DataArray>" << std::endl;
		ofs << "</PointData>" << std::endl;

		ofs << "<CellData>" << std::endl;
		ofs << "</CellData>" << std::endl;

		ofs << "</Piece>" << std::endl;

		ofs << "</PolyData>" << std::endl;
		ofs << "</VTKFile>" << std::endl;
		ofs.close();
	}

private:



private:
	int PNX;
	int PNY;
	int PNZ;

	float *pPointData;
	std::vector<Vertex> vVertexList;
};

#endif 

