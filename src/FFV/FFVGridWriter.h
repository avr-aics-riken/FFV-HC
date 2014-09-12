#ifndef GRID_WRITER_H
#define GRID_WRITER_H

#include <string>
#include <sstream>
#include <fstream>
#include <cstring>
#include <sys/stat.h>
#include <mpi.h>

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

	class GridWriter {
		BlockManager& blockManager;
		const MPI::Intracomm& comm;

		public:
		GridWriter()
			: blockManager(BlockManager::getInstance()),
			comm(blockManager.getCommunicator()) {
				int myrank = comm.Get_rank();

				for (int id = 0; id < blockManager.getNumBlock(); ++id) {
					BlockBase* block = blockManager.getBlock(id);
				}

				if (myrank == 0) {
				} else {
				}
			}

		~GridWriter() {
		}

		void writeGrid(
				RootGrid* rootGrid,
				BCMOctree* tree,
				Partition* partition,
				Vec3r rootOrigin,
				double rootLength) {
			int myrank = comm.Get_rank();
			if (myrank == 0) {
				ofstream ofs;
				ofs.open("data-grid.obj", ios::out);
				ofs << "o grid" << std::endl;
				ofs << "g grid" << std::endl;

				std::vector<Node*>& leafNodeArray = tree->getLeafNodeArray();
				for (int iRank = 0; iRank < comm.Get_size(); iRank++) {
					for (int id = partition->getStart(iRank); id < partition->getEnd(iRank); id++) {
						Node* node = leafNodeArray[id];
						Vec3r origin = tree->getOrigin(node) * rootLength;
						Vec3r blockSize = node->getBlockSize() * rootLength;

						float x0 = origin[0];
						float y0 = origin[1];
						float z0 = origin[2];
						float x1 = origin[0] + blockSize[0];
						float y1 = origin[1] + blockSize[1];
						float z1 = origin[2] + blockSize[2];
						float x[2] = {x0, x1};
						float y[2] = {y0, y1};
						float z[2] = {z0, z1};
						float v[9];

						ofs << "v " << x[0] << " " << y[0] << " " << z[0] << std::endl;
						ofs << "v " << x[1] << " " << y[0] << " " << z[0] << std::endl;
						ofs << "v " << x[1] << " " << y[1] << " " << z[0] << std::endl;
						ofs << "v " << x[0] << " " << y[1] << " " << z[0] << std::endl;
						ofs << "v " << x[0] << " " << y[0] << " " << z[1] << std::endl;
						ofs << "v " << x[1] << " " << y[0] << " " << z[1] << std::endl;
						ofs << "v " << x[1] << " " << y[1] << " " << z[1] << std::endl;
						ofs << "v " << x[0] << " " << y[1] << " " << z[1] << std::endl;
					}
				}
				int m = 0;
				for (int iRank = 0; iRank < comm.Get_size(); iRank++) {
					for (int id = partition->getStart(iRank); id < partition->getEnd(iRank); id++) {
						Node* node = leafNodeArray[id];
						Vec3r origin = tree->getOrigin(node) * rootLength;
						Vec3r blockSize = node->getBlockSize() * rootLength;
						int i0 = 8*m + 1;
						int i1 = 8*m + 2;
						int i2 = 8*m + 3;
						int i3 = 8*m + 4;
						int i4 = 8*m + 5;
						int i5 = 8*m + 6;
						int i6 = 8*m + 7;
						int i7 = 8*m + 8;
						ofs << "f " << i0 << " " << i1 << " " << i2 << " " << i3 << std::endl;
						ofs << "f " << i4 << " " << i7 << " " << i6 << " " << i5 << std::endl;
						ofs << "f " << i1 << " " << i5 << " " << i6 << " " << i2 << std::endl;
						ofs << "f " << i0 << " " << i3 << " " << i7 << " " << i4 << std::endl;
						ofs << "f " << i0 << " " << i4 << " " << i5 << " " << i1 << std::endl;
						ofs << "f " << i2 << " " << i6 << " " << i7 << " " << i3 << std::endl;
						m++;
					}
				}
				ofs.close();
			}
		}

		void WritePolygon(std::ofstream& ofs, float* pv) {
			ofs << "facet normal 0 0 0" << std::endl;
			ofs << "outer loop" << std::endl;
			ofs << "vertex " << pv[0] << " " << pv[1] << " " << pv[2] << std::endl;
			ofs << "vertex " << pv[3] << " " << pv[4] << " " << pv[5] << std::endl;
			ofs << "vertex " << pv[6] << " " << pv[7] << " " << pv[8] << std::endl;
			ofs << "endloop" << std::endl;
			ofs << "endfacet" << std::endl;
		}

	};

#ifdef BCMT_NAMESPACE
} // namespace BCMT_NAMESPACE
#endif

#endif

