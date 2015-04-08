#ifndef SOLVER_H
#define SOLVER_H

#include "BCMTools.h"
#include "BlockManager.h"

#include "LocalScalar3D.h"
#include "FFVVersion.h"
#include "FFVILS.h"
#include "FFVVTKWriter.h"
#include "FFVPlot3DWriter.h"
#include "FFVMC.h"
#include "FFVGridWriter.h"
#include "BCMFileSaver.h"

#include "real.h"

class Solver {
	public:
		Solver();
		~Solver();

	private:
		BlockManager& blockManager;
		int myrank;
		int vc;
		int diffLevel;
		int maxLevel;
		int minLevel;
		std::string updateMethod;

		RootGrid* rootGrid;
		BCMOctree* tree;
		Partition* partition;

		bool bRestart;

	private:
		real dt;

		real rhof;
		real rhos;
		real cpf;
		real cps;
		real kf;
		real ks;
		real mu;
		real csf;

	private:
		LocalScalar3D<real> *plsUX0;
		LocalScalar3D<real> *plsUX1;
		LocalScalar3D<real> *plsUXC;
		LocalScalar3D<real> *plsUXCP;
		LocalScalar3D<real> *plsUXD;
		LocalScalar3D<real> *plsUXDP;
		LocalScalar3D<real> *plsUY0;
		LocalScalar3D<real> *plsUY1;
		LocalScalar3D<real> *plsUYC;
		LocalScalar3D<real> *plsUYCP;
		LocalScalar3D<real> *plsUYD;
		LocalScalar3D<real> *plsUYDP;
		LocalScalar3D<real> *plsUZ0;
		LocalScalar3D<real> *plsUZ1;
		LocalScalar3D<real> *plsUZC;
		LocalScalar3D<real> *plsUZCP;
		LocalScalar3D<real> *plsUZD;
		LocalScalar3D<real> *plsUZDP;

		LocalScalar3D<real> *plsFX;
		LocalScalar3D<real> *plsFY;
		LocalScalar3D<real> *plsFZ;

		LocalScalar3D<real> *plsVw;
		LocalScalar3D<real> *plsVe;
		LocalScalar3D<real> *plsVs;
		LocalScalar3D<real> *plsVn;
		LocalScalar3D<real> *plsVb;
		LocalScalar3D<real> *plsVt;

		LocalScalar3D<real> *plsP0;
		LocalScalar3D<real> *plsP1;

		LocalScalar3D<real> *plsT0;
		LocalScalar3D<real> *plsT1;
		LocalScalar3D<real> *plsTC;
		LocalScalar3D<real> *plsTCP;
		LocalScalar3D<real> *plsTD;
		LocalScalar3D<real> *plsTDP;

		LocalScalar3D<int> *plsPhaseId;
		LocalScalar3D<int> *plsRegionId;
		LocalScalar3D<real> *plsCut0;
		LocalScalar3D<real> *plsCut1;
		LocalScalar3D<real> *plsCut2;
		LocalScalar3D<real> *plsCut3;
		LocalScalar3D<real> *plsCut4;
		LocalScalar3D<real> *plsCut5;
		LocalScalar3D<int> *plsCutId0;
		LocalScalar3D<int> *plsCutId1;
		LocalScalar3D<int> *plsCutId2;
		LocalScalar3D<int> *plsCutId3;
		LocalScalar3D<int> *plsCutId4;
		LocalScalar3D<int> *plsCutId5;

		int *pNormalN;
		real **pNormalX;
		real **pNormalY;
		real **pNormalZ;
		LocalScalar3D<int> *plsNormalIndex0;
		LocalScalar3D<int> *plsNormalIndex1;
		LocalScalar3D<int> *plsNormalIndex2;
		LocalScalar3D<int> *plsNormalIndex3;
		LocalScalar3D<int> *plsNormalIndex4;
		LocalScalar3D<int> *plsNormalIndex5;

		LocalScalar3D<real> *plsAp;
		LocalScalar3D<real> *plsAw;
		LocalScalar3D<real> *plsAe;
		LocalScalar3D<real> *plsAs;
		LocalScalar3D<real> *plsAn;
		LocalScalar3D<real> *plsAb;
		LocalScalar3D<real> *plsAt;
		LocalScalar3D<real> *plsb;

		LocalScalar3D<real> *plsr;
		LocalScalar3D<real> *plsr0;
		LocalScalar3D<real> *plsp;
		LocalScalar3D<real> *plsp_;
		LocalScalar3D<real> *plsq_;
		LocalScalar3D<real> *plss;
		LocalScalar3D<real> *plss_;
		LocalScalar3D<real> *plst_;

		LocalScalar3D<real> *plsLapP;

		LocalScalar3D<real> *plsFspx;
		LocalScalar3D<real> *plsFspy;
		LocalScalar3D<real> *plsFspz;
		LocalScalar3D<real> *plsFsvx;
		LocalScalar3D<real> *plsFsvy;
		LocalScalar3D<real> *plsFsvz;

		LocalScalar3D<real> *plsQx;
		LocalScalar3D<real> *plsQy;
		LocalScalar3D<real> *plsQz;

		LocalScalar3D<real> *plsM;
		LocalScalar3D<int> *plsMaskId;

	private:
		FFVILS* pils;

		real omegaU;
		int countMaxU;
		real epsilonU;
		int countPreConditionerU;
		int countUX;
		real residualUX;
		int countUY;
		real residualUY;
		int countUZ;
		real residualUZ;

		real omegaP;
		int countMaxP;
		real epsilonP;
		int countPreConditionerP;
		int countP;
		real residualP;

		real omegaT;
		int countMaxT;
		real epsilonT;
		int countPreConditionerT;
		int countT;
		real residualT;

	public:
		int Init(int argc, char** argv);
		int Loop();
		int Post();

	private:
		int Update(int step);
		void UpdateUX(int step);
		void UpdateUY(int step);
		void UpdateUZ(int step);
		void UpdateP(int step);
		void UpdateU(int step);
		void UpdateT(int step);
		void UpdateUXe(int step);
		void UpdateUYe(int step);
		void UpdateUZe(int step);
		void UpdateTe(int step);
		void UpdateF(int step);
		int FillRegion(LocalScalar3D<int> *plsId, int value, real xs, real ys, real zs);

		int Print(int step);
		double times[32];
		void PrintTime(int step);
		void PrintILS(int step);
		void PrintStats(int step);
		void PrintForceCID(int step, int cid_target);
		void PrintForce(int step);
		void PrintHeatFluxCID(int step, int cid_target);
		void PrintHeatFlux(int step);
		void PrintLog(int level, const char* format, ...);
		void PrintCS(int step);

		void WritePolygon(std::ofstream& ofs, float *pv);
		void PrintCut(int id);
		void PrintHole(int id);

		void PrintData(int step);
		void PrintBasicVariablesVTK(int step);
		void PrintBasicVariablesPLOT3D(int step);
		void PrintBasicVariablesBCM(int step);
		void PrintBasicVariablesSILO(int step);
		void PrintDerivedVariablesVTK(int step);
		void PrintDerivedVariablesPLOT3D(int step);
		void PrintDerivedVariablesBCM(int step);
		void PrintDerivedVariablesSILO(int step);
		void PrintContourQcriterion(int step);

	private:
		void Dump(const int step);
		void Load(const int step);

	private:
		void WriteGrid(
				RootGrid* rootGrid,
				BCMOctree* tree,
				Partition* partition) {
			GridWriter writer;
			writer.writeGrid(
					rootGrid,
					tree,
					partition,
					g_pFFVConfig->RootBlockOrigin,
					g_pFFVConfig->RootBlockLength);
		}

		void WriteBasicVariablesInVTKFormat(
				int step,
				int difflevel,
				RootGrid* rootGrid,
				BCMOctree* tree,
				Partition* partition) {
			VtkWriter writer;

			writer.writePUT<real>(
					this->plsP0->GetID(),
					this->plsUX0->GetID(),
					this->plsUY0->GetID(),
					this->plsUZ0->GetID(),
					this->plsT0->GetID(),
					this->vc,
					g_pFFVConfig->OutputDataFormatOptionVTKPath,
					g_pFFVConfig->OutputDataFormatOptionVTKPrefix,
					string("flow"),
					step,
					difflevel,
					rootGrid,
					tree,
					partition,
					g_pFFVConfig->RootBlockOrigin,
					g_pFFVConfig->RootBlockLength);

			writer.writeVtkOverlappingAMR_PUT<real>(
					this->plsP0->GetID(),
					this->plsUX0->GetID(),
					this->plsUY0->GetID(),
					this->plsUZ0->GetID(),
					this->plsT0->GetID(),
					this->vc,
					g_pFFVConfig->OutputDataFormatOptionVTKPath,
					g_pFFVConfig->OutputDataFormatOptionVTKPrefix,
					string("flow"),
					step,
					maxLevel,
					minLevel,
					rootGrid,
					tree,
					partition,
					g_pFFVConfig->RootBlockOrigin,
					g_pFFVConfig->RootBlockLength);
		}

		void WriteDerivedVariablesInVTKFormat(
				int step,
				int difflevel,
				RootGrid* rootGrid,
				BCMOctree* tree,
				Partition* partition) {
			VtkWriter writer;

			if( g_pFFVConfig->OutputDataDerivedVariablesHeatFlux ) {
				writer.writePU<real>(
						this->plsT0->GetID(),
						this->plsQx->GetID(),
						this->plsQy->GetID(),
						this->plsQz->GetID(),
						this->vc,
						g_pFFVConfig->OutputDataFormatOptionVTKPath,
						g_pFFVConfig->OutputDataFormatOptionVTKPrefix,
						string("q"),
						step,
						difflevel,
						rootGrid,
						tree,
						partition,
						g_pFFVConfig->RootBlockOrigin,
						g_pFFVConfig->RootBlockLength);
			}

			if( g_pFFVConfig->OutputDataDerivedVariablesQcriterion ) {
				writer.writeScalar<real>(
						this->plsLapP->GetID(),
						this->vc,
						g_pFFVConfig->OutputDataFormatOptionVTKPath,
						g_pFFVConfig->OutputDataFormatOptionVTKPrefix,
						string("lapp"),
						step,
						difflevel,
						rootGrid,
						tree,
						partition,
						g_pFFVConfig->RootBlockOrigin,
						g_pFFVConfig->RootBlockLength);
			}

			if( g_pFFVConfig->OutputDataDerivedVariablesForce ) {
				writer.writePU<real>(
						this->plsP0->GetID(),
						this->plsFspx->GetID(),
						this->plsFspy->GetID(),
						this->plsFspz->GetID(),
						this->vc,
						g_pFFVConfig->OutputDataFormatOptionVTKPath,
						g_pFFVConfig->OutputDataFormatOptionVTKPrefix,
						string("fp"),
						step,
						difflevel,
						rootGrid,
						tree,
						partition,
						g_pFFVConfig->RootBlockOrigin,
						g_pFFVConfig->RootBlockLength);

				writer.writePU<real>(
						this->plsP0->GetID(),
						this->plsFsvx->GetID(),
						this->plsFsvy->GetID(),
						this->plsFsvz->GetID(),
						this->vc,
						g_pFFVConfig->OutputDataFormatOptionVTKPath,
						g_pFFVConfig->OutputDataFormatOptionVTKPrefix,
						string("fv"),
						step,
						difflevel,
						rootGrid,
						tree,
						partition,
						g_pFFVConfig->RootBlockOrigin,
						g_pFFVConfig->RootBlockLength);
			}
		}

		void WriteXYZInPlot3DFormat(
				const char* dataname,
				int difflevel,
				RootGrid* rootGrid,
				BCMOctree* tree,
				Partition* partition) {
			Plot3DWriter writer;
			writer.writeXYZ(
					this->plsPhaseId->GetID(),
					this->vc,
					string(dataname),
					difflevel,
					rootGrid,
					tree,
					partition,
					g_pFFVConfig->RootBlockOrigin,
					g_pFFVConfig->RootBlockLength);
		}

		void WriteBasicVariablesInPlot3DFormat(
				const char* dataname,
				int step,
				int difflevel,
				RootGrid* rootGrid,
				BCMOctree* tree,
				Partition* partition) {
			Plot3DWriter writer;
			writer.writeData<real>(
					this->plsP0->GetID(),
					this->plsUX0->GetID(),
					this->plsUY0->GetID(),
					this->plsUZ0->GetID(),
					this->vc,
					string(dataname),
					step,
					difflevel,
					rootGrid,
					tree,
					partition,
					g_pFFVConfig->RootBlockOrigin,
					g_pFFVConfig->RootBlockLength);
		}

		void WriteContourQriterion(
				int step,
				int maxLevel,
				int minLevel,
				RootGrid* rootGrid,
				BCMOctree* tree,
				Partition* partition) {
			FFVMC mc;
			mc.writeContour<real>(
					this->plsLapP->GetID(),
					this->vc,
					g_pFFVConfig->OutputDataFormatOptionVTPPath,
					g_pFFVConfig->OutputDataFormatOptionVTPPrefix,
					string("Q"),
					step,
					maxLevel,
					minLevel,
					rootGrid,
					tree,
					partition,
					g_pFFVConfig->RootBlockOrigin,
					g_pFFVConfig->RootBlockLength,
					g_pFFVConfig->OutputDataContourQcriterionValue);
		}

		BCMFileIO::BCMFileSaver *psaver;

		void BCMFileSaverInit(
				RootGrid* rootGrid,
				BCMOctree* tree,
				Partition* partition) {
			Vec3r origin = g_pFFVConfig->RootBlockOrigin;
			double rootLength = g_pFFVConfig->RootBlockLength;
			Vec3i rootN = g_pFFVConfig->RootBlockGrid;
			Vec3r region(rootLength*rootN.x,
					rootLength*rootN.y,
					rootLength*rootN.z);
			psaver = new BCMFileIO::BCMFileSaver(origin, region, tree, "BCM_OUT");

			BCMFileIO::IdxUnit unit;
			unit.length = std::string("m");
			unit.L0_scale = 1.0;
			unit.velocity = std::string("Dimensional");
			unit.V0_scale = 1.0;

			BCMFileIO::IdxStep timeStep(
					g_pFFVConfig->TimeControlSessionStartI,
					g_pFFVConfig->TimeControlSessionEndI,
					g_pFFVConfig->OutputDataBasicVariablesIntervalI);

			int vc = g_pFFVConfig->LeafBlockNumberOfVirtualCells;
			int id_cellId = blockManager.setDataClass< Scalar3D<unsigned char> >(vc);
			for (int id = 0; id < blockManager.getNumBlock(); ++id){
				BlockBase* block = blockManager.getBlock(id);
				Vec3i sz = block->getSize();
				Scalar3D<unsigned char>* mesh = dynamic_cast< Scalar3D<unsigned char>* >(block->getDataClass(id_cellId));
				unsigned char* data = mesh->getData();
				Index3DS idx = mesh->getIndex();
				int* pPhaseId = plsPhaseId->GetBlockData(block);
				for(int z = -vc; z < sz.z + vc; z++){
					for(int y = -vc; y < sz.y + vc; y++){
						for(int x = -vc; x < sz.x + vc; x++){
							data[idx(x, y, z)] = pPhaseId[idx(x, y, z)];
						}
					}
				}
			}
			psaver->RegisterCellIDInformation(
					id_cellId,
					5,
					vc,
					"CellID",
					"cid",
					"lb",
					"cid");

			int id_p[1] = {0};
			int id_u[3] = {0, 0, 0};
			id_p[0] = plsP0->GetID();
			id_u[0] = plsUX0->GetID();
			id_u[1] = plsUY0->GetID();
			id_u[2] = plsUZ0->GetID();

#ifdef _REAL_IS_DOUBLE_
			psaver->RegisterDataInformation(
					id_p,
					BCMFileIO::LB_SCALAR,
					BCMFileIO::LB_FLOAT64,
					vc,
					"P64",
					"P64",
					"lb",
					timeStep,
					"P",
					true);
			psaver->RegisterDataInformation(
					id_u,
					BCMFileIO::LB_VECTOR3,
					BCMFileIO::LB_FLOAT64,
					vc,
					"Vel64",
					"Vel64",
					"lb",
					timeStep,
					"VEL",
					true);
#else
			psaver->RegisterDataInformation(
					id_p,
					BCMFileIO::LB_SCALAR,
					BCMFileIO::LB_FLOAT32,
					vc,
					"P32",
					"P32",
					"lb",
					timeStep,
					"P",
					true);
			psaver->RegisterDataInformation(
					id_u,
					BCMFileIO::LB_VECTOR3,
					BCMFileIO::LB_FLOAT32,
					vc,
					"Vel32",
					"Vel32",
					"lb",
					timeStep,
					"VEL",
					true);
#endif

			psaver->SetUnit(unit);
			psaver->Save();
			psaver->SaveLeafBlock("CellID");
		}

		void BCMFileSaverPrint(int step) {
#ifdef _REAL_IS_DOUBLE_
			psaver->SaveLeafBlock("P64", step);
			psaver->SaveLeafBlock("Vel64", step);
#else
			psaver->SaveLeafBlock("P32", step);
			psaver->SaveLeafBlock("Vel32", step);
#endif
		}

		void WriteBasicVariablesInSILOFormat(
				int step,
				int difflevel,
				RootGrid* rootGrid,
				BCMOctree* tree,
				Partition* partition) {
			/*
				 SiloWriter writer("data-flow", "mesh", false);
				 VtkWriter writer;
				 writer.writePUT<real>(
				 this->plsP0->GetID(),
				 this->plsUX0->GetID(),
				 this->plsUY0->GetID(),
				 this->plsUZ0->GetID(),
				 this->plsT0->GetID(),
				 this->vc,
				 g_pFFVConfig->OutputDataFormatOptionVTKPath,
				 g_pFFVConfig->OutputDataFormatOptionVTKPrefix,
				 string("flow"),
				 step,
				 difflevel,
				 rootGrid,
				 tree,
				 partition,
				 g_pFFVConfig->RootBlockOrigin,
				 g_pFFVConfig->RootBlockLength);

				 writer.writeVtkOverlappingAMR_PUT<real>(
				 this->plsP0->GetID(),
				 this->plsUX0->GetID(),
				 this->plsUY0->GetID(),
				 this->plsUZ0->GetID(),
				 this->plsT0->GetID(),
				 this->vc,
				 g_pFFVConfig->OutputDataFormatOptionVTKPath,
				 g_pFFVConfig->OutputDataFormatOptionVTKPrefix,
				 string("flow"),
				 step,
				 difflevel,
				 rootGrid,
				 tree,
				 partition,
				 g_pFFVConfig->RootBlockOrigin,
				 g_pFFVConfig->RootBlockLength);
			 */
		}

	private:
		double GetTime();
};

#endif

