#ifndef FFVCONFIG_H
#define FFVCONFIG_H

#include <string>

#include <BCMTools.h>
#include <PolygonBBoxDivider.h>
#include <TextParser.h>
#include <Vec3.h>

#include "real.h"
#include "FFVBC.h"

//3d vector value of double
typedef ::Vec3<double> Vec3d;

//PhysicalParametersForFluid
typedef struct _PPF {
	double rho;
	double cp;
	double k;
	double mu;
	std::string color;
}PPF;

//PhysicalParametersForSolid
typedef struct _PPS {
	double rho;
	double cp;
	double k;
	std::string color;
}PPS;

class FFVConfig {
public:
	FFVConfig();
	~FFVConfig();

private:
	TextParser* tp;

public:
//ApplicationContorl
	bool					CheckParameter;
	std::string   OperatorName;
	Vec3d					FillingOrigin;
	std::string		FillingMedium;
	std::string		OperationMode;

//ConvectionTerm
	std::string		ConvectionTermScheme;

//DomainInfo
	Vec3d         RootBlockOrigin;
	Vec3i         RootBlockGrid;
	double        RootBlockLength;
	bool          RootBlockPeriodicX;
	bool          RootBlockPeriodicY;
	bool          RootBlockPeriodicZ;
	int           LeafBlockNumberOfCells;
	int           LeafBlockNumberOfVirtualCells;
	int           LeafBlockNumberOfMarginalCells;
	std::string		TreeType;
	int           TreeMinLevel;
	int           TreeMaxLevel;
	std::vector<PolygonGroupSpec>
								PolygonGroupList;
	std::vector<BoundingBoxSpec>
								BoundingBoxList;
	std::vector<BoundingBoxSpec>
								SphericalBoxList;
	Vec3d         TreeDividerCenter;
	double        TreeDividerRadius;
	double        TreeDividerDeltaR;
	Vec3d         TreeDividerBBOrigin;
	Vec3d         TreeDividerBBEnd;
	bool          TreeDividerHollow;

//GeometryModel
	std::string		PolylibConfig;

//Iteration
	std::string		IterationSolverP;
	int						IterationMaxCountP;
	double				IterationEpsilonP;
	double				IterationOmegaP;
	int						IterationPreCountP;

	std::string		IterationSolverU;
	int						IterationMaxCountU;
	double				IterationEpsilonU;
	double				IterationOmegaU;
	int						IterationPreCountU;

	std::string		IterationSolverT;
	int						IterationMaxCountT;
	double				IterationEpsilonT;
	double				IterationOmegaT;
	int						IterationPreCountT;

	bool					IterationReferencePressureActive;
	Vec3d					IterationReferencePressurePoint;
	double				IterationReferencePressureValue;

//MediumTable
	std::vector<PPF> MediumTableFluid;
	std::vector<PPS> MediumTableSolid;

//Output
	bool					OutputLogBase;
	bool					OutputLogBlock;
	bool					OutputLogLaptime;
	bool					OutputLogIteration;
	bool					OutputLogProfiling;
	bool					OutputLogStatistics;
	bool					OutputLogForce;
	std::string		OutputLogFilenameBase;
	std::string		OutputLogFilenameProfiling;
	std::string		OutputLogFilenameLaptime;
	std::string		OutputLogFilenameIteration;
	std::string		OutputLogFilenameBlock;
	std::string		OutputLogFilenameStatistics;
	std::string		OutputLogFileIntervalType;
	int						OutputLogFileIntervalI;
	double				OutputLogFileIntervalD;

	int						OutputLogFileIntervalBase;
	int						OutputLogFileIntervalBlock;
	int						OutputLogFileIntervalLaptime;
	int						OutputLogFileIntervalIteration;
	int						OutputLogFileIntervalProfiling;
	int						OutputLogFileIntervalStatistics;
	int						OutputLogFileIntervalForce;

	std::string		OutputLogConsoleIntervalType;
	int						OutputLogConsoleIntervalI;
	double				OutputLogConsoleIntervalD;
	std::string		OutputDataBasicVariablesFormat;
	std::string		OutputDataBasicVariablesTemporalType;
	int						OutputDataBasicVariablesIntervalI;
	double				OutputDataBasicVariablesIntervalD;
	std::string		OutputDataDerivedVariablesFormat;
	std::string		OutputDataDerivedVariablesTemporalType;
	int						OutputDataDerivedVariablesIntervalI;
	double				OutputDataDerivedVariablesIntervalD;
	bool					OutputDataDerivedVariablesVorticity;
	bool					OutputDataDerivedVariablesHelicity;
	bool					OutputDataDerivedVariablesQcriterion;
	bool					OutputDataDerivedVariablesForce;
	std::string		OutputDataFormatOptionPLOT3DPath;
	std::string		OutputDataFormatOptionPLOT3DPrefix;
	std::string		OutputDataFormatOptionVTKPath;
	std::string		OutputDataFormatOptionVTKPrefix;

	bool					OutputDataBasicVariablesFormatVTK;
	bool					OutputDataDerivedVariablesFormatVTK;
	bool					OutputDataBasicVariablesFormatPLOT3D;
	bool					OutputDataDerivedVariablesFormatPLOT3D;
	bool					OutputDataBasicVariablesFormatBCM;
	bool					OutputDataDerivedVariablesFormatBCM;

//PhysicalParameter
	double				GravityX;
	double				GravityY;
	double				GravityZ;

//ShapeApproximation
	std::string		ShapeApproximationMethod;
	double				ShapeApproximationCutoff;
	bool					ShapeApproximationVoxelization;
	bool					ShapeApproximationSymmetrization;

//SolvingMethod
	std::string		TimeIntegrationMethodForFlow;

//StartCondition
	std::string		RestartInputPath;
	std::string		RestartOutputPath;
	std::string		RestartPrefix;
	int						RestartInterval;
	double				InitialValueP;
	Vec3d					InitialValueU;
	double				InitialValueT;
	double				InitialValueDP;
	Vec3d					InitialValueDU;
	double				InitialValueDT;

//TimeControl
	std::string		TimeControlAccelerationTemporalType;
	int						TimeControlAccelerationAcceleratingTimeI;
	double				TimeControlAccelerationAcceleratingTimeD;
	std::string		TimeControlTimeStepMode;
	double				TimeControlTimeStepDeltaT;
	std::string		TimeControlSessionTemporalType;
	int						TimeControlSessionStartI;
	int						TimeControlSessionEndI;
	double				TimeControlSessionStartD;
	double				TimeControlSessionEndD;

//GridGeneration
	bool					GridGenerationHoleFilling;
	bool					GridGenerationHoleFilling2;
	bool					GridGenerationOutputSTL;

//Tuning
	bool					TuningMasking;
	std::string   TuningBlockOrdering;
	std::string   TuningVCUpdate;

//BC
	std::string		BCOuterBoundaryFaceBCXminus;
	std::string		BCOuterBoundaryFaceBCXplus;
	std::string		BCOuterBoundaryFaceBCYminus;
	std::string		BCOuterBoundaryFaceBCYplus;
	std::string		BCOuterBoundaryFaceBCZminus;
	std::string		BCOuterBoundaryFaceBCZplus;

	std::vector<OBC> OuterBCP;
	std::vector<OBC> OuterBCUX;
	std::vector<OBC> OuterBCUY;
	std::vector<OBC> OuterBCUZ;
	std::vector<OBC> OuterBCT;

	void GetOuterBoundary(std::string FaceId, OBC& obcP, OBC& obcUX, OBC& obcUY, OBC& obcUZ, OBC& obcT);

public:
	void Load(std::string filename);
	void Check();

private:
	int GetNumOfNodes(const std::string& label) {
		tp->changeNode(label);

		std::vector<std::string> nodes;
		tp->getNodes(nodes);

		return nodes.size();
	}

	int GetNumOfLabels(const std::string& label) {
		tp->changeNode(label);

		std::vector<std::string> labels;
		tp->getLabels(labels);

		return labels.size();
	}

	template<class T>
	T convertToT(const std::string& value);

	template<class T>
	T Read(const std::string& label) {
		if( tp->chkLabel(label) != true ) {
			std::cout << "Label not found: ";
			std::cout << label << std::endl;
			Exit(EX_FAILURE);
		}
		std::string value;
		TextParserError error = tp->getValue(label, value);
		if( error != 0 ) {
			Exit(EX_FAILURE);
		}
		return convertToT<T>(value);
	}

	template<class T>
	T Read(const std::string& label, const T& default_value) {
		if( tp->chkLabel(label) != true ) {
			return default_value;
		}
		std::string value;
		TextParserError error = tp->getValue(label, value);
		if( error != 0 ) {
			Exit(EX_FAILURE);
		}
		return convertToT<T>(value);
	}
};

template<class T>
T FFVConfig::convertToT(const std::string& value) {
	T t;
	std::istringstream ist(value);
	ist >> t;
	return t;
}

template<>
inline std::string FFVConfig::convertToT<std::string>(const std::string& value) {
	return value;
}

template<>
inline bool FFVConfig::convertToT<bool>(const std::string& value) {
	int error = 0;
	bool b = tp->convertBool(value, &error);
	return b;
}

template<>
inline Vec3i FFVConfig::convertToT<Vec3i>(const std::string& value) {
	std::vector<std::string> svec;
	tp->splitVector(value, svec);
	int error = 0;
	Vec3i v;
	v.x = tp->convertInt(svec[0], &error);
	v.y = tp->convertInt(svec[1], &error);
	v.z = tp->convertInt(svec[2], &error);
	return v;
}

template<>
inline Vec3d FFVConfig::convertToT<Vec3d>(const std::string& value) {
	std::vector<std::string> svec;
	tp->splitVector(value, svec);
	int error = 0;
	Vec3d v;
	v.x = tp->convertDouble(svec[0], &error);
	v.y = tp->convertDouble(svec[1], &error);
	v.z = tp->convertDouble(svec[2], &error);
	return v;
}

#endif

