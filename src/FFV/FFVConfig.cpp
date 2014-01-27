#include "FFVConfig.h"

#include <BCMTools.h>
#include <TextParser.h>

#include <string.h>

FFVConfig::FFVConfig() :
	tp( TextParser::get_instance_singleton() ) {
}

FFVConfig::~FFVConfig() {
}

void FFVConfig::Load(std::string filename) {
	TextParserError TPError = tp->read(filename.c_str());
	if( TPError != 0 ) {
		Exit(EX_READ_CONFIG);
	}

//ApplicationContorl
	CheckParameter					= Read<bool>					("/ApplicationControl/CheckParameter");
	OperatorName						= Read<std::string>		("/ApplicationControl/Operator");
	FillingMedium						= Read<std::string>		("/ApplicationControl/Filling/Medium");
	FillingOrigin						= Read<Vec3d>					("/ApplicationControl/Filling/Origin");
	OperationMode						= "Normal";
	OperationMode						= Read<std::string>		("/ApplicationControl/OperationMode", "Normal");
	if( CheckParameter == true ) {
		OperationMode					= "CheckParameter";
	}

//ConvectionTerm
	ConvectionTermScheme		= Read<std::string>		("/ConvectionTerm/Scheme");

//DomainInfo
	RootBlockOrigin					= Read<Vec3d>					("/DomainInfo/RootBlock/Origin");
	RootBlockGrid						= Read<Vec3i>					("/DomainInfo/RootBlock/Grid");
	RootBlockLength					= Read<double>				("/DomainInfo/RootBlock/Length");
	RootBlockPeriodicX			= Read<bool>					("/DomainInfo/RootBlock/PeriodicX");
	RootBlockPeriodicY			= Read<bool>					("/DomainInfo/RootBlock/PeriodicY");
	RootBlockPeriodicZ			= Read<bool>					("/DomainInfo/RootBlock/PeriodicZ");
	LeafBlockNumberOfCells	= Read<int>						("/DomainInfo/LeafBlock/NumberOfCells");
	LeafBlockNumberOfVirtualCells
													= Read<int>						("/DomainInfo/LeafBlock/NumberOfVirtualCells", 2);
	LeafBlockNumberOfMarginalCells
													= Read<int>						("/DomainInfo/LeafBlock/NumberOfMarginalCells", LeafBlockNumberOfVirtualCells);
	TreeType								= Read<std::string>		("/DomainInfo/Tree/Type");
	TreeMinLevel						= Read<int>						("/DomainInfo/Tree/MinLevel");
	TreeMaxLevel						= Read<int>						("/DomainInfo/Tree/MaxLevel");

	if( !strcasecmp(TreeType.c_str(), "flat") ) {
	} else if( !strcasecmp(TreeType.c_str(), "simple") ) {
	} else if( !strcasecmp(TreeType.c_str(), "polygon") ) {
		std::vector<std::string> pglist;
		tp->getArrayLabels													("/DomainInfo/Tree/PolygonGroupList/PolygonGroup[@]", pglist);	
		PolygonGroupList.clear();
		for(int n=0; n<pglist.size(); n++) {
			std::string name		= Read<std::string>		(pglist[n] + "/Name");
			int level						= Read<int>						(pglist[n] + "/Level");
			PolygonGroupSpec pgs = {name, level};
			PolygonGroupList.push_back(pgs);
		}

		std::vector<std::string> bblist;
		tp->getArrayLabels													("/DomainInfo/Tree/BoundingBoxList/BoundingBox[@]", bblist);	
		BoundingBoxList.clear();
		for(int n=0; n<bblist.size(); n++) {
			Vec3d	origin				= Read<Vec3d>					(bblist[n] + "/Origin");
			Vec3d	end						= Read<Vec3d>					(bblist[n] + "/End");
			int level						= Read<int>						(bblist[n] + "/Level");
			BoundingBox bb			= BoundingBox(origin, end);
			BoundingBoxSpec bbs = {bb, level};
			BoundingBoxList.push_back(bbs);
		}
	} else if( !strcasecmp(TreeType.c_str(), "sphere_old") ) {
		TreeDividerCenter			= Read<Vec3d>					("/DomainInfo/Tree/DividerParams/Center");
		TreeDividerRadius			= Read<double>				("/DomainInfo/Tree/DividerParams/Radius");
		TreeDividerDeltaR			= Read<double>				("/DomainInfo/Tree/DividerParams/DeltaR");
		TreeDividerBBOrigin		= Read<Vec3d>					("/DomainInfo/Tree/DividerParams/BoundingBox/Origin");
		TreeDividerBBEnd			= Read<Vec3d>					("/DomainInfo/Tree/DividerParams/BoundingBox/End");
		TreeDividerHollow			= Read<bool>					("/DomainInfo/Tree/DividerParams/Hollow");

		std::vector<std::string> pglist;
		tp->getArrayLabels													("/DomainInfo/Tree/PolygonGroupList/PolygonGroup[@]", pglist);	
		PolygonGroupList.clear();
		for(int n=0; n<pglist.size(); n++) {
			std::string name		= Read<std::string>		(pglist[n] + "/Name");
			int level						= Read<int>						(pglist[n] + "/Level");
			PolygonGroupSpec pgs = {name, level};
			PolygonGroupList.push_back(pgs);
		}
	} else if( !strcasecmp(TreeType.c_str(), "sphere") ) {
		std::vector<std::string> pglist;
		tp->getArrayLabels													("/DomainInfo/Tree/PolygonGroupList/PolygonGroup[@]", pglist);	
		PolygonGroupList.clear();
		for(int n=0; n<pglist.size(); n++) {
			std::string name		= Read<std::string>		(pglist[n] + "/Name");
			int level						= Read<int>						(pglist[n] + "/Level");
			PolygonGroupSpec pgs = {name, level};
			PolygonGroupList.push_back(pgs);
		}

		std::vector<std::string> bblist;
		tp->getArrayLabels													("/DomainInfo/Tree/BoundingBoxList/BoundingBox[@]", bblist);	
		BoundingBoxList.clear();
		for(int n=0; n<bblist.size(); n++) {
			Vec3d	origin				= Read<Vec3d>					(bblist[n] + "/Origin");
			Vec3d	end						= Read<Vec3d>					(bblist[n] + "/End");
			int level						= Read<int>						(bblist[n] + "/Level");
			BoundingBox bb			= BoundingBox(origin, end);
			BoundingBoxSpec bbs = {bb, level};
			BoundingBoxList.push_back(bbs);
		}

		std::vector<std::string> sblist;
		tp->getArrayLabels													("/DomainInfo/Tree/SphericalBoxList/SphericalBox[@]", sblist);	
		SphericalBoxList.clear();
		for(int n=0; n<sblist.size(); n++) {
			Vec3d	origin				= Read<Vec3d>					(sblist[n] + "/Origin");
			Vec3d	end						= Read<Vec3d>					(sblist[n] + "/End");
			int level						= Read<int>						(sblist[n] + "/Level");
			BoundingBox bb			= BoundingBox(origin, end);
			BoundingBoxSpec bbs = {bb, level};
			SphericalBoxList.push_back(bbs);
		}
	} else {
	}

//GeometryModel
	PolylibConfig						= Read<std::string>		("/GeometryModel/PolylibFile");

//Iteration
	std::string lsp					= Read<std::string>		("/Iteration/Pressure");
	std::string lsu					= Read<std::string>		("/Iteration/Velocity");
	std::string lst					= Read<std::string>		("/Iteration/Temperature");

	std::vector<std::string> lslist;
	tp->getArrayLabels														("/Iteration/LinearSolver[@]", lslist);	
	for(int n=0; n<lslist.size(); n++) {
		std::string alias			= Read<std::string>		(lslist[n] + "/Alias");
		if( !strcasecmp(alias.c_str(), lsp.c_str()) ) {
			IterationSolverP		= Read<std::string>		(lslist[n] + "/Class");
			IterationMaxCountP	= Read<int>						(lslist[n] + "/MaxIteration");
			IterationEpsilonP		= Read<double>				(lslist[n] + "/ConvergenceCriterion");
			IterationOmegaP			= Read<double>				(lslist[n] + "/Omega", 1.0);
			IterationPreCountP	= Read<int>						(lslist[n] + "/PreCount", 0);
		}
		if( !strcasecmp(alias.c_str(), lsu.c_str() ) ) {
			IterationSolverU		= Read<std::string>		(lslist[n] + "/Class");
			IterationMaxCountU	= Read<int>						(lslist[n] + "/MaxIteration");
			IterationEpsilonU		= Read<double>				(lslist[n] + "/ConvergenceCriterion");
			IterationOmegaU			= Read<double>				(lslist[n] + "/Omega", 1.0);
			IterationPreCountU	= Read<int>						(lslist[n] + "/PreCount", 0);
		}
		if( !strcasecmp(alias.c_str(), lst.c_str() ) ) {
			IterationSolverT		= Read<std::string>		(lslist[n] + "/Class");
			IterationMaxCountT	= Read<int>						(lslist[n] + "/MaxIteration");
			IterationEpsilonT		= Read<double>				(lslist[n] + "/ConvergenceCriterion");
			IterationOmegaT			= Read<double>				(lslist[n] + "/Omega", 1.0);
			IterationPreCountT	= Read<int>						(lslist[n] + "/PreCount", 0);
		}
	}

	IterationReferencePressureActive
													= Read<bool>					("/Iteration/ReferencePressure/Active", false);
	if( IterationReferencePressureActive == true ) {
		IterationReferencePressurePoint
													= Read<Vec3d>					("/Iteration/ReferencePressure/Point");
		IterationReferencePressureValue
													= Read<double>				("/Iteration/ReferencePressure/Value");
	}

//MediumTable
	std::string FillingMediumState
													= Read<std::string>		("/MediumTable/" + FillingMedium + "/State");
	if( !strcasecmp(FillingMediumState.c_str(), "fluid") ) {
		PPF ppf0;
		ppf0.rho							= Read<double>				("/MediumTable/" + FillingMedium + "/MassDensity");
		ppf0.cp								= Read<double>				("/MediumTable/" + FillingMedium + "/SpecificHeat");
		ppf0.k								= Read<double>				("/MediumTable/" + FillingMedium + "/ThermalConductivity");
		ppf0.mu								= Read<double>				("/MediumTable/" + FillingMedium + "/Viscosity");
		ppf0.color						= Read<std::string>		("/MediumTable/" + FillingMedium + "/Color");
		MediumTableFluid.clear();
		MediumTableFluid.push_back(ppf0);
	}

//Output
	OutputLogBase						= Read<bool>					("/Output/Log/Base");
	OutputLogLaptime				= Read<bool>					("/Output/Log/Laptime");
	OutputLogIteration			= Read<bool>					("/Output/Log/Iteration");
	OutputLogStatistics			= Read<bool>					("/Output/Log/Statistics");
	OutputLogForce					= Read<bool>					("/Output/Log/Force", false);
	OutputLogHeatFlux				= Read<bool>					("/Output/Log/HeatFlux", false);
	OutputLogHeatFluxTargetID
													= Read<int>						("/Output/Log/HeatFluxTargetID", -1);
	OutputLogBlock					= Read<bool>					("/Output/Log/Block");
	OutputLogProfiling			= Read<bool>					("/Output/Log/Profiling");

	OutputLogFilenameBase		= Read<std::string>		("/Output/Log/FilenameBase");
	OutputLogFilenameProfiling
													= Read<std::string>		("/Output/Log/FilenameProfiling");
	OutputLogFilenameLaptime
													= Read<std::string>		("/Output/Log/FilenameLaptime");
	OutputLogFilenameIteration
													= Read<std::string>		("/Output/Log/FilenameIteration");
	OutputLogFilenameBlock	= Read<std::string>		("/Output/Log/FilenameBlock");
	OutputLogFilenameStatistics
													= Read<std::string>		("/Output/Log/FilenameStatistics");
	OutputLogFileIntervalType
													= Read<std::string>		("/Output/Log/History/TemporalType");
	if( !strcasecmp(OutputLogFileIntervalType.c_str(), "time") ) {
		OutputLogFileIntervalD
													= Read<double>				("/Output/Log/History/Interval");
	} else {
		OutputLogFileIntervalI
													= Read<int>						("/Output/Log/History/Interval");
	}
	OutputLogConsoleIntervalType
													= Read<std::string>		("/Output/Log/Console/TemporalType");
	if( !strcasecmp(OutputLogConsoleIntervalType.c_str(), "time") ) {
		OutputLogConsoleIntervalD
													= Read<double>				("/Output/Log/Console/Interval");
	} else {
		OutputLogConsoleIntervalI
													= Read<int>						("/Output/Log/Console/Interval");
	}

	OutputLogFileIntervalBase
													= Read<int>						("/Output/Log/Interval/Base", OutputLogFileIntervalI);
	OutputLogFileIntervalProfiling
													= Read<int>						("/Output/Log/Interval/Profiling", OutputLogFileIntervalI);
	OutputLogFileIntervalLaptime
													= Read<int>						("/Output/Log/Interval/Laptime", OutputLogFileIntervalI);
	OutputLogFileIntervalIteration
													= Read<int>						("/Output/Log/Interval/Iteration", OutputLogFileIntervalI);
	OutputLogFileIntervalBlock
													= Read<int>						("/Output/Log/Interval/Block", OutputLogFileIntervalI);
	OutputLogFileIntervalStatistics
													= Read<int>						("/Output/Log/Interval/Statistics", OutputLogFileIntervalI);
	OutputLogFileIntervalForce
													= Read<int>						("/Output/Log/Interval/Force", OutputLogFileIntervalI);
	OutputLogFileIntervalHeatFlux
													= Read<int>						("/Output/Log/Interval/HeatFlux", OutputLogFileIntervalI);

	OutputDataBasicVariablesFormat
													= Read<std::string>		("/Output/Data/BasicVariables/Format");
	OutputDataBasicVariablesTemporalType
													= Read<std::string>		("/Output/Data/BasicVariables/TemporalType");
	if( !strcasecmp(OutputDataBasicVariablesTemporalType.c_str(), "time") ) {
		OutputDataBasicVariablesIntervalD
													= Read<double>				("/Output/Data/BasicVariables/Interval");
	} else {
		OutputDataBasicVariablesIntervalI
													= Read<int>						("/Output/Data/BasicVariables/Interval");
	}

	OutputDataDerivedVariablesFormat
													= Read<std::string>		("/Output/Data/DerivedVariables/Format");
	OutputDataDerivedVariablesTemporalType
													= Read<std::string>		("/Output/Data/DerivedVariables/TemporalType");
	if( !strcasecmp(OutputDataDerivedVariablesTemporalType.c_str(), "time") ) {
		OutputDataDerivedVariablesIntervalD
													= Read<double>				("/Output/Data/DerivedVariables/Interval");
	} else {
		OutputDataDerivedVariablesIntervalI
													= Read<int>						("/Output/Data/DerivedVariables/Interval");
	}
	OutputDataDerivedVariablesVorticity
													= Read<bool>					("/Output/Data/DerivedVariables/Vorticity");
	OutputDataDerivedVariablesHelicity
													= Read<bool>					("/Output/Data/DerivedVariables/Helicity");
	OutputDataDerivedVariablesQcriterion
													= Read<bool>					("/Output/Data/DerivedVariables/Qcriterion");
	OutputDataDerivedVariablesForce
													= Read<bool>					("/Output/Data/DerivedVariables/Force");
	OutputDataDerivedVariablesHeatFlux
													= Read<bool>					("/Output/Data/DerivedVariables/HeatFlux", false);

	OutputDataBasicVariablesFormatVTK
													= false;
	OutputDataBasicVariablesFormatPLOT3D
													= false;
	OutputDataBasicVariablesFormatBCM
													= false;
	OutputDataBasicVariablesFormatSILO
													= false;
	if( !strcasecmp(OutputDataBasicVariablesFormat.c_str()  , "VTK") ) {
		OutputDataBasicVariablesFormatVTK
													= true;
	} else if( !strcasecmp(OutputDataBasicVariablesFormat.c_str()  , "PLOT3D") ) {
		OutputDataBasicVariablesFormatPLOT3D
													= true;
	} else if( !strcasecmp(OutputDataBasicVariablesFormat.c_str()  , "BCM") ) {
		OutputDataBasicVariablesFormatBCM
													= true;
	} else if( !strcasecmp(OutputDataBasicVariablesFormat.c_str()  , "SILO") ) {
		OutputDataBasicVariablesFormatSILO
													= true;
	}

	OutputDataDerivedVariablesFormatVTK
													= false;
	OutputDataDerivedVariablesFormatPLOT3D
													= false;
	OutputDataDerivedVariablesFormatBCM
													= false;
	OutputDataDerivedVariablesFormatSILO
													= false;
	if( !strcasecmp(OutputDataDerivedVariablesFormat.c_str()  , "VTK") ) {
		OutputDataDerivedVariablesFormatVTK
													= true;
	} else if( !strcasecmp(OutputDataDerivedVariablesFormat.c_str()  , "PLOT3D") ) {
		OutputDataDerivedVariablesFormatPLOT3D
													= true;
	} else if( !strcasecmp(OutputDataDerivedVariablesFormat.c_str()  , "BCM") ) {
		OutputDataDerivedVariablesFormatBCM
													= true;
	} else if( !strcasecmp(OutputDataDerivedVariablesFormat.c_str()  , "SILO") ) {
		OutputDataDerivedVariablesFormatSILO
													= true;
	}

	if( OutputDataBasicVariablesFormatVTK || OutputDataDerivedVariablesFormatVTK ) {
		OutputDataFormatOptionVTKPath
														= Read<std::string>		("/Output/Data/FormatOption/VTK/Path");
		OutputDataFormatOptionVTKPrefix
														= Read<std::string>		("/Output/Data/FormatOption/VTK/Prefix");
	} else if( OutputDataBasicVariablesFormatPLOT3D || OutputDataDerivedVariablesFormatPLOT3D ) {
		OutputDataFormatOptionPLOT3DPath
														= Read<std::string>		("/Output/Data/FormatOption/PLOT3D/Path");
		OutputDataFormatOptionPLOT3DPrefix
														= Read<std::string>		("/Output/Data/FormatOption/PLOT3D/Prefix");
	} else if( OutputDataBasicVariablesFormatBCM || OutputDataDerivedVariablesFormatBCM ) {
	} else if( OutputDataBasicVariablesFormatSILO || OutputDataDerivedVariablesFormatSILO ) {
	}

//PhysicalParameter
	GravityX								= Read<double>					("/PhysicalParameter/GravityX", 0.0);
	GravityY								= Read<double>					("/PhysicalParameter/GravityY", 0.0);
	GravityZ								= Read<double>					("/PhysicalParameter/GravityZ", 0.0);

//ShapeApproximation
	ShapeApproximationMethod
													= Read<std::string>		("/ShapeApproximation/Method");
	if( !strcasecmp(ShapeApproximationMethod.c_str(), "voxel") ) {
		ShapeApproximationVoxelization
													= true;
		ShapeApproximationSymmetrization
													= false;
		ShapeApproximationCutoff
													= 0.0;
	} else {
		ShapeApproximationVoxelization
													= Read<bool>					("/ShapeApproximation/Voxelization", false);
		ShapeApproximationSymmetrization
													= Read<bool>					("/ShapeApproximation/Symmetrization", false);
		ShapeApproximationCutoff
													= Read<double>				("/ShapeApproximation/Cutoff");
	}

//SolvingMethod
	TimeIntegrationMethodForFlow
													= Read<std::string>		("/SolvingMethod/Flow");

//StartCondition
	InitialValueP						= Read<double>				("/StartCondition/InitialState/Pressure");
	InitialValueU						= Read<Vec3d>					("/StartCondition/InitialState/Velocity");
	InitialValueT						= Read<double>				("/StartCondition/InitialState/Temperature");
	InitialValueDP					= Read<double>				("/StartCondition/InitialState/DPressure", 0.0);
	InitialValueDU					= Read<Vec3d>					("/StartCondition/InitialState/DVelocity", 0.0);
	InitialValueDT					= Read<double>				("/StartCondition/InitialState/DTemperature", 0.0);
	RestartInputPath				= Read<std::string>		("/StartCondition/Restart/InputPath");
	RestartOutputPath				= Read<std::string>		("/StartCondition/Restart/OutputPath");
	RestartPrefix						= Read<std::string>		("/StartCondition/Restart/Prefix");
	RestartInterval					= Read<int>						("/StartCondition/Restart/Interval");

//TimeControl
	TimeControlAccelerationTemporalType
													= Read<std::string>		("/TimeControl/Acceleration/TemporalType");
	if( !strcasecmp(TimeControlAccelerationTemporalType.c_str(), "time") ) {
		TimeControlAccelerationAcceleratingTimeD
													= Read<double>				("/TimeControl/Acceleration/AcceleratingTime");
	} else {
		TimeControlAccelerationAcceleratingTimeI
													= Read<int>						("/TimeControl/Acceleration/AcceleratingTime");
	}
	TimeControlTimeStepMode	= Read<std::string>		("/TimeControl/TimeStep/Mode");
	TimeControlTimeStepDeltaT
													= Read<double>				("/TimeControl/TimeStep/DeltaT");
	TimeControlSessionTemporalType
													= Read<std::string>		("/TimeControl/Session/TemporalType");
	if( !strcasecmp(TimeControlSessionTemporalType.c_str(), "time") ) {
		TimeControlSessionStartD
													= Read<double>				("/TimeControl/Session/Start");
		TimeControlSessionEndD
													= Read<double>				("/TimeControl/Session/End");
	} else {
		TimeControlSessionStartI
													= Read<int>						("/TimeControl/Session/Start");
		TimeControlSessionEndI
													= Read<int>						("/TimeControl/Session/End");
	}

//GridGeneration
	GridGenerationHoleFilling
													= Read<bool>					("/GridGeneration/HoleFilling", true);
	GridGenerationHoleFilling2
													= Read<bool>					("/GridGeneration/HoleFilling2", true);
	GridGenerationOutputSTL	= Read<bool>					("/GridGeneration/OutputSTL", false);


//Tuning
	TuningMasking						= Read<bool>					("/Tuning/Masking", true);
	TuningBlockOrdering			= Read<std::string>		("/Tuning/BlockOrdering", "Hilbert");
	TuningVCUpdate					= Read<std::string>		("/Tuning/VCUpdate", "AtOnce");

//BC
	OuterBCP.clear();
	OuterBCUX.clear();
	OuterBCUY.clear();
	OuterBCUZ.clear();
	OuterBCT.clear();
	OBC obcP;
	OBC obcUX;
	OBC obcUY;
	OBC obcUZ;
	OBC obcT;

	GetOuterBoundary("Xminus", obcP, obcUX, obcUY, obcUZ, obcT);
	OuterBCP.push_back(obcP);
	OuterBCUX.push_back(obcUX);
	OuterBCUY.push_back(obcUY);
	OuterBCUZ.push_back(obcUZ);
	OuterBCT.push_back(obcT);

	GetOuterBoundary("Xplus" , obcP, obcUX, obcUY, obcUZ, obcT);
	OuterBCP.push_back(obcP);
	OuterBCUX.push_back(obcUX);
	OuterBCUY.push_back(obcUY);
	OuterBCUZ.push_back(obcUZ);
	OuterBCT.push_back(obcT);

	GetOuterBoundary("Yminus", obcP, obcUX, obcUY, obcUZ, obcT);
	OuterBCP.push_back(obcP);
	OuterBCUX.push_back(obcUX);
	OuterBCUY.push_back(obcUY);
	OuterBCUZ.push_back(obcUZ);
	OuterBCT.push_back(obcT);

	GetOuterBoundary("Yplus" , obcP, obcUX, obcUY, obcUZ, obcT);
	OuterBCP.push_back(obcP);
	OuterBCUX.push_back(obcUX);
	OuterBCUY.push_back(obcUY);
	OuterBCUZ.push_back(obcUZ);
	OuterBCT.push_back(obcT);

	GetOuterBoundary("Zminus", obcP, obcUX, obcUY, obcUZ, obcT);
	OuterBCP.push_back(obcP);
	OuterBCUX.push_back(obcUX);
	OuterBCUY.push_back(obcUY);
	OuterBCUZ.push_back(obcUZ);
	OuterBCT.push_back(obcT);

	GetOuterBoundary("Zplus" , obcP, obcUX, obcUY, obcUZ, obcT);
	OuterBCP.push_back(obcP);
	OuterBCUX.push_back(obcUX);
	OuterBCUY.push_back(obcUY);
	OuterBCUZ.push_back(obcUZ);
	OuterBCT.push_back(obcT);
}

void FFVConfig::GetOuterBoundary(std::string FaceId, OBC& obcP, OBC& obcUX, OBC& obcUY, OBC& obcUZ, OBC& obcT) {
	std::string	BCId				= Read<std::string>		("/BCTable/OuterBoundary/FaceBC/" + FaceId);
	std::string BCClass			= Read<std::string>		("/BCTable/OuterBoundary/" + BCId + "/Class");

	std::string typeP				= "Neumann";
	std::string typeUX			= "Neumann";
	std::string typeUY			= "Neumann";
	std::string typeUZ			= "Neumann";
	std::string typeT				= "Neumann";
	double			valueP			= 0.0;
	double			valueUX			= 0.0;
	double			valueUY			= 0.0;
	double			valueUZ			= 0.0;
	double			valueT			= 0.0;
	double			xcP					= 0.0;
	double			xcUX				= 0.0;
	double			xcUY				= 0.0;
	double			xcUZ				= 0.0;
	double			xcT					= 0.0;
	double			ycP					= 0.0;
	double			ycUX				= 0.0;
	double			ycUY				= 0.0;
	double			ycUZ				= 0.0;
	double			ycT					= 0.0;
	double			zcP					= 0.0;
	double			zcUX				= 0.0;
	double			zcUY				= 0.0;
	double			zcUZ				= 0.0;
	double			zcT					= 0.0;
	double			rcP					= 0.0;
	double			rcUX				= 0.0;
	double			rcUY				= 0.0;
	double			rcUZ				= 0.0;
	double			rcT					= 0.0;
	double			lcP					= 0.0;
	double			lcUX				= 0.0;
	double			lcUY				= 0.0;
	double			lcUZ				= 0.0;
	double			lcT					= 0.0;
	double			hcP					= 0.0;
	double			hcUX				= 0.0;
	double			hcUY				= 0.0;
	double			hcUZ				= 0.0;
	double			hcT					= 0.0;
	double			wcP					= 0.0;
	double			wcUX				= 0.0;
	double			wcUY				= 0.0;
	double			wcUZ				= 0.0;
	double			wcT					= 0.0;
	if( !strcasecmp(BCClass.c_str(), "Direct") ) {
		typeP									= Read<std::string>		("/BCTable/OuterBoundary/" + BCId + "/TypeP");
		typeUX								= Read<std::string>		("/BCTable/OuterBoundary/" + BCId + "/TypeUX");
		typeUY								= Read<std::string>		("/BCTable/OuterBoundary/" + BCId + "/TypeUY");
		typeUZ								= Read<std::string>		("/BCTable/OuterBoundary/" + BCId + "/TypeUZ");
		typeT									= Read<std::string>		("/BCTable/OuterBoundary/" + BCId + "/TypeT");
		valueP								= Read<double>				("/BCTable/OuterBoundary/" + BCId + "/ValueP");
		valueUX								= Read<double>				("/BCTable/OuterBoundary/" + BCId + "/ValueUX");
		valueUY								= Read<double>				("/BCTable/OuterBoundary/" + BCId + "/ValueUY");
		valueUZ								= Read<double>				("/BCTable/OuterBoundary/" + BCId + "/ValueUZ");
		valueT								= Read<double>				("/BCTable/OuterBoundary/" + BCId + "/ValueT");
	} else if( !strcasecmp(BCClass.c_str(), "Inlet") ) {
		typeP									= "Neumann";
		typeUX								= "Dirichlet";
		typeUY								= "Dirichlet";
		typeUZ								= "Dirichlet";
		typeT									= "Neumann";
		valueP								= 0.0;
		valueUX								= Read<double>				("/BCTable/OuterBoundary/" + BCId + "/Velocity");
		valueUY								= 0.0;
		valueUZ								= 0.0;
		valueT								= 0.0;
	} else if( !strcasecmp(BCClass.c_str(), "Outlet") ) {
		typeP									= "Dirichlet";
		typeUX								= "Neumann";
		typeUY								= "Neumann";
		typeUZ								= "Neumann";
		typeT									= "Neumann";
		valueP								= Read<double>				("/BCTable/OuterBoundary/" + BCId + "/Pressure");
		valueUX								= 0.0;
		valueUY								= 0.0;
		valueUZ								= 0.0;
		valueT								= 0.0;
	} else if( !strcasecmp(BCClass.c_str(), "Slip") ) {
		typeP									= "Neumann";
		typeUX								= "Neumann";
		typeUY								= "Neumann";
		typeUZ								= "Neumann";
		typeT									= "Neumann";
		valueP								= 0.0;
		valueUX								= 0.0;
		valueUY								= 0.0;
		valueUZ								= 0.0;
		valueT								= 0.0;
	} else {
	}

	obcP.type								= GetBoundaryType(typeP);
	obcP.value							= valueP;
	obcP.xc									= xcP;
	obcP.yc									= ycP;
	obcP.zc									= zcP;
	obcP.rc									= rcP;
	obcP.lc									= lcP;
	obcP.hc									= hcP;
	obcP.wc									= wcP;

	obcUX.type							= GetBoundaryType(typeUX);
	obcUX.value							= valueUX;
	obcUX.xc								= xcUX;
	obcUX.yc								= ycUX;
	obcUX.zc								= zcUX;
	obcUX.rc								= rcUX;
	obcUX.lc								= lcUX;
	obcUX.hc								= hcUX;
	obcUX.wc								= wcUX;

	obcUY.type							= GetBoundaryType(typeUY);
	obcUY.value							= valueUY;
	obcUY.xc								= xcUY;
	obcUY.yc								= ycUY;
	obcUY.zc								= zcUY;
	obcUY.rc								= rcUY;
	obcUY.lc								= lcUY;
	obcUY.hc								= hcUY;
	obcUY.wc								= wcUY;

	obcUZ.type							= GetBoundaryType(typeUZ);
	obcUZ.value							= valueUZ;
	obcUZ.xc								= xcUZ;
	obcUZ.yc								= ycUZ;
	obcUZ.zc								= zcUZ;
	obcUZ.rc								= rcUZ;
	obcUZ.lc								= lcUZ;
	obcUZ.hc								= hcUZ;
	obcUZ.wc								= wcUZ;

	obcT.type								= GetBoundaryType(typeT);
	obcT.value							= valueT;
	obcT.xc									= xcT;
	obcT.yc									= ycT;
	obcT.zc									= zcT;
	obcT.rc									= rcT;
	obcT.lc									= lcT;
	obcT.hc									= hcT;
	obcT.wc									= wcT;

	for(int n=0; n<32; n++) {
		std::ostringstream bcid;
		bcid.width(2);
		bcid.setf(std::ios::fixed);
		bcid.fill('0');
		bcid << n;
		BCInternalBoundaryType[n]  = Read<int>   ("/BCTable/LocalBoundary/ID" + bcid.str() + "/Type", -1);
		BCInternalBoundaryValue[n] = Read<double>("/BCTable/LocalBoundary/ID" + bcid.str() + "/Value", 0.0);
	}
}

void FFVConfig::Check() {
	MPI::Comm& comm = MPI::COMM_WORLD;
	if( comm.Get_rank() != 0 ) {
		return;
	}

}

