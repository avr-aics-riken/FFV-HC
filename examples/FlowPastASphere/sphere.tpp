GeometryModel
{
	PolylibFile						= "sphere.tpl"
}

DomainInfo
{
	RootBlock
	{
		Origin							= (-7.5, -7.5, -7.5)
		Grid								= (1, 1, 1)
		Length							= 15.0
		PeriodicX						= "false"
		PeriodicY						= "false"
		PeriodicZ						= "false"
	}
	LeafBlock
	{
		NumberOfCells				= 8
	}
	Tree
	{
		Type								= "polygon"
		MinLevel						= 0
		MaxLevel						= 0
		PolygonGroupList
		{
			PolygonGroup[@]
			{
				Name						= "sphere"
				Level						= 5
			}
		}
		BoundingBoxList
		{
/*
			BoundingBox[@]
			{
				Origin					= (0.0, -0.55, -0.55)
				End							= (1.0, 0.55, 0.55)
				Level						= 5
			}
			BoundingBox[@]
			{
				Origin					= (-0.6, -0.6, -0.6)
				End							= (4.0, 0.6, 0.6)
				Level						= 4 
			}
			BoundingBox[@]
			{
				Origin					= (-0.6, -0.6, -0.6)
				End							= (8.0, 0.6, 0.6)
				Level						= 3
			}
			BoundingBox[@]
			{
				Origin					= (-0.6, -0.6, -0.6)
				End							= (16.0, 0.6, 0.6)
				Level						= 2
			}
*/
		}
	}
}

BCTable
{
	OuterBoundary
	{
		FaceBC
		{
			Xminus						= "Inlet"
			Xplus							= "Outlet"
			Yminus						= "Neumann"
			Yplus							= "Neumann"
			Zminus						= "Neumann"
			Zplus							= "Neumann"
		}

		Inlet
		{
			Class							= "Direct"
			TypeP							= "Neumann"
			TypeUX						= "Dirichlet"
			TypeUY						= "Dirichlet"
			TypeUZ						= "Dirichlet"
			TypeT							= "Dirichlet"
			ValueP						= 0.0
			ValueUX						= 1.0
			ValueUY						= 0.0
			ValueUZ						= 0.0
			ValueT						= 0.0
		}
		Outlet
		{
			Class							= "Direct"
			TypeP							= "Dirichlet"
			TypeUX						= "Neumann"
			TypeUY						= "Neumann"
			TypeUZ						= "Neumann"
			TypeT							= "Neumann"
			ValueP						= 0.0
			ValueUX						= 0.0
			ValueUY						= 0.0
			ValueUZ						= 0.0
			ValueT						= 0.0
		}
		Neumann
		{
			Class							= "Direct"
			TypeP							= "Neumann"
			TypeUX						= "Neumann"
			TypeUY						= "Neumann"
			TypeUZ						= "Neumann"
			TypeT							= "Neumann"
			ValueP						= 0.0
			ValueUX						= 0.0
			ValueUY						= 0.0
			ValueUZ						= 0.0
			ValueT						= 0.0
		}

	}
}

TimeControl
{
	Acceleration
	{
		TemporalType				= "step"
		AcceleratingTime		= 0
	}
	TimeStep
	{
		Mode								= "direct"
		DeltaT							= 0.005859375
	}
	Session
	{
		TemporalType				= "step"
		Start								= 0
		End									= 1000
	}
}

ApplicationControl
{
	CheckParameter				= "false"
	Operator							= "Junya_Onishi"	
	Filling
	{
		Origin								= (-7.5, -7.5, -7.5)
		Medium								= "fluid0"
	}
}

MediumTable
{
	fluid0
	{
		State								= "fluid"
		MassDensity					= 1.0
		SpecificHeat				= 1.0
		ThermalConductivity	= 0.01
		Viscosity						= 0.01
		Color								= "0000FF"
	}
}



ConvectionTerm
{
	Scheme								= "C2"
}

Iteration
{
	LinearSolver[@]
	{
		Alias								= "pbicgstabU"
		Class								= "pbicgstab"
		MaxIteration				= 1000
		ConvergenceCriterion
												= 1.0e-5
		Omega								= 1.0
		PreCount						= 2
	}
	LinearSolver[@]
	{
		Alias								= "pbicgstabP"
		Class								= "pbicgstab"
		MaxIteration				= 1000
		ConvergenceCriterion
												= 1.0e-5
		Omega								= 1.0
		PreCount						= 2
	}
	Pressure							= "pbicgstabP"
	Velocity							= "pbicgstabU"
	Temperature						= "pbicgstabU"
}

Output
{
	Log
	{
		Base								= "true"
		Iteration						= "true"
		Profiling						= "true"
		Block								= "true"
		Laptime							= "true"
		Statistics					= "true"
		Force								= "true"
		FilenameBase				= "log-base.txt"
		FilenameIteration		= "log-iter.txt"
		FilenameProfiling		= "log-prof.txt"
		FilenameBlock				= "log-block.txt"
		FilenameLaptime			= "log-lap.txt"
		FilenameStatistics	= "log-stats.txt"
		History	
		{
			TemporalType			= "step"
			Interval					= 1
		}
		Console
		{
			TemporalType			= "step"
			Interval					= 1
		}
	}
	Data
	{
		BasicVariables
		{
			Format						= "VTK"
			TemporalType			= "step"
			Interval					= 100
		}
		DerivedVariables
		{
			Format						= "VTK"
			TemporalType			= "step"
			Interval					= 100
			Vorticity					= "false"
			Helicity					= "false"
			Qcriterion				= "true"
			Force							= "false"
		}
		FormatOption
		{
			PLOT3D {
				Path						= "PLOT3D"
				Prefix					= "data-"
			}
			VTK {
				Path						= "VTK"
				Prefix					= "data-"
			}
		}
	}
}

ShapeApproximation
{
	Method								= "cut"
	Cutoff								= 0.01
}

SolvingMethod {
//	Flow									= "implicit"
	Flow									= "explicit"
}

StartCondition
{
	InitialState
	{
		Pressure						= 0.0
		Velocity						= (0.0, 0.0, 0.0)
		Temperature					= 0.0

		DPressure						= 0.0
		DVelocity						= (0.0, 0.01, 0.0)
		DTemperature					= 0.0
	}
	Restart
	{
		InputPath						= "BIN"
		OutputPath					= "BIN"
		Prefix							= "dump-"
		Interval						= "1000"
	}
}

GridGeneration
{
	OutputSTL							= "true"
}

