include(Visualization/CMakeLists.txt)
add_definitions(-DPBD_DATA_PATH="../data")

# search all demos
set(PBD_DEMOS 
	# BarDemo
	# ClothDemo
	# CosseratRodsDemo
	# CouplingDemos
	# DistanceFieldDemos
	FluidDemo
	# GenericConstraintsDemos
	# PositionBasedElasticRodsDemo
	# RigidBodyDemos
	# SceneLoaderDemo
)

add_subdirectory(Simulation)
add_subdirectory(Utils)

if (NOT PBD_NO_DEMOS)
	foreach (_demo_name ${PBD_DEMOS})
		option(Build_${_demo_name} "Build ${_demo_name}"	ON)
		if (Build_${_demo_name})
			add_subdirectory (${CMAKE_CURRENT_SOURCE_DIR}/${_demo_name})
		endif (Build_${_demo_name})
	endforeach ()
endif()

install(DIRECTORY ./Common
DESTINATION include/Demos
FILES_MATCHING PATTERN "*.h")

install(DIRECTORY ./Visualization
DESTINATION include/Demos
FILES_MATCHING PATTERN "*.h")