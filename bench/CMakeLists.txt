# Specify the (private) source for the -bench executables

target_sources(PDE-dense-bench
    PRIVATE
	# Headers
	dense.h

        # Sources
        PDE-dense-bench.cpp
)

target_sources(PDE-sparse-man-bench
    PRIVATE
	# Headers
	sparse-man.h

        # Sources
        PDE-sparse-man-bench.cpp
)