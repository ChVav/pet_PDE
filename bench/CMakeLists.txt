# Specify the (private) source for the -bench executables

target_sources(PDE-dense-bench
    PRIVATE
	# Headers
	dense.h

        # Sources
        PDE-dense-bench.cpp
)