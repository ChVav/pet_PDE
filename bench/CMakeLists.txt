# Specify the (private) source for the -bench executables

target_sources(PDE-dense-bench
    PRIVATE
	# Headers
	common.h
	dense.h

        # Sources
	common.cpp
        PDE-dense-bench.cpp
)

target_sources(PDE-sparse-man-bench
    PRIVATE
	# Headers
	common.h
	sparse-man.h

        # Sources
	common.cpp
        PDE-sparse-man-bench.cpp
)

target_sources(PDE-sparse-csr-bench
    PRIVATE
	# Headers
	common.h
	sparse-csr.h

        # Sources
	common.cpp
        PDE-sparse-csr-bench.cpp
)

target_sources(PDE-sparse-eigen-bench
    PRIVATE
	# Headers
	common.h
	sparse-eigen.h

        # Sources
	common.cpp
        PDE-sparse-eigen-bench.cpp
)