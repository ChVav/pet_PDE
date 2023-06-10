# Specify the (private) sources of the PDEsolver target
target_sources(PDEsolver-dense
    PRIVATE
        # Headers
        "common.h"
	"dense.h"

        # Sources
        "common.cpp"
        "dense.cpp"
)

target_sources(PDEsolver-sparse-man
    PRIVATE
        # Headers
	"common.h"
        "sparse-man.h"

        # Sources
        "common.cpp"
	"sparse-man.cpp"
)