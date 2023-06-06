# Specify the (private) sources of the PDEsolver target
target_sources(PDEsolver-dense
    PRIVATE
        # Headers
        "common.h"

        # Sources
        "common.cpp"
        "dense.cpp"
)