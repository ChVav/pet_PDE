# Specify the (private) sources of the PDEsolver target
target_sources(PDEsolver
    PRIVATE
        # Headers
        "dense.h"

        # Sources
        "dense.cpp"
        "main.cpp"
)