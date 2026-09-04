# Distributed arrays module (DDI-compatible)
target_sources(${main_lib} PRIVATE darrays_types.F90 darrays_distrib.F90
                                   darrays_core.F90 darrays_dlb.F90 darrays.F90)
