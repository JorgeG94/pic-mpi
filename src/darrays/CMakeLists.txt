# Distributed arrays module (DDI-compatible)
target_sources(${main_lib} PRIVATE darrays_types.f90 darrays_distrib.f90
                                   darrays_core.f90 darrays_dlb.f90 darrays.f90)
