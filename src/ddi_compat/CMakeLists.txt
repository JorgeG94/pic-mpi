# DDI compatibility wrapper for GAMESS Files must be listed in dependency order:
# 1. ddi_impl.f90 - implementation module (no dependencies on other ddi files)
# 2. ddi_compat.f90 - backward-compatible module interface (depends on ddi_impl)
# 3. ddi_api.f90 - external bind(c) wrappers for GAMESS linkage (depends on
# ddi_impl)
target_sources(${main_lib} PRIVATE ddi_impl.f90 ddi_compat.f90 ddi_api.f90)
