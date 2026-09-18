set(_lib "vapaa")
set(_pkg "VAPAA")
set(_url "https://github.com/JorgeG94/vapaa")
# TEMPORARY -- testing branch, do not merge.
#
# Points at the branch of vapaa#3, which replaces vapaa's hardcoded -fPIE with
# the POSITION_INDEPENDENT_CODE property. Without it a consumer that links a
# shared core (metalquicha with MQC_SHARED_LIB, roundabout) fails at link time
# with "relocation R_X86_64_PC32 ... can not be used when making a shared
# object". Revert to "dev" once vapaa#3 has merged.
set(_rev "fix/position-independent-code")

include("${CMAKE_CURRENT_LIST_DIR}/sample_utils.cmake")

my_fetch_package("${_lib}" "${_url}" "${_rev}")

unset(_lib)
unset(_pkg)
unset(_url)
unset(_rev)
