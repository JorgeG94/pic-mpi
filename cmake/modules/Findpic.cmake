set(_lib "pic")
set(_url "https://github.com/JorgeG94/pic/")

# Use a specific tag set(_rev "v0.5.0") my_fetch_package("${_lib}" "${_url}"
# "${_rev}")

include("${CMAKE_CURRENT_LIST_DIR}/sample_utils.cmake")

# Or use a branch
set(_rev "feat/move_blas_out")
my_fetch_package("${_lib}" "${_url}" "${_rev}")

# Or use a commit hash set(_rev "abc123def456") my_fetch_package("${_lib}"
# "${_url}" "${_rev}")

# Or omit to use HEAD (default branch) my_fetch_package("${_lib}" "${_url}")

unset(_lib)
unset(_url)
unset(_rev)
