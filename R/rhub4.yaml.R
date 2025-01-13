# rhubv2 configuration for testing on r-devel-linux-x86_64-debian-clang

version: 2

tasks:
  - platform: r-devel-linux-x86_64-debian-clang
timeout: 3600
env:
  _R_CHECK_FORCE_SUGGESTS_: false
_R_CHECK_CRAN_INCOMING_: true
skip:
  examples: false
tests: false
options:
  manual: false
build_vignettes: true
