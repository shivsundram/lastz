include make-include.mak

default: build_lastz

lastz_32: build_lastz_32

lastz_40: build_lastz_40

#---------
# builds/installation
#---------

build: build_lastz

build_lastz:
	cd src && ${MAKE} lastz lastz_D

build_lastz_32:
	cd src && ${MAKE} lastz_32

build_lastz_40:
	cd src && ${MAKE} lastz_40

# Stage-timer-instrumented variant. Built side-by-side with the un-instrumented
# `lastz`. Prints a per-stage timing report on completion (to stderr by default;
# set the LASTZ_STAGE_REPORT env var to redirect to a file).
build_lastz_timed:
	cd src && ${MAKE} lastz_T
lastz_T: build_lastz_timed

# Substage variant: same as lastz_T plus per-iteration rdtsc timers around the
# seed_hit_search inner loop (chain walk vs. processor callback) and inside
# process_for_simple_hit (dedup short-circuit vs. xdrop_extend_seed_hit vs.
# HSP reporter). Adds a couple percent overhead via per-iteration rdtsc — use
# lastz_T for absolute wall-time, lastz_TS for substage breakdown.
build_lastz_substages:
	cd src && ${MAKE} lastz_TS
lastz_TS: build_lastz_substages

build_test_version:
	cd src && ${MAKE} lastz-test lastz_D-test

install: install_lastz

install_lastz:
	cd src && ${MAKE} install

install_32:
	cd src && ${MAKE} install_32

install_40:
	cd src && ${MAKE} install_40

install_test_version:
	cd src && ${MAKE} install_test_version

# cleanup

clean:
	cd src && ${MAKE} clean

cleano:
	cd src && ${MAKE} cleano

#---------
# testing
#
# test:
#	A small test to give some comfort level that the program has built properly,
#	or that changes you've made to the source code haven't broken it. If the
#	test succeeds, there will be no output from the diff.
# base_tests:
#	More extensive tests (but still small). The results should be of this form:
#	SUCCESS: ../test_data/xxx and ../test_results/yyy are equivalent
#---------

test:
	cd src && ${MAKE} test

base_tests:
	cd src && ${MAKE} base_tests

clean_test: clean_tests

clean_tests:
	cd src && ${MAKE} clean_tests

