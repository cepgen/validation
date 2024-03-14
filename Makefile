ROOT_BIN=root -l -b -q

.PHONY: all

all: chi2_test

chi2_test: chi2_test.C
	${ROOT_BIN} '$^(0, 0, 1, "$@")'
