OUTPUT_DIR=output
ROOT_BIN=root -l -b -q

.PHONY: all

all: chi2_test

chi2_test: chi2_test.C | output
	${ROOT_BIN} '$^(0, 0, 1, "${OUTPUT_DIR}/$@")'

compare_datasets: compare_datasets.C | output
	${ROOT_BIN} '$^(1, 1, 0, 0, "'${OUTPUT_DIR}'/datasets_comparison", {"pdf", "png"})'
	${ROOT_BIN} '$^(1, 1, 1, 0, "'${OUTPUT_DIR}'/datasets_comparison_w", {"pdf", "png"})'

output:
	mkdir -p ${OUTPUT_DIR}
