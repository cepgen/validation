OUTPUT_DIR=output
ROOT_BIN=root -l -b -q

.PHONY: all

all: chi2_test compare_datasets compare_sf_vs_xbj compare_sf_vs_w

chi2_test: chi2_test.C | output
	${ROOT_BIN} '$^(0, 0, 1, "${OUTPUT_DIR}/$@")' > /dev/null

compare_datasets: compare_datasets.C | output
	${ROOT_BIN} '$^(1, 1, 0, 0, "'${OUTPUT_DIR}'/datasets_comparison", {"pdf", "png"})' > /dev/null
	${ROOT_BIN} '$^(1, 1, 1, 0, "'${OUTPUT_DIR}'/datasets_comparison_w", {"pdf", "png"})' > /dev/null

compare_sf_vs_xbj: compare_sf_vs_xbj.C | output
	${ROOT_BIN} '$^(0.11, 1, 0, 0, 1, 0, "${OUTPUT_DIR}/$@_xbj_q2-0p11gev2", {"pdf", "png"})' > /dev/null
	${ROOT_BIN} '$^(0.425, 0, 0, 0, 1, 0, "${OUTPUT_DIR}/$@_xbj_q2-0p425gev2", {"pdf", "png"})' > /dev/null
	${ROOT_BIN} '$^(0.625, 0, 1, 0, 1, 0, "${OUTPUT_DIR}/$@_xbj_q2-0p625gev2", {"pdf", "png"})' > /dev/null
	${ROOT_BIN} '$^(0.75, 1, 0, 0, 1, 0, "${OUTPUT_DIR}/$@_xbj_q2-0p75gev2", {"pdf", "png"})' > /dev/null
	${ROOT_BIN} '$^(0.975, 0, 0, 0, 1, 0, "${OUTPUT_DIR}/$@_xbj_q2-0p975gev2", {"pdf", "png"})' > /dev/null
	${ROOT_BIN} '$^(1.225, 0, 0, 0, 1, 0, "${OUTPUT_DIR}/$@_xbj_q2-1p225gev2", {"pdf", "png"})' > /dev/null
	${ROOT_BIN} '$^(2.5, 1, 0, 0, 1, 0, "${OUTPUT_DIR}/$@_xbj_q2-2p5gev2", {"pdf", "png"})' > /dev/null
	${ROOT_BIN} '$^(5., 1, 0, 0, 1, 0, "${OUTPUT_DIR}/$@_xbj_q2-5gev2", {"pdf", "png"})' > /dev/null
	${ROOT_BIN} '$^(20., 1, 0, 0, 1, 0, "${OUTPUT_DIR}/$@_xbj_q2-20gev2", {"pdf", "png"})' > /dev/null
	${ROOT_BIN} '$^(50., 1, 0, 0, 1, 0, "${OUTPUT_DIR}/$@_xbj_q2-50gev2", {"pdf", "png"})' > /dev/null
	${ROOT_BIN} '$^(100., 0, 0, 0, 1, 0, "${OUTPUT_DIR}/$@_xbj_q2-100gev2", {"pdf", "png"})' > /dev/null
	${ROOT_BIN} '$^(1000., 1, 0, 0, 1, 0, "${OUTPUT_DIR}/$@_xbj_q2-1000gev2", {"pdf", "png"})' > /dev/null
	${ROOT_BIN} '$^(12000., 0, 0, 0, 1, 0, "${OUTPUT_DIR}/$@_xbj_q2-12000gev2", {"pdf", "png"})' > /dev/null
	${ROOT_BIN} '$^(20000., 0, 0, 0, 1, 0, "${OUTPUT_DIR}/$@_xbj_q2-20000gev2", {"pdf", "png"})' > /dev/null

compare_sf_vs_w: compare_sf_vs_w.C | output
	${ROOT_BIN} '$^(0.11, 0, 0, 0, 1, 0, "${OUTPUT_DIR}/$@_w_q2-0p11gev2", {"pdf", "png"})' > /dev/null
	${ROOT_BIN} '$^(0.425, 0, 0, 0, 1, 0, "${OUTPUT_DIR}/$@_w_q2-0p425gev2", {"pdf", "png"})' > /dev/null
	${ROOT_BIN} '$^(0.625, 1, 1, 0, 1, 0, "${OUTPUT_DIR}/$@_w_q2-0p625gev2", {"pdf", "png"})' > /dev/null
	${ROOT_BIN} '$^(0.75, 1, 0, 0, 1, 0, "${OUTPUT_DIR}/$@_w_q2-0p75gev2", {"pdf", "png"})' > /dev/null
	${ROOT_BIN} '$^(0.975, 1, 0, 0, 1, 0, "${OUTPUT_DIR}/$@_w_q2-0p975gev2", {"pdf", "png"})' > /dev/null
	${ROOT_BIN} '$^(1.225, 1, 0, 0, 1, 0, "${OUTPUT_DIR}/$@_w_q2-1p225gev2", {"pdf", "png"})' > /dev/null
	${ROOT_BIN} '$^(2.5, 1, 0, 0, 1, 0, "${OUTPUT_DIR}/$@_w_q2-2p5gev2", {"pdf", "png"})' > /dev/null
	${ROOT_BIN} '$^(5., 0, 0, 0, 1, 0, "${OUTPUT_DIR}/$@_w_q2-5gev2", {"pdf", "png"})' > /dev/null
	${ROOT_BIN} '$^(20., 1, 0, 0, 1, 0, "${OUTPUT_DIR}/$@_w_q2-20gev2", {"pdf", "png"})' > /dev/null
	${ROOT_BIN} '$^(50., 0, 0, 0, 1, 0, "${OUTPUT_DIR}/$@_w_q2-50gev2", {"pdf", "png"})' > /dev/null
	${ROOT_BIN} '$^(100., 0, 0, 0, 1, 0, "${OUTPUT_DIR}/$@_w_q2-100gev2", {"pdf", "png"})' > /dev/null
	${ROOT_BIN} '$^(1000., 0, 0, 0, 1, 0, "${OUTPUT_DIR}/$@_w_q2-1000gev2", {"pdf", "png"})' > /dev/null
	${ROOT_BIN} '$^(12000., 0, 1, 0, 1, 0, "${OUTPUT_DIR}/$@_w_q2-12000gev2", {"pdf", "png"})' > /dev/null
	${ROOT_BIN} '$^(20000., 0, 1, 0, 1, 0, "${OUTPUT_DIR}/$@_w_q2-20000gev2", {"pdf", "png"})' > /dev/null

output:
	mkdir -p ${OUTPUT_DIR}
