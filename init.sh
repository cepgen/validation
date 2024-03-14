if [ -z $CEPGEN_PATH ]; then
  echo "Variable \$CEPGEN_PATH is unset.";
  return;
else
  echo "Variable \$CEPGEN_PATH is set to '$CEPGEN_PATH'."
fi

cat << EOF > CepGenEnvironment.h
#ifndef CepGenEnvironment_h
#define CepGenEnvironment_h

#define CEPGEN_PATH "$CEPGEN_PATH"

// clang-format off
R__ADD_INCLUDE_PATH($CEPGEN_PATH)
R__ADD_LIBRARY_PATH($CEPGEN_PATH/build)
R__LOAD_LIBRARY(libCepGen.so)
// clang-format on

#endif
EOF

if [ ! -f a08tmc.dat ]; then
  ln -s $CEPGEN_PATH/External/a08tmc.dat a08tmc.dat
fi
if [ ! -f mstw_sf_scan_nnlo.dat ]; then
  ln -s $CEPGEN_PATH/External/mstw_sf_scan_nnlo.dat mstw_sf_scan_nnlo.dat
fi
