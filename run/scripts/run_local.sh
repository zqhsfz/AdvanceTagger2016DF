Reco_tf.py \
--preExec 'rec.doApplyAODFix.set_Value_and_Lock(True);from BTagging.BTaggingFlags import BTaggingFlags;BTaggingFlags.CalibrationTag = "BTagCalibRUN12-08-18"' \
--inputAODFile \
/afs/cern.ch/user/n/nwhallon/work/public/xAOD_samples/mc15_13TeV.301523.MadGraphPythia8EvtGen_A14NNPDF23LO_RS_G_hh_bbbb_c20_M2000.merge.AOD.e3820_s2608_s2183_r7772_r7676/AOD.08176783._000001.pool.root.1 \
--outputDAODFile FTAG5.pool.root \
--maxEvents 10 \
--reductionConf FTAG5
# --reductionConf FTAG5 | tee log_FTAG5.dat
