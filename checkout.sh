echo " "
echo "Checking out packages ..."
echo " "
pkgco.py DerivationFrameworkFlavourTag-00-01-55
pkgco.py JetTagTools-01-00-96
pkgco.py BTagging-00-07-64
pkgco.py ParticleJetTools
pkgco.py AODFix

echo " "
echo "Applying patches ..."
echo " "
cp patches/FTAG5.py PhysicsAnalysis/DerivationFramework/DerivationFrameworkFlavourTag/share
cp patches/FlavourTagCommon.py PhysicsAnalysis/DerivationFramework/DerivationFrameworkFlavourTag/python

cp patches/KmeansbbTagTool.h PhysicsAnalysis/JetTagging/JetTagTools/JetTagTools
cp patches/KmeansbbTagTool.cxx PhysicsAnalysis/JetTagging/JetTagTools/src
cp patches/ExKtbbTagTool.h PhysicsAnalysis/JetTagging/JetTagTools/JetTagTools
cp patches/ExKtbbTagTool.cxx PhysicsAnalysis/JetTagging/JetTagTools/src
cp patches/ExKtbbTag.cxx PhysicsAnalysis/JetTagging/JetTagTools/src
cp patches/JetTagTools_entries.cxx PhysicsAnalysis/JetTagging/JetTagTools/src/components
cp patches/requirements PhysicsAnalysis/JetTagging/JetTagTools/cmt

cp patches/BTaggingConfiguration_LoadTools.py PhysicsAnalysis/JetTagging/JetTagAlgs/BTagging/python

cp patches/ParticleToJetAssociator.cxx PhysicsAnalysis/AnalysisCommon/ParticleJetTools/Root
cp patches/ParticleToJetAssociator.h PhysicsAnalysis/AnalysisCommon/ParticleJetTools/ParticleJetTools

cp patches/AODFix_r207.py Reconstruction/AODFix/python

