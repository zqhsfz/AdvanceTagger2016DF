#====================================================================
# FTAG5.py
# It requires the reductionConf flag FTAG5 in Reco_tf.py
#====================================================================

# Set up common services and job objects
# This should appear in ALL derivation job options
from DerivationFrameworkCore.DerivationFrameworkMaster import *
from DerivationFrameworkInDet.InDetCommon import *
from DerivationFrameworkJetEtMiss.JetCommon import *
from DerivationFrameworkJetEtMiss.ExtendedJetCommon import *
from DerivationFrameworkJetEtMiss.METCommon import *
from DerivationFrameworkEGamma.EGammaCommon import *
from DerivationFrameworkMuons.MuonsCommon import *

from AthenaCommon.AthenaCommonFlags import jobproperties as jp
jp.AthenaCommonFlags.EvtMax.set_Value_and_Lock(-1)

#===================================================================
# Variable R track jets
#===================================================================
from DerivationFrameworkExotics.JetDefinitions import *
from JetRec.JetRecStandard import jtm

FTAG5Seq = CfgMgr.AthSequencer("FTAG5Sequence")

jtm.modifiersMap["smallvr_track_modifiers"] = jtm.modifiersMap["pv0track"]
jtm.modifiersMap["largevr_track_modifiers"] = [jtm.ktsplitter] # (nikola: what is this used for?)

jfind_smallvr_track = jtm.addJetFinder("AntiKtVR50Rmax4Rmin0TrackJets", "AntiKt", 0.4, "pv0track", "smallvr_track_modifiers",
                                        ghostArea = 0 , ptmin = 2000, ptminFilter = 7000,
                                        variableRMinRadius = 0, variableRMassScale = 50000, calibOpt = "none")

from JetRec.JetRecConf import JetAlgorithm
jetalg_smallvr_track= JetAlgorithm("jfind_smallvr_track", Tools = [jfind_smallvr_track])

FTAG5Seq += jetalg_smallvr_track

from JetRec.JetRecConf import PseudoJetGetter
jtm += PseudoJetGetter(
     "gvr50rmax4rmin0trackget", # give a unique name
     InputContainer = jetFlags.containerNamePrefix() + "AntiKtVR50Rmax4Rmin0TrackJets", # SG key
     Label = "GhostVR50Rmax4Rmin0TrackJet", # this is the name you'll use to retrieve ghost associated VR track jets
     OutputContainer = "PseudoJetGhostVR50Rmax4Rmin0TrackJet",
     SkipNegativeEnergy = True,
     GhostScale = 1.e-20, # this makes the PseudoJet Ghosts, and thus the reco flow will treat them as such
   )

jtm.gettersMap["lctopo"]+= [jtm.gvr50rmax4rmin0trackget] # has to happen before the trimmed jet gets clustered for the VR track jets to be ghost associated

#===================================================================
# Build Trimmed large-R jet
#===================================================================

# trim the large-R jets - the VR track jets will become ghost associated to the large-R jets and stored in the trimmed large-R jet collection (nikola: confirm)
addDefaultTrimmedJets(FTAG5Seq, "FTAG5")
applyJetCalibration_CustomColl("AntiKt10LCTopoTrimmedPtFrac5SmallR20", FTAG5Seq)

#===================================================================
# Build AntiKt R=1.0 TrackJet
#===================================================================

jfind_akt10trackjet = jtm.addJetFinder("AntiKt10TrackJets", "AntiKt", 1.0, "pv0track", ghostArea=0.00, ptmin=2000, ptminFilter=7000, calibOpt="none")
jetalg_akt10trackjet = JetAlgorithm("jfind_akt10trackjet", Tools = [jfind_akt10trackjet])

FTAG5Seq += jetalg_akt10trackjet

#===================================================================
# Build ExKt Subjets
#===================================================================

# make exkt subjet finding tool
def buildExclusiveSubjets(JetCollectionName, nsubjet, ToolSvc = ToolSvc):
    from JetSubStructureMomentTools.JetSubStructureMomentToolsConf import SubjetFinderTool
    from JetSubStructureMomentTools.JetSubStructureMomentToolsConf import SubjetRecorderTool

    SubjetContainerName = "%sExKt%iSubJets" % (JetCollectionName.replace("Jets", ""), nsubjet)

    subjetrecorder = SubjetRecorderTool("subjetrecorder%i_%s" % (nsubjet, JetCollectionName))
    ToolSvc += subjetrecorder

    subjetlabel = "ExKt%iSubJets" % (nsubjet)

    subjetrecorder.SubjetLabel = subjetlabel
    subjetrecorder.SubjetContainerName = SubjetContainerName

    from JetTagTools.JetTagToolsConf import Analysis__ExKtbbTagTool
    ExKtbbTagToolInstance = Analysis__ExKtbbTagTool(
      name = "ExKtbbTagTool%i_%s" % (nsubjet, JetCollectionName),
      JetAlgorithm = "Kt",
      JetRadius = 10.0,
      PtMin = 5000,
      ExclusiveNJets = 2,
      InputJetContainerName = JetCollectionName,   # touch
      SubjetContainerName = SubjetContainerName,
      # SubjetFinder = subjetfinder,
      SubjetRecorder = subjetrecorder,
      SubjetLabel = subjetlabel,
      SubjetAlgorithm_BTAG = "AntiKt",
      SubjetRadius_BTAG = 0.4
    )
    ToolSvc += ExKtbbTagToolInstance

    return (ExKtbbTagToolInstance, SubjetContainerName)

ExKtJetCollection__FatJet = ["AntiKt10LCTopoTrimmedPtFrac5SmallR20Jets", "AntiKt10TrackJets"]
# ExKtJetCollection__FatJet = ["AntiKt10LCTopoTrimmedPtFrac5SmallR20Jets"]
ExKtJetCollection__SubJet = []
for JetCollectionExKt in ExKtJetCollection__FatJet:
    # build ExKtbbTagTool instance
    (ExKtbbTagToolInstance, SubjetContainerName) = buildExclusiveSubjets(JetCollectionExKt, 2)
    ExKtJetCollection__SubJet += [SubjetContainerName]

    # build subjet collection through JetRecTool
    from JetRec.JetRecConf import JetRecTool
    jetrec = JetRecTool(
                         name = "JetRecTool_ExKtbb_%s" % (JetCollectionExKt),
                         OutputContainer = JetCollectionExKt,
                         InputContainer = JetCollectionExKt,
                         JetModifiers = [ExKtbbTagToolInstance],
                       )

    ToolSvc += jetrec
    FTAG5Seq += JetAlgorithm(
                             name = "JetAlgorithm_ExKtbb_%s" % (JetCollectionExKt),
                             Tools = [jetrec],
                            )

#===================================================================
# Run b-tagging
#===================================================================

defaultTaggers = ['IP2D', 'IP3D', 'SV0', 'MultiSVbb1', 'MultiSVbb2', 'SV1', 'BasicJetFitter', 'JetFitterTag', 'JetFitterNN', 'GbbNNTag', 'MV2c00', 'MV2c10', 'MV2c20', 'MV2c100', 'MV2m']
specialTaggers = ['ExKtbb_Hbb_MV2Only', 'ExKtbb_Hbb_MV2andJFDRSig', 'ExKtbb_Hbb_MV2andTopos']

from BTagging.BTaggingFlags import BTaggingFlags

# alias for VR
BTaggingFlags.CalibrationChannelAliases += ["AntiKtVR50Rmax4Rmin0Track->AntiKt4EMTopo"]

# alias for ExKt
BTaggingFlags.CalibrationChannelAliases += [
                                             "AntiKt10LCTopoTrimmedPtFrac5SmallR20->AntiKt10LCTopo,AntiKt6LCTopo,AntiKt6TopoEM,AntiKt4LCTopo,AntiKt4TopoEM,AntiKt4EMTopo",
                                             "AntiKt10LCTopoTrimmedPtFrac5SmallR20ExKt2Sub->AntiKt4LCTopo,AntiKt4TopoEM,AntiKt4EMTopo",

                                             "AntiKt10Track->AntiKt4EMTopo",
                                             "AntiKt10TrackExKt2Sub->AntiKt4EMTopo",
                                           ]

from DerivationFrameworkFlavourTag.FlavourTagCommon import FlavorTagInit
# must re-tag AntiKt4LCTopoJets and AntiKt4PV0TrackJets to make JetFitterNN work with corresponding VR jets (nikola: why?)
# also, re-tag R=0.2 track jets
FlavorTagInit( myTaggers = defaultTaggers, JetCollections = ["AntiKt4PV0TrackJets", "AntiKtVR50Rmax4Rmin0TrackJets", "AntiKt2PV0TrackJets"] + ExKtJetCollection__SubJet, Sequencer = FTAG5Seq )
FlavorTagInit( myTaggers = defaultTaggers + specialTaggers, JetCollections = ExKtJetCollection__FatJet, Sequencer = FTAG5Seq )

#====================================================================
# Schedule k-means clustering
#====================================================================

def buildKmeansSubjets(JetCollectionName, nsubjet, recordNameAppendix = "", ToolSvc = ToolSvc):
    RecordSubjetContainerName = "%sKmeans%i%sSubJets" % (JetCollectionName.replace("Jets", ""), nsubjet, recordNameAppendix)
    RecordSubjetLabel = "Kmeans%i%sSubJets" % (nsubjet, recordNameAppendix)

    from JetTagTools.JetTagToolsConf import Analysis__KmeansbbTagTool
    KmeansbbTagToolInstance = Analysis__KmeansbbTagTool(
      name = "KmeansbbTagTool%i%s_%s" % (nsubjet, recordNameAppendix, JetCollectionName),
      Debug = False, # True,
      PrimaryVertexContainerName = "PrimaryVertices",
      nAxis = nsubjet,
      SubjetLabel = "ExKt2SubJets", # -- hard-coded, since we want to make sure inputs are the same for all variaitons #  ( "ExKt%iSubJets" % (nsubjet) if nsubjet != 1 else "" ),
      UseJetSeed = False,
      JetRadiusCut = -1,
      PVDistance = 1.,
      InheritSeed = False,
      MaxChi2 = 25.,
      InputJetContainerName = JetCollectionName,
      RecordSubjetLabel = RecordSubjetLabel,
      RecordSubjetContainerName = RecordSubjetContainerName,
      RecordSubjetAlgorithm_BTAG = "AntiKt",
      RecordSubjetRadius_BTAG = 0.4
    )
    ToolSvc += KmeansbbTagToolInstance

    return (KmeansbbTagToolInstance, RecordSubjetContainerName)

KmeansJetCollection__FatJet = ["AntiKt10LCTopoTrimmedPtFrac5SmallR20Jets", "AntiKt10TrackJets"]
KmeansJetCollection__SubJet = []

def addKmeansCollection(nsubjet, recordNameAppendix = "", ToolSvc = ToolSvc, Sequencer = FTAG5Seq):
  global KmeansJetCollection__SubJet

  output = []
  for JetCollectionKmeans in KmeansJetCollection__FatJet:
      # build KmeansbbTagTool instance
      (KmeansbbTagToolInstance, SubjetContainerName) = buildKmeansSubjets(JetCollectionKmeans, nsubjet, recordNameAppendix, ToolSvc)
      KmeansJetCollection__SubJet += [SubjetContainerName]
      output += [SubjetContainerName]

      # build subjet collection through JetRecTool
      from JetRec.JetRecConf import JetRecTool
      jetrec = JetRecTool(
                           name = "JetRecTool_Kmeansbb%i%s_%s" % (nsubjet, recordNameAppendix, JetCollectionKmeans),
                           OutputContainer = JetCollectionKmeans,
                           InputContainer = JetCollectionKmeans,
                           JetModifiers = [KmeansbbTagToolInstance],
                         )

      ToolSvc += jetrec
      Sequencer += JetAlgorithm(
                                name = "JetAlgorithm_Kmeansbb%i%s_%s" % (nsubjet, recordNameAppendix, JetCollectionKmeans),
                                Tools = [jetrec],
                               )

  return output

KmeansSubJetCollection__StandardAssocTracks = []
KmeansSubJetCollection__StandardAssocTracks += addKmeansCollection(1, "StandardAssocTracks", ToolSvc, FTAG5Seq)
KmeansSubJetCollection__StandardAssocTracks += addKmeansCollection(2, "StandardAssocTracks", ToolSvc, FTAG5Seq)

KmeansSubJetCollection__MSVAssocTracks = []
KmeansSubJetCollection__MSVAssocTracks += addKmeansCollection(1, "MSVAssocTracks", ToolSvc, FTAG5Seq)
KmeansSubJetCollection__MSVAssocTracks += addKmeansCollection(2, "MSVAssocTracks", ToolSvc, FTAG5Seq)

#====================================================================
# b-Tagging on k-means subjets
#====================================================================

BTaggingFlags.CalibrationChannelAliases += [item[:-4]+"->AntiKt4LCTopo,AntiKt4TopoEM,AntiKt4EMTopo" for item in KmeansJetCollection__SubJet]

FlavorTagInit( myTaggers = defaultTaggers, JetCollections = KmeansSubJetCollection__StandardAssocTracks, Sequencer = FTAG5Seq )
FlavorTagInit( myTaggers = defaultTaggers, JetCollections = KmeansSubJetCollection__MSVAssocTracks, DoFullRetag="Kmeans", Sequencer = FTAG5Seq )

#====================================================================
# ATTEMPT TO FIX TRUTH LABELLING BUG
#====================================================================
#from JetRec.JetRecConf import JetRecTool
#JetRecTool_AntiKt10LCTopo = JetRecTool("AODFix_AntiKt10LCTopoJets", InputContainer = "AntiKt10LCTopoJets", OutputContainer = "AntiKt10LCTopoJets", JetModifiers = [ToolSvc.AODFix_jetdrlabeler]) # does not use the right DRMax for large-R jets!!
#ToolSvc += JetRecTool_AntiKt10LCTopo  

#====================================================================
# SKIMMING TOOLS
#====================================================================
# this is a basic cut for Higgs tagging studies
# with the mass cut
#offlineExpression = '((count (AntiKt10LCTopoTrimmedPtFrac5SmallR20Jets.m > 50 * GeV && AntiKt10LCTopoTrimmedPtFrac5SmallR20Jets.pt > 250 * GeV && (abs(AntiKt10LCTopoTrimmedPtFrac5SmallR20Jets.eta) < 2.0)) > 0 ))'
# without the mass cut
offlineExpression = '((count (AntiKt10LCTopoTrimmedPtFrac5SmallR20Jets.pt > 250 * GeV && (abs(AntiKt10LCTopoTrimmedPtFrac5SmallR20Jets.eta) < 2.0)) > 0 ))'

if globalflags.DataSource()=='data':
    triggers=[
        # primay triggers suggested for periodC
        "HLT_e24_lhmedium_L1EM18VH",
        "HLT_e60_lhmedium",
        "HLT_e120_lhloose",
        "HLT_e24_lhmedium_iloose_L1EM18VH",
        "HLT_e60_lhmedium",
        "HLT_e120_lhloose",
        "HLT_mu20_iloose_L1MU15",
        "HLT_mu40",
        # other single lepton triggers
        "HLT_e24_lhtight_iloose_L1EM20VH",
        "HLT_e24_tight_iloose_L1EM20VH",
        "HLT_mu14_iloose",
        "HLT_e17_loose"
        ]

    ORStr=" || "
    triggerStr=ORStr.join(triggers)
    triggerExpression = "((EventInfo.eventTypeBitmask==1) || (" + triggerStr +" ))"
    expression = offlineExpression+' && '+triggerExpression

else:
    triggers = []
    expression = offlineExpression

#====================================================================
# CREATE THE DERIVATION KERNEL ALGORITHM AND PASS THE ABOVE TOOLS
#====================================================================

# The name of the kernel (LooseSkimKernel in this case) must be unique to this derivation
from DerivationFrameworkCore.DerivationFrameworkCoreConf import DerivationFramework__DerivationKernel
DerivationFrameworkJob += FTAG5Seq
FTAG5Seq += CfgMgr.DerivationFramework__DerivationKernel("FTAG5Kernel")

#====================================================================
# SET UP STREAM
#====================================================================

# The base name (DAOD_FTAG5 here) must match the string in
streamName = derivationFlags.WriteDAOD_FTAG5Stream.StreamName
fileName   = buildFileName( derivationFlags.WriteDAOD_FTAG5Stream )
FTAG5Stream = MSMgr.NewPoolRootStream( streamName, fileName )
# Only events that pass the filters listed below are written out.
# Name must match that of the kernel above
# AcceptAlgs  = logical OR of filters
# RequireAlgs = logical AND of filters
FTAG5Stream.AcceptAlgs(["FTAG5Kernel"])

from DerivationFrameworkCore.SlimmingHelper import SlimmingHelper
FTAG5SlimmingHelper = SlimmingHelper("FTAG5SlimmingHelper")

# NB: the BTagging_AntiKt4EMTopo smart collection includes both AntiKt4EMTopoJets and BTagging_AntiKt4EMTopo
# container variables. Thus BTagging_AntiKt4EMTopo is needed in SmartCollections as well as AllVariables
FTAG5SlimmingHelper.AppendToDictionary = {
"AntiKtVR50Rmax4Rmin0TrackJets"               :   "xAOD::JetContainer"        ,
"AntiKtVR50Rmax4Rmin0TrackJetsAux"            :   "xAOD::JetAuxContainer"     ,
"BTagging_AntiKtVR50Rmax4Rmin0Track"          :   "xAOD::BTaggingContainer"   ,
"BTagging_AntiKtVR50Rmax4Rmin0TrackAux"       :   "xAOD::BTaggingAuxContainer",
"AntiKt10LCTopoTrimmedPtFrac5SmallR20Jets"    :   "xAOD::JetContainer"        ,
"AntiKt10LCTopoTrimmedPtFrac5SmallR20JetsAux" :   "xAOD::JetAuxContainer"     ,
"AntiKt10TrackJets"                           :   "xAOD::JetContainer"        ,
"AntiKt10TrackJetsAux"                        :   "xAOD::JetAuxContainer"     ,
}

# Write ExKt contents to output stream
# Unfortunately it seems it has to be done in this way ... 
for JetCollectionExKt in ExKtJetCollection__SubJet + ExKtJetCollection__FatJet + KmeansJetCollection__SubJet:
    JetName = JetCollectionExKt[:-4]

    # Jets (sub-jets only) #

    if JetCollectionExKt in ExKtJetCollection__SubJet + KmeansJetCollection__SubJet:
        FTAG5SlimmingHelper.StaticContent.append("xAOD::JetContainer#"+JetCollectionExKt)
        FTAG5SlimmingHelper.StaticContent.append("xAOD::JetAuxContainer#"+JetCollectionExKt+"Aux.-Parent")    # "Parent" link is broken after deep copy of parent jet in b-tagging module

    # b-tagging #

    FTAG5SlimmingHelper.StaticContent.append("xAOD::BTaggingContainer#BTagging_"+JetName)
    FTAG5SlimmingHelper.StaticContent.append("xAOD::BTaggingAuxContainer#BTagging_" + JetName + "Aux.")

    FTAG5SlimmingHelper.StaticContent.append("xAOD::VertexContainer#BTagging_" + JetName + "SecVtx")
    FTAG5SlimmingHelper.StaticContent.append("xAOD::VertexAuxContainer#BTagging_" + JetName + "SecVtx" + "Aux.")

    FTAG5SlimmingHelper.StaticContent.append("xAOD::BTagVertexContainer#BTagging_" + JetName + "JFVtx")
    FTAG5SlimmingHelper.StaticContent.append("xAOD::BTagVertexAuxContainer#BTagging_" + JetName + "JFVtx" + "Aux.-vxTrackAtVertex")


FTAG5SlimmingHelper.SmartCollections = ["Electrons","Muons",
                                        "MET_Reference_AntiKt4EMTopo",
                                        "AntiKt4EMTopoJets",
                                        "BTagging_AntiKt4EMTopo"]

FTAG5SlimmingHelper.AllVariables = ["AntiKt3PV0TrackJets",
                                    "AntiKt2PV0TrackJets",
                                    "AntiKt4PV0TrackJets",
                                    "AntiKt4TruthJets",
                                    "BTagging_AntiKt4EMTopo",
                                    "BTagging_AntiKt2Track",
                                    "BTagging_AntiKt3Track",
                                    "BTagging_AntiKt4EMTopoJFVtx",
                                    "BTagging_AntiKt2TrackJFVtx",
                                    "BTagging_AntiKt3TrackJFVtx",
                                    "BTagging_AntiKt4EMTopoSecVtx",
                                    "BTagging_AntiKt2TrackSecVtx",
                                    "BTagging_AntiKt3TrackSecVtx",
                                    "TruthVertices",
                                    "TruthParticles",
                                    "TruthEvents",
                                    "MET_Truth",
                                    "MET_TruthRegions",
                                    "InDetTrackParticles",
                                    "PrimaryVertices",
                                    "AntiKtVR50Rmax4Rmin0TrackJets", 
                                    "BTagging_AntiKtVR50Rmax4Rmin0Track",
                                    "AntiKt10LCTopoJets",
                                    "AntiKt10LCTopoTrimmedPtFrac5SmallR20Jets",
                                    "AntiKt10TrackJets",
                                    "CaloCalTopoClusters",
                                    ]

from DerivationFrameworkJetEtMiss.AntiKt4EMTopoJetsCPContent import AntiKt4EMTopoJetsCPContent

FTAG5SlimmingHelper.ExtraVariables.append(AntiKt4EMTopoJetsCPContent[1].replace("AntiKt4EMTopoJetsAux","AntiKt10LCTopoJets"))
FTAG5SlimmingHelper.ExtraVariables.append("AntiKt10LCTopoJets.GhostAntiKt2TrackJet")
FTAG5SlimmingHelper.ExtraVariables.append("AntiKt10LCTopoJets.GhostAntiKt2TrackJetPt")
FTAG5SlimmingHelper.ExtraVariables.append("AntiKt10LCTopoJets.GhostAntiKt2TrackJetCount")
FTAG5SlimmingHelper.ExtraVariables.append("AntiKt10LCTopoJets.ConeExclBHadronsFinal")
FTAG5SlimmingHelper.ExtraVariables.append("AntiKt10LCTopoJets.ConeExclCHadronsFinal")

FTAG5SlimmingHelper.IncludeMuonTriggerContent = True
FTAG5SlimmingHelper.IncludeEGammaTriggerContent = True
FTAG5SlimmingHelper.IncludeJetTriggerContent = True
FTAG5SlimmingHelper.IncludeEtMissTriggerContent = True
FTAG5SlimmingHelper.IncludeBJetTriggerContent = True

FTAG5SlimmingHelper.AppendContentToStream(FTAG5Stream)
