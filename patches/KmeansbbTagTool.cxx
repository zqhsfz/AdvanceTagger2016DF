#include "JetTagTools/KmeansbbTagTool.h"

#include "GaudiKernel/IToolSvc.h"
#include "JetTagTools/JetTagUtils.h"     // get jet author
#include "xAODBTagging/SecVtxHelper.h"   // extract SV information

#include <algorithm>
#include "TList.h"

#include "TMVA/Reader.h"


namespace Analysis {

/** 
    @class KmeansbbTagTool 
    Identification of double b-hadron jets using k-means clusterization of seed vertices
    This class only serves to do clusterization, and build "subjet" based on it.
    @author Qi. Zeng
*/

KmeansbbTagTool::KmeansbbTagTool(std::string myname):
  JetModifierBase(myname),
  m_debug(false),
  m_primaryVertxContainerName("PrimaryVertices"),
  m_nAxis(2),
  m_subjetLabel(""),
  m_useJetSeed(false),
  m_JetRadiusCut(-1.),
  m_PVDistance(1.),
  m_InheritSeed(false),
  m_MaxChi2(9.),
  m_inputJetContainerName(""),
  m_recordSubjetLabel(""),
  m_recordSubjetContainerName(""),
  m_recordSubjetAlgorithm_BTAG(""),
  m_recordSubjetRadius_BTAG(0.)
{
  // options
  declareProperty("Debug",                      m_debug);

  declareProperty("PrimaryVertexContainerName", m_primaryVertxContainerName);
  declareProperty("nAxis",                      m_nAxis);
  declareProperty("SubjetLabel",                m_subjetLabel);
  declareProperty("UseJetSeed",                 m_useJetSeed);
  declareProperty("JetRadiusCut",               m_JetRadiusCut);
  declareProperty("PVDistance",                 m_PVDistance);
  declareProperty("InheritSeed",                m_InheritSeed);
  declareProperty("MaxChi2",                    m_MaxChi2);

  declareProperty("InputJetContainerName",      m_inputJetContainerName);
  declareProperty("RecordSubjetLabel",          m_recordSubjetLabel);
  declareProperty("RecordSubjetContainerName",  m_recordSubjetContainerName);
  declareProperty("RecordSubjetAlgorithm_BTAG", m_recordSubjetAlgorithm_BTAG);
  declareProperty("RecordSubjetRadius_BTAG",    m_recordSubjetRadius_BTAG);

  // other initialization
}


StatusCode KmeansbbTagTool::initialize() {

  return StatusCode::SUCCESS;
}


int KmeansbbTagTool::modifyJet(xAOD::Jet& jetToTag) const{

  /////////////////////////////////////
  // Initialize the subjet container //
  /////////////////////////////////////

  if( (m_recordSubjetLabel == "") || (m_recordSubjetContainerName == "") ){
    ATH_MSG_ERROR("Invalid subjet label or container name");
    return 0;
  }

  auto subjet_container = evtStore()->tryRetrieve<xAOD::JetContainer>(m_recordSubjetContainerName);

  if(subjet_container == 0){
    // first time building the subjet
    subjet_container = new xAOD::JetContainer;
    subjet_container->setStore(new xAOD::JetAuxContainer);

    // record xAOD
    StatusCode sc = evtStore()->record(subjet_container, m_recordSubjetContainerName);
    if(sc.isFailure()){
      ATH_MSG_ERROR("Error recording subjet container " << m_recordSubjetContainerName);
      return 0;
    }

    // record xAODAux
    sc = evtStore()->record(dynamic_cast<xAOD::JetAuxContainer*>(subjet_container->getStore()), m_recordSubjetContainerName+"Aux.");
    if(sc.isFailure()){
      ATH_MSG_ERROR("Error recording subjet aux container " << m_recordSubjetContainerName << "Aux.");
      return 0;
    }
  }

  // set dummy association to subjets
  std::vector<const xAOD::Jet*> v_subjet_const__dummy;
  std::vector<int> v_subjet_index__dummy;

  jetToTag.setAssociatedObjects(m_recordSubjetLabel, v_subjet_const__dummy);
  jetToTag.auxdata<std::string>(m_recordSubjetLabel + "_ContainerName") = m_recordSubjetContainerName;
  jetToTag.auxdata<std::vector<int> >(m_recordSubjetLabel + "_IndexList") = v_subjet_index__dummy;

  // set dummy association to input vertices
  std::vector<const xAOD::Vertex*> v_inputVtxList__dummy;
  jetToTag.setAssociatedObjects("Kmeans_inputVtxList", v_inputVtxList__dummy);

  ////////////////////////
  // Get Primary Vertex //
  ////////////////////////

  const xAOD::Vertex* m_priVtx = RetrievePrimaryVertex();
  if(m_priVtx == 0){
    ATH_MSG_ERROR("Failure in retrieving primary vertex!");
    return 0;
  }

  ///////////////////////////////////
  // Prepare raw input vertex list //
  ///////////////////////////////////

  std::vector<const xAOD::Vertex*> rawInputVtxList = RetrieveMSVInJet(jetToTag);

  std::vector<const xAOD::Jet*> subjetList;
  if(m_subjetLabel != ""){
    // check existence
    std::string subjetContainerName = jetToTag.auxdata<std::string>(m_subjetLabel+"_ContainerName");
    auto subjetContainer = evtStore()->tryRetrieve<xAOD::JetContainer>(subjetContainerName);
    if(subjetContainer == 0){
      ATH_MSG_ERROR("Unable to get subjet container " << subjetContainerName << "!");
      return 0;
    }

    std::vector<int> subjetIndexList = jetToTag.auxdata<std::vector<int> >(m_subjetLabel+"_IndexList");
    for(auto index : subjetIndexList){
      auto subjet = subjetContainer->at(index);
      subjetList.push_back(subjet);

      std::vector<const xAOD::Vertex*> inputVtxListToBeInsert = RetrieveMSVInJet(*subjet);
      rawInputVtxList.insert(rawInputVtxList.end(), inputVtxListToBeInsert.begin(), inputVtxListToBeInsert.end());
    }

    // In case number of subjets less than number of requested clusters, we still use vertices found by any existing subjet.
    // However, subjet information will no longer be used any further (seed direction etc.)
    if(int(subjetList.size()) < m_nAxis){
      ATH_MSG_WARNING("Number of subjets less than number of requested clusters. No subjet information will be further used.");
      subjetList.clear();
    }
  }

  ///////////////////////////////////
  // Clean up the input vertx list //
  ///////////////////////////////////

  std::vector<const xAOD::Vertex*> selectInputVtxList = InputVtxCleaning(m_priVtx, jetToTag, rawInputVtxList);
  jetToTag.setAssociatedObjects("Kmeans_inputVtxList", selectInputVtxList);

  // sanity check
  if(int(selectInputVtxList.size()) < m_nAxis){
    if(m_debug){
      ATH_MSG_WARNING("Impossible to run k-means: Number of input vertices less than number of requested clusters. Before cleaning: " << rawInputVtxList.size() << " ; After cleaning: " << selectInputVtxList.size());
    }

    return 0;
  }

  ////////////////////////////////
  // Prepare the seed direction //
  ////////////////////////////////

  std::vector<Amg::Vector3D> m_JetSeeds;

  if(m_useJetSeed){
    if(m_nAxis == 1){
      m_JetSeeds.push_back( FromTVector3ToVector3D(jetToTag.p4().Vect()) );
    }
    else if(int(subjetList.size()) >= m_nAxis){
      std::vector<const xAOD::Jet*> subjetList_sorted = subjetList;
      std::sort(subjetList_sorted.begin(), subjetList_sorted.end(), KmeansbbTagTool::SortPt);

      for(int iclass = 0; iclass < m_nAxis; iclass++) m_JetSeeds.push_back( FromTVector3ToVector3D(subjetList_sorted[iclass]->p4().Vect()) );
    }
    else{
      ATH_MSG_WARNING("Insufficient subjets for seeding. Random initialization will be used.");
      m_JetSeeds.clear();
    }
  }

  ////////////////////////////////
  // Now the k-means clustering //
  ////////////////////////////////

  std::vector<KmeansElement>                history;
  std::vector<std::vector<KmeansElement> >  detailHistory;

  KmeansElement result_kmeans = AdaptiveKmeans(m_priVtx, selectInputVtxList, history, detailHistory, m_JetSeeds);

  /////////////////////
  // Build "subjets" //
  /////////////////////

  // sanity check
  if(int(result_kmeans.class_axis.size()) != m_nAxis){
    ATH_MSG_WARNING("Failure in finding " << m_nAxis << " clusters through k-means. No further subjet will be built.");
    return 0;
  }

  std::vector<const xAOD::Jet*> v_subjet_const;
  std::vector<int> v_subjet_index;
  for(int iclass = 0; iclass < m_nAxis; iclass++){
    Amg::Vector3D axis = result_kmeans.class_axis[iclass];
    TVector3 axis_HEP = FromVector3DToTVector3(axis);
    std::vector<const xAOD::Vertex*> vtxList = result_kmeans.class_vtxList[iclass];

    xAOD::Jet* subjet = new xAOD::Jet();
    subjet_container->push_back(subjet);
    v_subjet_const.push_back(subjet);
    v_subjet_index.push_back(subjet->index());

    // set 4-mom
    subjet->setJetP4(xAOD::JetFourMom_t(jetToTag.pt()/m_nAxis, axis_HEP.Eta(), axis_HEP.Phi(), 0.));      // TODO: pT and mass are dummy for now. While it is OK to leave mass to be dummy, we can assign pT based on some heuristic rule
    subjet->setJetP4(xAOD::JetConstitScaleMomentum, subjet->jetP4());

    // set author/radius/inputType etc. (all alias)
    // Just b-tagging alias. User should make sure this is consistent with the b-tagging calibration alias.
    subjet->setAlgorithmType(xAOD::JetAlgorithmType::algId(m_recordSubjetAlgorithm_BTAG));
    subjet->setSizeParameter(m_recordSubjetRadius_BTAG);
    subjet->setInputType(jetToTag.getInputType());                                         // TODO: should we set this to tracks? Or, does it matter?
    subjet->setConstituentsSignalState(jetToTag.getConstituentsSignalState());             // TODO: same as above

    // set association to parent jet
    const xAOD::JetContainer* parent_container = dynamic_cast<const xAOD::JetContainer*>(jetToTag.container());
    ElementLink<xAOD::JetContainer> el_parent(*parent_container, jetToTag.index());
    subjet->setAttribute("Parent", el_parent);
    subjet->auxdata<std::string>("Parent_ContainerName") = m_inputJetContainerName;
    subjet->auxdata<int>("Parent_Index") = jetToTag.index();

    // set association to vertices / tracks
    std::vector<ElementLink<xAOD::VertexContainer> > v_el_vtx;
    std::vector<ElementLink<xAOD::TrackParticleContainer> > v_el_tracks;
    for(auto vtx : vtxList){
      const xAOD::VertexContainer* vtx_container = dynamic_cast<const xAOD::VertexContainer*>(vtx->container());
      ElementLink<xAOD::VertexContainer> el_vtx(*vtx_container, vtx->index());
      v_el_vtx.push_back(el_vtx);

      std::vector<const xAOD::TrackParticle*> v_track_checkduplicate;   // only for checking dupication purpose
      for(auto track : VertexTrack(vtx)){
        const xAOD::TrackParticleContainer* track_container = dynamic_cast<const xAOD::TrackParticleContainer*>(track->container());
        ElementLink<xAOD::TrackParticleContainer> el_track(*track_container, track->index());

        if(CheckDuplicate<const xAOD::TrackParticle*>(v_track_checkduplicate, track)) continue;

        v_el_tracks.push_back(el_track);
        v_track_checkduplicate.push_back(track);
      }
    }
    subjet->setAttribute("Kmeans_AssocVertices", v_el_vtx);
    subjet->setAttribute("Kmeans_AssocTracks", v_el_tracks);
  }

  // set association to subjets
  jetToTag.setAssociatedObjects(m_recordSubjetLabel, v_subjet_const);
  jetToTag.auxdata<std::string>(m_recordSubjetLabel + "_ContainerName") = m_recordSubjetContainerName;
  jetToTag.auxdata<std::vector<int> >(m_recordSubjetLabel + "_IndexList") = v_subjet_index;

  if(m_debug){
   DumpKmeansElement(m_priVtx, jetToTag, history[0], "Initial");
   DumpKmeansElement(m_priVtx, jetToTag, result_kmeans, "Final");
  }

  return StatusCode::SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////////////////////

std::vector<const xAOD::Vertex*> KmeansbbTagTool::InputVtxCleaning(const xAOD::Vertex* primaryVertex, const xAOD::Jet & jet, const std::vector<const xAOD::Vertex*> & inputVtxList) const{
  std::vector<const xAOD::Vertex*> output;

  for(auto vtx : inputVtxList){
    // PV-SV should point to roughly same direction as jet
    Amg::Vector3D PV_SV = vtx->position() - primaryVertex->position();
    TVector3      PV_SV_HEP = FromVector3DToTVector3(PV_SV);
    Amg::Vector3D jetDirection = FromTVector3ToVector3D(jet.p4().Vect());

    if(m_JetRadiusCut >= 0.){
      if(PV_SV_HEP.DeltaR(jet.p4().Vect()) > m_JetRadiusCut) continue;
    }
    else if(m_JetRadiusCut == -1){
      // semi-sphere defined by jet direction
      if(PV_SV.dot(jetDirection) <= 0.) continue;
    }
    else{
      ATH_MSG_WARNING("Undefined JetRadiusCut: " << m_JetRadiusCut << ". -1 scheme will be used.");
      if(PV_SV.dot(jetDirection) <= 0.) continue;
    }
    
    // SV cannot be too close to PV
    double dR3D = 0;
    double dR3D_error = 0;
    double dR3D_sig = GeneralDRSignificance3D(primaryVertex->position(), primaryVertex->covariancePosition(), vtx->position(), vtx->covariancePosition(), dR3D, dR3D_error);

    // debug -- calcualte this in another way
    // --> It has been shown that two method gives almost the same answer.
    // Don't remove this for future reference.
    // std::cout << "======== DEBUG ON DR SIGNIFICANCE ========" << std::endl;
    // std::cout << "From Method 1: " << dR3D << " " << dR3D_error << std::endl;
    // double dR3D_new = PV_SV_HEP.Mag();
    // auto PV_SV_CovMatrix = primaryVertex->covariancePosition() + vtx->covariancePosition();
    // Amg::Vector3D J = PV_SV/dR3D_new;
    // double dR3D_error_new = TMath::Sqrt( J.transpose() * PV_SV_CovMatrix * J );
    // std::cout << "From Method 2: " << dR3D_new << " " << dR3D_error_new << std::endl;
    // std::cout << "==========================================" << std::endl;


    if(dR3D_sig < m_PVDistance) continue;

    // remove redundant vertex
    bool f_redundant = false;

    for(auto vtx_select : output){
      auto tracks_current = VertexTrack(vtx);
      auto tracks_select  = VertexTrack(vtx_select);
      auto tracks_common  = FindCommon(tracks_current, tracks_select);

      if( (tracks_common.size() == tracks_current.size()) && (tracks_common.size() == tracks_select.size()) ){
        f_redundant = true;
        break;
      }
    }

    if(f_redundant) continue;

    // store
    output.push_back(vtx);
  }

  return output;
}

KmeansElement KmeansbbTagTool::AdaptiveKmeans(const xAOD::Vertex* primaryVertex, std::vector<const xAOD::Vertex*> & inputVtxList, std::vector<KmeansElement> & history, std::vector<std::vector<KmeansElement> > & detailHistory, std::vector<Amg::Vector3D> ClusterSeeds) const{

  KmeansElement empty;
  int count = 0;

  while(true){
    // InheritSeed
    if( (history.size() > 0) && (m_InheritSeed) ){
      ClusterSeeds.clear();

      KmeansElement LastResult = history.back();
      for(auto axis : LastResult.class_axis){
        ClusterSeeds.push_back(axis);
      }
    }

    // kernel k-means
    std::vector<KmeansElement> history_onekmeans;
    KmeansElement result_onekmeans = SimpleKmeans(primaryVertex, inputVtxList, history_onekmeans, ClusterSeeds);

    // store history
    history.push_back(result_onekmeans);
    detailHistory.push_back(history_onekmeans);
    count++;

    // check if the k-means is successful
    if((int)(result_onekmeans.class_axis.size()) != m_nAxis){
      ATH_MSG_WARNING("K-means failure in SimpleKmeans module.");
      return empty;
    }

    if(m_debug){
      std::cout << "Checking the chi-square of cluster results ..." << std::endl;
    }

    // remove the vtx with largest chi2 constribution that is above threshold
    const xAOD::Vertex* vtx_outlier = 0;
    double largest_chisquare = -1.;

    for(unsigned int iclass = 0; iclass < result_onekmeans.class_vtxList.size(); iclass++){
      Amg::Vector3D                     axis    = result_onekmeans.class_axis[iclass];
      std::vector<const xAOD::Vertex*>  vtxList = result_onekmeans.class_vtxList[iclass];

      for(auto vtx : vtxList){
        std::vector<const xAOD::Vertex*> temp; temp.push_back(vtx);

        double chisquare = GetChisquare(primaryVertex, axis, temp);

        if(m_debug){
          std::cout << "vtx " << "(" << vtx->position()(0) << "," << vtx->position()(1) << "," << vtx->position()(2) << ")" << " has chi-square " << chisquare << " w.r.t axis " << "(" << axis(0) << "," << axis(1) << "," << axis(2) << ")" << std::endl;
        }

        if( (chisquare > m_MaxChi2) && (chisquare > largest_chisquare) ){
          vtx_outlier = vtx;
          largest_chisquare = chisquare;
        }
      }
    }

    if(m_debug){
      std::cout << "Checking DONE" << std::endl;
    }

    if(largest_chisquare == -1){
      // converge
      return result_onekmeans;
    }
    else{
      if(m_debug){
        ATH_MSG_WARNING("Vertex with chi-square: " << largest_chisquare << " is removed.");
      }

      // remove it
      std::vector<const xAOD::Vertex*>::iterator FindIter = std::find(inputVtxList.begin(), inputVtxList.end(), vtx_outlier);
      if(FindIter == inputVtxList.end()){
        ATH_MSG_ERROR("Cannot find the outlier vertex in original input list. There must be something crazy going on ...");
        return empty;
      }
      else{
        inputVtxList.erase(FindIter);
      }
    }
  }

  ATH_MSG_ERROR("Wait! You are not supposed to see this!");
  return empty;
}

// assuming input size is not smaller than k
// The protection should be done in wherever that calls this function
// We also do not check the sanity of output here (if it is same as k). This should be done at where calls this function
// inputVtxList must NOT pass in by reference, since there is erasing operation on it.
std::vector<Amg::Vector3D> KmeansbbTagTool::InitializeKmeans(const xAOD::Vertex* primaryVertex, std::vector<const xAOD::Vertex*> inputVtxList, std::vector<Amg::Vector3D> ClusterSeeds, bool doKmeansPP) const{
  if(ClusterSeeds.size() != 0){
    return ClusterSeeds;
  }

  std::vector<Amg::Vector3D> output;

  if(doKmeansPP){
    // Kmeans++ initialization
    std::default_random_engine generator;

    std::vector<const xAOD::Vertex*> chosenVertices;

    std::vector<double> weights(inputVtxList.size());
    for(unsigned int index = 0; index < weights.size(); index++) weights[index] = 1.;

    int counter = 0;
    while(int(output.size()) < m_nAxis){
      // protection against infinite loop
      counter++;
      if(counter > (m_nAxis * 2)){
        ATH_MSG_ERROR("Potential infinite loop during k-means++ initialization! This should not happen. Program will still proceed, but you should check what is going on!");
        break;
      }

      // randomly choose a vtx according to weight
      std::discrete_distribution<int> distribution(weights.begin(), weights.end());
      int chosenIndex = distribution(generator);
      const xAOD::Vertex* chosenVtx = inputVtxList[chosenIndex];

      // push back
      if(std::find(chosenVertices.begin(), chosenVertices.end(), chosenVtx) == chosenVertices.end()){
        chosenVertices.push_back(chosenVtx);
        output.push_back(chosenVtx->position() - primaryVertex->position());
      }
      else{
        ATH_MSG_ERROR("Duplicate vertex chosen during k-means++. This should not happen. Program will still proceed, but you should check what is going on!");
        continue;
      }

      // erase chosen one
      // By definition of k-means++, this is not necessary, since the weight for chosen one will be 0.
      // However, to prevent numerical errors, we would simply remove chosen items from the input list
      inputVtxList.erase(inputVtxList.begin() + chosenIndex);
      weights.erase(weights.begin() + chosenIndex);

      // calculate new weights
      for(unsigned int index = 0; index < inputVtxList.size(); index++){
        auto vtx = inputVtxList[index];
        
        int assignIndex__temporary = HardAssignment(primaryVertex, vtx, output);
        std::vector<const xAOD::Vertex*> inputVtxList__temporary { vtx };

        double newWeight = GetChisquare(primaryVertex, output[assignIndex__temporary], inputVtxList__temporary);
        if(newWeight < 0.){
          ATH_MSG_ERROR("Negative chi-square obtained during k-means++. This should not happend. Program will still proceed by setting it to 0, but you should check what is going on!");
          newWeight = 0.;
        }
        weights[index] = newWeight;
      }
    }
  }
  else{
    // Random initialization
    for(int iclass = 0; iclass < m_nAxis; iclass++){
      int random_index = random() % (inputVtxList.size());

      const xAOD::Vertex* random_vtx = inputVtxList[random_index];
      Amg::Vector3D random_axis = random_vtx->position() - primaryVertex->position();
      output.push_back(random_axis);

      inputVtxList.erase(inputVtxList.begin() + random_index);
    }
  }

  return output;
}

KmeansElement KmeansbbTagTool::SimpleKmeans(const xAOD::Vertex* primaryVertex, const std::vector<const xAOD::Vertex*> & inputVtxList, std::vector<KmeansElement> & ClusHistory, std::vector<Amg::Vector3D> ClusterSeeds) const{
  KmeansElement empty;

  // should have enough input for clustering //
  if((int)(inputVtxList.size()) < m_nAxis){
    if(m_debug){
      ATH_MSG_WARNING("Number of input vertices less than number of requested clusters: " << inputVtxList.size() << " / " << m_nAxis);
    }
    return empty;
  }

  // initialization //
  KmeansElement output;
  output.class_vtxList.clear();
  output.class_axis.clear();
  output.class_chisquare.clear();

  // Initialization //
  // k-means++ is turned on automatically, since this is already a pretty standard thing nowadays
  output.class_axis = InitializeKmeans(primaryVertex, inputVtxList, ClusterSeeds, true);
  if(int(output.class_axis.size()) != m_nAxis){
    ATH_MSG_ERROR("Size of seeds not equal to number of requested clusters: " << output.class_axis.size() << " / " << m_nAxis);
    return empty;
  }
  for(int iclass = 0; iclass < m_nAxis; iclass++){
    output.class_chisquare.push_back(-1.);

    std::vector<const xAOD::Vertex*> dummyVtxList; 
    output.class_vtxList.push_back(dummyVtxList);
  }

  // start k-means clustering //
  double precision = 0.001;
  double last_sum_chisquare = 0.;
  int count = 0;

  while(true){
    count++;

    // store thie history (after axis update, before re-assignment)
    ClusHistory.push_back(output);

    // reset
    for(int iclass = 0; iclass < m_nAxis; iclass++){
      output.class_chisquare[iclass] = -1.;
      output.class_vtxList[iclass].clear();
    }

    // new (hard) assignment
    for(unsigned int ivtx = 0; ivtx < inputVtxList.size(); ivtx++){
      auto vtx_current = inputVtxList[ivtx];

      int index_assignment = HardAssignment(primaryVertex, vtx_current, output.class_axis);
      if(index_assignment == -1){
        ATH_MSG_ERROR("Failure in assignment step within SimpleKmeans!");
        return empty;
      }

      output.class_vtxList[index_assignment].push_back(vtx_current);
    }

    // check if any cluster is empty, which is disallowed (otherwise it would carsh)
    for(int iclass = 0; iclass < m_nAxis; iclass++){
      if(output.class_vtxList[iclass].size() == 0){
        ATH_MSG_ERROR("Empty cluster found! You will crash very soon ...");
      }
    }

    // store the history (after re-assignment, before axis update)
    ClusHistory.push_back(output);

    // update axis
    double sum_chisquare = 0.;
    for(unsigned int iclass = 0; iclass < output.class_vtxList.size(); iclass++){
      Amg::Vector3D new_axis = UpdateAxis(primaryVertex, output.class_vtxList[iclass]);
      output.class_axis[iclass] = new_axis;
      output.class_chisquare[iclass] = GetChisquare(primaryVertex, new_axis, output.class_vtxList[iclass]);
      sum_chisquare += output.class_chisquare[iclass];
    }

    // exit condition
    if(fabs(sum_chisquare - last_sum_chisquare) < precision){
      return output;
    }

    // emergency -- infinite loop
    if(count > 1000){
      ATH_MSG_WARNING("Too many iteration in SimpleKmeans. Whatever current result would be returned.");
      return output;
    }

    // update chisquare
    last_sum_chisquare = sum_chisquare;
  }

  ATH_MSG_ERROR("Wait! You are not supposed to see this!");
  return empty;
}

double KmeansbbTagTool::GetChisquare(const xAOD::Vertex* primaryVertex, const Amg::Vector3D & axis, const std::vector<const xAOD::Vertex*> & inputVtxList) const{
  double sum = 0.;

  for(auto vtx : inputVtxList){
    Amg::Vector3D  PV_SV = vtx->position() - primaryVertex->position();
    AmgSymMatrix3D PV_SV_CovMatrix = vtx->covariancePosition() + primaryVertex->covariancePosition();   // OK we are assuming correlation are null ... (which might not be exacty true)

    TVector3 PV_SV_HEP = FromVector3DToTVector3(PV_SV);
    double theta = PV_SV_HEP.Theta();
    double phi = PV_SV_HEP.Phi();

    auto directionCovMatrix = GetDirectionCovMatrix(PV_SV, PV_SV_CovMatrix);

    // take inverse
    bool invertible;
    AmgSymMatrix2D inverse = Inverse<AmgSymMatrix2D, AmgSymMatrix2D>(directionCovMatrix, invertible);

    if(!invertible){
      ATH_MSG_WARNING("Direction covariant matrix is unexpectedly non-invertible (in GetChisquare)! Identity matrix is used as placeholder.");

      if(m_debug){
        std::cout << "Here is the matrix to be inversed: " << std::endl << directionCovMatrix << std::endl;
        std::cout << "Here is the inversed matrix: " << std::endl << inverse << std::endl;
        std::cout << "Here is the PV_SV: " << std::endl << PV_SV << std::endl;
        std::cout << "Here is the covariant matrix in XYZ: " << std::endl << PV_SV_CovMatrix << std::endl;
        std::cout << "Here is the Jacobian: " << std::endl << GetJacobian_FromXYZToSphere(PV_SV) << std::endl;
      }

      inverse(0, 0) = 1; inverse(0, 1) = 0;
      inverse(1, 0) = 0; inverse(1, 1) = 1;
    }

    // // temporary debug
    // std::cout << "==============================" << std::endl;
    // std::cout << "Here is the matrix to be inversed: " << std::endl << directionCovMatrix << std::endl;
    // std::cout << "Here is the inversed matrix: " << std::endl << inverse << std::endl;
    // std::cout << "Here is the PV_SV: " << std::endl << PV_SV << std::endl;
    // std::cout << "Here is the covariant matrix in XYZ: " << std::endl << PV_SV_CovMatrix << std::endl;
    // std::cout << "Here is the Jacobian: " << std::endl << GetJacobian_FromXYZToSphere(PV_SV) << std::endl;
    // std::cout << "==============================" << std::endl;

    // Gaussian chi-square
    TVector3 axis_HEP = FromVector3DToTVector3(axis);
    double axis_theta = axis_HEP.Theta();
    double axis_phi = axis_HEP.Phi();

    Amg::Vector2D diff;
    diff(0) = axis_theta - theta;
    diff(1) = axis_phi - phi;

    sum += ( diff.transpose() * inverse * diff );
  }

  return sum;
}

Amg::Vector3D KmeansbbTagTool::UpdateAxis(const xAOD::Vertex* primaryVertex, const std::vector<const xAOD::Vertex*> & inputVtxList) const{

  AmgSymMatrix2D sum_InverseCov;
  sum_InverseCov(0, 0) = 0.; sum_InverseCov(0, 1) = 0.;
  sum_InverseCov(1, 0) = 0.; sum_InverseCov(1, 1) = 0.;

  Amg::Vector2D sum_InverseCovX;
  sum_InverseCovX(0) = 0.;
  sum_InverseCovX(1) = 0.;

  double sum_theta_backup = 0.;
  double sum_phi_backup = 0.;

  for(auto vtx : inputVtxList){
    Amg::Vector3D  PV_SV = vtx->position() - primaryVertex->position();
    AmgSymMatrix3D PV_SV_CovMatrix = vtx->covariancePosition() + primaryVertex->covariancePosition();   // OK we are assuming correlation are null ... (which might not be exacty true)

    TVector3 PV_SV_HEP = FromVector3DToTVector3(PV_SV);
    double theta = PV_SV_HEP.Theta();
    double phi = PV_SV_HEP.Phi();

    auto directionCovMatrix = GetDirectionCovMatrix(PV_SV, PV_SV_CovMatrix);

    // take inverse
    bool invertible;
    AmgSymMatrix2D inverse = Inverse<AmgSymMatrix2D, AmgSymMatrix2D>(directionCovMatrix, invertible);

    if(!invertible){
      ATH_MSG_WARNING("Direction covariant matrix is unexpectedly non-invertible (in UpdateAxis)! Identity matrix is used as placeholder.");

      if(m_debug){
        std::cout << "Here is the matrix to be inversed: " << std::endl << directionCovMatrix << std::endl;
        std::cout << "Here is the inversed matrix: " << std::endl << inverse << std::endl;
        std::cout << "Here is the PV_SV: " << std::endl << PV_SV << std::endl;
        std::cout << "Here is the covariant matrix in XYZ: " << std::endl << PV_SV_CovMatrix << std::endl;
        std::cout << "Here is the Jacobian: " << std::endl << GetJacobian_FromXYZToSphere(PV_SV) << std::endl;
      }

      inverse(0, 0) = 1; inverse(0, 1) = 0;
      inverse(1, 0) = 0; inverse(1, 1) = 1;
    }

    // gather necessary factors
    Amg::Vector2D x;
    x(0) = theta;
    x(1) = phi;

    sum_InverseCov += inverse;
    sum_InverseCovX += (inverse * x);

    sum_theta_backup += theta;
    sum_phi_backup += phi;
  }

  // take inverse
  bool invertible;
  AmgSymMatrix2D sum_InverseCov_Inverse = Inverse<AmgSymMatrix2D, AmgSymMatrix2D>(sum_InverseCov, invertible);

  // prepare output
  Amg::Vector2D output_2D;

  if(!invertible){
    ATH_MSG_WARNING("Sum-of-Inverse-CovMatrix is unexpectedly non-invertible! Simply average is used as placeholder.");

    if(m_debug){
      std::cout << "Here is the matrix to be inversed: " << std::endl << sum_InverseCov << std::endl;
      std::cout << "Here is the inversed matrix: " << std::endl << sum_InverseCov_Inverse << std::endl;
    }

    output_2D(0) = sum_theta_backup / inputVtxList.size();
    output_2D(1) = sum_phi_backup / inputVtxList.size();
  }
  else{
    output_2D = sum_InverseCov_Inverse * sum_InverseCovX;
  }

  // convert to 3D vector
  TVector3 output_3D;
  output_3D.SetMagThetaPhi(1., output_2D(0), output_2D(1));

  return FromTVector3ToVector3D(output_3D);
}

int KmeansbbTagTool::HardAssignment(const xAOD::Vertex* primaryVertex, const xAOD::Vertex* vtx, const std::vector<Amg::Vector3D> & class_axis) const{
  int assign_index = -1;
  double assign_chisquare = 9.E+10;

  std::vector<const xAOD::Vertex*> vtxList; 
  vtxList.push_back(vtx);

  for(unsigned int iclass = 0; iclass < class_axis.size(); iclass++){
    double current_chisquare = GetChisquare(primaryVertex, class_axis[iclass], vtxList);
    if(current_chisquare < assign_chisquare){
      assign_chisquare = current_chisquare;
      assign_index = iclass;
    }
  }

  return assign_index;
}

///////////////////////////////////////////////////////////////////////////////////////////////

// following JetTagAlgs/BTagTool.cxx
const xAOD::Vertex* KmeansbbTagTool::RetrievePrimaryVertex() const{
  const xAOD::VertexContainer* vxContainer(0);
  StatusCode sc = evtStore()->retrieve(vxContainer, m_primaryVertxContainerName);

  if(sc.isFailure()){
    ATH_MSG_ERROR("Unable to get primary vertex collection " << m_primaryVertxContainerName);
    return 0;
  }

  // Hmm, it seems we should use retrieve(), instead of tryRetrieve() ... 
  // Updates: OK, for things in AOD, one should use retrieve. But for things in transient, seems tryRetrieve is better ... 
  // auto vxContainer = evtStore()->tryRetrieve<xAOD::VertexContainer>(m_primaryVertxContainerName);
  // if(vxContainer == 0){
  //   ATH_MSG_ERROR("Unable to get primary vertex collection " << m_primaryVertxContainerName);
  //   return 0;
  // }

  if(vxContainer->size() == 0){
    ATH_MSG_ERROR("Empty primary vertex collection " << m_primaryVertxContainerName << "!");
    return 0;
  }

  const xAOD::Vertex* output = 0;

  for(auto vtx : *vxContainer){
    if(vtx->vertexType() == xAOD::VxType::PriVtx){
      output = vtx;
      break;
    }
  }

  if(!output){
    ATH_MSG_WARNING("No vertex labeled as VxType::PriVtx!");
    output = vxContainer->at(0);
    if(output->nTrackParticles() == 0){
      ATH_MSG_WARNING("PV==BeamSpot: probably poor tagging");
    }
  }

  return output;
}

std::vector<const xAOD::Vertex*> KmeansbbTagTool::RetrieveMSVInJet(const xAOD::Jet& jet) const{
  std::vector<ElementLink<xAOD::VertexContainer> > msvVertices;
  try{
    jet.btagging()->variable<std::vector<ElementLink<xAOD::VertexContainer> > >("MSV", "vertices", msvVertices);
  }
  catch(...){
    ATH_MSG_WARNING("Unable to get MSV vertices");
  }

  std::vector<const xAOD::Vertex*> output;
  for(auto el_vtx : msvVertices){
    if(!el_vtx.isValid()){
      ATH_MSG_WARNING("Invalid element link to MSV vtx. It will be skipped.");
      continue;
    }

    const xAOD::Vertex* vtx = (*el_vtx);
    if(vtx == 0){
      ATH_MSG_WARNING("Null ptr returned for MSV vtx. It will be skipped.");
      continue;
    }

    output.push_back(vtx);
  }

  return output;
}

std::vector<const xAOD::TrackParticle*> KmeansbbTagTool::VertexTrack(const xAOD::Vertex* vertex) const{
  std::vector<const xAOD::TrackParticle*> output;

  for(auto el_track : vertex->trackParticleLinks()){
    if(!el_track.isValid()){
      ATH_MSG_WARNING("Invalid element link to tracks from MSV vtx. It will be skipped.");
      continue;
    }

    const xAOD::TrackParticle* track = (*el_track);
    if(track == 0){
      ATH_MSG_WARNING("Null ptr returned for track from MSV vtx. It will be skipped.");
      continue;
    }

    output.push_back(track);
  }

  return output;
}

double KmeansbbTagTool::MaxDeltaR(const xAOD::Jet& jet) const{
  double max_dR = -1;

  auto v_constituents = jet.getConstituents();
  for(xAOD::JetConstituentVector::iterator iter = v_constituents.begin(); iter != v_constituents.end(); iter++){
    TLorentzVector p4_constituent;
    p4_constituent.SetPtEtaPhiE((*iter)->pt(), (*iter)->eta(), (*iter)->phi(), (*iter)->e());

    double dR = jet.p4().DeltaR(p4_constituent);
    if(dR > max_dR) max_dR = dR;
  }

  if(max_dR < 0){
    ATH_MSG_WARNING("Negative maximal dR! This means jet has non constituent ... Jet radius parameter would be used as placeholder.");
    max_dR = jet.getSizeParameter();
  }

  return max_dR;
}

template<class T>
std::vector<T> KmeansbbTagTool::FindCommon(const std::vector<T> & collection1, const std::vector<T> & collection2) const{
  std::vector<T> output;

  for(auto item1 : collection1){
    if(!CheckDuplicate(collection2, item1)) continue;
    if(CheckDuplicate(output, item1)) continue;

    output.push_back(item1);
  }

  return output;
}

template<class T>
bool KmeansbbTagTool::CheckDuplicate(const std::vector<T> & collection, T & item) const{
  return (std::find(collection.begin(), collection.end(), item) != collection.end());
}

Amg::Vector3D KmeansbbTagTool::FromTVector3ToVector3D(TVector3 input) const{
  Amg::Vector3D output;
  output(0) = input(0);
  output(1) = input(1);
  output(2) = input(2);

  return output;
}

TVector3 KmeansbbTagTool::FromVector3DToTVector3(Amg::Vector3D input) const{
  TVector3 output;
  output(0) = input(0);
  output(1) = input(1);
  output(2) = input(2);

  return output;
}

// calculate all eigen variation of vtx position
std::vector<Amg::Vector3D> KmeansbbTagTool::GetVtxEigenVariation(Amg::Vector3D vtxPosition, AmgSymMatrix3D vtxCovMatrix) const{
  std::vector<Amg::Vector3D> output;

  // solve eigen system
  Eigen::SelfAdjointEigenSolver<AmgSymMatrix3D> EigenSolver(vtxCovMatrix);
  if(EigenSolver.info() != Eigen::Success){
    ATH_MSG_WARNING("Input matrix is not Hermitian, which should be impossible for covariant matrix!");
    return output;
  }

  auto EigenValues = EigenSolver.eigenvalues();
  auto EigenVectors = EigenSolver.eigenvectors();

  // get variation
  for(unsigned int index = 0; index < 3; index++){
    Amg::Vector3D vtxPositionVariation_Up;
    vtxPositionVariation_Up(0) = vtxPosition(0) + TMath::Sqrt(EigenValues(index)) * EigenVectors(0, index);
    vtxPositionVariation_Up(1) = vtxPosition(1) + TMath::Sqrt(EigenValues(index)) * EigenVectors(1, index);
    vtxPositionVariation_Up(2) = vtxPosition(2) + TMath::Sqrt(EigenValues(index)) * EigenVectors(2, index);

    Amg::Vector3D vtxPositionVariation_Down;
    vtxPositionVariation_Down(0) = vtxPosition(0) - TMath::Sqrt(EigenValues(index)) * EigenVectors(0, index);
    vtxPositionVariation_Down(1) = vtxPosition(1) - TMath::Sqrt(EigenValues(index)) * EigenVectors(1, index);
    vtxPositionVariation_Down(2) = vtxPosition(2) - TMath::Sqrt(EigenValues(index)) * EigenVectors(2, index);

    output.push_back(vtxPositionVariation_Up);
    output.push_back(vtxPositionVariation_Down);
  }          

  return output;
}

// 
double KmeansbbTagTool::GeneralDRSignificance3D(Amg::Vector3D vtx1, AmgSymMatrix3D vtx1_CovMatrix, Amg::Vector3D vtx2, AmgSymMatrix3D vtx2_CovMatrix, double& dR3D, double& dR3D_error) const{
  dR3D = -1.;
  dR3D_error = -1.;

  std::vector<Amg::Vector3D> vtx1_variations = GetVtxEigenVariation(vtx1, vtx1_CovMatrix);
  std::vector<Amg::Vector3D> vtx2_variations = GetVtxEigenVariation(vtx2, vtx2_CovMatrix);

  auto function_GetdR3D = [](Amg::Vector3D input1, Amg::Vector3D input2)->double{
    return TMath::Sqrt((input1 - input2).dot(input1 - input2));
  };

  dR3D = function_GetdR3D(vtx1, vtx2);

  double sumerror_vtx1 = 0.;
  double sumerror_vtx2 = 0.;
  for(unsigned int index = 0; index < 3; index++){
    double error_vtx1 = std::max(fabs(function_GetdR3D(vtx1_variations[2*index], vtx2) - dR3D), fabs(function_GetdR3D(vtx1_variations[2*index+1], vtx2) - dR3D));
    sumerror_vtx1 += (error_vtx1 * error_vtx1);

    double error_vtx2 = std::max(fabs(function_GetdR3D(vtx1, vtx2_variations[2*index]) - dR3D), fabs(function_GetdR3D(vtx1, vtx2_variations[2*index+1]) - dR3D));
    sumerror_vtx2 += (error_vtx2 * error_vtx2);
  }

  double sumerror = sumerror_vtx1 + sumerror_vtx2;
  dR3D_error = TMath::Sqrt(sumerror);

  double significance = 0.;
  if(dR3D_error == 0.){
    ATH_MSG_WARNING("Zero error obtained in GeneralDRSignificance3D, which is very unlikely.");
    significance = 999.;
  }
  else{
    significance = dR3D/dR3D_error;
  }

  significance = LimitUpperBound(significance, 50, 100);
  return significance;
}

// ATLAS coordination:
// x = r*sin(theta)*cos(phi)
// y = r*sin(theta)*sin(phi)
// z = r*cos(theta)
AmgJacobMatrix3b2 KmeansbbTagTool::GetJacobian_FromSphereToXYZ(Amg::Vector3D vtxPosition) const{
  AmgJacobMatrix3b2 output;

  TVector3 vtxPosition_hep = FromVector3DToTVector3(vtxPosition);
  double r = vtxPosition_hep.Mag();
  double theta = vtxPosition_hep.Theta();
  double phi = vtxPosition_hep.Phi();

  output(0, 0) = r * TMath::Cos(theta) * TMath::Cos(phi); output(0, 1) = -r * TMath::Sin(theta) * TMath::Sin(phi);
  output(1, 0) = r * TMath::Cos(theta) * TMath::Sin(phi); output(1, 1) = r * TMath::Sin(theta) * TMath::Cos(phi);
  output(2, 0) = -r * TMath::Sin(theta);                  output(2, 1) = 0.;

  return output;
}

// ATLAS coordination:
// theta = arctan(sqrt(x^2+y^2)/z)
// phi = arctan(y/x)
AmgJacobMatrix2b3 KmeansbbTagTool::GetJacobian_FromXYZToSphere(Amg::Vector3D vtxPosition) const{
  AmgJacobMatrix2b3 output;

  TVector3 vtxPosition_hep = FromVector3DToTVector3(vtxPosition);
  double r = vtxPosition_hep.Mag();
  double perp = vtxPosition_hep.Perp();
  double x = vtxPosition_hep.x();
  double y = vtxPosition_hep.y();
  double z = vtxPosition_hep.z();

  if( (r == 0.) || (perp == 0.) ){
    ATH_MSG_WARNING("All zero in GetJacobian_FromXYZToSphere.");
    return output;
  }
  else{
    output(0, 0) = (x * z) / (perp * r * r); output(0, 1) = (y * z) / (perp * r *r); output(0, 2) = -perp / (r * r);
    output(1, 0) = -y / (perp * perp);       output(1, 1) = x / (perp * perp);       output(1, 2) = 0;

    return output;
  }
}

AmgSymMatrix2D KmeansbbTagTool::GetDirectionCovMatrix(Amg::Vector3D vtxPosition, AmgSymMatrix3D vtxCovMatrix) const{
  auto J = GetJacobian_FromXYZToSphere(vtxPosition);
  return J * vtxCovMatrix * J.transpose();
}

template<class T, class TI>
TI KmeansbbTagTool::Inverse(T inputMatrix, bool& is_invertible) const{
  // dummy output
  TI empty;
  empty.setZero();

  // zero matrix
  if(inputMatrix.isZero(1e-10)){
    ATH_MSG_WARNING("Matrix for inversion is a zero matrix! Zero matrix would be returned.");
    is_invertible = false;
    return empty;
  }

  // scale
  double scale = 1.;
  auto ref_number = std::max(std::abs(inputMatrix.maxCoeff()), std::abs(inputMatrix.minCoeff()));
  if(ref_number < 0.01) scale = 10./ref_number;
  T inputMatrix_scaled = scale * inputMatrix;

  // get inverse
  TI inverse;
  inputMatrix_scaled.computeInverseWithCheck(inverse, is_invertible);

  // scale back
  return scale * inverse;
}

void KmeansbbTagTool::DumpKmeansElement(const xAOD::Vertex* primaryVertex, xAOD::Jet& jetToTag, KmeansElement& result, std::string description) const{

  std::cout << "=========== Dumping Kmeans Info =============" << std::endl;

  ////////////////////////////////////////
  // Making Demo Plot for visualization //
  ////////////////////////////////////////

  const xAOD::EventInfo* evtInfo(0);
  StatusCode sc = evtStore()->retrieve(evtInfo, "EventInfo");
  if(sc.isFailure()){
    ATH_MSG_ERROR("Unable to get EventInfo!");
    return;
  }

  std::string plotName( TString::Format("VisualKmeansElement_EventNumber%i_JetIndex%i_%s", int(evtInfo->eventNumber()), int(jetToTag.index()), name().data()).Data() );
  if(description != ""){
    plotName = plotName + "_" + description;
  }
  DumpKmeansElement__Visual(plotName, primaryVertex, jetToTag, result);

  ////////////////////////////////
  // Raw information for record //
  ////////////////////////////////

  DumpKmeansElement__Raw(primaryVertex, jetToTag, result);

  std::cout << "=============== End of Dump =================" << std::endl;

}

void KmeansbbTagTool::DumpKmeansElement__Visual(std::string name, const xAOD::Vertex* primaryVertex, xAOD::Jet& jetToTag, KmeansElement& result) const{

  // determine the boundary of graph

  std::vector<double> x;
  std::vector<double> y;

  double x_low = 99999.;
  double x_up  = -99999.;
  double y_low = 99999.;
  double y_up  = -99999.;
  double radius_max = -1.;

  for(unsigned int iclass = 0; iclass < result.class_vtxList.size(); iclass++){
    for(auto vtx : result.class_vtxList[iclass]){
      Amg::Vector3D SV_PV = vtx->position() - primaryVertex->position();

      x.push_back(SV_PV(0));
      y.push_back(SV_PV(1));

      if(x.back() < x_low) x_low = x.back();
      if(x.back() > x_up)  x_up  = x.back();
      if(y.back() < y_low) y_low = y.back();
      if(y.back() > y_up)  y_up  = y.back();

      double radius = TMath::Sqrt(x.back()*x.back() + y.back()*y.back());
      if(radius > radius_max) radius_max = radius;
    }
  }

  double x_boundary = ((std::abs(x_low) > std::abs(x_up)) ? std::abs(x_low) : std::abs(x_up));
  double y_boundary = ((std::abs(y_low) > std::abs(y_up)) ? std::abs(y_low) : std::abs(y_up));

  // initialize canvas

  TCanvas c(name.data(), name.data(), 800, 600);
  c.cd();

  // prepare multigraph

  TMultiGraph mg("multigraph", "multigraph");

  // draw input points

  std::vector<TGraph*> v_graph_inputs;
  for(unsigned int iclass = 0; iclass < result.class_vtxList.size(); iclass++){
    std::vector<double> class_x;
    std::vector<double> class_y;
    for(auto vtx : result.class_vtxList[iclass]){
      Amg::Vector3D SV_PV = vtx->position() - primaryVertex->position();

      class_x.push_back(SV_PV(0));
      class_y.push_back(SV_PV(1));
    }

    TGraph* graph_inputs = new TGraph(class_x.size(), &class_x[0], &class_y[0]);
    graph_inputs->SetMarkerColor(2 + iclass);
    graph_inputs->SetMarkerStyle(22);
    graph_inputs->SetMarkerSize(1.5);
    graph_inputs->SetName(TString::Format("Inputs_class%i", iclass).Data());

    v_graph_inputs.push_back(graph_inputs);
    mg.Add(graph_inputs);
  }

  // draw truth bhadron 

  std::vector<const xAOD::TruthParticle*> Bhadrons_truth;
  try{
    Bhadrons_truth = jetToTag.getAssociatedObjects<xAOD::TruthParticle>("GhostBHadronsFinal");

    // Only works by using getAssociatedObjects. Not work by trying to access vector of EL ... 
    // Keep this for future reference
    // auto v_el_Bhadrons = jetToTag.auxdata<std::vector<ElementLink<xAOD::TruthParticleContainer> > >("GhostBHadronsFinal");

    // for(auto el_Bhadron : v_el_Bhadrons){
    //   if(!el_Bhadron.isValid()){
    //     ATH_MSG_WARNING("Invalid element link to GhostBHadronsFinal. It will be skipped");
    //     continue;
    //   }
    //   const xAOD::TruthParticle* Bhadron_truth = (*el_Bhadron);
    //   if(Bhadron_truth == 0){
    //     ATH_MSG_WARNING("Null ptr returned for truth Bhadron. It will be skipped");
    //     continue;
    //   }
    //   Bhadrons_truth.push_back(Bhadron_truth);
    // }

  }
  catch(...){
    ATH_MSG_WARNING("Failure in retrieving truth Bhadrons.");
  }

  std::vector<double> mc_vx_x;
  std::vector<double> mc_vx_y;
  std::vector<TLine>  truth_bhadrons_axis;
  for(auto Bhadron : Bhadrons_truth){
    // axis
    TLine arrow(0., 0., radius_max*TMath::Cos(Bhadron->p4().Phi()), radius_max*TMath::Sin(Bhadron->p4().Phi()));
    arrow.SetLineStyle(1);
    arrow.SetLineWidth(2);
    arrow.SetLineColor(4);
    truth_bhadrons_axis.push_back(arrow);

    // bhadron decay vtx
    const xAOD::TruthVertex* decayVtx = Bhadron->decayVtx();
    Amg::Vector3D BhadronVtx;
    if(decayVtx == 0){
      ATH_MSG_WARNING("Unable to get MC truth decay vertex of Bhadron. A dummy at origin would be set");
      BhadronVtx = primaryVertex->position() - primaryVertex->position();
    }
    else{
      BhadronVtx(0) = decayVtx->x() - primaryVertex->x();
      BhadronVtx(1) = decayVtx->y() - primaryVertex->y();
      BhadronVtx(2) = decayVtx->z() - primaryVertex->z();
    }

    mc_vx_x.push_back(BhadronVtx(0));
    mc_vx_y.push_back(BhadronVtx(1));
  }

  TGraph* mc_vx_graph = 0;
  if(mc_vx_x.size() > 0){
    mc_vx_graph = new TGraph(mc_vx_x.size(), &mc_vx_x[0], &mc_vx_y[0]);
    mc_vx_graph->SetMarkerColor(4);
    mc_vx_graph->SetMarkerStyle(29);
    mc_vx_graph->SetMarkerSize(2.5);
    mc_vx_graph->SetName("mc_vx");
    mg.Add(mc_vx_graph);
  }

  // we draw multigraph

  mg.Draw("AP");

  // mg.GetXaxis()->SetLimits(-x_boundary-5, x_boundary+5);
  // mg.GetYaxis()->SetLimits(-y_boundary-5, y_boundary+5);
  // mg.GetXaxis()->SetRangeUser(-x_boundary-5, x_boundary+5);
  // mg.GetYaxis()->SetRangeUser(-y_boundary-5, y_boundary+5);

  x_boundary *= 1.5;
  y_boundary *= 1.5;
  mg.GetXaxis()->SetLimits(-x_boundary, x_boundary);
  mg.GetYaxis()->SetLimits(-y_boundary, y_boundary);
  mg.GetXaxis()->SetRangeUser(-x_boundary, x_boundary);
  mg.GetYaxis()->SetRangeUser(-y_boundary, y_boundary);

  // then we actually draw truth bhadron axis here
  for(unsigned int index = 0; index < truth_bhadrons_axis.size(); index++){
    truth_bhadrons_axis[index].Draw();
  }

  // draw fitted axis
  // NOTE: This is NOT necessarily the final reconstructed vertex direction

  for(unsigned int iclass = 0; iclass < result.class_axis.size(); iclass++){
    TVector3 axis_HEP = FromVector3DToTVector3(result.class_axis[iclass]);

    TLine arrow;
    arrow.SetLineStyle(2);
    arrow.SetLineWidth(2);
    arrow.SetLineColor(2 + iclass);

    arrow.DrawLine(0., 0., radius_max*TMath::Cos(axis_HEP.Phi()), radius_max*TMath::Sin(axis_HEP.Phi()));
  }

  // jet cone

  // jet-axis first

  TLine arrow_JetAxis;
  arrow_JetAxis.SetLineStyle(2);
  arrow_JetAxis.SetLineColor(6);
  arrow_JetAxis.SetLineWidth(1);
  arrow_JetAxis.DrawLine(0., 0., radius_max*TMath::Cos(jetToTag.p4().Phi()), radius_max*TMath::Sin(jetToTag.p4().Phi()));

  // jet-boundary
  double JetRadius = MaxDeltaR(jetToTag);

  TLine edge1; edge1.SetLineStyle(2); edge1.SetLineColor(6); edge1.SetLineWidth(1);
  edge1.DrawLine(0., 0., radius_max*TMath::Cos(jetToTag.p4().Phi() - JetRadius), radius_max*TMath::Sin(jetToTag.p4().Phi() - JetRadius));

  TLine edge2; edge2.SetLineStyle(2); edge2.SetLineColor(6); edge2.SetLineWidth(1);
  edge2.DrawLine(0., 0., radius_max*TMath::Cos(jetToTag.p4().Phi() + JetRadius), radius_max*TMath::Sin(jetToTag.p4().Phi() + JetRadius));

  TLine edge3; edge3.SetLineStyle(2); edge3.SetLineColor(6); edge3.SetLineWidth(1);
  edge3.DrawLine(radius_max*TMath::Cos(jetToTag.p4().Phi() - JetRadius), radius_max*TMath::Sin(jetToTag.p4().Phi() - JetRadius), radius_max*TMath::Cos(jetToTag.p4().Phi() + JetRadius), radius_max*TMath::Sin(jetToTag.p4().Phi() + JetRadius));

  c.Update();

  c.SaveAs((name+".root").data());
  c.SaveAs((name+".pdf").data());

  // delete
  for(auto g : v_graph_inputs){
    delete g;
    g = 0;
  }

  delete mc_vx_graph;
  mc_vx_graph = 0;

}

// TODO
void KmeansbbTagTool::DumpKmeansElement__Raw(const xAOD::Vertex* primaryVertex, xAOD::Jet& jetToTag, KmeansElement& result) const{

  std::cout << "================================================" << std::endl
            << "Jet Index: " << jetToTag.index() << std::endl
            << "Primary vertex position: " << primaryVertex->position() << std::endl
            << "- Number of clusters: " << result.class_axis.size() << std::endl;

  for(unsigned int iclass = 0; iclass < result.class_axis.size(); iclass++){
    std::cout << "- Print out of cluster " << iclass << std::endl
              << "- - Axis: " << result.class_axis[iclass] << std::endl
              << "- - Vtx: " << std::endl;

    for(auto vtx : result.class_vtxList[iclass]){
      std::cout << "- - - " << vtx->position() << std::endl;
    }
  }

  std::cout << "================================================" << std::endl;

}


///////////////////////////////////////////////////////////////////////////////////////////////

}//end namespace



