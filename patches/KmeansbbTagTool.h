#ifndef JETTAGTOOLS_KMEANSBBTAGTOOL_H
#define JETTAGTOOLS_KMEANSBBTAGTOOL_H

/******************************************************
    @class KmeansbbTagTool
    Identification of double b-hadron jets using k-means clusterization of seed vertices
    This class only serves to do clusterization, and build "subjet" based on it. 
    A separate tagger, KmenasbbTag, is used to collect all information for final discrimination
    @ author Qi Zeng
********************************************************/

// base class
#include "GaudiKernel/ToolHandle.h"
#include "JetRec/JetModifierBase.h"

// std
#include <string>

// xAOD EDM
#include "xAODTracking/VertexContainer.h"
#include "xAODTracking/TrackParticleContainer.h"
#include "xAODJet/JetContainer.h"
#include "xAODJet/JetAuxContainer.h"
#include "xAODBTagging/BTaggingContainer.h"
#include "xAODBTagging/BTagVertexContainer.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthVertex.h"
#include "xAODEventInfo/EventInfoContainer.h"

// ROOT
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLine.h"
#include "TMultiGraph.h"
#include "TString.h"
#include "TAxis.h"

typedef Eigen::Matrix<double, 3, 3, 0, 3, 3> AmgSymMatrix3D;
typedef Eigen::Matrix<double, 2, 2, 0, 2, 2> AmgSymMatrix2D;
typedef Eigen::Matrix<double, 3, 2, 0, 3, 2> AmgJacobMatrix3b2;
typedef Eigen::Matrix<double, 2, 3, 0, 2, 3> AmgJacobMatrix2b3;

namespace Analysis { 

  struct KmeansElement{
    std::vector<std::vector<const xAOD::Vertex*> >       class_vtxList;      // vtx to be clustered
    std::vector<Amg::Vector3D>                           class_axis;         // w.r.t. ref vtx
    std::vector<double>                                  class_chisquare;    // w.r.t. ref vtx
  };

  class KmeansbbTagTool : public JetModifierBase {
    ASG_TOOL_CLASS(KmeansbbTagTool, IJetModifier)
   
  public:
    // Constructor from tool name.
    KmeansbbTagTool(std::string myname);

    virtual StatusCode initialize();

    // Inherited method to modify a jet.
    virtual int modifyJet(xAOD::Jet& jet) const;


  private:

    // main sequence //

    // clean up of input vertex list
    std::vector<const xAOD::Vertex*>  InputVtxCleaning(const xAOD::Vertex* primaryVertex, const xAOD::Jet & jet, const std::vector<const xAOD::Vertex*> & inputVtxList) const;
    // repeating k-means adaptively to remove outliers
    KmeansElement                     AdaptiveKmeans(const xAOD::Vertex* primaryVertex, std::vector<const xAOD::Vertex*> & inputVtxList, std::vector<KmeansElement> & history, std::vector<std::vector<KmeansElement> > & detailHistory, std::vector<Amg::Vector3D> ClusterSeeds) const;
    // core k-means part
    std::vector<Amg::Vector3D>        InitializeKmeans(const xAOD::Vertex* primaryVertex, std::vector<const xAOD::Vertex*> inputVtxList, std::vector<Amg::Vector3D> ClusterSeeds, bool doKmeansPP = true) const;
    KmeansElement                     SimpleKmeans(const xAOD::Vertex* primaryVertex, const std::vector<const xAOD::Vertex*> & inputVtxList, std::vector<KmeansElement> & ClusHistory, std::vector<Amg::Vector3D> ClusterSeeds) const;

    // kernel clustering calculation //

    double         GetChisquare(const xAOD::Vertex* primaryVertex, const Amg::Vector3D & axis, const std::vector<const xAOD::Vertex*> & inputVtxList) const;
    Amg::Vector3D  UpdateAxis(const xAOD::Vertex* primaryVertex, const std::vector<const xAOD::Vertex*> & inputVtxList) const;
    int            HardAssignment(const xAOD::Vertex* primaryVertex, const xAOD::Vertex* vtx, const std::vector<Amg::Vector3D> & class_axis) const;

    // utils //

    const xAOD::Vertex*                       RetrievePrimaryVertex() const;
    std::vector<const xAOD::Vertex*>          RetrieveMSVInJet(const xAOD::Jet & jet) const;
    std::vector<const xAOD::TrackParticle*>   VertexTrack(const xAOD::Vertex* vertex) const;
    double                                    MaxDeltaR(const xAOD::Jet & jet) const;

    template<class T> 
    std::vector<T>                            FindCommon(const std::vector<T> & collection1, const std::vector<T> & collection2) const;

    template<class T>
    bool                                      CheckDuplicate(const std::vector<T> & collection, T & item) const;

    std::vector<const xAOD::TrackParticle*>   TrackCollectionRemoveVtxs(const std::vector<const xAOD::TrackParticle*> & inputTracks, const std::vector<const xAOD::Vertex*> & removeVtxs) const;
    Amg::Vector3D                             FromTVector3ToVector3D(TVector3 input) const;
    TVector3                                  FromVector3DToTVector3(Amg::Vector3D input) const;
    static bool                               SortPt(const xAOD::Jet* j1, const xAOD::Jet* j2) {return j1->pt() > j2->pt();}

    double LimitUpperBound(double input, double start, double limit) const{
      if(input > start) return start + (limit-start)/TMath::Pi()*2.*TMath::ATan(TMath::Pi()/2./(limit-start)*(input-start));
      else              return input;
    }
    std::vector<Amg::Vector3D>                GetVtxEigenVariation(Amg::Vector3D vtxPosition, AmgSymMatrix3D vtxCovMatrix) const;
    double                                    GeneralDRSignificance3D(Amg::Vector3D vtx1, AmgSymMatrix3D vtx1_CovMatrix, Amg::Vector3D vtx2, AmgSymMatrix3D vtx2_CovMatrix, double& dR3D, double& dR3D_error) const;

    AmgJacobMatrix3b2                         GetJacobian_FromSphereToXYZ(Amg::Vector3D vtxPosition) const;
    AmgJacobMatrix2b3                         GetJacobian_FromXYZToSphere(Amg::Vector3D vtxPosition) const;
    AmgSymMatrix2D                            GetDirectionCovMatrix(Amg::Vector3D vtxPosition, AmgSymMatrix3D vtxCovMatrix) const;

    // a wrapper of Eigen inverse function, that provides:
    // 1. scaling the covaraint matrix for numerical purpose
    // caveat: 
    // 1. always assuming this is a two-dimensional matrix
    // 2. Up to 4*4
    template<class T, class TI>
    TI                                        Inverse(T inputMatrix, bool& is_invertible) const;

    void                                      DumpKmeansElement(const xAOD::Vertex* primaryVertex, xAOD::Jet& jetToTag, KmeansElement& result, std::string description = "") const;
    void                                      DumpKmeansElement__Visual(std::string name, const xAOD::Vertex* primaryVertex, xAOD::Jet& jetToTag, KmeansElement& result) const;
    void                                      DumpKmeansElement__Raw(const xAOD::Vertex* primaryVertex, xAOD::Jet& jetToTag, KmeansElement& result) const;


    // options //

    bool        m_debug;

    std::string m_primaryVertxContainerName;
    int         m_nAxis;
    std::string m_subjetLabel;
    bool        m_useJetSeed;
    double      m_JetRadiusCut;
    double      m_PVDistance;
    bool        m_InheritSeed;
    double      m_MaxChi2;

    std::string m_inputJetContainerName;
    std::string m_recordSubjetLabel;
    std::string m_recordSubjetContainerName;
    std::string m_recordSubjetAlgorithm_BTAG;
    float       m_recordSubjetRadius_BTAG;


    // class members //


  }; // End class

} // End namespace 

#endif
