/* **************************************************************************
   ParticleToJetAssociator.h  -  Description
   -------------------
begin   : 10.09.04
authors : Kyle Cranmer
email   : kyle.cranmer@cern.ch

comments: A. Wildauer: 27.01.2005 rewritten into an AlgTool

 ***************************************************************************/

#ifndef PARTICLEJETTOOLS_PARTICLETOJETASSOCIATOR_H
#define PARTICLEJETTOOLS_PARTICLETOJETASSOCIATOR_H

#include "AsgTools/AsgTool.h"
#include "AsgTools/ToolHandle.h"
//#include "AnalysisTools/AnalysisTools.h"
// #include "JetEvent/JetCollection.h"
// #include "JetEvent/JetAssociationBase.h"

#include "xAODJet/Jet.h"

#include <vector>
#include <list>
#include <string>

namespace Analysis
{

typedef std::string             NameType;
typedef std::vector<xAOD::Jet*> jetcollection_t;

#ifndef ROOTCORE
  static const InterfaceID IID_ParticleToJetAssociator("Analysis::ParticleToJetAssociator", 1, 0);
#endif
typedef std::vector<std::string> StringVector;

/**  Class ParticleToJetAssociator.  This athena AlgTool
  reads the Jet collection previously built by a jet finder
  (ex. KtJet) and a collection of particles.  It then makes new
  jets by adding particles in a cone or kt region.*/

class ParticleToJetAssociator : public asg::AsgTool {
    ASG_TOOL_CLASS0(ParticleToJetAssociator)
    public:

        /** Constructors and destructors */
        ParticleToJetAssociator(const std::string& name);
        virtual ~ParticleToJetAssociator();

        /** AlgTool interface methods (only in Athena) */
#ifndef ROOTCORE
        static const InterfaceID& interfaceID() { return IID_ParticleToJetAssociator; };
#endif

        /** AlgTool initailize method */
        StatusCode initialize();
        /** AlgTool finalize method */
        StatusCode finalize();

        // template<typename ConstituentType, typename coll> void associateParticlesToJets(
        //                       JetCollection* SGParticleJetContainer,
        //                       const coll* particleContainer,
        //                       const std::string& constituentName
        // 			    );

        // For the new way
        // FF: return the vector of associations 
        template<typename ConstituentType, typename coll>
            std::vector<ConstituentType*>
            associateParticlesToJets(jetcollection_t* theJets,
                    const coll* particleContainer,
                    const std::string& constituentName);

        inline double getConeSize(double pt) {	    //!< helper method to get the cone size for a given jet pT
            return (m_coneSizeFitPar1 + exp(m_coneSizeFitPar2 + m_coneSizeFitPar3*pt));
        }

    private:
//        ToolHandle<AnalysisTools> m_analysisTools;
        double m_trackCone;		    //!< size in deltaR of cone for association
        bool   m_useVariableSizedTrackCone;   //!< if true uses a deltaR=f(pT)
        std::string m_existLabel;            //!< Directly use an existing label for association
        double m_coneSizeFitPar1;		    //!< parameters defining f function
        double m_coneSizeFitPar2;
        double m_coneSizeFitPar3;
        bool   m_shareTracks;		    //!< if true allows a track to belong to many jets
}; // End class


//  template<typename ConstituentType, typename coll> inline 
//    void ParticleToJetAssociator::associateParticlesToJets( JetCollection* SGParticleJetContainer,
// 							    const coll* particleContainer,
// 							    const std::string& constituentName
// 							    ) {

//  ATH_MSG_VERBOSE("Number of Jets      : " << SGParticleJetContainer->size());
//  if(SGParticleJetContainer->size()<1) return;
//  ATH_MSG_VERBOSE("Number of Particles of Type " << constituentName << " : " << particleContainer->size());

//  // every jet needs TrackConstiuents (do it like this to avoid const clashes
//  // in ParticleJet::constituent(const std::string& key ) CONST
//  std::vector<ConstituentType*> vecTrackConst;
//  for (JetCollection::iterator itr = SGParticleJetContainer->begin();
//       itr != SGParticleJetContainer->end(); ++itr) {
//    NameType name(constituentName);
//    vecTrackConst.push_back(new ConstituentType(name));
//  }

//  int index(0);
//  if(!m_shareTracks) {
//    // Case 1: associate a track to only one jet:
//    // loop through track particles and add to closest jet
//    typename coll::const_iterator tpItr = particleContainer->begin();
//    double deltaRMatch;
//    for (; tpItr!=particleContainer->end() ; tpItr++) {
// 	double pt = (*tpItr)->pt();
// 	ATH_MSG_VERBOSE("Particle-Jet association (no sharing). Jet pT= " << pt);
// 	bool findAMatch = m_analysisTools->matchR((*tpItr), SGParticleJetContainer, index, deltaRMatch);
// 	ATH_MSG_VERBOSE(" matching: find= " << findAMatch << " dR=" << deltaRMatch);
// 	double trackCone = m_trackCone;
// 	if( m_useVariableSizedTrackCone ) {
// 	  trackCone = this->getConeSize(SGParticleJetContainer->at(index)->pt());
// 	  ATH_MSG_VERBOSE("Using Variable Sized Cone of " << trackCone); 
// 	} else {
// 	  ATH_MSG_VERBOSE("Using Fixed Sized Cone of " << trackCone);
// 	}
// 	if (findAMatch && deltaRMatch < trackCone ) {
// 	  vecTrackConst[index]->set_association(particleContainer, (*tpItr), 0.);
// 	  ATH_MSG_VERBOSE(" added constituent to Jet " << index);
// 	}
//    }

//  } else {
//    // Case 2: can associate a track to many jets:
//    // loop on jets and tracks
//    JetCollection::const_iterator jetItr = SGParticleJetContainer->begin();
//    JetCollection::const_iterator jetEnd = SGParticleJetContainer->end();
//    for(; jetItr!=jetEnd; jetItr++) {
//      typename coll::const_iterator tpItr = particleContainer->begin();
//      typename coll::const_iterator tpEnd = particleContainer->end();
//      for(; tpItr!=tpEnd; tpItr++) {
// 	double dR = m_analysisTools->deltaR( (*jetItr), (*tpItr) );
// 	double trackCone = m_trackCone;
// 	if( m_useVariableSizedTrackCone ) {
// 	  trackCone = this->getConeSize((*jetItr)->pt());
// 	  ATH_MSG_VERBOSE("Using Variable Sized Cone of " << trackCone); 
// 	} else {
// 	  ATH_MSG_VERBOSE("Using Fixed Sized Cone of " << trackCone);
// 	}
// 	if(dR<trackCone) {
// 	  ATH_MSG_VERBOSE( "Particle-jet association (with sharing): "
// 	       << " jet pt,phi,eta= " << (*jetItr)->pt() << " " << (*jetItr)->phi() << " " << (*jetItr)->eta()
// 	       << " trk pt,phi,eta= " << (*tpItr)->pt() << " " << (*tpItr)->phi() << " " << (*tpItr)->eta() );
// 	  vecTrackConst[index]->set_association( particleContainer, (*tpItr), 0);
// 	}
//      }
//      index++;
//    }
//  }
//  index=0;
//  for (JetCollection::iterator itr = SGParticleJetContainer->begin();
//       itr != SGParticleJetContainer->end(); ++itr) {
//    (*itr)->setAssociation(vecTrackConst[index]);
//    index++;
//  }
//  return;
// }

// For the new way
template<typename ConstituentType, typename coll>
    inline std::vector<ConstituentType*>
    ParticleToJetAssociator::associateParticlesToJets(jetcollection_t* theJets,
            const coll* particleContainer,
            const std::string& constituentName) {

        std::vector<ConstituentType*> vecTrackConst;
        ATH_MSG_VERBOSE("Number of Jets      : " << theJets->size());
        if (theJets->size()<1) return vecTrackConst;
        ATH_MSG_VERBOSE("Number of Particles of Type " << constituentName << " : " << particleContainer->size());

        for (jetcollection_t::iterator itr = theJets->begin(); itr != theJets->end(); ++itr) {
            //    NameType name(constituentName);
            vecTrackConst.push_back(new ConstituentType);
        }

        if(m_existLabel != ""){ // use existing label directly
            
            // deduce the value_type without const and pointer ... 
            // It is a bit dangerous here, since that means we are expecting a vector of const pointer for ConstituentType, which is not necessarily the case ... (though it is the case at this moment)
            typedef typename ConstituentType::value_type const_T_ptr;
            typedef typename std::remove_pointer<const_T_ptr>::type const_T;
            typedef typename std::remove_const<const_T>::type T;

            int index(0);
            jetcollection_t::iterator jetItr = theJets->begin();
            jetcollection_t::iterator jetEnd = theJets->end();

            for(; jetItr!=jetEnd; jetItr++){
                auto jet = (*jetItr);

                ConstituentType v_assoc_particles;
                v_assoc_particles.clear();

                try{ // getAssociatedObjects
                    ATH_MSG_DEBUG("Trying getAssociatedObjects ...");
                    v_assoc_particles = jet->getAssociatedObjects<T>(m_existLabel);
                }
                catch(...){
                    v_assoc_particles.clear();
                    try{ // direct access vector of element link
                        ATH_MSG_DEBUG("Trying getAttribute ...");
                        auto v_el_obj = jet->getAttribute<std::vector<ElementLink<DataVector<T> > > >(m_existLabel);
                        for(auto el_obj : v_el_obj){
                            if(!el_obj.isValid()){
                                ATH_MSG_WARNING("Invalid element link to object. Skipped.");
                                continue;
                            }

                            const T* obj = (*el_obj);

                            if(obj == 0){
                                ATH_MSG_WARNING("Null ptr returned for object. Skipped.");
                                continue;
                            }

                            if(obj->container() != particleContainer){
                                ATH_MSG_WARNING("Object not from expected container. Skipped");
                                continue;
                            }

                            v_assoc_particles.push_back(obj);
                        }
                    }
                    catch(...){
                        v_assoc_particles.clear();
                        ATH_MSG_WARNING("Unable to find associated object " << m_existLabel << ". Empty collection will be returned.");
                    }
                }

                for(auto track : v_assoc_particles){
                    // making sure there is no duplicated objects in the vector.
                    // Otherwise, some b-tagger (especially JetFitter) will crash
                    if(std::find(vecTrackConst[index]->begin(), vecTrackConst[index]->end(), track) != vecTrackConst[index]->end()){
                        ATH_MSG_DEBUG("Duplicate obj found. It would be skipped.");
                        continue;
                    }

                    vecTrackConst[index]->push_back(track);
                }

                index++;
            }

        }
        else{ // cone-based search

            int index(0);
            if(!m_shareTracks) {
                // Case 1: associate a track to only one jet:
                // loop through track particles and add to closest jet
                typename coll::const_iterator tpItr = particleContainer->begin();
                for (; tpItr!=particleContainer->end() ; tpItr++) {
                    double pt = (*tpItr)->pt();
                    double deltaRMatch(1000);
                    ATH_MSG_VERBOSE("Particle-Jet association (no sharing). Jet pT= " << pt);
                    bool findAMatch = false;
                    int i(0);
                    xAOD::Jet *theJet = 0;
                    for (jetcollection_t::iterator itr = theJets->begin(); itr != theJets->end(); ++itr) {
                        double r = (*itr)->p4().DeltaR((*tpItr)->p4());
                        if ( r < deltaRMatch ) {
                            deltaRMatch = r;
                            index = i;
                            findAMatch = true;
                            theJet = (*itr);
                        }
                        i++;
                    }

                    ATH_MSG_VERBOSE(" matching: find= " << findAMatch << " dR=" << deltaRMatch);
                    double trackCone = m_trackCone;
                    if( m_useVariableSizedTrackCone ) {
                        if(0 != theJet) trackCone = getConeSize(theJet->pt());
                        ATH_MSG_VERBOSE("Using Variable Sized Cone of " << trackCone); 
                    } else {
                        ATH_MSG_VERBOSE("Using Fixed Sized Cone of " << trackCone);
                    }
                    if (findAMatch && deltaRMatch < trackCone ) {
                        // FF: changed this to assume an STL container or at least push_back() functionality
                        //     note that this does not have any mechanism for updating particle weights etc.
                        // vecTrackConst[index]->set_association(particleContainer, (*tpItr), 0.);
                        vecTrackConst[index]->push_back(*tpItr);
                        ATH_MSG_VERBOSE(" added constituent to Jet " << index);
                    }
                }
            } else {
                // Case 2: can associate a track to many jets:
                // loop on jets and tracks
                jetcollection_t::iterator jetItr = theJets->begin();
                jetcollection_t::iterator jetEnd = theJets->end();
                for(; jetItr!=jetEnd; jetItr++) {
                    typename coll::const_iterator tpItr = particleContainer->begin();
                    typename coll::const_iterator tpEnd = particleContainer->end();
                    for(; tpItr!=tpEnd; tpItr++) {
                        double dR = (*jetItr)->p4().DeltaR((*tpItr)->p4());
                        double trackCone = m_trackCone;
                        if( m_useVariableSizedTrackCone ) {
                            trackCone = getConeSize((*jetItr)->pt());
                            ATH_MSG_VERBOSE("Using Variable Sized Cone of " << trackCone); 
                        } else {
                            ATH_MSG_VERBOSE("Using Fixed Sized Cone of " << trackCone);
                        }
                        if(dR<trackCone) {
                            ATH_MSG_VERBOSE("Particle-jet association (with sharing): "
                                    << " jet pt,phi,eta= " << (*jetItr)->pt() << " " << (*jetItr)->phi() << " " << (*jetItr)->eta()
                                    << " trk pt,phi,eta= " << (*tpItr)->pt() << " " << (*tpItr)->phi() << " " << (*tpItr)->eta());
                            // FF: changed this to assume an STL container or at least push_back() functionality
                            //     note that this does not have any mechanism for updating particle weights etc.
                            // vecTrackConst[index]->set_association( particleContainer, (*tpItr), 0);
                            vecTrackConst[index]->push_back(*tpItr);
                        }
                    }
                    index++;
                }
            }

        }
        // FF: the code below cannot be used for xAOD::Jet
        // index=0;
        // for (jetcollection_t::iterator itr = theJets->begin(); itr != theJets->end(); ++itr) {
        //     (*itr)->setAssociation( vecTrackConst[index] );
        //     index++;
        //   } 

        return vecTrackConst;
    }


} // End namespace
#endif
