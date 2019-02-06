//
// Created by Michael Heuer on 2019-02-02.
//

#include "MotifAnalysis.h"


namespace MotifAnalysis {

    std::string toString(MotifType type) {
        switch(type) {
            case MotifType::Core :
                return "Core";
            case MotifType::Valence :
                return "Valence";
            default:
                return "unassigned";
        }
    }

    MotifType fromString(const std::string& string) {
        if(string == "Core")
            return MotifType::Core;
        else if (string == "Valence")
            return MotifType::Valence;
        else
            return MotifType::unassigned;
    }


    Motif::Motif(const std::list<Eigen::Index> &electronIndices, MotifType type)
            : type_(type), electronIndices_(electronIndices) {};

    bool Motif::containsQ(Eigen::Index i) const {
        return std::find(electronIndices_.begin(), electronIndices_.end(), i) != electronIndices_.end();
    }

    bool Motif::operator<(const Motif &rhs) const {
        return electronIndices_ < rhs.electronIndices_;
    }

    bool Motif::operator>(const Motif &rhs) const {
        return rhs < *this;
    }

    bool Motif::operator<=(const Motif &rhs) const {
        return !(rhs < *this);
    }

    bool Motif::operator>=(const Motif &rhs) const {
        return !(*this < rhs);
    }

    MotifType Motif::type() const {
        return type_;
    }

    void Motif::setType(MotifType type_) {
        Motif::type_ = type_;
    }

    const std::list<Eigen::Index> &Motif::electronIndices() const {
        return electronIndices_;
    }

    void Motif::setElectronIndices(const std::list<Eigen::Index> &electronIndices_) {
        Motif::electronIndices_ = electronIndices_;
    }

    bool coreElectronQ(const Electron& e, const AtomsVector& atoms, double threshold) {
        for (Eigen::Index k = 0; k < atoms.numberOfEntities(); ++k)
            if(Metrics::distance(e.position(), atoms[k].position()) <= threshold)
                return true;

        return false;
    }

    Motifs::Motifs(const Eigen::MatrixXb &adjacencyMatrix)
            : motifVector() {
        auto electronIndicesLists = GraphAnalysis::findGraphClusters(adjacencyMatrix);

        for (const auto &list : electronIndicesLists)
            motifVector.emplace_back(list);
    };

    Motifs::Motifs(std::vector<Motif> motifs)
            : motifVector(std::move(motifs)) {
    };

    void Motifs::classifyMotifs(const MolecularGeometry& molecule) {
        for(auto& motif : motifVector) {
            std::vector<bool> atCore;

            for(auto i : motif.electronIndices())
                atCore.emplace_back(coreElectronQ(molecule.electrons()[i], molecule.atoms()));

            if(std::all_of(atCore.begin(), atCore.end(), [](bool b){ return b; }))
                motif.setType(MotifType::Core);
            else if(std::none_of(atCore.begin(), atCore.end(), [](bool b){ return b; }))
                motif.setType(MotifType::Valence);
            else {
                YAML::Emitter out; out << molecule.electrons() <<  motif;
                spdlog::warn("Found a motif that is not clearly separable into Valence and Core. "
                             "ElectronsVector:\n{0} Motif:\n{1}", out.c_str());
            }
        }
    }

    void Motifs::splitCoreMotifs(const MolecularGeometry& molecule) {
        std::vector<Motif> newMotifVector{};

        for(const auto& m : motifVector) {
            if(m.type() == MotifType::Core){
                for (Eigen::Index k = 0; k < molecule.atoms().numberOfEntities(); ++k) {
                    auto atCoreK = molecule.coreElectronsIndices(k);

                    std::list<Eigen::Index> intersection;
                    std::set_intersection(
                            m.electronIndices().begin(), m.electronIndices().end(),
                            atCoreK.begin(), atCoreK.end(), std::back_inserter(intersection));
                    assert(intersection.size() <= 2
                    && "The number of electrons at a nuclear cusp must be less than 2 due to the antisymmetry principle.");

                    newMotifVector.emplace_back(Motif(intersection, MotifType::Core));
                }
            } else
                newMotifVector.emplace_back(m);
        }
        motifVector = newMotifVector;
    }

    void Motifs::sort(){
        std::sort(std::begin(motifVector), std::end(motifVector),
                  [] (const auto& lhs, const auto& rhs) {
                      return lhs.type() < rhs.type();
                  });
    }

}


namespace YAML {
    Node convert<MotifAnalysis::Motif>::encode(const MotifAnalysis::Motif &rhs) {
        Node node;
        node["Type"] = MotifAnalysis::toString(rhs.type());
        node["Indices"] = rhs.electronIndices();
        return node;
    }

    bool convert<MotifAnalysis::Motif>::decode(const Node &node, MotifAnalysis::Motif &rhs) {
        if (!node.IsSequence())
            return false;

        rhs.setType(MotifAnalysis::fromString(node["Type"].as<std::string>()));
        rhs.setElectronIndices(node["Indices"].as<std::list<Eigen::Index>>());
        return true;
    }

    Emitter &operator<<(Emitter &out, const MotifAnalysis::Motif &rhs) {
        out << YAML::Flow << BeginMap
        << Key << "Type" << Value << MotifAnalysis::toString(rhs.type())
        << Key << "Indices" << Value << BeginSeq;
        for(auto i : rhs.electronIndices())
            out <<  i;
        auto copy = rhs;
        out << EndSeq << EndMap;
        return out;
    };
}