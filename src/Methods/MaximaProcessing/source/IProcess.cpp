// Copyright (C) 2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include "IProcess.h"

IProcess::ProcessType IProcess::typeFromString(const std::string& processName) {
    // Clusterer Processes
    if(processName == "IdentityClusterer")
        return ProcessType::IdentityClusterer;
    else if(processName == "PreClusterer")
        return ProcessType::DistanceClusterer;
    else if(processName == "DensityBasedClusterer")
        return ProcessType::DensityBasedClusterer;
    else if(processName == "ReferencePositionsClusterer")
        return ProcessType::ReferencePositionsClusterer;
    else if(processName == "SOAPClusterer")
        return ProcessType::SOAPClusterer;

    // Unknown Process name
    else
        return ProcessType::invalid;
};

IClusterer::IClusterer(std::vector<Sample> & samples)
: samples_(samples) {}
