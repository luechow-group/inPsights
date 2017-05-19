//
// Created by Moria on 19.05.2017.
//

#include "SpinDeterminer.h"

SpinDeterminer::SpinDeterminer(int highestAlphaIndex): highestAlphaIndex(highestAlphaIndex) {

}

SpinDeterminer::~SpinDeterminer() {

}

spintype SpinDeterminer::determineSpin(int index) {
    return index<=highestAlphaIndex?spintype::SPIN_ALPHA:spintype::SPIN_BETA;
}
