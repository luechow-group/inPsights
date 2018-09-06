//
// Created by Michael Heuer on 30.08.18.
//

#ifndef AMOLQCPP_RETURNANDRESET_H
#define AMOLQCPP_RETURNANDRESET_H

#include <utility>


/*
 * Automatic: waits for the next resetStrategy() call and resets to all
 * OnFinished: on the next resetStrategy() call, the ref is reset to all if usage is Finished
 * Manual: the user manually needs to call resetRef()
 */


enum class Reset{Automatic, Manual, OnFinished};
enum class Usage{Standard, NotFinished, Finished};


template <typename ObjectType, typename ReturnType>
struct RETURN_AND_RESET{
    RETURN_AND_RESET(ObjectType& obj, ReturnType vec) : obj_(obj), vec_(std::move(vec)){};

    ReturnType returnAndReset(){ return vec_; };

    ~RETURN_AND_RESET(){
        obj_.resetRef();
    }
    ObjectType& obj_;
    ReturnType vec_;
};

#endif //AMOLQCPP_RETURNANDRESET_H
