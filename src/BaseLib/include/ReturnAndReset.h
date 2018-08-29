//
// Created by Michael Heuer on 30.08.18.
//

#ifndef AMOLQCPP_RETURNANDRESET_H
#define AMOLQCPP_RETURNANDRESET_H

template <typename ObjectType, typename ReturnType>
struct RETURN_AND_RESET{
    RETURN_AND_RESET(ObjectType& obj, ReturnType vec) : obj_(obj), vec_(std::move(vec)){};

    ReturnType returnAndReset(){ return vec_; };

    ~RETURN_AND_RESET(){
        obj_.resetRefToAll();
    }
    ObjectType& obj_;
    ReturnType vec_;
};

#endif //AMOLQCPP_RETURNANDRESET_H
