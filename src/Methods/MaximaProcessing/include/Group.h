//
// Created by heuer on 11.04.19.
//

#ifndef INPSIGHTS_GROUP_H
#define INPSIGHTS_GROUP_H

#include <vector>
#include <Eigen/Core>
#include <memory>
#include <Enumerate.h>
#include <Sample.h>

// In order to use Group class,
// <Reference.h> has to be included as well, due to forward declaration

class Group;
class Reference;

class Group : public std::vector<Group> {
public:
    Group() = default;
    Group(const Group& group) = default;

    explicit Group(std::vector<Group>::size_type size);
    Group(std::initializer_list<Group> group);
    explicit Group(Reference reference);

    bool isLeaf() const;

    Group::size_type numberOfLeaves() const;

    void sort();
    void sortAll();

    Group& operator+= (const Group& other);

    void makeSubgroup(std::vector<Group::iterator> its);

    void permuteAll(const Eigen::PermutationMatrix<Eigen::Dynamic>& perm, std::vector<Sample>& samples);

    std::shared_ptr<const Reference> representative() const;
    std::shared_ptr<Reference> representative();

    bool operator<(const Group& other) const;

    std::vector<size_t > allSampleIds() const;

    friend std::ostream& operator<<(std::ostream& os, const Group & g);


private:
    virtual void updateRepresentative();
    std::shared_ptr<Reference> representative_;
};

inline Group operator+ (Group lhs, const Group& rhs) {
    lhs += rhs;
    return lhs;
}


#endif //INPSIGHTS_GROUP_H
