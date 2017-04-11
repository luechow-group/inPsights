//
// Created by heuer on 10.04.17.
//

#ifndef AMOLQCGUI_PROBLEMOBSERVER_H
#define AMOLQCGUI_PROBLEMOBSERVER_H

class ProblemObserver{
public:
    ProblemObserver(){};
    virtual ~ProblemObserver(){};
    virtual void stepPerformed(){};
};

class ObservableProblem{
public:
    void addObserver(ProblemObserver *observer) {
        std::cout << "Obersver added: " << observer << std::endl;
        observers_.push_back(observer);
    }

    void removeObserver(ProblemObserver *observer) {
        observers_.erase(std::remove(observers_.begin(), observers_.end(), observer), observers_.end());
    }

    void notifyObserversAboutPerformedStep() {
        std::cout << "step performed" << std::endl;
        for (auto &observer: observers_) observer->stepPerformed();
    }

    unsigned long getObserverCount(){
        return observers_.size();
    }

private:
    std::vector<ProblemObserver*> observers_;
};

#endif //AMOLQCGUI_PROBLEMOBSERVER_H
