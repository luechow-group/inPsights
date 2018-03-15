//
// Created by Michael Heuer on 13.03.18.
//

#ifndef AMOLQCPP_INTEGRANDEXAMPLEFUNCTION_H
#define AMOLQCPP_INTEGRANDEXAMPLEFUNCTION_H

class IntegrandExampleFunctor
{
public:
    IntegrandExampleFunctor(const double alpha)
            : m_alpha(alpha) {
        assert(alpha>0);
    }

    double operator()(const double x) const {
        assert(x>0);
        return log(m_alpha*x) / sqrt(x);
    }

    void setAlpha(const double alpha) {
        m_alpha = alpha;
    }

private:
    double m_alpha;
};


#endif //AMOLQCPP_INTEGRANDEXAMPLEFUNCTION_H
