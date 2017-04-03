//
// Created by Michael Heuer on 17.06.16.
//

#ifndef RTQC_BSPLINEKNOTINSERTER_H
#define RTQC_BSPLINEKNOTINSERTER_H

class BSpline;

/*! Implementation of the knot insertion algorithm by Böhm.
 * doi: 10.1016/0010-4485(80)90154-2
 * */
class BSplineKnotInserter{
public:
  BSplineKnotInserter();

  void insertKnotByReference(const double uInsert, BSpline &bs) const;
  BSpline insertKnotByCopy(const double u, const BSpline &bs) const;
};

#endif //RTQC_BSPLINEKNOTINSERTER_H
