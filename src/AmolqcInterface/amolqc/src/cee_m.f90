module primitive_cee
    use statistics
    implicit none
    type cee
        real*8 :: kin
        real*8 :: ne
        real*8 :: ee
        real*8 :: nn
    end type cee
contains
   real*8 function Esum(p)
     type(cee) :: p
     Esum = p%kin + p%ne + p%ee + p%nn
   end function Esum
   real*8 function Esigma(p)
    type(cee) :: p
     Esigma = sqrt((p%kin)**2 + (p%ne)**2 + (p%ee)**2 )
   end function Esigma

end module primitive_cee

module weighted_cee
    use statistics
    use primitive_cee
    implicit none

    type w_cee
        private
        type(weightStat) :: kin
        type(weightStat) :: ne
        type(weightStat) :: ee
        real*8           :: nn
    end type w_cee
contains
    !--------------------!
    subroutine w_reset(w)
        !--------------------!
        type(w_cee) :: w
        call reset(w%kin)
        call reset(w%ne)
        call reset(w%ee)
        w%nn = 0
    end subroutine w_reset
    !-------------------------!
    subroutine w_add(wce,pce,w)
        !-------------------------!
        type(cee) :: pce
        type(w_cee) :: wce
        real*8      :: w

        call addData(wce%kin,pce%kin,w)
        call addData(wce%ne,pce%ne,w)
        call addData(wce%ee,pce%ee,w)
        wce%nn = pce%nn

    end subroutine w_add
    !----------------!
    function w_mean(w)
        !----------------!
        type(w_cee) :: w
        type(cee) :: w_mean

        w_mean%kin = mean(w%kin)
        w_mean%ne = mean(w%ne)
        w_mean%ee = mean(w%ee)
        w_mean%nn = w%nn

    end function w_mean
    !------------------------!
    function w_meanAllNodes(w)
        !------------------------!
        type(w_cee) :: w
        type(cee) :: w_meanAllNodes

        w_meanAllNodes%kin = meanAllNodes(w%kin)
        w_meanAllNodes%ne  = meanAllNodes(w%ne)
        w_meanAllNodes%ee  = meanAllNodes(w%ee)
        w_meanAllNodes%nn  = w%nn
    end function w_meanAllNodes

    !----------------!
    function w_stdDevmean(w)
        !----------------!
        type(w_cee) :: w
        type(cee) :: w_stdDevmean

        w_stdDevmean%kin = stdDevmean(w%kin)
        w_stdDevmean%ne = stdDevmean(w%ne)
        w_stdDevmean%ee = stdDevmean(w%ee)
        w_stdDevmean%nn = 0
    end function w_stdDevmean
    !------------------------!
    function w_stdDevmeanAllNodes(w)
        !------------------------!
        type(w_cee) :: w
        type(cee) :: w_stdDevmeanAllNodes

        w_stdDevmeanAllNodes%kin = stdDevmeanAllNodes(w%kin)
        w_stdDevmeanAllNodes%ne = stdDevmeanAllNodes(w%ne)
        w_stdDevmeanAllNodes%ee = stdDevmeanAllNodes(w%ee)
        w_stdDevmeanAllNodes%nn = 0
    end function w_stdDevmeanAllNodes
    !----------------!
    function w_variance(w)
        !----------------!
        type(w_cee) :: w
        type(cee) :: w_variance

        w_variance%kin = variance(w%kin)
        w_variance%ne = variance(w%ne)
        w_variance%ee = variance(w%ee)
        w_variance%nn = 0
    end function w_variance
    !------------------------!
    function w_varianceAllNodes(w)
        !------------------------!
        type(w_cee) :: w
        type(cee) :: w_varianceAllNodes

        w_varianceAllNodes%kin = varianceAllNodes(w%kin)
        w_varianceAllNodes%ne = varianceAllNodes(w%ne)
        w_varianceAllNodes%ee = varianceAllNodes(w%ee)
        w_varianceAllNodes%nn = 0
    end function w_varianceAllNodes

end module weighted_cee

module simple_cee
    use statistics
    use primitive_cee
    implicit none
    type s_cee
        private
        type(simpleStat) :: kin
        type(simpleStat) :: ne
        type(simpleStat) :: ee
        real*8           :: nn
    end type s_cee
contains
    !---------------------!
    subroutine s_reset(c)
        !---------------------!
        type(s_cee) :: c

        call reset(c%kin)
        call reset(c%ne)
        call reset(c%ee)
        c%nn = 0

    end subroutine s_reset
    !-------------------------!
    subroutine s_add(sce,pce)
        !-------------------------!
        type(cee) :: pce
        type(s_cee) :: sce

        call addData(sce%kin,pce%kin)
        call addData(sce%ne,pce%ne)
        call addData(sce%ee,pce%ee)
        sce%nn = pce%nn
    end subroutine s_add
    !------------------!
    function s_mean(s)
        !------------------!
        type(s_cee) :: s
        type(cee) :: s_mean

        s_mean%kin = mean(s%kin)
        s_mean%ne = mean(s%ne)
        s_mean%ee = mean(s%ee)
        s_mean%nn = s%nn
    end function s_mean
    !--------------------------!
    function s_meanAllNodes(s)
        !--------------------------!
        type(s_cee) :: s
        type(cee) :: s_meanAllNodes

        s_meanAllNodes%kin = meanAllNodes(s%kin)
        s_meanAllNodes%ne  = meanAllNodes(s%ne)
        s_meanAllNodes%ee  = meanAllNodes(s%ee)
        s_meanAllNodes%nn  = s%nn
    end function s_meanAllNodes

    !----------------!
    function s_stdDevmean(s)
        !----------------!
        type(s_cee) :: s
        type(cee) :: s_stdDevmean

        s_stdDevmean%kin = stdDevmean(s%kin)
        s_stdDevmean%ne = stdDevmean(s%ne)
        s_stdDevmean%ee = stdDevmean(s%ee)
        s_stdDevmean%nn = 0

    end function s_stdDevmean
    !------------------------!
    function s_stdDevmeanAllNodes(s)
        !------------------------!
        type(s_cee) :: s
        type(cee) :: s_stdDevmeanAllNodes

        s_stdDevmeanAllNodes%kin = stdDevmeanAllNodes(s%kin)
        s_stdDevmeanAllNodes%ne = stdDevmeanAllNodes(s%ne)
        s_stdDevmeanAllNodes%ee = stdDevmeanAllNodes(s%ee)
        s_stdDevmeanAllNodes%nn = 0
    end function s_stdDevmeanAllNodes
    !----------------!
    function s_variance(s)
        !----------------!
        type(s_cee) :: s
        type(cee) :: s_variance

        s_variance%kin = variance(s%kin)
        s_variance%ne = variance(s%ne)
        s_variance%ee = variance(s%ee)
        s_variance%nn = 0
    end function s_variance
    !------------------------!
    function s_varianceAllNodes(s)
        !------------------------!
        type(s_cee) :: s
        type(cee) :: s_varianceAllNodes

        s_varianceAllNodes%kin = varianceAllNodes(s%kin)
        s_varianceAllNodes%ne = varianceAllNodes(s%ne)
        s_varianceAllNodes%ee = varianceAllNodes(s%ee)
        s_varianceAllNodes%nn = 0
    end function s_varianceAllNodes

end module simple_cee

module ceestat
    use simple_cee
    use weighted_cee
    implicit none

    interface addData
        module procedure s_add, w_add
    end interface

    interface reset
        module procedure s_reset, w_reset
    end interface

    interface mean
        module procedure s_mean, w_mean
    end interface

    interface meanAllNodes
        module procedure s_meanAllNodes, w_meanAllNodes
    end interface

    interface stdDevmean
        module procedure s_stdDevmean, w_stdDevmean
    end interface

    interface stdDevmeanAllNodes
        module procedure s_stdDevmeanAllNodes, w_stdDevmeanAllNodes
    end interface

    interface variance
        module procedure s_variance, w_variance
    end interface

    interface varianceAllNodes
        module procedure s_varianceAllNodes, w_varianceAllNodes
    end interface

end module ceestat


