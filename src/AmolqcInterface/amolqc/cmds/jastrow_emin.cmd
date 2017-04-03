x$begin_loop(count=5)
$sample(remove_outliers,tol=5.0,no_replace)
x$optimize_parameters(jastrow,method=newton,eq_iter=0,optmode=1)
$sample(remove_outliers,no_replace)
$sample(change_size,new_size=init_size)
$qmc(vmc,move=umr,steps=10,blocks=3,persist=9,discard_all)
x$qmc(vmc,move=umr,steps=10,blocks=20,discard_all)
$end_loop()
!

