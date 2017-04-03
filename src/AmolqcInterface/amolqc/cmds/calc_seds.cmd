$sample(create,start=density,generate=random,size=200)
$sample(remove_outliers)
$qmc(vmc,move=umr,steps=10,blocks=5,discard_all,tau=0.020,accept_ratio=0.5,persist=9)
$qmc(vmc,move=umr,steps=200,discard=3,blocks=200,accept_ratio=0.5,stddev=0.0002,sed,epart,readref,ref=1,grid=4d0)

