$begin_subroutine(name=eq1)
$sample(change_size,new_size=init_size,last)
$qmc(vmc,steps=1500,block_len=500,accept_ratio=0.5,discard_all)
x$qmc(vmc,accumulate,step_stride=10,acc_size=1k,block_len=200,
accept_ratio=0.5,discard=1,move=umr)
$sample(remove_outliers,no_replace)
$end_subroutine()


$begin_subroutine(name=eq2)
$sample(change_size,new_size=init_size,last)
$qmc(vmc,steps=1500,block_len=500,accept_ratio=0.5,discard_all)
x$qmc(vmc,accumulate,step_stride=10,acc_size=1k,block_len=200,
accept_ratio=0.5,discard=1,move=umr)
$save_result()
$sample(remove_outliers,no_replace)
$end_subroutine()


$sample(change_size,new_size=1,last)
x$qmc(vmc,accumulate,step_stride=10,acc_size=1k,block_len=200,
accept_ratio=0.5,discard=1,move=umr)
$sample(remove_outliers,no_replace)
x$optimize_parameters(params=jastrow,variance_min,method=varmin,E_ref=-133.98,E_ref_adp,
NL2SOL_D_mode=1,max_iter=1,eq_iter=5,eq_call=eq1)
$sample(change_size,new_size=1,last)
x$qmc(vmc,accumulate,step_stride=10,acc_size=1k,block_len=200,
accept_ratio=0.5,discard=1,move=umr)
x$optimize_parameters(params=jastrow,variance_min,method=varmin,E_ref=-133.98,E_ref_adp,
NL2SOL_D_mode=0,max_iter=1,eq_iter=5,write_wf,eq_call=eq2)
$sample(change_size,new_size=1,last)
$qmc(vmc,steps=1500,block_len=500,accept_ratio=0.5,discard_all)
x$qmc(vmc,accumulate,step_stride=10,acc_size=1k,block_len=200,
accept_ratio=0.5,discard=1,move=umr)
$save_result()
$print_results()
