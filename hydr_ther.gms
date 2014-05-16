$ontext
一个包含4水电一火电的系统，考虑24时段下的火电耗煤最小
使用表格为hydr_ther.xlsx
$offtext
set
i        hydro plant     /1*4/
j        thermal plant   /1/
pa_hy    hydro parameter /c1*c6/
pa_hys                   /up_pl,dlay/
pa_th    thermal parameter /d1*d3/
st_hy    st of hydro      /vmin,vmax,vini,vend,qmin,qmax,pmin,pmax/
t        time horizon    /1*24/;

$include 1.inc
$ontext
table par_hydro(i,pa_hy) parameter of hydro plant
table par_ther(j,pa_th)  parameter of thermal plant
table st_hydro(i,st_hy)  STs of hydro plant
table inflows(t,i)       inflow for 24 hours
parameter load(t)        load for 24 hours
table flow_delay(i,pa_hys)
$offtext
positive variable
         q(i,t),p_f_t(j,t);
variable
         v(i,t),
         in_delay(i,t)
         p_hy_t(i,t),P_hys_t(t),p_hys_ts
         p_fs_t(t),p_f_ts(j)
         o_f_t(j,t),o_f_ts(j),o_fs_ts,q_a;

q.up(i,t)=st_hydro(i,'qmax');
q.lo(i,t)=st_hydro(i,'qmin');
v.up(i,t)=st_hydro(i,'vmax');
v.lo(i,t)=st_hydro(i,'vmin');
v.fx(i,'24')=st_hydro(i,'vend');
p_hy_t.lo(i,t)=st_hydro(i,'pmin');
p_hy_t.up(i,t)=st_hydro(i,'pmax');
p_f_t.lo(j,t)=500;
p_f_t.up(j,t)=2500;
binary variable
         bi(i,t);
integer variable
         b_a;
equation
Cal_p_hy_t_1(i,t)
Cal_p_hy_t(i,t)
Cal_p_hys_t(t)  per hydro plant per hour output
Cal_p_hys_ts    all hydro plant 24 hour output
Cal_delay(i,t)  calculate the inflow from upstream considering the time delay
Cal_delay3(i,t)
Cal_delay4(i,t)
Cal_v_ini(i,t)
Cal_v(i,t)      balance of v considering the delay of water flow

Cal_p_fs_t(t)   all thermal plant per hour output
Cal_P_f_ts(j)   for load balance
St_p(t)         balance of p
Cal_o_f_t(j,t)   per thermal plant 24 hours output
Obj_o_fs_ts     objs
q_all
bin_all
;

*Cal_p_hy_t_1(i,t)$(ord(t) eq 1).. p_hy_t(i,t)=e=par_hydro(i,'c1')*power(st_hydro(i,'vini'),2)+par_hydro(i,'c2')*power(q(i,t),2)+par_hydro(i,'c3')*st_hydro(i,'vini')*q(i,t)+par_hydro(i,'c4')*st_hydro(i,'vini')+par_hydro(i,'c5')*q(i,t)+par_hydro(i,'c6');;
*Cal_p_hy_t(i,t)$(ord(t) ne 1)..   p_hy_t(i,t)=e=par_hydro(i,'c1')*power(v(i,t-1),2)+par_hydro(i,'c2')*power(q(i,t),2)+par_hydro(i,'c3')*v(i,t-1)*q(i,t)+par_hydro(i,'c4')*v(i,t-1)+par_hydro(i,'c5')*q(i,t)+par_hydro(i,'c6');;


Cal_p_hy_t_1(i,t)$(ord(t) eq 1).. p_hy_t(i,t)=e=bi(i,t)*(par_hydro(i,'c1')*power(st_hydro(i,'vini'),2)+par_hydro(i,'c2')*power(q(i,t),2)+par_hydro(i,'c3')*st_hydro(i,'vini')*q(i,t)+par_hydro(i,'c4')*st_hydro(i,'vini')+par_hydro(i,'c5')*q(i,t)+par_hydro(i,'c6'));
Cal_p_hy_t(i,t)..   p_hy_t(i,t)=e=bi(i,t)*(par_hydro(i,'c1')*power(v(i,t),2)+par_hydro(i,'c2')*power(q(i,t),2)+par_hydro(i,'c3')*v(i,t)*q(i,t)+par_hydro(i,'c4')*v(i,t)+par_hydro(i,'c5')*q(i,t)+par_hydro(i,'c6'));

Cal_p_hys_t(t).. p_hys_t(t)=e=sum(i,p_hy_t(i,t));
Cal_p_hys_ts..   p_hys_ts=e=sum(t,p_hys_t(t));

Cal_p_fs_t(t)..  p_fs_t(t)=e=sum(j,p_f_t(j,t));
Cal_P_f_ts(j)..  p_f_ts(j)=e=sum(t,p_f_t(j,t));

*Cal_o_f_t(j,t)..  o_f_t(j,t)=e=par_ther(j,'d1')*power(p_f_t(j,t),2)+par_ther(j,'d2')*p_f_t(j,t)+par_ther(j,'d3');
*Cal_o_f_t(j,t)..  o_f_t(j,t)=e=sum(pa_th,par_ther(j,pa_th)*power(p_f_t(j,t),3-ord(pa_th)))+700*sin(0.085*(500-p_f_t(j,t)));
Cal_o_f_t(j,t)..  o_f_t(j,t)=e=sum(pa_th,par_ther(j,pa_th)*power(p_f_t(j,t),3-ord(pa_th)));
Obj_o_fs_ts..    o_fs_ts=e=sum((j,t),o_f_t(j,t));

Cal_delay(i,t)$(ord(i) <3)..    in_delay(i,t)=e=0;
Cal_delay3(i,t)$(ord(i) = 3)..  in_delay(i,t)=e=bi(i,t)*(q('2',t-3)$(ord(t) >3)+bi(i,t)*q('1',t-2))$(ord(t) >2)+0$(ord(t) = 1 or ord(t) =2);
Cal_delay4(i,t)$(ord(i) = 4)..  in_delay(i,t)=e=bi(i,t)*q('3',t-4)$(ord(t) >4)+0$(ord(t) <= 4);
Cal_v_ini(i,t)$(ord(t) eq 1)..  v(i,t)=e=st_hydro(i,'vini')-bi(i,t)*q(i,t)+inflows(t,i)+in_delay(i,t);
Cal_v(i,t)$(ord(t) ne 1)..      v(i,t)=e=v(i,t-1)-bi(i,t)*q(i,t)+inflows(t,i)+in_delay(i,t);
St_P(t)..                       load(t)=e=p_fs_t(t)+p_hys_t(t);

q_all..  q_a=e=sum((i,t),bi(i,t)*q(i,t));
bin_all.. b_a=e=sum((i,t),bi(i,t));
model mincost /bin_all,q_all,Cal_p_hy_t_1,Cal_p_hy_t,Cal_p_hys_t,Cal_p_hys_ts,Cal_delay,Cal_delay3,Cal_delay4,Cal_v_ini,Cal_v,St_p,Cal_p_fs_t,Cal_P_f_ts,Cal_o_f_t,Obj_o_fs_ts/;
$call gams hydr_ther.gms lo=%GAMS.lo% pf4=0
option minlp=bonmin;
solve mincost using minlp minimizing o_fs_ts;
display load,q.l,v.l,p_hy_t.l,p_hys_t.l,p_hys_ts.l,p_f_t.l,p_f_ts.l,p_fs_t.l,o_fs_ts.l,q_a.l,bi.l,b_a.l;
