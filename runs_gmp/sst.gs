'reinit'
pat='/rhome/soo/NEMO/couple/GRIMs/coupled_GRIMs/gmp_v3.2/runs_gmp'
'open 'pat'/test_2020_yng/fgb.ctl'


lon1=160
lon2=200 
lat1=-10
lat2=10
time1='00Z20MAY2020'
time2='00Z31AUG2020'

'set background 1'
'c'
'set grads off'
'set hershey off'

'set x 1';'set y 1';'set z 1';'set time 'time1' 'time2''
'define new = aave(tmpsfc,lon='lon1',lon='lon2',lat='lat1',lat='lat2')'

'set time 'time1' 'time2''
'set vrange 29 31.0'
'set ylint 0.2'
'd new-273.15'
