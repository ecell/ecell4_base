#!/usr/bin/env python


from matplotlib.pylab import *


D=[0.25,0.5,1,2,4]

data_0=[
0.862862179761,
0.753491410993,
0.60815952851,
0.445889601095,
0.305622223357
]

data_em6=[
0.788132485945,
0.663137430836,
0.507006568658,
0.335227272727,
0.208129584352
]

data_em5=[
0.662728062554,
0.521941198262,
0.360172138544,
0.218285548167,
0.12251119726
]

data_em4=[
0.44779842834,
0.31161462002,
0.19177344261,
0.111646217941,
0.0709871388007
]


data_em3=[
0.211627566698,
0.135715851787,
0.0858701895462,
0.0565463070573,
0.0380785002929
]

data_em2=[
0.0804020100503,
0.0543882484195,
0.0390695283383,
0.0297150280101,
0.0308418302857
]


data_em1=[
0.0304764199241,
0.0269461077844,
0.0247045664543,
0.0233066617045,
0.0225048923679
]


axes([.13,.13,.8,.8])

semilogx(D, data_0, 'o:', label=r'$\tau_{\rm rel}=0$')
semilogx(D, data_em6, 'o:', label=r'$\tau_{\rm rel}=1 \ {\rm \mu s}$')
semilogx(D, data_em5, 'o:', label=r'$\tau_{\rm rel}=10 \ {\rm \mu s}$')
semilogx(D, data_em4, 'o:', label=r'$\tau_{\rm rel}=100 \ {\rm \mu s}$')
semilogx(D, data_em3, 'o:', label=r'$\tau_{\rm rel}=1 \ {\rm ms}$')
semilogx(D, data_em2, 'o:', label=r'$\tau_{\rm rel}=10 \ {\rm ms}$')
semilogx(D, data_em1, 'o:', label=r'$\tau_{\rm rel}=100 \  {\rm ms}$')



xlim(0.2,5)
xticks(D,[str(i) for i in D],fontsize=20)
yticks(fontsize=20)

ylim(0,1)
xlabel(r'Diffusion speed [${\rm \mu m^2 / s }$]',fontsize=24)
#ylabel('ratio',fontsize=20)
#legend()


leg = legend(loc=1,
              shadow=True,
              pad=0.05
              )
#for l in leg.get_lines():
#    l.set_linewidth(1.5)  # the legend line width


show()
